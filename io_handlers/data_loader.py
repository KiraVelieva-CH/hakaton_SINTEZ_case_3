from __future__ import annotations
from pathlib import Path
import csv
import re
import pandas as pd
from core.models import Phenotype, Variant
from core.exceptions import InvalidInputError, InvalidVCFError
import json 

def load_phenotypes(path: str | Path) -> list[Phenotype]:
    p = Path(path)
    if not p.exists():
        raise InvalidInputError(f"Phenotype file not found: {p}")
    
    if p.suffix.lower() == ".json":
        try:
            data = json.loads(p.read_text(encoding="utf-8"))
        except Exception as e:
            raise InvalidInputError(f"Invalid JSON in phenotype file: {e}")

        out = []
        
        features = data.get("phenotypicFeatures", [])
        
        
        if not features and isinstance(data, list):
            features = data

        for item in features:
            
            hpo_id = None
            term = ""
            
            if isinstance(item, dict):
                type_info = item.get("type", {})
                if isinstance(type_info, dict):
                    hpo_id = type_info.get("id")
                    term = type_info.get("label", "")
                elif isinstance(type_info, str): 
                    hpo_id = type_info
                
                
                if not hpo_id:
                    hpo_id = item.get("id") or item.get("hpo_id")
                if not term:
                    term = item.get("label") or item.get("term") or ""
                    
                
                excluded = item.get("excluded", False)
                present = not excluded
            elif isinstance(item, str): 
                hpo_id = item.strip() if item.startswith("HP:") else None
                term = "" if hpo_id else item.strip()
                present = True
            else:
                continue

            if hpo_id or term:
                out.append(Phenotype(
                    term=term,
                    hpo_id=hpo_id,
                    present=present,
                    source_text=str(item)
                ))
                
        if not out:
            raise InvalidInputError("No usable phenotypes found in JSON file.")
        return out

    text = p.read_text(encoding="utf-8")
    if p.suffix.lower() in {".csv", ".tsv"}:
        sep = "\t" if p.suffix.lower() == ".tsv" else ","
        df = pd.read_csv(p, sep=sep).fillna("")
        cols = {c.lower(): c for c in df.columns}
        term_col = cols.get("term") or cols.get("phenotype") or cols.get("symptom")
        if term_col:
            out = []
            for _, row in df.iterrows():
                term = str(row[term_col]).strip()
                if term:
                    hpo_col = cols.get("hpo_id") or cols.get("hpo")
                    out.append(Phenotype(term=term, hpo_id=str(row[hpo_col]).strip() or None if hpo_col else None))
            return out
            
    lines = [x.strip(" \t-•") for x in text.splitlines() if x.strip()]
    if len(lines) == 1 and len(lines[0]) > 180:
        return [Phenotype(term=lines[0], source_text=lines[0])]
    return [Phenotype(term=x, source_text=x) for x in lines]


def _normalize_chrom(chrom: str) -> str:
    c = str(chrom).strip()
    return c if c.lower().startswith("chr") else f"chr{c}"


def load_variant_table(path: str | Path) -> list[Variant]:
    p = Path(path)
    if not p.exists():
        raise InvalidInputError(f"Variant table not found: {p}")
    df = pd.read_csv(p, sep=None, engine="python", comment="#").fillna("")
    cols = {str(c).strip().lower(): c for c in df.columns}
    def get(row, *names, default=""):
        for n in names:
            if n.lower() in cols:
                return row[cols[n.lower()]]
        return default
    result = []
    for i, row in df.iterrows():
        chrom = get(row, "chrom", "chr", "chromosome")
        pos = get(row, "pos", "position")
        ref, alt = get(row, "ref"), get(row, "alt")
        if not chrom or not pos or not ref or not alt:
            continue
        chrom = _normalize_chrom(chrom)
        pos = int(float(pos))
        gene = str(get(row, "gene", "gene_symbol", "symbol")).strip() or None
        vid = str(get(row, "variant_id", "id")).strip() or f"{chrom}:{pos}:{ref}:{alt}"

        def num(*names):
            v = get(row, *names)
            if str(v).strip().upper() in {"NA", "N/A", "", "NONE", "."}:
                return None
            try: return float(v) 
            except (ValueError, TypeError): return None

        result.append(Variant(
            variant_id=vid, chrom=chrom, pos=pos, ref=str(ref), alt=str(alt),
            gene=gene, consequence=str(get(row, "consequence")).strip() or None,
            genotype=str(get(row, "genotype", "gt")).strip() or None,
            qual=num("qual", "quality"), depth=int(num("depth", "dp")) if num("depth", "dp") is not None else None,
            allele_balance=num("allele_balance", "ab"), population_af=num("population_af", "af", "gnomad_af"),
            clinvar_significance=str(get(row, "clinvar_significance", "clinical_significance")).strip() or None,
            clinvar_review_status=str(get(row, "clinvar_review_status", "review_status")).strip() or None
        ))
    if not result:
        raise InvalidInputError("No usable variants found in variant table")
    return result


def load_vcf(path: str | Path) -> list[Variant]:
    p = Path(path)
    if not p.exists():
        raise InvalidVCFError(f"VCF file not found: {p}")
    try:
        from cyvcf2 import VCF
        reader = VCF(str(p))
        out = []
        for rec in reader:
            for alt in rec.ALT:
                if alt is None or alt.startswith("<") or alt == "*":
                    continue
                gt = rec.genotypes[0] if rec.genotypes else None
                gt_text = None
                if gt:
                    gt_text = f"{gt[0]}{'|' if gt[2] else '/'}{gt[1]}"
                dp = rec.INFO.get("DP")
                out.append(Variant(
                    variant_id=f"{_normalize_chrom(rec.CHROM)}:{rec.POS}:{rec.REF}:{alt}",
                    chrom=_normalize_chrom(rec.CHROM), pos=rec.POS, ref=rec.REF, alt=alt,
                    genotype=gt_text, qual=float(rec.QUAL) if rec.QUAL is not None else None,
                    depth=int(dp) if dp is not None else None,
                    annotations={"filter": rec.FILTER}
                ))
        reader.close()
        if not out:
            raise InvalidVCFError("VCF contains no supported small variants")
        return out
    except ImportError:
        return _load_vcf_minimal(p)


def _load_vcf_minimal(path: Path) -> list[Variant]:
    out = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line or line.startswith("#"): continue
        parts = line.split("\t")
        if len(parts) < 8: raise InvalidVCFError("Malformed VCF row: expected >=8 columns")
        chrom, pos, _id, ref, alts, qual, filt, info = parts[:8]
        inf = dict(x.split("=", 1) for x in info.split(";") if "=" in x)
        sample = parts[9] if len(parts) > 9 else ""
        fmt = parts[8].split(":") if len(parts) > 8 else []
        vals = sample.split(":") if sample else []
        sample_map = dict(zip(fmt, vals))
        for alt in alts.split(","):
            if alt.startswith("<") or alt == "*": continue
            out.append(Variant(
                variant_id=f"{_normalize_chrom(chrom)}:{int(pos)}:{ref}:{alt}", chrom=_normalize_chrom(chrom),
                pos=int(pos), ref=ref, alt=alt, genotype=sample_map.get("GT"),
                qual=float(qual) if qual not in {".", ""} else None,
                depth=int(sample_map.get("DP", inf.get("DP", 0))) or None,
                annotations={"filter": filt}
            ))
    if not out: raise InvalidVCFError("VCF contains no supported small variants")
    return out
