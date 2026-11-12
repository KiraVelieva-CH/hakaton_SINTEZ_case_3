import requests
from core.exceptions import ResourceUnavailableError
import time

class VEPAnnotator:
    BASE = "https://rest.ensembl.org/vep/human/region"
    
    def __init__(self, config):
        self.config = config
        self.timeout = int(config.get("analysis", {}).get("vep_timeout_seconds", 60))
        self.batch_size = 100  

    def annotate(self, variants):
        if not self.config.get("analysis", {}).get("use_vep", False):
            return variants
        if not self.config.get("analysis", {}).get("allow_network", False):
            raise ResourceUnavailableError("VEP requested but analysis.allow_network is false")

        for i in range(0, len(variants), self.batch_size):
            batch = variants[i:i+self.batch_size]
            
            
            vcf_strings = [f"{v.chrom.removeprefix('chr')} {v.pos} . {v.ref} {v.alt}" for v in batch]
            
            payload = {"variants": vcf_strings}
            headers = {"Content-Type": "application/json", "Accept": "application/json"}
            
            try:
                r = requests.post(self.BASE, params={"content-type": "application/json"}, 
                                  json=payload, headers=headers, timeout=self.timeout)
                r.raise_for_status()
                results = r.json()
                
                
                res_dict = {f"{res['seq_region_name']}:{res['start']}:{res['allele_string'].split('/')[0]}:{res['allele_string'].split('/')[1]}": res 
                            for res in results if 'seq_region_name' in res}
                
                for v in batch:
                    
                    key = f"{v.chrom.removeprefix('chr')}:{v.pos}:{v.ref}:{v.alt}"
                    res = res_dict.get(key)
                    
                    if res:
                        v.annotations["vep"] = res
                        
                        
                        for tx in res.get("transcript_consequences", []):
                            if tx.get("gene_symbol"):
                                if not v.gene or tx.get("canonical") == 1:
                                    v.gene = tx["gene_symbol"]
                                    v.consequence = tx.get("consequence_terms", [None])[0]
                        
                        
                        af = None
                        for colo in res.get("colocated_variants", []):
                            if "frequencies" in colo:
                                gnomad = colo["frequencies"].get("gnomADe") or colo["frequencies"].get("gnomADg")
                                if gnomad and "af" in gnomad:
                                    af = gnomad["af"]
                                    break
                        
                        if af is not None:
                            v.population_af = float(af)
                            
            except requests.RequestException as e:
                print(f"VEP Batch Error: {e}")
            time.sleep(1) 

        return variants