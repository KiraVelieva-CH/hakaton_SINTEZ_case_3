from __future__ import annotations
from typing import Any, Optional
from pydantic import BaseModel, Field, ConfigDict


class Phenotype(BaseModel):
    model_config = ConfigDict(extra="allow")
    term: str
    hpo_id: Optional[str] = None
    present: bool = True
    source_text: Optional[str] = None


class Variant(BaseModel):
    model_config = ConfigDict(extra="allow")
    variant_id: str
    chrom: str
    pos: int = Field(gt=0)
    ref: str
    alt: str
    gene: Optional[str] = None
    consequence: Optional[str] = None
    genotype: Optional[str] = None
    qual: Optional[float] = None
    depth: Optional[int] = None
    allele_balance: Optional[float] = None
    population_af: Optional[float] = None
    clinvar_significance: Optional[str] = None
    clinvar_review_status: Optional[str] = None
    scores: dict[str, Optional[float]] = Field(default_factory=dict)
    filters: list[str] = Field(default_factory=list)
    annotations: dict[str, Any] = Field(default_factory=dict)


class Gene(BaseModel):
    model_config = ConfigDict(extra="allow")
    gene_id: str
    symbol: Optional[str] = None
    score: float = 0.0
    evidence: list[str] = Field(default_factory=list)
    variants: list[str] = Field(default_factory=list)
    diseases: list[str] = Field(default_factory=list)


class EvidenceItem(BaseModel):
    source: str
    statement: str
    strength: str = "supporting"
    details: dict[str, Any] = Field(default_factory=dict)


class Candidate(BaseModel):
    gene: str
    variant_id: Optional[str] = None
    score: float
    component_scores: dict[str, Optional[float]] = Field(default_factory=dict)
    evidence: list[EvidenceItem] = Field(default_factory=list)
    caveats: list[str] = Field(default_factory=list)


class AnalysisResult(BaseModel):
    case_id: str
    genome_build: str = "GRCh38"
    phenotypes: list[Phenotype] = Field(default_factory=list)
    variants: list[Variant] = Field(default_factory=list)
    candidates: list[Candidate] = Field(default_factory=list)
    gene_ranking: list[Gene] = Field(default_factory=list)
    disease_hypotheses: list[dict[str, Any]] = Field(default_factory=list)
    block_results: dict[str, Any] = Field(default_factory=dict)
    consistency: list[dict[str, Any]] = Field(default_factory=list)
    explanation: str = ""
    status: str = "completed"
    warnings: list[str] = Field(default_factory=list)
    metadata: dict[str, Any] = Field(default_factory=dict)
