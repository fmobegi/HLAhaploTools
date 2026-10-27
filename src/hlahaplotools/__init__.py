"""HLA haplotype frequency estimation: a Rust port of Hapl-o-Mat with
handling of untyped and half-typed loci, and of family data.

Typical use:

    import hlahaplotools as h

    ds = h.read_mac("population.dat", mac_codes="data/mac/alpha.v3.txt")
    nom = h.Nomenclature("data/derived/3.66.0")
    est = h.estimate(ds, nom, {"A": "g", "B": "g", "DRB1": "g"})
    est.frequencies()[:5]

Families (FAMILY_ID, Family_Member = F/M/C1.., A_1, A_2, ...):

    fam = h.read_family_table("families.csv")
    h.check_mendelian(fam, nom, "g")
    est, phases = h.estimate_families(fam, nom, "g")
"""

from hlahaplotools._core import (
    Dataset,
    Estimate,
    FamilyData,
    FamilyPhase,
    MacCodes,
    Nomenclature,
    audit_families,
    check_mendelian,
    classify_typing,
    compare_to_families,
    estimate,
    estimate_families,
    phase_families,
    read_family_table,
    read_gls,
    read_glsc,
    read_hml,
    read_mac,
    read_ngsengine,
    sniff_xml,
)
from hlahaplotools.cli import main
from hlahaplotools.multistart import estimate_multistart

__all__ = [
    "Dataset",
    "Estimate",
    "FamilyData",
    "FamilyPhase",
    "MacCodes",
    "Nomenclature",
    "audit_families",
    "check_mendelian",
    "classify_typing",
    "compare_to_families",
    "estimate",
    "estimate_families",
    "estimate_multistart",
    "main",
    "phase_families",
    "read_family_table",
    "read_glsc",
    "read_gls",
    "read_hml",
    "read_mac",
    "read_ngsengine",
    "sniff_xml",
]
