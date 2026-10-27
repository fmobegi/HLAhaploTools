"""Hardy–Weinberg tests and linkage disequilibrium."""

import json
import subprocess
import sys

from conftest import DERIVED, ROOT, needs_derived

import hlahaplotools as h


def write(tmp_path, text, name="pop.dat"):
    p = tmp_path / name
    p.write_text(text)
    return p


@needs_derived
def test_hardy_weinberg_flags_missing_heterozygotes(tmp_path):
    rows = ["id\tA\tA\tB\tB"]
    # A: only homozygotes (far from HWE). B: 25/50/25 (exactly HWE).
    for i in range(100):
        a = "01:01\t01:01" if i < 50 else "02:01\t02:01"
        b = ["07:02\t07:02", "07:02\t08:01", "07:02\t08:01", "08:01\t08:01"][i % 4]
        rows.append(f"s{i}\t{a}\t{b}")
    ds = h.read_mac(str(write(tmp_path, "\n".join(rows) + "\n")))
    res = {r["locus"]: r for r in h.hardy_weinberg(ds, h.Nomenclature(str(DERIVED)), "2f")}
    assert res["A"]["observed_heterozygosity"] == 0 and res["A"]["p_value"] < 1e-3 and res["A"]["f"] > 0.99
    assert abs(res["B"]["f"]) < 0.02 and res["B"]["p_value"] > 0.5
    assert res["A"]["n"] == 100 and res["A"]["alleles"] == 2


@needs_derived
def test_untyped_half_typed_and_blank_homozygotes_are_counted_apart(tmp_path):
    p = write(tmp_path, "id\tA\tA\n1\t01:01\t02:01\n2\t01:01\t\n3\t\t\n4\t02:01\t02:01\n")
    nom = h.Nomenclature(str(DERIVED))
    r = h.hardy_weinberg(h.read_mac(str(p)), nom, "2f")[0]
    assert (r["n"], r["untyped"], r["homozygous_filled"]) == (3, 1, 1)
    r = h.hardy_weinberg(h.read_mac(str(p), half_typed="unknown"), nom, "2f")[0]
    assert (r["n"], r["half_typed"], r["homozygous_filled"]) == (2, 1, 0)


@needs_derived
def test_linkage_disequilibrium_on_linked_and_unlinked_loci(tmp_path):
    rows = ["id\tA\tA\tB\tB\tC\tC"]
    for i in range(400):
        h1 = ("01:01", "08:01") if i % 2 == 0 else ("02:01", "07:02")
        h2 = ("01:01", "08:01") if (i // 2) % 2 == 0 else ("02:01", "07:02")
        c = ("03:04", "07:01")[(i // 4) % 2], ("03:04", "07:01")[(i // 8) % 2]
        rows.append(f"s{i}\t{h1[0]}\t{h2[0]}\t{h1[1]}\t{h2[1]}\t{c[0]}\t{c[1]}")
    ds = h.read_mac(str(write(tmp_path, "\n".join(rows) + "\n")))
    nom = h.Nomenclature(str(DERIVED))
    ld = {
        (r["locus_a"], r["locus_b"]): r for r in h.linkage_disequilibrium(ds, nom, "2f", epsilon=1e-10, permutations=99)
    }
    assert set(ld) == {("A", "B"), ("A", "C"), ("B", "C")}
    ab, ac = ld[("A", "B")], ld[("A", "C")]
    assert ab["d_prime"] > 0.999 and ab["wn"] > 0.999 and ab["p_value"] == 0.01
    assert ac["d_prime"] < 0.1 and ac["p_value"] > 0.05
    top = ab["pairs"][0]
    assert {top["allele_a"], top["allele_b"]} <= {"A*01:01", "A*02:01", "B*08:01", "B*07:02"}
    only = h.linkage_disequilibrium(ds, nom, "2f", pairs=[("B", "C")])
    assert [(r["locus_a"], r["locus_b"]) for r in only] == [("B", "C")] and only[0]["p_value"] is None


@needs_derived
def test_stats_command_writes_tables(tmp_path):
    out = tmp_path / "st"
    run = subprocess.run(
        [
            sys.executable, "-m", "hlahaplotools", "stats", str(ROOT / "examples" / "population" / "population.dat"),
            "--loci", "A:g,C:g,B:g", "-o", str(out), "--nomenclature", str(DERIVED), "--permutations", "500",
            "--pairs", "adjacent", "-q",
        ],
        capture_output=True,
        text=True,
    )  # fmt: skip
    assert run.returncode == 0, run.stderr
    hwe = (out / "hwe.tsv").read_text().splitlines()
    assert hwe[0].startswith("locus\tn\talleles") and len(hwe) == 4
    ld = (out / "ld.tsv").read_text().splitlines()
    assert len(ld) == 3 and ld[1].startswith("A\tC\t") and ld[2].startswith("C\tB\t")
    assert (out / "ld_alleles.tsv").read_text().count("\n") > 10
    assert hwe[0].split("\t")[7] == "q_value" and ld[0].split("\t")[-1] == "q_value"
    assert json.loads((out / "run.json").read_text())["parameters"]["permutations"] == 500


def test_adjust_fdr_matches_benjamini_hochberg():
    q = h.adjust_fdr([0.01, 0.04, None, 0.03, 0.2])
    # Sorted 0.01, 0.03, 0.04, 0.2 over m = 4: 0.04, 0.0533, 0.0533, 0.2 (monotone).
    assert q[2] is None and abs(q[0] - 0.04) < 1e-12 and abs(q[4] - 0.2) < 1e-12
    assert abs(q[1] - 0.16 / 3) < 1e-12 and abs(q[3] - 0.16 / 3) < 1e-12
    assert h.adjust_fdr([0.001] * 55) == [0.001] * 55
