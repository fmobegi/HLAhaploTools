"""Renaming alleles that IMGT/HLA has deleted or renamed."""

import subprocess
import sys

from hlahaplotools.rename import load_renames, rename_file


def imgt(tmp_path):
    d = tmp_path / "imgt"
    d.mkdir()
    (d / "Allelelist.txt").write_text(
        "# file: Allelelist.txt\nAlleleID,Allele\nHLA1,C*04:09L\nHLA2,B*47:01:01:03\nHLA3,A*01:01:01:01\n"
    )
    (d / "Deleted_alleles.txt").write_text(
        "# file: Deleted_alleles.txt\nAlleleID,Allele,Description\n"
        "HLA9,C*04:09N,Suffix changed becoming C*04:09L (December 2024)\n"
        "HLA8,B*47:01:01:01,Sequence shown to contain errors and be identical to B*47:01:01:03 (April 2018)\n"
        "HLA7,A*01:99,Sequence renamed A*01:98 (2010)\n"  # A*01:98 later renamed again
        "HLA6,A*01:98,Sequence named in error and renamed A*01:01:01:01 (2012)\n"
        "HLA5,B*99:01,Name never assigned\n"
    )
    return d


def test_load_renames_follows_chains_and_skips_unknown_targets(tmp_path):
    r = load_renames(imgt(tmp_path))
    assert r == {
        "C*04:09N": "C*04:09L",
        "B*47:01:01:01": "B*47:01:01:03",
        "A*01:99": "A*01:01:01:01",
        "A*01:98": "A*01:01:01:01",
    }


def test_rename_keeps_each_format(tmp_path):
    r = load_renames(imgt(tmp_path))
    mac = tmp_path / "x.dat"
    mac.write_text("id\tB\tB\tC\tC\ns1\t47:01:01:01\t08:01\t04:09N\t\ns2\t07:02\t08:01\t04:09N/04:01\t07:01\n")
    seen = rename_file(mac, tmp_path / "y.dat", "mac", r)
    assert (tmp_path / "y.dat").read_text() == (
        "id\tB\tB\tC\tC\ns1\t47:01:01:03\t08:01\t04:09L\t\ns2\t07:02\t08:01\t04:09L/04:01\t07:01\n"
    )
    assert seen[("C*04:09N", "C*04:09L")] == 2

    glc = tmp_path / "x.glc"
    glc.write_text("1 C*04:09N+C*07:01|C*04:01+C*07:01 B*08:01+\n")
    rename_file(glc, tmp_path / "y.glc", "glsc", r)
    assert (tmp_path / "y.glc").read_text() == "1 C*04:09L+C*07:01|C*04:01+C*07:01 B*08:01+\n"

    fam = tmp_path / "x.csv"
    fam.write_text('FAMILY_ID,Family_Member,C_1,C_2,B_1,B_2\nf1,F,C*04:09N,*07:01,"47:01:01:01",08:01\n')
    rename_file(fam, tmp_path / "y.csv", "families", r)
    assert (tmp_path / "y.csv").read_text() == (
        'FAMILY_ID,Family_Member,C_1,C_2,B_1,B_2\nf1,F,C*04:09L,*07:01,"47:01:01:03",08:01\n'
    )


def test_rename_command(tmp_path):
    mac = tmp_path / "x.dat"
    mac.write_text("id\tC\tC\ns1\t04:09N\t07:01\n")
    out = tmp_path / "y.dat"
    run = subprocess.run(
        [sys.executable, "-m", "hlahaplotools", "rename", str(mac), str(out), "--imgt-dir", str(imgt(tmp_path))],
        capture_output=True,
        text=True,
    )
    assert run.returncode == 0, run.stderr
    assert "C*04:09N -> C*04:09L: 1" in run.stdout
    assert out.read_text() == "id\tC\tC\ns1\t04:09L\t07:01\n"
    assert (tmp_path / "y.dat.renames.tsv").read_text() == "old\tnew\ttypings\nC*04:09N\tC*04:09L\t1\n"
