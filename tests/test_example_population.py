"""The example population must read whole: 1,000 people, none unreadable, with DQB1 untyped in some.

Rows with an untyped last locus end in tab characters; a whitespace "fixer" that strips them
makes those rows too short to read (this happened once, so it is tested).
"""

from conftest import ROOT

import hlahaplotools as h

POPULATION = ROOT / "examples" / "population" / "population.dat"


def test_every_row_has_every_column():
    lines = POPULATION.read_text(encoding="utf-8").splitlines()
    widths = {len(line.split("\t")) for line in lines}
    assert widths == {11}, widths


def test_example_population_reads_whole():
    ds = h.read_mac(str(POPULATION))
    assert len(ds) == 1000
    assert not ds.issues
    assert ds.missing_counts()["DQB1"] > 0
