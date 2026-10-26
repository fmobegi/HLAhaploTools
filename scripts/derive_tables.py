"""Build the lookup tables used to resolve typings to g, G, P or n-field codes.

Input: one IMGT/HLA release directory (hla_nom_g.txt, hla_nom_p.txt).
Output, in Hapl-o-Mat's own file formats so both tools can share them:

    AllAllelesExpanded.txt  any typing prefix -> full allele names
    LargeG.txt              G group -> member alleles
    P.txt                   P group -> member alleles
    Smallg.txt              g group -> member alleles

The rules, including their odd corners, follow Hapl-o-Mat's prepareData
scripts, so that results can be compared with the original tool on equal
terms. `scripts/build_nomenclature.py` calls `derive_all`; the comparison
test is tests/test_derived_tables.py.
"""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path

EXPRESSION_LETTERS = ("N", "L", "S", "Q")


def _data_lines(path: Path) -> list[str]:
    with path.open(encoding="utf-8-sig") as f:
        return [line.rstrip("\r\n") for line in f if not line.startswith("#")]


def all_alleles(hla_nom_g: Path) -> list[str]:
    """Every allele named in hla_nom_g.txt, in file order (e.g. 'A*01:01:01:01')."""
    out = []
    for line in _data_lines(hla_nom_g):
        locus, alleles, _group = line.split(";")[:3]
        out.extend(locus + a for a in alleles.split("/"))
    return out


def expanded(alleles: list[str]) -> list[tuple[str, list[str]]]:
    """Map every field prefix of every allele to the alleles it covers.

    'A*01:01:01:02N' is listed under A*01, A*01N, A*01:01, A*01:01N, ...,
    A*01:01:01:02N and A*01:01:01:02.
    """
    table: dict[str, list[str]] = defaultdict(list)
    for code in alleles:
        prefix = ""
        for field in code.split(":"):
            prefix += field
            table[prefix].append(code)
            if code.endswith(EXPRESSION_LETTERS):
                if prefix.endswith(EXPRESSION_LETTERS):
                    table[prefix[:-1]].append(code)
                else:
                    table[prefix + code[-1]].append(code)
            prefix += ":"
    rows = sorted(table.items())
    # Hapl-o-Mat's placeholders for absent DRB3/4/5 genes.
    rows += [(f"DRB{i}*NNNN", [f"DRB{i}*NNNN"]) for i in (3, 4, 5)]
    return rows


def groups(hla_nom: Path, suffix: str) -> list[tuple[str, list[str]]]:
    """G groups (suffix 'G') or P groups (suffix 'P') with their member alleles."""
    rows = []
    with hla_nom.open(encoding="utf-8-sig") as f:
        lines = sorted(f)  # Hapl-o-Mat sorts whole raw lines
    for line in lines:
        if line.startswith("#"):
            continue
        locus, alleles, group = line.rstrip("\r\n").split(";")[:3]
        if group.endswith(suffix):
            rows.append((locus + group, [locus + a for a in alleles.split("/")]))
    return rows


def _smallest_allele(stripped: list[str], previous: str | None) -> str | None:
    """Pick the 'smallest' allele field by field, as BuildSmallg.py does.

    Returns `previous` when no single allele wins, matching the original,
    where the variable silently kept its value from the last group.
    """
    columns: list[list[str]] = [[], [], [], []]
    for a in stripped:
        parts = a.split(":")
        for i in range(4):
            columns[i].append(parts[i] if i < len(parts) else "")
    last_index: list[int] = []
    for column in columns:
        try:
            values: list = [int(v) for v in column]
        except ValueError:
            values = column
        lowest = min(values)
        index = [i for i, v in enumerate(values) if v == lowest]
        if not last_index:
            last_index = index
        else:
            if len(set(index) & set(last_index)) == 1:
                return stripped[index[0]]
            last_index = index
    return previous


def small_g(p_rows: list[tuple[str, list[str]]], g_rows: list[tuple[str, list[str]]]) -> dict[str, list[str]]:
    """g groups: P groups widened with the null/low alleles of matching G groups."""
    p_groups: list[tuple[str, set[str]]] = [(name, set(members)) for name, members in p_rows]
    min_code: str | None = None
    for _g_name, members in g_rows:
        with_letter = {c for c in members if c.endswith(EXPRESSION_LETTERS)}
        plain = {c for c in members if not c.endswith(EXPRESSION_LETTERS)}
        matched = False
        for _p_name, p_members in p_groups:
            if plain and plain <= p_members:
                p_members |= with_letter
                matched = True
        if not matched:
            # No P group covers this G group (e.g. a G group of one expressed
            # and one null allele): make a new group named after its smallest allele.
            stripped = [s.rstrip("NLSQ") for s in members]
            min_code = _smallest_allele(stripped, min_code)
            if min_code is None:
                raise ValueError(f"cannot name a g group for {members}")
            f1, f2 = min_code.split(":")[:2]
            p_groups.append((f"{f1}:{f2}P", set(members)))

    table: dict[str, list[str]] = {}
    for p_name, p_members in p_groups:
        name = p_name[:-1]
        if not all(":".join(c.split(":")[:2]) == name for c in p_members):
            name += "g"
        table[name] = sorted(p_members)
    return table


def add_alleles_missing_from_g(
    g_table: dict[str, list[str]], expanded_rows: list[tuple[str, list[str]]]
) -> list[tuple[str, list[str]]]:
    """Put alleles that no g group covers into the g group of their 2-field stem."""
    allele_to_g: dict[str, str] = {}
    stem_to_g: dict[str, str] = {}
    for g, alleles in g_table.items():
        for a in alleles:
            allele_to_g[a] = g
            stem_to_g[_stem(a)] = g

    missing: dict[str, set[str]] = defaultdict(set)
    for key, alleles in expanded_rows:
        if len(alleles) > 1 and key.count(":") == 1:
            for a in alleles:
                if a in allele_to_g:
                    continue
                stem = _stem(a)
                if stem in stem_to_g:
                    g = stem_to_g[stem]
                    g_table[g].append(a)
                    g_table[g].sort()
                    allele_to_g[a] = g
                else:
                    missing[stem].add(a)
    for stem, alleles in missing.items():
        if stem not in g_table:
            g_table[stem] = sorted(alleles)
    return sorted(g_table.items())


def _stem(allele: str) -> str:
    f = allele.split(":", 2)
    return f[0] + ":" + f[1]


def _write(path: Path, rows: list[tuple[str, list[str]]]) -> None:
    with path.open("w", newline="\n") as f:
        for key, values in rows:
            f.write(key + "\t" + "\t".join(values) + "\n")


def derive_all(release_dir: Path, out_dir: Path) -> dict[str, int]:
    """Write the four tables into out_dir; return the row count of each."""
    out_dir.mkdir(parents=True, exist_ok=True)
    nom_g, nom_p = release_dir / "hla_nom_g.txt", release_dir / "hla_nom_p.txt"

    expanded_rows = expanded(all_alleles(nom_g))
    p_rows = groups(nom_p, "P")
    g_rows = groups(nom_g, "G")
    g_small = add_alleles_missing_from_g(small_g(p_rows, g_rows), expanded_rows)

    tables = {
        "AllAllelesExpanded.txt": expanded_rows,
        "LargeG.txt": g_rows,
        "P.txt": p_rows,
        "Smallg.txt": g_small,
    }
    for name, rows in tables.items():
        _write(out_dir / name, rows)
    return {name: len(rows) for name, rows in tables.items()}


if __name__ == "__main__":
    import sys

    counts = derive_all(Path(sys.argv[1]), Path(sys.argv[2]))
    for name, n in counts.items():
        print(f"{name}: {n} rows")
