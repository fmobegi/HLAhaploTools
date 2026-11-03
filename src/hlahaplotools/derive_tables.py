"""Build the lookup tables used to resolve typings to g, G, P or n-field codes.

Input: one IMGT/HLA release directory (hla_nom_g.txt, hla_nom_p.txt, and
Allelelist.txt for the loci hla_nom_g.txt does not cover).
Output, in Hapl-o-Mat's own file formats so both tools can share them:

    AllAllelesExpanded.txt  any typing prefix -> full allele names
    LargeG.txt              G group -> member alleles
    P.txt                   P group -> member alleles
    Smallg.txt              g group -> member alleles

The rules, including their odd corners, follow Hapl-o-Mat's prepareData
scripts, so that results can be compared with the original tool on equal
terms. `hlahaplotools download-data` (reference_data.py) calls `derive_all`;
the comparison test is tests/test_derived_tables.py.

hla_nom_g.txt only covers loci with G groups. With `extra_loci=True` (the
default), alleles of every other locus in Allelelist.txt (MICA, MICB, H, J, K,
L, DQA2, TAP1, ...) are added too, so they can be resolved to 1-4 fields and
to g (their 2-field stem); IMGT/HLA defines no G or P groups for them. Rows
for the loci Hapl-o-Mat covers are unchanged.
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


def extra_alleles(allelelist: Path, covered: set[str]) -> list[str]:
    """Alleles in Allelelist.txt whose locus has no entry in hla_nom_g.txt."""
    out = []
    for line in _data_lines(allelelist):
        parts = line.split(",")
        if len(parts) < 2 or "*" not in parts[1]:
            continue
        name = parts[1].strip()
        if name.split("*", 1)[0] not in covered:
            out.append(name)
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
    with path.open("w", encoding="utf-8", newline="\n") as f:
        for key, values in rows:
            f.write(key + "\t" + "\t".join(values) + "\n")


# ---------------------------------------------------------------- serology and V2 names

# Serological locus (rel_ser_ser.txt) -> the allele loci whose alleles carry its antigens.
SEROLOGY_LOCI = {"A": ["A"], "B": ["B"], "C": ["C"], "Cw": ["C"], "DR": ["DRB1", "DRB3", "DRB4", "DRB5"],
                 "DQ": ["DQB1"], "DP": ["DPB1"], "DPw": ["DPB1"]}  # fmt: skip


def serology(rel_dna_ser: Path, rel_ser_ser: Path) -> list[tuple[str, list[str]]]:
    """Serological antigen -> full allele names, from the WMDA tables.

    rel_dna_ser.txt gives each allele its antigen (the unambiguous, assumed and expert
    assignments; 'possible' ones are left out, as they are not assignments). rel_ser_ser.txt
    gives each broad antigen its splits and each antigen its associated antigens, so a
    broad antigen such as B15 covers the alleles of B62, B63, ... Keys are written like
    typings, `B*62`, one per allele locus (DR52 gives `DRB3*52`).
    """
    by_antigen: dict[tuple[str, str], list[str]] = defaultdict(list)
    for line in _data_lines(rel_dna_ser):
        c = line.split(";")
        if len(c) < 6:
            continue
        locus, allele = c[0].rstrip("*"), c[1]
        for antigen in {c[2], c[4], c[5]}:
            if antigen and antigen not in ("0", "?"):
                by_antigen[(locus, antigen)].append(f"{locus}*{allele}")
    below: dict[tuple[str, str], set[str]] = defaultdict(set)  # (serological locus, antigen) -> splits/associated
    for line in _data_lines(rel_ser_ser):
        c = line.split(";")
        if len(c) < 4:
            continue
        for other in (c[2] + "/" + c[3]).split("/"):
            if other:
                below[(c[0], c[1])].add(other)
    rows: dict[str, list[str]] = {}
    sero_loci = {s for s, _ in below} | set(SEROLOGY_LOCI)
    antigens = {(s, a) for s, a in below} | {(s, a) for s in sero_loci for (l, a) in by_antigen
                                                if l in SEROLOGY_LOCI.get(s, [])}  # fmt: skip
    for s, antigen in sorted(antigens):
        seen, todo = set(), [antigen]
        while todo:
            x = todo.pop()
            if x not in seen:
                seen.add(x)
                todo.extend(below.get((s, x), ()))
        for locus in SEROLOGY_LOCI.get(s, []):
            alleles = sorted({a for x in seen for a in by_antigen.get((locus, x), [])})
            if alleles:
                key = f"{locus}*{antigen.zfill(2)}"
                rows.setdefault(key, alleles)
    return sorted(rows.items())


def v2_names(nomenclature_2009: Path) -> list[tuple[str, list[str]]]:
    """Pre-2010 allele names -> current names, for every field level.

    Nomenclature_2009.txt lists full names (A*01010101 -> A*01:01:01:01); the shorter
    typings are derived from where the colons fall (A*0101 -> A*01:01, A*010101 ->
    A*01:01:01). A short V2 name that would map to two current names is left out.
    """
    out: dict[str, set[str]] = defaultdict(set)
    with nomenclature_2009.open(encoding="utf-8-sig") as f:
        for line in f:
            parts = line.split()
            if len(parts) != 2 or "*" not in parts[0] or ":" not in parts[1]:
                continue
            old, new = parts
            locus, fields = new.split("*", 1)
            suffix = fields[-1] if fields[-1].isalpha() else ""
            fields = fields.rstrip("NLSQCAE")
            parts_new = fields.split(":")
            for n in range(2, len(parts_new) + 1):
                digits = "".join(parts_new[:n])
                end = suffix if n == len(parts_new) else ""
                out[f"{locus}*{digits}{end}"].add(f"{locus}*{':'.join(parts_new[:n])}{end}")
            out[old].add(new)
    return sorted((k, sorted(v)) for k, v in out.items() if len(v) == 1)


def derive_all(release_dir: Path, out_dir: Path, extra_loci: bool = True) -> dict[str, int]:
    """Write the tables into out_dir; return the row count of each.

    ``extra_loci=False`` reproduces Hapl-o-Mat's prepareData output exactly.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    nom_g, nom_p = release_dir / "hla_nom_g.txt", release_dir / "hla_nom_p.txt"

    alleles = all_alleles(nom_g)
    if extra_loci and (release_dir / "Allelelist.txt").exists():
        covered = {a.split("*", 1)[0] for a in alleles}
        alleles += extra_alleles(release_dir / "Allelelist.txt", covered)
    expanded_rows = expanded(alleles)
    p_rows = groups(nom_p, "P")
    g_rows = groups(nom_g, "G")
    g_small = add_alleles_missing_from_g(small_g(p_rows, g_rows), expanded_rows)

    tables = {
        "AllAllelesExpanded.txt": expanded_rows,
        "LargeG.txt": g_rows,
        "P.txt": p_rows,
        "Smallg.txt": g_small,
    }
    # Optional tables (serological antigens, pre-2010 names), when their sources are present.
    if (release_dir / "rel_dna_ser.txt").exists() and (release_dir / "rel_ser_ser.txt").exists():
        tables["Serology.txt"] = serology(release_dir / "rel_dna_ser.txt", release_dir / "rel_ser_ser.txt")
    if (release_dir / "Nomenclature_2009.txt").exists():
        tables["V2.txt"] = v2_names(release_dir / "Nomenclature_2009.txt")
    for name, rows in tables.items():
        _write(out_dir / name, rows)
    return {name: len(rows) for name, rows in tables.items()}


if __name__ == "__main__":
    import sys

    counts = derive_all(Path(sys.argv[1]), Path(sys.argv[2]))
    for name, n in counts.items():
        print(f"{name}: {n} rows")
