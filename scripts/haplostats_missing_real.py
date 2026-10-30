"""Real cohort with its own untyped DQA1/DPA1, for the haplo.stats comparison (hlahaplotools, 10 starts, 8 threads).

    uv run python scripts/haplostats_missing_real.py results/10July2025/haplomap_input_current_names.dat bench/real_missing
"""

import json
import sys
from pathlib import Path

import hlahaplotools as h

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
from benchmark import prepare_real, run_ours, save_freqs, write_mac  # noqa: E402

SETS = {  # name: (loci, loci that may be untyped)
    "DQA1_DQB1": (["DQA1", "DQB1"], {"DQA1"}),
    "DRB1_DQA1_DQB1": (["DRB1", "DQA1", "DQB1"], {"DQA1"}),
    "class2": (["DRB1", "DQA1", "DQB1", "DPA1", "DPB1"], {"DQA1", "DPA1"}),
    "classical8": (["A", "B", "C", "DRB1", "DQB1", "DPB1", "DQA1", "DPA1"], {"DQA1", "DPA1"}),
}

src, out = Path(sys.argv[1]).resolve(), Path(sys.argv[2]).resolve()
manifest = json.loads((ROOT / "data" / "manifest.json").read_text())
nom = h.Nomenclature(str(ROOT / "data" / "derived" / manifest["imgthla"]["release"]))
for name, (loci, untyped) in SETS.items():
    rows, info = prepare_real(src, loci, "2f", nom, 5000, 1, may_be_untyped=untyped)
    sim = out / name
    sim.mkdir(parents=True, exist_ok=True)
    data = sim / "population.dat"
    write_mac(data, loci, rows)
    info["incompletely_typed"] = sum(1 for _, per in rows if not all(a for a, _ in per))
    r = run_ours(sim, data, loci, 1e-6, 8, 3600, starts=10)
    freqs = r.pop("freqs", None)
    r.pop("stdout", None)
    if freqs is not None:
        save_freqs(sim / "freqs_hlahaplotools.tsv", freqs)
        r.pop("stderr", None)
    r["info"] = info
    (sim / "ours.json").write_text(json.dumps(r, indent=2) + "\n")
    print(name, "~".join(loci), info["used"], f"({info['incompletely_typed']} incomplete)", r.get("status"),
          f"{r.get('seconds', 0):.1f}s {r.get('peak_mb', 0):.0f} MB", flush=True)  # fmt: skip
