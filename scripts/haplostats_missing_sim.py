"""Missing-typing sets for the haplo.stats comparison: simulate, then hlahaplotools (10 starts, 8 threads).

    uv run python scripts/haplostats_missing_sim.py bench/missing_typing
"""

import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
from benchmark import compare, read_hfs, run_ours, save_freqs  # noqa: E402

C8 = ["A", "B", "C", "DRB1", "DQB1", "DPB1", "DQA1", "DPA1"]
SETS = {  # name: (loci, the two loci left untyped, rate)
    "l4_r30": (C8[:4], C8[2:4], 0.3),
    "l6_r30": (C8[:6], C8[4:6], 0.3),
    "l8_r30": (C8, C8[6:8], 0.3),
    "l10_r30": (C8 + ["MICA", "MICB"], ["MICA", "MICB"], 0.3),
    "l6_r10": (C8[:6], C8[4:6], 0.1),
    "l8_r10": (C8, C8[6:8], 0.1),
    "l10_r10": (C8 + ["MICA", "MICB"], ["MICA", "MICB"], 0.1),
}

out = Path(sys.argv[1]).resolve()
for name, (loci, miss, rate) in SETS.items():
    sim = out / name
    sim.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [sys.executable, str(ROOT / "scripts" / "simulate.py"), str(sim), "--samples", "2000",
         "--loci", ",".join(loci), "--haplotypes", "300", "--alleles", "25", "--seed", "1",
         "--missing", ",".join(f"{l}:{rate}" for l in miss)],
        check=True, capture_output=True,
    )  # fmt: skip
    r = run_ours(sim, sim / "population.dat", loci, 1e-6, 8, 3600, starts=10)
    freqs = r.pop("freqs", None)
    r.pop("stdout", None)
    if freqs is not None:
        save_freqs(sim / "freqs_hlahaplotools.tsv", freqs)
        r["tvd_truth"] = compare(freqs, read_hfs(sim / "truth.tsv"))["tvd"]
        r.pop("stderr", None)
    (sim / "ours.json").write_text(json.dumps(r, indent=2) + "\n")
    print(name, "~".join(loci), r.get("status"), f"{r.get('seconds', 0):.1f}s", f"TVD {r.get('tvd_truth')}", flush=True)
