"""Reproducible synthetic benchmarking and auditable plain-text outputs."""

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
import time

import numpy as np

from .model import Config, PreparedModel, simulate


def json_write(path, data):
    path.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def write_csv(path, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def wilson_interval(successes, total):
    z = 1.959963984540054
    p = successes / total
    centre = (p + z*z/(2*total)) / (1 + z*z/total)
    half = z * np.sqrt(p*(1-p)/total + z*z/(4*total*total)) / (1 + z*z/total)
    return [float(centre-half), float(centre+half)]


def summarize(rows):
    summaries, calibration = [], []
    for scenario in dict.fromkeys(row["scenario"] for row in rows):
        for mode in dict.fromkeys(row["mode"] for row in rows):
            selected = [r for r in rows if r["scenario"] == scenario and r["mode"] == mode]
            p = np.array([r["probability"] for r in selected])
            truth = np.array([r["truth"] for r in selected])
            clipped = np.clip(p, 1e-12, 1-1e-12)
            replicate_rows = [r for r in selected if r["gene"] == "SYN_CORE_01"]
            n = len(replicate_rows)
            covered = sum(r["joint_90_covered"] for r in replicate_rows)
            # Genes are dependent within a replicate. SE uses replicate mean scores.
            replicate_brier = ((p-truth)**2).reshape(n, 4).mean(axis=1)
            summaries.append({"scenario": scenario, "mode": mode, "replicates": n,
                "gene_decisions": len(selected), "brier": float(np.mean((p-truth)**2)),
                "brier_mcse_by_replicate": float(np.std(replicate_brier, ddof=1)/np.sqrt(n)),
                "log_loss": float(np.mean(-truth*np.log(clipped)-(1-truth)*np.log1p(-clipped))),
                "accuracy_at_0_5": float(np.mean((p >= .5) == truth)),
                "observed_core_prevalence": float(truth.mean()),
                "majority_noncore_accuracy": float(1-truth.mean()),
                "joint_90_coverage": covered/n,
                "joint_90_coverage_wilson_95": wilson_interval(covered, n),
                "mean_joint_credible_mass": float(np.mean([r["joint_credible_mass"] for r in replicate_rows])),
                "mean_joint_set_size": float(np.mean([r["joint_set_size"] for r in replicate_rows]))})
            for low in np.arange(0, 1, .1):
                high = min(1.0, round(float(low + .1), 10))
                mask = (p >= low) & ((p < high) if high < 1 else (p <= high))
                calibration.append({"scenario": scenario, "mode": mode,
                    "bin_low": round(float(low), 1), "bin_high": high, "count": int(mask.sum()),
                    "mean_probability": float(p[mask].mean()) if mask.any() else "",
                    "observed_fraction": float(truth[mask].mean()) if mask.any() else ""})
    return summaries, calibration


def run(args):
    if args.replicates < 2 or args.stress_replicates < 2:
        raise ValueError("At least two replicates per scenario are required")
    started = time.perf_counter()
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    config = Config()
    models = {mode: PreparedModel(config, mode, args.network_order, args.mapping_order)
              for mode in PreparedModel.MODES}
    rng = np.random.default_rng(args.seed)
    rows, datasets = [], []
    for scenario, n, multiplier in (("matched", args.replicates, 1.),
                                     ("noise_misspecification", args.stress_replicates, 1.8)):
        for replicate in range(n):
            sample = simulate(config, rng, multiplier, args.missing_rate)
            datasets.append({"scenario": scenario, "replicate": replicate,
                **{key: [None if not np.isfinite(x) else float(x) for x in value]
                   for key, value in sample.items()}})
            for mode, model in models.items():
                inference = model.infer(sample["variant"], sample["direct"])
                selected, mass = model.credible_set(inference["state_probability"])
                true_state = np.flatnonzero(np.all(model.states == sample["truth"], axis=1))[0]
                for gene, probability in enumerate(inference["core_probability"]):
                    rows.append({"scenario": scenario, "replicate": replicate, "mode": mode,
                        "gene": f"SYN_CORE_{gene+1:02}", "truth": int(sample["truth"][gene]),
                        "probability": float(probability), "joint_90_covered": int(true_state in selected),
                        "joint_credible_mass": mass, "joint_set_size": len(selected)})
            if (replicate+1) % 50 == 0:
                print(f"{scenario}: {replicate+1}/{n} completed", flush=True)
    summaries, calibration = summarize(rows)
    # Check convergence on the first ten matched datasets, without changing parameters.
    higher = PreparedModel(config, "full", args.network_order+4, args.mapping_order+8)
    convergence = []
    for sample in datasets[:min(10, args.replicates)]:
        variant = np.array([np.nan if x is None else x for x in sample["variant"]])
        direct = np.array([np.nan if x is None else x for x in sample["direct"]])
        a, b = models["full"].infer(variant, direct), higher.infer(variant, direct)
        convergence.append({"replicate": sample["replicate"],
            "max_absolute_marginal_difference": float(np.max(np.abs(a["core_probability"]-b["core_probability"]))),
            "joint_total_variation_distance": float(.5*np.abs(a["state_probability"]-b["state_probability"]).sum())})
    write_csv(output / "predictions.csv", rows)
    write_csv(output / "calibration.csv", calibration)
    json_write(output / "summary.json", summaries)
    json_write(output / "datasets.json", datasets)
    json_write(output / "convergence.json", convergence)
    source_hashes = {str(p.relative_to(Path(__file__).parent.parent)).replace("\\", "/"):
                     hashlib.sha256(p.read_bytes()).hexdigest()
                     for p in sorted(Path(__file__).parent.glob("*.py"))}
    metadata = {"created_utc": datetime.now(timezone.utc).isoformat(), "seed": args.seed,
        "status": "new synthetic educational prototype; not validated on genomic data",
        "python": platform.python_version(), "numpy": np.__version__, "config": config.to_dict(),
        "replicates_matched": args.replicates, "replicates_noise_stress": args.stress_replicates,
        "network_quadrature_order": args.network_order, "mapping_quadrature_order": args.mapping_order,
        "missing_rate": args.missing_rate, "duration_seconds": round(time.perf_counter()-started, 3),
        "source_sha256": source_hashes,
        "no_training": True, "source_data": "All observations and truth are synthetic; no biological identifiers"}
    json_write(output / "run_metadata.json", metadata)
    lines = ["# Executed synthetic benchmark", "", f"Seed: {args.seed}. Created: {metadata['created_utc']}.", "",
        "All numbers below concern generated toy data. They are not evidence of biological validity.", "",
        "| Scenario | Model | Brier (lower better) | Log loss | Joint 90% set coverage | Mean set mass |",
        "|---|---|---:|---:|---:|---:|"]
    for row in summaries:
        lines.append(f"| {row['scenario']} | {row['mode']} | {row['brier']:.4f} | {row['log_loss']:.4f} | {row['joint_90_coverage']:.3f} | {row['mean_joint_credible_mass']:.3f} |")
    lines += ["", "Coverage uses a joint credible set over the 16 possible core-membership states, not a confidence interval for a gene.",
        "Discrete sets often contain more than 90% posterior mass; compare coverage with mean set mass as well as 0.90.",
        "The stress scenario increases observation noise by 1.8 without informing inference.",
        "Calibration bins are descriptive: gene outcomes within each replicate are dependent.",
        "Monte Carlo standard errors and Wilson coverage intervals are recorded in summary.json.", "",
        f"Largest marginal change in the ten-dataset quadrature check: {max(x['max_absolute_marginal_difference'] for x in convergence):.3g}.",
        "This is a local convergence diagnostic, not a proof for every possible input.", "",
        "There is no train/test split: no model parameter was fitted or tuned. Independent datasets are prior-predictive draws."]
    (output / "RESULTS.md").write_text("\n".join(lines)+"\n", encoding="utf-8")
    print(json.dumps({"output": str(output), "prediction_rows": len(rows),
                      "duration_seconds": metadata["duration_seconds"]}, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default="results")
    parser.add_argument("--seed", type=int, default=20260925)
    parser.add_argument("--replicates", type=int, default=300)
    parser.add_argument("--stress-replicates", type=int, default=150)
    parser.add_argument("--network-order", type=int, default=8)
    parser.add_argument("--mapping-order", type=int, default=16)
    parser.add_argument("--missing-rate", type=float, default=.15)
    run(parser.parse_args())


if __name__ == "__main__":
    main()
