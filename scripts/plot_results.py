"""Plot recorded results; does not rerun, fit or modify the statistical model."""

import argparse
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", default="results")
    args = parser.parse_args()
    folder = Path(args.results)
    with (folder / "calibration.csv").open(encoding="utf-8", newline="") as handle:
        bins = list(csv.DictReader(handle))
    summary = json.loads((folder / "summary.json").read_text(encoding="utf-8"))
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                        "axes.spines.top": False, "axes.spines.right": False,
                        "svg.fonttype": "none", "svg.hashsalt": "coregene-synthetic-20260925"})
    colours = {"full": "#086e82", "fixed_both": "#c66a21"}
    labels = {"full": "Full uncertainty", "fixed_both": "Fixed mapping + network"}
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5), gridspec_kw={"width_ratios": [1, 1, 1.15]})
    fig.patch.set_facecolor("#fcfcfa")
    for axis, scenario, title in zip(axes[:2], ("matched", "noise_misspecification"),
                                     ("Matched model · 300 datasets", "Noise stress · 150 datasets")):
        axis.plot([0, 1], [0, 1], linestyle="--", color="#94a3a8", lw=1)
        for mode in ("full", "fixed_both"):
            rows = [r for r in bins if r["scenario"] == scenario and r["mode"] == mode and int(r["count"]) > 0]
            x, y = [float(r["mean_probability"]) for r in rows], [float(r["observed_fraction"]) for r in rows]
            sizes = [16 + 2*np.sqrt(int(r["count"])) for r in rows]
            axis.plot(x, y, color=colours[mode], alpha=.7, lw=1.1)
            axis.scatter(x, y, s=sizes, color=colours[mode], alpha=.85, label=labels[mode], zorder=3)
        axis.set(xlim=(-.02, 1.02), ylim=(-.02, 1.02), xlabel="Mean predicted core probability",
                 ylabel="Observed fraction of synthetic cores", title=title)
        axis.grid(alpha=.15)
        axis.set_aspect("equal")
    axes[0].legend(frameon=False, fontsize=8, loc="upper left")
    axis = axes[2]
    y_positions = [3, 2, 1, 0]
    ticklabels = []
    for y_pos, (scenario, mode) in zip(y_positions, (("matched", "full"), ("matched", "fixed_both"),
                                                    ("noise_misspecification", "full"), ("noise_misspecification", "fixed_both"))):
        row = next(r for r in summary if r["scenario"] == scenario and r["mode"] == mode)
        coverage = row["joint_90_coverage"]
        lo, hi = row["joint_90_coverage_wilson_95"]
        axis.errorbar(coverage, y_pos, xerr=np.array([[coverage-lo], [hi-coverage]]), fmt="o",
                      color=colours[mode], capsize=3, markersize=6)
        axis.scatter(row["mean_joint_credible_mass"], y_pos, marker="|", s=150, color=colours[mode], zorder=4)
        ticklabels.append(("Matched" if scenario == "matched" else "Noise stress") + " / " + ("full" if mode == "full" else "fixed"))
    axis.axvline(.9, ls="--", color="#94a3a8", lw=1)
    axis.set(yticks=y_positions, yticklabels=ticklabels, xlim=(.45, 1.02), ylim=(-.5, 3.5),
             xlabel="Empirical coverage", title="Joint 90% credible-set coverage")
    axis.grid(axis="x", alpha=.15)
    fig.suptitle("CoreGene Bayes Benchmark | Synthetic calibration audit", x=.05, ha="left", fontsize=15, fontweight="bold")
    fig.text(.05, .065, "Reliability points: ten probability bins; marker size reflects count. Bins are descriptive; within-dataset genes are dependent.", fontsize=9, color="#444444")
    fig.text(.05, .025, "Coverage: dots and bars = observed coverage and Wilson 95% interval; vertical ticks = actual mean posterior set mass. No biological validation.", fontsize=9, color="#444444")
    fig.tight_layout(rect=(0, .11, 1, .92), w_pad=2.5)
    fig.savefig(folder / "calibration.svg", metadata={"Date": None}, facecolor=fig.get_facecolor())
    fig.savefig(folder / "calibration.png", dpi=180, facecolor=fig.get_facecolor())
    plt.close(fig)
    (folder / "plot_metadata.json").write_text(json.dumps({"matplotlib": matplotlib.__version__, "numpy": np.__version__,
        "inputs": ["calibration.csv", "summary.json"], "outputs": ["calibration.svg", "calibration.png"],
        "note": "Plotting-only runtime; statistical execution versions are recorded in run_metadata.json"}, indent=2)+"\n", encoding="utf-8")
    print("Created calibration.svg and calibration.png from recorded results.")


if __name__ == "__main__":
    main()
