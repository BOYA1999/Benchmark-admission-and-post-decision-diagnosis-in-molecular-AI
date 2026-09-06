from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "additional_analyses" / "cross_case" / "ace_threshold_decision_map.csv"
OUT = ROOT / "supplementary_figures"
STYLE = OUT / "publication.mplstyle"

plt.style.use(STYLE)
plt.rcParams.update({
    "text.usetex": False,
})

frame = pd.read_csv(DATA)
curves = frame.groupby(["estimand", "threshold"], as_index=False)["positive_datasets_out_of_3"].first()
colors = {"endpoint_deduplicated": "#6B7280", "pair_equal": "#B55245"}
labels = {"endpoint_deduplicated": "Unique = inverse-multiplicity", "pair_equal": "Pair-equal"}

fig, ax = plt.subplots(figsize=(7.2, 4.4))
for estimand in ["endpoint_deduplicated", "pair_equal"]:
    rows = curves[curves["estimand"] == estimand]
    ax.step(rows["threshold"], rows["positive_datasets_out_of_3"], where="post", color=colors[estimand], linewidth=2.2, label=labels[estimand])

ax.axvline(1.25, color="#2A9D8F", linewidth=1.4, linestyle="--")
ax.text(1.2525, 0.28, "Recorded threshold 1.25", color="#237E72", fontsize=8.4, rotation=90, va="bottom")
ax.set_xlim(1.15, 1.35)
ax.set_ylim(-0.08, 3.12)
ax.set_xticks([1.15, 1.20, 1.25, 1.30, 1.35])
ax.set_yticks([0, 1, 2, 3])
ax.set_xlabel("RMSE-ratio threshold")
ax.set_ylabel("ECFP-positive data sets (of 3)")
ax.tick_params(direction="out")
ax.grid(False)
ax.legend(loc="upper right")
fig.tight_layout(pad=0.8)

stem = OUT / "supplementary_figure_s2_ace_threshold_map"
fig.savefig(stem.with_suffix(".png"), dpi=300, facecolor="white")
fig.savefig(stem.with_suffix(".svg"), facecolor="white")
fig.savefig(stem.with_suffix(".pdf"), facecolor="white")
fig.savefig(stem.with_suffix(".tiff"), dpi=600, facecolor="white", pil_kwargs={"compression": "tiff_lzw"})
plt.close(fig)
print(f"saved: {stem.with_suffix('.png')}")
