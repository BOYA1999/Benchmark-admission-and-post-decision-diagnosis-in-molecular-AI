from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "additional_analyses" / "cross_case" / "rso_identity_sensitivity.csv"
OUT = ROOT / "supplementary_figures"
STYLE = OUT / "publication.mplstyle"

plt.style.use(STYLE)
plt.rcParams.update({
    "text.usetex": False,
})

frame = pd.read_csv(DATA)
labels = ["Exact tuple", "Ignore reagent", "Product only"]
all_rows = 100 * frame["overlap_total_rate"].to_numpy()
valid_rows = 100 * frame["overlap_valid_rate"].to_numpy()
all_counts = frame["overlap_total_rows"].astype(int).to_numpy()
valid_counts = frame["overlap_valid_rows"].astype(int).to_numpy()
all_denominators = frame["test_total_rows"].astype(int).to_numpy()
valid_denominators = frame["test_valid_reactant_product_rows"].astype(int).to_numpy()
x = np.arange(len(labels))
width = 0.34

fig, ax = plt.subplots(figsize=(7.2, 4.4))
bars_all = ax.bar(x - width / 2, all_rows, width, color="#A9A9A9", edgecolor="white", linewidth=0.8, label="All test rows")
bars_valid = ax.bar(x + width / 2, valid_rows, width, color="#2A9D8F", hatch="//", edgecolor="white", linewidth=0.8, label="Parser-valid rows")
for bars, values, counts, denominators in [(bars_all, all_rows, all_counts, all_denominators), (bars_valid, valid_rows, valid_counts, valid_denominators)]:
    for bar, value, count, denominator in zip(bars, values, counts, denominators):
        ax.text(bar.get_x() + bar.get_width() / 2, value + 0.22, f"{count:,}/{denominator:,}\n{value:.3f}%", ha="center", va="bottom", fontsize=7.8, linespacing=1.05)

ax.set_xticks(x, labels)
ax.set_ylabel("Test rows overlapping train or validation (%)")
ax.set_xlabel("Reaction identity definition")
ax.set_ylim(0, 15.4)
ax.yaxis.grid(True, color="#E7E7E7", linewidth=0.7, linestyle="--")
ax.set_axisbelow(True)
ax.tick_params(length=0)
ax.legend(loc="upper left")
fig.tight_layout(pad=0.8)

stem = OUT / "supplementary_figure_s1_rso_identity_sensitivity"
fig.savefig(stem.with_suffix(".png"), dpi=300, facecolor="white")
fig.savefig(stem.with_suffix(".svg"), facecolor="white")
fig.savefig(stem.with_suffix(".pdf"), facecolor="white")
fig.savefig(stem.with_suffix(".tiff"), dpi=600, facecolor="white", pil_kwargs={"compression": "tiff_lzw"})
plt.close(fig)
print(f"saved: {stem.with_suffix('.png')}")
