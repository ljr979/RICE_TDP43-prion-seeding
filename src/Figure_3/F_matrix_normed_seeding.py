"""Plot particle-dilution effects by seed batch for Figure 3D.

Reads the prepared dilution-response table and exports a patterned
bar-and-point summary across batches.
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.patches as mpatches
from statannotations.Annotator import Annotator
from src.paperfig_utils import save_figure_svg_png, set_params_for_plotting
from matplotlib.colors import LogNorm

output_folder = 'output/Figure_3/D_E_F_particles_titration/'

summ = pd.read_csv('input_data/Figure_3/Fig3F_matrix_seeding_normed.csv')
#only plotting below here.
batches = [1,2,3]
norm = LogNorm(vmin=0.5, vmax=2)

fig, axes = plt.subplots(1, 2, figsize=(7, 3.2), sharey=True, layout="constrained")

for ax, cond, title in zip(axes, ["original_difference", "normed_difference"],
                        ["Original", "Particle # equalised"]):
    s = summ[summ.type_difference == cond]
    mean = s.pivot(index="div_by_batch", columns="original_batch", values="mean").reindex(index=batches, columns=batches)
    sd = s.pivot(index="div_by_batch", columns="original_batch", values="std").reindex(index=batches, columns=batches)
    annot = (mean.round(2).astype(str) + "\n± " + sd.round(2).astype(str)).where(mean.notna(), "")

    sns.heatmap(mean, ax=ax, annot=False, cmap="PRGn", norm=norm, cbar=False,
            linewidths=1, linecolor="white", square=True)
    ax.set_facecolor("#e6e6e6")   # greys out the undefined diagonal
        
    cmap = ax.collections[0].cmap
    for i, ref in enumerate(batches):
        for j, seed in enumerate(batches):
            m, sdev = mean.loc[ref, seed], sd.loc[ref, seed]
            if pd.notna(m):
                r, g, b, _ = cmap(norm(m))
                lum = 0.299 * r + 0.587 * g + 0.114 * b
                ax.text(j + 0.5, i + 0.5, f"{m:.2f}\n± {sdev:.2f}",
                        ha="center", va="center", fontsize=9,
                        color="white" if lum < 0.45 else "black")

    ax.set_title(title)
    ax.set_xlabel("Seed batch")
    ax.set_ylabel("Batch normalised to" if ax is axes[0] else "")
    ax.set_xticklabels([f"{b}" for b in batches], rotation=0)
    ax.set_yticklabels([f"{b}" for b in batches], rotation=0)

cb = fig.colorbar(axes[0].collections[0], ax=axes, shrink=0.8,
                label="Degree of change due to normalisation")
cb.set_ticks([0.5, 0.75, 1, 1.5, 2])
cb.set_ticklabels(["0.5", "0.75", "1", "1.5", "2"])


save_figure_svg_png(fig, f'{output_folder}3F_particle-norm_seeding', width=3.5, height=None, dpi=600, tight=False)