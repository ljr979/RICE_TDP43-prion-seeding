"""Shared plotting helpers used by the manuscript figure scripts."""

import os
import matplotlib.pyplot as plt

def save_figure_svg_png(fig, filename, width=3.5, height=None, dpi=300, tight=True):
    """
    Save a matplotlib figure as both SVG and PNG with standardized formatting.
    
    Parameters:
    - fig: matplotlib.figure.Figure
        The figure object to save.
    - filename: str
        Path (with or without extension). Saves as .svg and .png.
    - width: float
        Width in inches. Default is 3.5" (single column).
    - height: float or None
        Height in inches. If None, uses golden ratio.
    - dpi: int
        Dots per inch for PNG (raster export).
    - tight: bool
        Whether to use tight_layout before saving.
    """
    # Set height using golden ratio if not provided
    if height is None:
        height = width * 0.618  # Golden ratio

    fig.set_size_inches(width, height)

    if tight:
        plt.tight_layout()

    # Ensure filename has no extension so we can add both .svg and .png
    base, _ = os.path.splitext(filename)

    svg_path = f"{base}.svg"
    png_path = f"{base}.png"

    # Save SVG
    fig.savefig(svg_path, format='svg', dpi=dpi)
    # Save PNG
    fig.savefig(png_path, format='png', dpi=dpi)

    print(f"Saved SVG to: {svg_path}")
    print(f"Saved PNG to: {png_path}")

def set_params_for_plotting(fontsize=10):
    """Set the shared Matplotlib typography and SVG text settings."""
    new_rc_params = {'text.usetex': False,
    'font.family': 'Arial',
    'font.size': fontsize,            # Base font size
    'axes.titlesize': fontsize,
    'axes.labelsize': fontsize,
    'xtick.labelsize': fontsize,
    'ytick.labelsize': fontsize,      
    'legend.fontsize': fontsize,
    'figure.titlesize': fontsize,
    'svg.fonttype': 'none',     # Optional: keeps text as text in SVG
    }
    
    plt.rcParams.update(new_rc_params)

def make_palette():
    """Return the manuscript color mapping for the listed protein labels."""
    my_palette = {'Atto647 alone': (0.25, 0.25, 0.25),
    'WT': (0.617762399077278, 0.6021376393694733, 0.7834525182622069),
 'G335D': (0.9835755478662053, 0.4127950788158401, 0.28835063437139563),
 'A382P': (0.45176470588235296, 0.7670895809304115, 0.4612072279892349),
 'W385G': (0.41708573625528644, 0.6806305267204922, 0.8382314494425221)}
    return my_palette
