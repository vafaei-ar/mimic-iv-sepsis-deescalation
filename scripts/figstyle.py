#!/usr/bin/env python3
"""Shared Nature-specification figure style for the publication figures.

Nature's author guidelines fix the physical geometry of a figure: it is placed at
89 mm (single column), 120 mm (1.5 column) or 183 mm (double column), it may not
exceed 247 mm of depth, and type must stay legible at 5-7 pt after the figure is
reduced to that width. Sizing a figure in arbitrary inches and letting the
publisher scale it is what produces the mismatched type sizes seen across a
typical submission. Every builder here therefore starts from a column width.

`pdf.fonttype = 42` embeds TrueType outlines rather than Type 3, which is what
lets a production editor open and edit the vector text; Type 3 is rejected by
several journals outright.
"""
from __future__ import annotations

import matplotlib as mpl
import matplotlib.pyplot as plt

# Nature column widths, millimetres -> inches.
MM = 1.0 / 25.4
SINGLE = 89 * MM
ONE_HALF = 120 * MM
DOUBLE = 183 * MM
MAX_DEPTH = 247 * MM

# Categorical pair, validated for colour-vision deficiency: worst adjacent
# separation dE 21.9 (protan), 31.2 (normal vision), both well above the
# dE >= 8 target. Derived from the Okabe-Ito safe qualitative set.
BLUE = "#0072B2"      # MIMIC-IV
VERMILLION = "#D55E00"  # Penn State
GREEN = "#009E73"     # third slot, only if a third series is unavoidable

# Non-data ink. Reference lines and context series must never wear a data hue -
# the current figures draw the null line in the same blue as the estimates,
# which reads as a fifth data series.
INK = "#1A1A1A"
MUTED = "#6E6E6E"
FAINT = "#B8B8B8"
RULE = "#9A9A9A"


def use_nature_style() -> None:
    """Install the rcParams every figure in the submission shares."""
    mpl.rcParams.update({
        # Helvetica is the Nature house face; Nimbus Sans and Liberation Sans are
        # metrically compatible clones and are what is actually installed here.
        "font.family": "sans-serif",
        "font.sans-serif": ["Nimbus Sans", "Liberation Sans", "Helvetica", "Arial", "DejaVu Sans"],
        "font.size": 7,
        "axes.labelsize": 7,
        "axes.titlesize": 7,
        "xtick.labelsize": 6,
        "ytick.labelsize": 6,
        "legend.fontsize": 6,
        "figure.titlesize": 8,

        # Hairlines. Nature specifies 0.25-1 pt; matplotlib's 0.8 default axis
        # and 1.5 default line are both too heavy at column width.
        "axes.linewidth": 0.5,
        "lines.linewidth": 1.0,
        "xtick.major.width": 0.5,
        "ytick.major.width": 0.5,
        "xtick.major.size": 2.5,
        "ytick.major.size": 2.5,
        "xtick.direction": "out",
        "ytick.direction": "out",

        # Two spines, not a box. The top and right rules carry no information.
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.edgecolor": INK,
        "axes.labelcolor": INK,
        "text.color": INK,
        "xtick.color": INK,
        "ytick.color": INK,

        "legend.frameon": False,
        "axes.grid": False,

        # Vector output a production editor can edit.
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "savefig.dpi": 600,
        "figure.dpi": 200,
        # Deliberately NOT bbox="tight". Tight-cropping recomputes the canvas from
        # the rendered artists, so the PDF comes out whatever width the labels
        # happen to need - which defeats sizing to a column in the first place.
        # Margins are reserved explicitly instead, and the saved page is exactly
        # the figsize declared.
        "savefig.bbox": None,
        "savefig.transparent": False,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
    })


def panel_label(ax, letter: str, dx: float = -0.08, dy: float = 1.04) -> None:
    """Bold lowercase panel letter in the top-left, per Nature convention."""
    ax.text(
        dx, dy, letter,
        transform=ax.transAxes,
        fontsize=8, fontweight="bold", va="bottom", ha="left", color=INK,
    )


def null_line(ax, x: float = 0.0) -> None:
    """Neutral reference rule. Recessive, behind the data, never a data hue."""
    ax.axvline(x, color=RULE, linestyle=(0, (3, 2)), linewidth=0.5, zorder=0)


def strip_y_axis(ax) -> None:
    """Categorical y axes carry their meaning in the tick text, not a spine."""
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)


def savefig(fig, out_dir, stem: str) -> None:
    """Write the PDF a journal wants alongside a PNG for quick review."""
    out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir / f"{stem}.pdf")
    fig.savefig(out_dir / f"{stem}.png", dpi=600)
    plt.close(fig)
