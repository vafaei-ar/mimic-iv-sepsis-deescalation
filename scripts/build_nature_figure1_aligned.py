#!/usr/bin/env python3
"""Rebuild Figure 1 with aligned panel letters and tighter lower-panel spacing.

This is a presentation-only wrapper around the frozen Nature-style Figure 1 drawing
functions. It does not alter cohort definitions, counts, estimands, or source data.
It writes the same Figure 1 output stem used by ``build_nature_figures.py`` so the
corrected layout can replace the prior figure directly.
"""
from __future__ import annotations

import matplotlib.pyplot as plt
import pandas as pd

import figstyle as fs
from build_nature_figures import (
    FLOW,
    OUT,
    _draw_attrition,
    _draw_timeline,
    _draw_treatment_split,
)

PANEL_X = -0.08
PANEL_Y = 1.01


def _normalize_panel_letter(ax, letter: str) -> None:
    """Replace the helper-drawn panel letter with one shared anchor."""
    for text in list(ax.texts):
        if text.get_text() == letter and text.get_fontweight() == "bold":
            text.remove()
    fs.panel_label(ax, letter, dx=PANEL_X, dy=PANEL_Y)


def build_fig1() -> None:
    """Build the complete Figure 1 with a common a/b/c label position."""
    f = pd.read_csv(FLOW)
    fig = plt.figure(figsize=(fs.DOUBLE, 4.95))

    # Keep panel a visually separate while pulling b and c into a single lower
    # block. The three axes share the same left/right figure margins.
    gs = fig.add_gridspec(
        3,
        1,
        height_ratios=[1.15, 2.45, 0.62],
        hspace=0.18,
    )
    axes = [fig.add_subplot(gs[i]) for i in range(3)]

    _draw_timeline(axes[0])
    _draw_attrition(axes[1], f)
    _draw_treatment_split(axes[2], f)

    for ax, letter in zip(axes, "abc"):
        _normalize_panel_letter(ax, letter)

    fig.subplots_adjust(left=0.30, right=0.97, top=0.965, bottom=0.055)
    fs.savefig(fig, OUT, "Fig1_target_trial_and_cohort")


def main() -> None:
    fs.use_nature_style()
    OUT.mkdir(parents=True, exist_ok=True)
    if not FLOW.exists():
        raise FileNotFoundError(f"Missing required cohort-flow input: {FLOW}")
    build_fig1()
    print("Built aligned Figure 1")


if __name__ == "__main__":
    main()
