#!/usr/bin/env python3
"""
Generate a publication-quality figure of the 50 points on S^2 for counterexample 1.

Shows:
  - The unit sphere (transparent surface + faint grid)
  - 30 original icosidodecahedron points (filled circles, dark blue)
  - 20 new expansion points (diamonds, red)
  - The 40 great-circle triples drawn as arcs, faded for back-facing segments
  - A clean, TeX-friendly style (no title, tight layout, PDF output)
"""

from __future__ import annotations

import math
import sys
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from mpl_toolkits.mplot3d import Axes3D          # noqa: F401

# ── make the sibling module importable ──────────────────────────────────────
sys.path.insert(0, os.path.dirname(__file__))
from counterexample1_50pts import get_numeric_coords, TRIPLES, COORDS

# ── viewing direction (used for depth-based fading) ──────────────────────────
ELEV = 22
AZIM = 35

def view_dir(elev_deg: float, azim_deg: float) -> np.ndarray:
    """Unit vector pointing toward the viewer."""
    e = math.radians(elev_deg)
    a = math.radians(azim_deg)
    return np.array([math.cos(e) * math.cos(a),
                     math.cos(e) * math.sin(a),
                     math.sin(e)])

VIEW = view_dir(ELEV, AZIM)

# ── helpers ─────────────────────────────────────────────────────────────────

def is_new_point(idx: int) -> bool:
    return any(tok in "".join(COORDS[idx]) for tok in ("x", "y"))


def slerp(p: np.ndarray, q: np.ndarray, t: float) -> np.ndarray:
    dot = float(np.clip(np.dot(p, q), -1.0, 1.0))
    omega = math.acos(dot)
    if abs(omega) < 1e-10:
        return p.copy()
    return (math.sin((1 - t) * omega) * p + math.sin(t * omega) * q) / math.sin(omega)


def great_arc(p: np.ndarray, q: np.ndarray, n: int = 80) -> np.ndarray:
    return np.array([slerp(p, q, t / n) for t in range(n + 1)])


def depth(pt: np.ndarray) -> float:
    """Signed depth: positive = facing viewer, negative = behind."""
    return float(np.dot(pt, VIEW))


def sphere_surface(ax) -> None:
    """Draw a very faint sphere surface + light grid lines."""
    u = np.linspace(0, 2 * np.pi, 60)
    v = np.linspace(0, np.pi, 30)
    xs = np.outer(np.cos(u), np.sin(v))
    ys = np.outer(np.sin(u), np.sin(v))
    zs = np.outer(np.ones_like(u), np.cos(v))
    ax.plot_surface(xs, ys, zs, color="#d0dce8", alpha=0.18,
                    linewidth=0, antialiased=True, zorder=0)

    gc = "#9aafc0"
    lw_grid = 0.3
    # latitudes
    for lat in np.linspace(-60, 60, 5):
        th = np.linspace(0, 2 * np.pi, 200)
        r = math.cos(math.radians(lat))
        z = math.sin(math.radians(lat))
        ax.plot(r * np.cos(th), r * np.sin(th), z * np.ones_like(th),
                color=gc, lw=lw_grid, alpha=0.5, zorder=1)
    # longitudes
    for lon in np.linspace(0, 150, 6):
        th = np.linspace(0, 2 * np.pi, 200)
        lo = math.radians(lon)
        ax.plot(math.cos(lo) * np.cos(th), math.sin(lo) * np.cos(th), np.sin(th),
                color=gc, lw=lw_grid, alpha=0.5, zorder=1)


# ── main ────────────────────────────────────────────────────────────────────

def main() -> None:
    pts_list = get_numeric_coords()
    pts = np.array(pts_list)          # shape (50, 3)

    old_mask = np.array([not is_new_point(i) for i in range(50)])
    new_mask = ~old_mask

    # ── figure setup ────────────────────────────────────────────────────────
    fig = plt.figure(figsize=(5.0, 5.0))
    ax = fig.add_subplot(111, projection="3d")
    ax.set_box_aspect([1, 1, 1])
    ax.view_init(elev=ELEV, azim=AZIM)
    # zoom in to reduce white space around the sphere
    ax.set_xlim(-0.92, 0.92)
    ax.set_ylim(-0.92, 0.92)
    ax.set_zlim(-0.92, 0.92)
    # remove the extra padding matplotlib adds around 3-D axes
    ax.margins(0)

    # ── sphere ──────────────────────────────────────────────────────────────
    sphere_surface(ax)

    # ── great-circle arcs with depth-based fading ───────────────────────────
    ARC_COLOR_FRONT = np.array([0.18, 0.47, 0.71])   # steel blue
    ARC_COLOR_BACK  = np.array([0.65, 0.78, 0.90])   # pale blue

    N_SEG = 80
    for a, b, c in TRIPLES:
        for i, j in [(a, b), (b, c), (a, c)]:
            arc = great_arc(pts[i], pts[j], n=N_SEG)
            # draw segment by segment with depth-blended colour
            for k in range(len(arc) - 1):
                mid = (arc[k] + arc[k + 1]) / 2
                d = depth(mid)
                # map depth [-1,1] -> [0,1]
                t = (d + 1) / 2
                col = (1 - t) * ARC_COLOR_BACK + t * ARC_COLOR_FRONT
                alpha = 0.25 + 0.55 * t
                lw = 0.6 + 0.7 * t
                ax.plot([arc[k, 0], arc[k+1, 0]],
                        [arc[k, 1], arc[k+1, 1]],
                        [arc[k, 2], arc[k+1, 2]],
                        color=col, lw=lw, alpha=alpha,
                        solid_capstyle="round", zorder=2)

    # ── points ──────────────────────────────────────────────────────────────
    OLD_COLOR = "#154360"   # deep navy
    NEW_COLOR = "#922b21"   # deep red

    # depth-shade manually: front points larger & brighter
    def scatter_with_depth(mask, color, marker, base_s):
        idx = np.where(mask)[0]
        depths = np.array([depth(pts[k]) for k in idx])
        d_norm = (depths + 1) / 2
        sizes = base_s * (0.50 + 0.70 * d_norm)
        base_rgb = np.array(matplotlib.colors.to_rgb(color))
        for k, dn, s in zip(idx, d_norm, sizes):
            col = base_rgb + (1 - base_rgb) * (1 - dn) * 0.52
            ax.scatter([pts[k, 0]], [pts[k, 1]], [pts[k, 2]],
                       s=s, c=[col], marker=marker,
                       edgecolors="white", linewidths=0.5, zorder=6)

    scatter_with_depth(old_mask, OLD_COLOR, "o", 70)
    scatter_with_depth(new_mask, NEW_COLOR, "D", 48)

    # ── manual legend (proxy artists) ───────────────────────────────────────
    legend_elements = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor=OLD_COLOR,
               markeredgecolor="white", markersize=7,
               label="Icosidodecahedron (30)"),
        Line2D([0], [0], marker="D", color="w", markerfacecolor=NEW_COLOR,
               markeredgecolor="white", markersize=6,
               label="Expansion points (20)"),
    ]
    # place legend in the lower-right corner (inside the figure, not above)
    ax.legend(handles=legend_elements, loc="lower right", fontsize=7,
              framealpha=0.90, borderpad=0.5, handletextpad=0.4,
              labelspacing=0.35, edgecolor="#cccccc")

    ax.set_axis_off()
    fig.subplots_adjust(left=-0.05, right=1.05, bottom=-0.05, top=1.05)

    out_path = os.path.normpath(
        os.path.join(os.path.dirname(__file__), "..", "figures",
                     "counterexample1_sphere.pdf"))
    fig.savefig(out_path, format="pdf", bbox_inches="tight",
                pad_inches=0.01, dpi=300)
    print(f"Saved: {out_path}")

    png_path = out_path.replace(".pdf", ".png")
    fig.savefig(png_path, format="png", bbox_inches="tight",
                pad_inches=0.01, dpi=220)
    print(f"Saved: {png_path}")


if __name__ == "__main__":
    main()
