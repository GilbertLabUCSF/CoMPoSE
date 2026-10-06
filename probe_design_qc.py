"""QC figures for custom probe design.

plot_metrics() draws the QC figure for one target; qc_plot() calls it
from already-loaded dataframes. Metric functions and GC_TARGET are imported from
probe_design, so the figure cannot disagree with the pipeline about what passed
a filter.

Changes in this revision:
  * exclude_regions: sense-strand (start, end) pairs, 1-based inclusive, for
    windows removed by gap-fill masking. They are drawn on their own track and
    no longer counted as eligible / nominated.
  * close: the figure is closed before returning, so Jupyter displays it once
    and loops over many targets do not leak open figures.
  * detailed=True (opt-in): "Candidate, hit" / "Candidate, no hit" labels, dotted
    guides at selected probe centres, ringed selected hits. Default is the compact
    paper style. Selection penalises but does not exclude off-target hits, so a
    selected probe can sit on either the Off-target or Nominated row.
  * regions_from_gapfill(): builds exclude_regions from the gap-fill sheet.
"""

import math

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

from probe_design import (
    GC_TARGET,
    calculate_entropy as ent,
    calculate_gc_content as gc,
    max_homopolymer_length as max_hp,
    score_sequence_repetition as rep,
)

__all__ = ["plot_metrics", "qc_plot", "score_windows", "build_offtarget_lookup",
           "fill_fail", "draw_track", "regions_from_gapfill"]

plt.rcParams.update({
    "font.size":         8,
    "axes.labelsize":    8,
    "axes.titlesize":    9,
    "xtick.labelsize":   7,
    "ytick.labelsize":   7,
    "legend.fontsize":   7,
    "lines.linewidth":   0.8,
    "axes.linewidth":    0.5,
    "xtick.major.width": 0.5,
    "ytick.major.width": 0.5,
    "xtick.major.size":  2,
    "ytick.major.size":  2,
    "svg.fonttype":      "none",
})

GC_DIFF     = 14
HOMOPOLYMER = 4
ENTROPY_MIN = 0.8
REP_MAX     = 125

ARM_LEN   = 25
PROBE_LEN = 50
WINDOW_CENTER = 24.5


SWAP_ARM_LABELS = True

TRACK_FOOTPRINT = False

MAX_LABEL_CLUSTERS = 20

_KEYS = ["pos", "lhsGC", "rhsGC", "lhsHP", "rhsHP", "lhsEnt", "rhsEnt",
         "lhsRep", "rhsRep", "passes", "ct_pass", "gc_pass", "hp_pass",
         "en_pass", "rp_pass"]


def rc(s):
    """Reverse complement, uppercased."""
    return str(s).upper().translate(str.maketrans("ACGT", "TGCA"))[::-1]


def score_windows(seq):
    """Per-window arm metrics and filter outcomes, in reverse-complement space.

    'pos' is the 1-based window start, matching the annotation table's 'start'.
    Each 25mer is scored once and cached: the RHS arm of window i is the LHS arm
    of window i+25.
    """
    seq = rc(seq)
    n = len(seq) - PROBE_LEN + 1
    if n <= 0:
        return {k: np.array([]) for k in _KEYS}

    cache = {}

    def arm(i):
        if i not in cache:
            s = seq[i:i + ARM_LEN]
            cache[i] = (gc(s), max_hp(s), ent(s), rep(s))
        return cache[i]

    rows = []
    for i in range(n):
        first, second = arm(i), arm(i + ARM_LEN)
        (lgc, lhp, le, lr), (rgc, rhp, re_, rr) = (
            (second, first) if SWAP_ARM_LABELS else (first, second)
        )

        ct    = seq[i + ARM_LEN - 1] == "T"
        gc_ok = abs(lgc - GC_TARGET) <= GC_DIFF and abs(rgc - GC_TARGET) <= GC_DIFF
        hp_ok = lhp <= HOMOPOLYMER and rhp <= HOMOPOLYMER
        en_ok = le >= ENTROPY_MIN and re_ >= ENTROPY_MIN
        rp_ok = lr <= REP_MAX and rr <= REP_MAX

        rows.append((i + 1, lgc, rgc, lhp, rhp, le, re_, lr, rr,
                     ct and gc_ok and hp_ok and en_ok and rp_ok,
                     ct, gc_ok, hp_ok, en_ok, rp_ok))

    return {k: np.array(v) for k, v in zip(_KEYS, zip(*rows))}


# ── fill_fail ─────────────────────────────────────────────────────────────────
def fill_fail(ax, pos, y, threshold, direction):
    y = np.asarray(y, dtype=float)
    if direction == "max":
        ax.fill_between(pos, y, threshold, where=y > threshold,
                        color="#E24B4A", alpha=0.20, lw=0)
    elif direction == "min":
        ax.fill_between(pos, y, threshold, where=y < threshold,
                        color="#E24B4A", alpha=0.20, lw=0)
    elif direction == "band":
        lo, hi = threshold
        ax.fill_between(pos, y, lo, where=y < lo, color="#E24B4A", alpha=0.20, lw=0)
        ax.fill_between(pos, y, hi, where=y > hi, color="#E24B4A", alpha=0.20, lw=0)


# ── off-target lookup ─────────────────────────────────────────────────────────
def build_offtarget_lookup(blast_df):
    """Map window start -> [(side, gene, mismatches, matched_bases), ...].

    Keyed on bin start alone, so blast_df must be restricted to one target;
    plot_metrics() does that for you.

    Mismatch counts are cast to Python int because the plotting code tests
    isinstance(mm, int), which is False for numpy int64 and float64.
    """
    sides = [
        ("LHS", "lhs_off_target_gene", "lhs_mismatch", "lhs_max_off_target_binding"),
        ("RHS", "rhs_off_target_gene", "rhs_mismatch", "rhs_max_off_target_binding"),
    ]
    if SWAP_ARM_LABELS:
        sides = [("LHS", *sides[1][1:]), ("RHS", *sides[0][1:])]

    lookup = {}
    for _, row in blast_df.iterrows():
        try:
            start = int(str(row["id"]).rsplit("_bin_", 1)[-1])
        except (ValueError, KeyError):
            continue
        hits = []
        for side, gene_col, mm_col, len_col in sides:
            gene = row.get(gene_col)
            if not isinstance(gene, str):
                continue
            mm = row.get(mm_col)
            matched = row.get(len_col)
            hits.append((side, gene,
                         int(mm) if pd.notna(mm) else "?",
                         int(matched) if pd.notna(matched) else "?"))
        if hits:
            lookup[start] = hits
    return lookup


def draw_track(ax, pos, mask, color, alpha, y_center, height, to_sense=None,
               footprint=TRACK_FOOTPRINT):
    """Fill rectangles across positions where mask is True.

    Parameter order matches the call site: y_center before height. `to_sense`
    maps reverse-complement positions onto sense-strand coordinates; without it
    the row is drawn mirrored relative to the axis.
    """
    if to_sense is None:
        def to_sense(x):
            return x
    ylo, yhi = y_center - height / 2, y_center + height / 2

    idx = np.flatnonzero(np.asarray(mask))
    if idx.size == 0:
        return

    if footprint:
        per_window = alpha if idx.size == 1 else max(alpha / 6, 0.10)
        for p in np.asarray(pos)[idx]:
            ax.fill_betweenx([ylo, yhi], to_sense(p - 0.5),
                             to_sense(p + PROBE_LEN - 0.5),
                             color=color, alpha=per_window, lw=0, zorder=2)
        return

    for run in np.split(idx, np.flatnonzero(np.diff(idx) > 1) + 1):
        ax.fill_betweenx([ylo, yhi],
                         to_sense(pos[run[0]] + WINDOW_CENTER - 0.5),
                         to_sense(pos[run[-1]] + WINDOW_CENTER + 0.5),
                         color=color, alpha=alpha, lw=0, zorder=2)


def _subset(df, seq_id):
    """Restrict a table to one target, by 'target' column or probe id prefix."""
    if df is None or seq_id is None:
        return df
    if "target" in df.columns:
        sub = df[df["target"] == seq_id]
        if not sub.empty:
            return sub
    if "id" in df.columns:
        sub = df[df["id"].astype(str).str.rsplit("_bin_", n=1).str[0] == seq_id]
        if not sub.empty:
            return sub
    return df


def regions_from_gapfill(path, sheet, gene):
    """Exclusion regions for one gene from the gap-fill sheet, as a list of
    (start, end) sense-strand pairs (the sheet's own coordinates), for the
    exclude_regions argument of plot_metrics()."""
    g = pd.read_excel(path, sheet_name=sheet)
    g = g[g["gene"] == gene]
    return [(int(s), int(e)) for s, e in zip(g["start"], g["end"])]


# ── main plot — compact, 6.2in wide ──────────────────────────────────────────
def plot_metrics(seq, seq_id="sequence", blast_df=None, top_df=None,
                 save_path=None, svg_path=None, show_offtarget_labels=False,
                 show=False, exclude_regions=None, close=True, detailed=False):
    """QC figure for one target.

    blast_df and top_df may be the full multi-target tables; they are restricted
    to seq_id here. That matters because the off-target lookup and the selected
    mask key on bin start alone, so unfiltered tables silently mix targets.

    Parameters:
    - exclude_regions: list of (start, end) in sense-strand coordinates, 1-based
      inclusive (the gap-fill sheet's convention). Any passing window overlapping
      a region is drawn on an 'Excluded' track and removed from the eligible,
      off-target and nominated counts.
    - detailed: False (default) gives the compact paper style: tracks labelled
      Off-target / Nominated / Selected and the original title. True adds
      "Candidate, hit" / "Candidate, no hit" labels, dotted guides at selected
      probe centres, ringed selected hits, and a fuller title.
    - close: close the figure before returning (default True). The returned
      Figure still renders once when it is the last expression in a notebook cell.
    """
    blast_df = _subset(blast_df, seq_id)
    top_df = _subset(top_df, seq_id)

    d = score_windows(seq)
    pos, passes = d["pos"], d["passes"]
    if pos.size == 0:
        raise ValueError(
            f"{seq_id}: sequence is {len(seq)} bp, shorter than one {PROBE_LEN} bp probe."
        )

    L = len(seq)

    def to_sense(x):
        return (L + 1) - x

    excluded_mask = np.zeros(len(pos), dtype=bool)
    for gs, ge in (exclude_regions or []):
        excluded_mask |= (L - pos - 48 <= ge) & (gs <= L + 1 - pos)
    excluded_mask &= passes.astype(bool)     # only matters for windows that passed
    eligible = passes.astype(bool) & ~excluded_mask

    ot_lookup = build_offtarget_lookup(blast_df) if blast_df is not None else {}
    excluded_starts = set(pos[excluded_mask].tolist())
    ot_lookup = {s: h for s, h in ot_lookup.items() if s not in excluded_starts}
    ot_mask = np.array([p in ot_lookup for p in pos]) & eligible

    selected_starts = set(top_df["start"].values) if top_df is not None else set()
    selected_mask = np.array([p in selected_starts for p in pos])

    if (selected_mask & excluded_mask).any():
        print(f"WARNING {seq_id}: selected probe(s) overlap an excluded region.")

    if (selected_mask.sum() < len(selected_starts)) or (selected_mask & ~passes.astype(bool)).any():
        print(f"WARNING {seq_id}: some selected probes do not correspond to a window "
              f"that passes the filters in this sequence -- the table and sequence "
              f"may not match, or thresholds differ from nominate_probes().")

    OT_LABEL, NOM_LABEL = (("Candidate, hit", "Candidate, no hit") if detailed
                           else ("Off-target", "Nominated"))

    FAIL_C = "#A9525A"
    TRACK_COLORS = {
        "Ligation site":  FAIL_C,
        "GC content":     FAIL_C,
        "Homopolymer":    FAIL_C,
        "Entropy":        FAIL_C,
        "Repetitiveness": FAIL_C,
        "Excluded":       "#777777",
        OT_LABEL:         "#D4537E",
        NOM_LABEL:        "#185FA5",
        "Selected":       "#0C447C",
    }
    SEL_C = "#0C447C"

    def shade_passing(ax):
        idx = np.flatnonzero(eligible)
        if idx.size == 0:
            return
        for run in np.split(idx, np.flatnonzero(np.diff(idx) > 1) + 1):
            ax.axvspan(to_sense(pos[run[0]] + WINDOW_CENTER - .5),
                       to_sense(pos[run[-1]] + WINDOW_CENTER + .5),
                       color="#185FA5", alpha=0.10, lw=0)

    metrics = [
        ("GC (%)",         d["lhsGC"],  d["rhsGC"],
         [GC_TARGET - GC_DIFF, GC_TARGET + GC_DIFF], "band"),
        ("Homopolymer",    d["lhsHP"],  d["rhsHP"],  HOMOPOLYMER, "max"),
        ("Entropy",        d["lhsEnt"], d["rhsEnt"], ENTROPY_MIN, "min"),
        ("Repetitiveness", d["lhsRep"], d["rhsRep"], REP_MAX,     "max"),
    ]

    heights = [0.5, 0.5, 0.5, 0.5, 0.75, 1.2 if not exclude_regions else 1.35]
    fig, axes = plt.subplots(6, 1,
                             figsize=(6.2, 7.6),
                             sharex=True,
                             gridspec_kw={"hspace": 0.05,
                                          "height_ratios": heights,
                                          "top": 0.93,
                                          "bottom": 0.13,
                                          "left": 0.17,
                                          "right": 0.97})
    fig.patch.set_facecolor("white")

    x_center = to_sense(pos + WINDOW_CENTER)

    for ax, (ylabel, lhs, rhs, thresh, direction) in zip(axes[:4], metrics):
        shade_passing(ax)
        fill_fail(ax, x_center, lhs, thresh, direction)
        fill_fail(ax, x_center, rhs, thresh, direction)
        ax.plot(x_center, lhs, color="#1A1A1A", lw=0.8)
        ax.plot(x_center, rhs, color="#AAAAAA", lw=0.8, alpha=0.9)
        if direction == "band":
            ax.axhline(thresh[0], color="#E24B4A", lw=0.9, ls="--")
            ax.axhline(thresh[1], color="#E24B4A", lw=0.9, ls="--")
            ax.axhline(GC_TARGET, color="#888", lw=0.6, ls=":", alpha=0.6)
            bounds = [thresh[0], thresh[1]]
        else:
            ax.axhline(thresh, color="#E24B4A", lw=0.9, ls="--")
            bounds = [thresh]

        lo = min(bounds + [np.nanmin(lhs), np.nanmin(rhs)])
        hi = max(bounds + [np.nanmax(lhs), np.nanmax(rhs)])
        pad = max((hi - lo) * 0.12, 1e-6)
        ax.set_ylim(lo - pad, hi + pad)
        ax.set_ylabel(ylabel, labelpad=2)
        ax.yaxis.set_major_locator(
            ticker.MaxNLocator(nbins=3, integer=(direction != "band")))
        ax.grid(axis="y", color="#EEEEEE", lw=0.4)
        ax.spines[["top", "right"]].set_visible(False)

    ax_ot = axes[4]
    ax_ot.spines[["top", "right"]].set_visible(False)
    shade_passing(ax_ot)

    all_hits = [(to_sense(start + WINDOW_CENTER), side, gene, mm, matched / ARM_LEN)
                for start, hits in ot_lookup.items()
                for side, gene, mm, matched in hits
                if isinstance(matched, int)]

    if all_hits:
        h_min = min(h for *_, h in all_hits)
        y_data_bottom = max(0.0, math.floor((h_min - 0.03) * 20) / 20)
    else:
        y_data_bottom = 0.8
    y_data_top = 1.0

    x_span = abs(x_center[0] - x_center[-1]) or 1.0
    tier_h = 0.045
    label_y0 = y_data_top + 0.025
    placements = []

    if show_offtarget_labels and all_hits:
        cluster_gap = 0.02 * x_span
        by_key = {}
        for x, side, gene, mm, h in sorted(all_hits, key=lambda t: t[0]):
            by_key.setdefault((side, gene), []).append((x, mm, h))

        clustered = []
        for (side, gene), entries in by_key.items():
            group = [entries[0]]
            for e in entries[1:]:
                if e[0] - group[-1][0] <= cluster_gap:
                    group.append(e)
                else:
                    clustered.append((side, gene, group))
                    group = [e]
            clustered.append((side, gene, group))

        if len(clustered) > MAX_LABEL_CLUSTERS:
            print(f"{seq_id}: {len(clustered)} off-target clusters; labels suppressed "
                  f"(raise MAX_LABEL_CLUSTERS to override).")
        else:
            merged_hits = []
            for side, gene, group in clustered:
                xs = [g[0] for g in group]
                mms = [g[1] for g in group if isinstance(g[1], int)]
                label_text = (f"{side[0]}:{gene}"
                              + (f" \u00d7{len(group)}" if len(group) > 1 else ""))
                merged_hits.append((xs[len(xs) // 2], min(mms) if mms else None,
                                    max(g[2] for g in group), label_text))

            avg_char_w = 0.011 * x_span
            last_x_per_tier = []
            for rep_x, rep_mm, rep_h, label_text in sorted(merged_hits, key=lambda t: t[0]):
                needed_gap = max(0.025 * x_span, len(label_text) * avg_char_w)
                tier = 0
                while tier < len(last_x_per_tier) and rep_x - last_x_per_tier[tier] < needed_gap:
                    tier += 1
                if tier == len(last_x_per_tier):
                    last_x_per_tier.append(rep_x)
                else:
                    last_x_per_tier[tier] = rep_x
                placements.append((rep_x, rep_mm, rep_h, label_text, tier))

    y_top = (label_y0 + max(t for *_, t in placements) * tier_h + tier_h * 0.7
             if placements else y_data_top + 0.04)

    ax_ot.set_ylim(y_data_bottom, y_top)
    ax_ot.set_yticks([t for t in [0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0]
                      if t >= y_data_bottom - 0.01])
    ax_ot.set_ylabel("Max off-target\nhomology", labelpad=2)
    ax_ot.grid(axis="y", color="#EEEEEE", lw=0.4)

    def ot_color(mm):
        return "#A32D2D" if mm == 0 else "#D4537E" if mm == 1 else "#EF9F27"

    selected_x = {to_sense(s0 + WINDOW_CENTER) for s0 in selected_starts}
    for x, side, gene, mm, h in all_hits:
        chosen = detailed and x in selected_x
        ax_ot.scatter(x, h, color=ot_color(mm), s=22 if chosen else 13, zorder=4,
                      edgecolors="black" if chosen else "white",
                      linewidths=0.8 if chosen else 0.4)

    for rep_x, rep_mm, rep_h, label_text, tier in placements:
        color = ot_color(rep_mm)
        ly = label_y0 + tier * tier_h
        ax_ot.plot([rep_x, rep_x], [rep_h, ly - 0.008], color=color, lw=0.35, zorder=3)
        ax_ot.annotate(label_text, xy=(rep_x, ly), fontsize=6, color=color,
                       ha="center", va="bottom", zorder=5)

    ax_c = axes[5]
    x_min, x_max = to_sense(pos[-1] + 52), to_sense(pos[0] - 2)
    ax_c.set_xlim(x_min, x_max)
    ax_c.set_ylim(0, 1)
    ax_c.set_yticks([])
    ax_c.spines[["top", "right", "left", "bottom"]].set_visible(False)
    ax_c.set_xlabel("Position on sense strand (bp)", labelpad=2)

    ct, gcp, hp, en, rp = (d["ct_pass"].astype(bool), d["gc_pass"].astype(bool),
                           d["hp_pass"].astype(bool), d["en_pass"].astype(bool),
                           d["rp_pass"].astype(bool))
    tracks = [
        ("Ligation site",  ~ct,                      0.75),
        ("GC content",     ct & ~gcp,                0.75),
        ("Homopolymer",    ct & gcp & ~hp,           0.75),
        ("Entropy",        ct & gcp & hp & ~en,      0.75),
        ("Repetitiveness", ct & gcp & hp & en & ~rp, 0.75),
    ]
    if exclude_regions:
        tracks.append(("Excluded", excluded_mask, 0.75))
    tracks += [
        (OT_LABEL,         ot_mask,                  0.80),
        (NOM_LABEL,        eligible & ~ot_mask,      0.70),
        ("Selected",       selected_mask,            1.00),
    ]

    n_tracks = len(tracks)
    pad_top = pad_bot = 0.03
    gap = 0.008
    track_h = (1.0 - pad_top - pad_bot - gap * (n_tracks - 1)) / n_tracks
    label_x = x_min - 1

    for i, (label, mask, alpha) in enumerate(tracks):
        y_center = (1.0 - pad_top - i * (track_h + gap)) - track_h / 2
        color = TRACK_COLORS[label]

        if label == "Selected":
            bar_h = track_h * 0.88
            for p in pos[mask]:
                ax_c.fill_betweenx([y_center - bar_h / 2, y_center + bar_h / 2],
                                   to_sense(p - 0.5), to_sense(p + PROBE_LEN - 0.5),
                                   facecolor=SEL_C, edgecolor="white", lw=0.8,
                                   zorder=3, alpha=1.0)
        else:
            draw_track(ax_c, pos, mask, color, alpha, y_center, track_h * 0.88,
                       to_sense=to_sense)

        if not mask.any():
            ax_c.text(to_sense((pos[0] + pos[-1]) / 2), y_center, "none",
                      ha="center", va="center", fontsize=5,
                      color="#AAAAAA", style="italic", zorder=4)

        ax_c.text(label_x, y_center, label, ha="right", va="center",
                  fontsize=7, color=color, fontweight="normal")

    box_x_left = to_sense(pos[-1] + PROBE_LEN + 0.5)
    box_x_right = to_sense(pos[0] - 0.5)
    box_y_top = 1.0 - pad_top
    box_y_bottom = box_y_top - (n_tracks * track_h + (n_tracks - 1) * gap)
    ax_c.add_patch(mpatches.Rectangle(
        (box_x_left, box_y_bottom), box_x_right - box_x_left,
        box_y_top - box_y_bottom,
        fill=False, edgecolor="black", lw=0.8, zorder=5, clip_on=False))
    for i in range(1, n_tracks):
        y_sep = box_y_top - i * (track_h + gap) + gap / 2
        ax_c.plot([box_x_left, box_x_right], [y_sep, y_sep],
                  color="black", lw=0.6, ls=":", zorder=5, clip_on=False)

    title = (f"{seq_id}  |  {int(passes.sum())}/{len(pos)} pass filters  |  ")
    if exclude_regions:
        title += f"{int(excluded_mask.sum())} excluded  |  "
    if detailed:
        title += (f"{int((eligible & ~ot_mask).sum())} no-hit + {int(ot_mask.sum())} "
                  f"off-target-hit candidates  |  "
                  f"{len(selected_starts)} selected "
                  f"({int((selected_mask & ot_mask).sum())} with hit)")
    else:
        title += (f"{int(ot_mask.sum())} off-target  |  "
                  f"{int((eligible & ~ot_mask).sum())} nominated  |  "
                  f"{len(selected_starts)} selected")
    axes[0].set_title(title, fontsize=7, pad=4)

    legend_handles = [
        Line2D([0], [0], color="#1A1A1A", lw=1.2, label="LHS"),
        Line2D([0], [0], color="#AAAAAA", lw=1.2, label="RHS"),
        Line2D([0], [0], color="#E24B4A", lw=1.0, ls="--", label="Threshold"),
        mpatches.Patch(facecolor="#E24B4A", alpha=0.20, edgecolor="none",
                       label="Out of range"),
        mpatches.Patch(facecolor="#185FA5", alpha=0.18, edgecolor="none",
                       label="Eligible window" if exclude_regions else "Passing window"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#A32D2D",
               markersize=4.5, lw=0, label="0 mismatch"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#D4537E",
               markersize=4.5, lw=0, label="1 mismatch"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#EF9F27",
               markersize=4.5, lw=0, label="2+ mismatch"),
    ]
    fig.legend(handles=legend_handles, loc="lower center",
               bbox_to_anchor=(0.55, 0.015), ncol=5, fontsize=6.5, frameon=False,
               handlelength=1.2, handletextpad=0.3, columnspacing=0.8)

    # Dotted guide at the centre of every selected probe, through all panels, so
    # a selected footprint can be matched to its candidate tick and metric values.
    if detailed:
        for p0 in pos[selected_mask]:
            xc = to_sense(p0 + WINDOW_CENTER)
            for ax in axes:
                ax.axvline(xc, color=SEL_C, lw=0.5, ls=":", alpha=0.7, zorder=1)

    fig.align_ylabels(axes)

    if save_path:
        fig.savefig(save_path, dpi=300)
        print(f"Saved: {save_path}")
    if svg_path:
        fig.savefig(svg_path, format="svg")
        print(f"Saved: {svg_path}")
    if show:
        plt.show()
    if close:
        plt.close(fig)

    return fig


def qc_plot(target, fd_seq, fd_blast_scored, fd_top, svg_path=None, save_path=None,
            show_offtarget_labels=False, show=False, exclude_regions=None,
            close=True, detailed=False):
    """plot_metrics() for one target, from already-loaded dataframes.

    Parameters:
    - target: target name, matching fd_seq['header'] (or 'id') and the 'target'
      column of the probe tables.
    - fd_seq: reference table with 'sequence' and a 'header' or 'id' column.
    - fd_blast_scored: annotated probe table across all targets, after
      calculate_priority_score().
    - fd_top: selected probes across all targets.
    - exclude_regions, close, detailed: see plot_metrics().
    """
    id_col = next((c for c in ("header", "id") if c in fd_seq.columns), None)
    if id_col is None:
        raise KeyError(
            f"fd_seq has neither a 'header' nor an 'id' column -- found: "
            f"{list(fd_seq.columns)}"
        )

    rows = fd_seq.loc[fd_seq[id_col] == target, "sequence"]
    if rows.empty:
        near = [v for v in fd_seq[id_col].astype(str) if target in v]
        raise ValueError(
            f"'{target}' not found in fd_seq['{id_col}']."
            + (f" Close matches: {near[:5]}" if near else
               " Check for stray whitespace or a leading BOM character.")
        )

    nominated = _subset(fd_blast_scored, target)
    selected = _subset(fd_top, target)
    if nominated is None or nominated.empty:
        raise ValueError(f"No rows for target '{target}' in fd_blast_scored.")
    if selected is None or selected.empty:
        raise ValueError(f"No rows for target '{target}' in fd_top.")

    return plot_metrics(
        str(rows.values[0]).upper(), seq_id=target,
        blast_df=nominated, top_df=selected,
        svg_path=svg_path, save_path=save_path,
        show_offtarget_labels=show_offtarget_labels, show=show,
        exclude_regions=exclude_regions, close=close, detailed=detailed,
    )