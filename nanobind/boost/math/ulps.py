#  Copyright 2026 Nick Thompson
#  Distributed under the Boost Software License, Version 1.0.
#  https://www.boost.org/LICENSE_1_0.txt
"""ulps plots against an mpmath reference.

    import boost.math as bm, mpmath as mp
    from boost.math import ulps
    ulps.plot(bm.digamma, mp.digamma, np.linspace(-5, 5, 2000))
    ulps.heatmap(bm.lgamma_p, lambda a, x: mp.log(mp.gammainc(a, 0, x, regularized=True)),
                 np.logspace(0, 7, 300), np.linspace(0.01, 0.99, 300), y_of=lambda a, r: a * r)

Each point is classified as one of:
  ulps        |f - ref| / ulp(ref), both finite;
  overflow    f overflowed (raised or returned +-inf) where ref is finite: shown cyan;
  nan         f raised a domain error or returned NaN where ref is finite: shown green;
  handled     ref isn't a finite double and f agrees (inf, or underflow to 0): shown gray.
"""

import math

import mpmath as mp
import numpy as np

OVERFLOW, NAN, HANDLED = -1.0, -2.0, -3.0
SPECIAL = {OVERFLOW: ("#2ec4d6", "spurious overflow"), NAN: ("#4caf50", "NaN / domain error"),
           HANDLED: ("#8a8a85", "out of range, handled")}
SURFACE, INK, MUTED = "#1a1a19", "#f0efec", "#c3c2b7"


def error(f, ref, *args, dps=50):
    """The error of f(*args) against the mpmath reference ref(*args): ulps, or a SPECIAL code."""
    with mp.workdps(dps):
        try:
            r = ref(*[mp.mpf(a) for a in args])
        except Exception:
            return math.nan  # No reference here; left blank in the plots.
        try:
            got = f(*args)
        except OverflowError:
            got = math.inf
        except (ValueError, ArithmeticError):
            got = math.nan
        if isinstance(r, mp.mpc) or not mp.isfinite(r):
            return HANDLED if (math.isinf(got) or math.isnan(got)) else float("inf")
        rf = float(r)
        if math.isinf(rf) or (rf == 0 and r != 0):
            # The answer isn't representable: overflow to inf or underflow to 0 is right.
            return HANDLED if (got == rf or (rf == 0 and abs(got) <= np.finfo(float).tiny)) else float(
                abs(mp.mpf(got) - r) / np.spacing(np.finfo(float).tiny))
        if math.isinf(got):
            return OVERFLOW
        if math.isnan(got):
            return NAN
        return float(abs(mp.mpf(got) - r) / np.spacing(abs(rf)) if rf != 0 else abs(got) / np.spacing(0.0))


def _summary(err):
    e = err[err >= 0]
    missing = int(np.sum(np.isnan(err)))
    counts = {name: int(np.sum(err == code)) for code, (_, name) in SPECIAL.items() if np.any(err == code)}
    parts = [f"max {e.max():.3g} ulps, median {np.median(e):.2g}" if e.size else "no finite errors"]
    parts += [f"{n} {name}" for name, n in counts.items()]
    if missing:
        parts.append(f"{missing} without a reference")
    return "; ".join(parts)


def _style(ax, fig):
    fig.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)
    ax.tick_params(colors=MUTED)
    for s in ax.spines.values():
        s.set_color(MUTED)


def _decades(cb):
    from matplotlib.ticker import FuncFormatter
    cb.ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: "$10^{%d}$" % round(np.log10(v))))
    cb.ax.tick_params(colors=MUTED)
    cb.outline.set_edgecolor(MUTED)
    cb.set_label("error in ulps", color=INK)


def plot(f, ref, x, title=None, xlabel="x", xscale="linear", out=None, dps=50, reference=None):
    """1D: error in ulps of f(x) against ref(x). Returns (fig, err)."""
    import matplotlib.pyplot as plt
    x = np.asarray(x, dtype=float)
    err = np.array([error(f, ref, xi, dps=dps) for xi in x])
    fig, ax = plt.subplots(figsize=(8, 4.5))
    _style(ax, fig)
    ok = err >= 0
    ax.scatter(x[ok], np.maximum(err[ok], 0.01), s=3, color="#f5a524", rasterized=True)
    top = np.max(err[ok]) if ok.any() else 1
    for code, (color, name) in SPECIAL.items():
        m = err == code
        if m.any():
            ax.scatter(x[m], np.full(m.sum(), max(top, 1) * 2), s=6, color=color, label=name, rasterized=True)
    ax.set_xscale(xscale)
    ax.set_yscale("log")
    ax.set_xlabel(xlabel, color=INK)
    ax.set_ylabel("error in ulps", color=INK)
    ax.set_title(title or getattr(f, "__name__", ""), color=INK)
    if any(np.any(err == c) for c in SPECIAL):
        ax.legend(frameon=False, labelcolor=INK)
    fig.text(0.5, 0.01, f"Reference: {reference or f'mpmath at {dps} digits'}; {_summary(err)}.", ha="center",
             color=MUTED, fontsize=8, wrap=True)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    if out:
        fig.savefig(out, dpi=300, facecolor=SURFACE)
    return fig, err


def heatmap(f, ref, x, y, title=None, xlabel="x", ylabel="y", xscale="log", yscale="log",
            y_of=None, out=None, dps=50, reference=None):
    """2D: error in ulps of f(x, y) over the grid x by y. With y_of, the second argument is
    y_of(x, y) instead of y (e.g. y_of=lambda a, r: a * r to plot against x/a). Returns (fig, err)."""
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap, LogNorm
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    err = np.empty((y.size, x.size))
    for i, yi in enumerate(y):
        for j, xj in enumerate(x):
            err[i, j] = error(f, ref, xj, y_of(xj, yi) if y_of else yi, dps=dps)
    fig, ax = plt.subplots(figsize=(7.2, 6))
    _style(ax, fig)
    finite = np.where(err >= 0, np.maximum(err, 0.1), np.nan)
    vmax = max(np.nanmax(finite) if np.isfinite(finite).any() else 1, 1)
    im = ax.pcolormesh(x, y, finite, norm=LogNorm(vmin=0.1, vmax=vmax), cmap="inferno", shading="nearest",
                       rasterized=True)
    handles = []
    for code, (color, name) in SPECIAL.items():
        m = err == code
        if m.any():
            ax.pcolormesh(x, y, np.where(m, 1.0, np.nan), cmap=ListedColormap([color]), shading="nearest",
                          rasterized=True)
            handles.append(plt.Rectangle((0, 0), 1, 1, color=color, label=name))
    ax.set_xscale(xscale)
    ax.set_yscale(yscale)
    ax.set_xlabel(xlabel, color=INK, fontsize=12)
    ax.set_ylabel(ylabel, color=INK, fontsize=12)
    ax.set_title(title or getattr(f, "__name__", ""), color=INK, fontsize=12)
    if handles:
        ax.legend(handles=handles, frameon=False, labelcolor=INK, loc="upper left", fontsize=8)
    _decades(fig.colorbar(im, ax=ax))
    fig.text(0.5, 0.01, f"Reference: {reference or f'mpmath at {dps} digits'}; {y.size} x {x.size} grid.\n"
             f"{_summary(err)}; errors below 0.1 ulp shown as 0.1.", ha="center", color=MUTED, fontsize=8)
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    if out:
        fig.savefig(out, dpi=300, facecolor=SURFACE)
    return fig, err
