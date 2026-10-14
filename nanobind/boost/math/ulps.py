#  Copyright 2026 Nick Thompson
#  Distributed under the Boost Software License, Version 1.0.
#  https://www.boost.org/LICENSE_1_0.txt
"""ulps plots and heatmaps from boost::math::tools::ulps_plot and ulps_heatmap; the reference is the same
function at cpp_bin_float_50, or reference= any Python callable (an mpmath function, say). For example:

    import boost.math as bm; from boost.math import ulps
    ulps.plot(bm.digamma, -5, 5).save("digamma.svg")
    ulps.render_to_png(ulps.heatmap(bm.lgamma_p, (1, 1e7), (0.01, 1e7), log_x=True, log_y=True)[0], "lgamma_p.png")
"""

import math
import os
import re
import shutil
import subprocess
import tempfile

from . import _ext


class Svg:
    """An SVG: str() is its text, summary the C++ summary(); shown inline in Jupyter."""

    def __init__(self, text, summary):
        self.text = text
        self.summary = summary

    def __str__(self):
        return self.text

    def _repr_svg_(self):
        return self.text

    def save(self, path):
        with open(path, "w", encoding="utf-8") as f:
            f.write(self.text)
        return path


def _name(f):
    if not hasattr(_ext, "_ulps_plot"):
        raise RuntimeError("boost.math.ulps was built without Boost.Multiprecision or zlib; see nanobind/CMakeLists.txt")
    if isinstance(f, str):
        return f
    if getattr(_ext, getattr(f, "__name__", ""), None) is not f:
        raise TypeError("ulps plots need a boost.math function or its name; arbitrary Python callables are not supported")
    return f.__name__


def _python_reference(reference):
    """Wraps reference for the C++ side: arguments come as floats, or as 55-digit strings (made mpmath numbers here)
    where they aren't doubles, and the result goes back as a float or a decimal string; NaN if it raises or isn't real."""
    try:
        import mpmath as mp
    except ImportError:
        raise ImportError("reference= needs mpmath: pip install mpmath") from None

    def call(*args):
        try:
            r = reference(*[mp.mpf(a) if isinstance(a, str) else a for a in args])
            if isinstance(r, mp.mpf):
                return mp.nstr(r, 55) if mp.isfinite(r) else float(r)
            return float(r)
        except Exception:
            return math.nan

    return call


def _run(fn, reference, threads):
    if reference is None:
        return fn(None, threads)
    import mpmath as mp
    call = _python_reference(reference)
    # The Python calls are serialized anyway; more threads only add GIL churn:
    with mp.workdps(60):
        return fn(call, 1)


def plot(f, a, b, *, samples=10000, log_x=False, title=None, reference=None, clip=None, show_function=False,
         threads=0):
    """The ulps plot of f on [a, b]; threads=0 uses them all. reference= an mpmath-based callable replaces
    the cpp_bin_float_50 reference, and runs on one thread."""
    name = _name(f)
    svg, summary = _run(lambda ref, th: _ext._ulps_plot(name, a, b, samples, log_x, name if title is None else title, ref,
                                                        clip or -1, show_function, th), reference, threads)
    return Svg(svg, summary)


def heatmap(f, x_range, y_range, *, width=800, height=500, log_x=False, log_y=False, title=None, x_label="x",
            y_label="y", log_scale=False, clip=None, reference=None, show_function=False, threads=0):
    """The ulps heatmap of f(x, y) over x_range by y_range: (ulps, ulps/envelope), two Svgs. reference= as for plot."""
    name = _name(f)
    ulps, envelope, summary = _run(
        lambda ref, th: _ext._ulps_heatmap(name, *x_range, *y_range, width, height, log_x, log_y,
                                           name if title is None else title, x_label, y_label, log_scale, clip or -1,
                                           ref, show_function, th), reference, threads)
    return Svg(ulps, summary), Svg(envelope, summary)


def _chrome():
    candidates = [os.environ.get("CHROME"), "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"]
    candidates += [shutil.which(c) for c in ("google-chrome", "google-chrome-stable", "chromium", "chromium-browser")]
    return next((c for c in candidates if c and os.path.isfile(c)), None)


def render_to_png(svg, path, *, width=None, height=None):
    """Renders an Svg (or SVG text) to a PNG with headless Chrome, or cairosvg. Returns path."""
    text = str(svg)
    size = lambda attr: int(float(re.search(r"<svg[^>]*?\b" + attr + r"=['\"]([0-9.]+)", text).group(1)))
    width = width or size("width")
    height = height or size("height")
    path = os.path.abspath(path)
    chrome = _chrome()
    if chrome:
        with tempfile.TemporaryDirectory() as tmp:
            src = os.path.join(tmp, "plot.svg")
            with open(src, "w", encoding="utf-8") as f:
                f.write(text)
            subprocess.run([chrome, "--headless", "--disable-gpu", "--hide-scrollbars", f"--screenshot={path}",
                            f"--window-size={width},{height}", "file://" + src],
                           check=True, capture_output=True, timeout=120)
        return path
    try:
        import cairosvg
    except ImportError:
        raise RuntimeError("render_to_png needs Google Chrome (or set $CHROME) or `pip install cairosvg`") from None
    cairosvg.svg2png(bytestring=text.encode(), write_to=path, output_width=width, output_height=height)
    return path
