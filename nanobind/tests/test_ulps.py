#  Copyright 2026 Nick Thompson
#  Distributed under the Boost Software License, Version 1.0.
#  https://www.boost.org/LICENSE_1_0.txt
import importlib.util

import mpmath as mp
import pytest

import boost.math as bm
from boost.math import ulps

pytestmark = pytest.mark.skipif(not hasattr(bm._ext, "_ulps_plot"), reason="built without boost.math.ulps")


def check(svg):
    assert str(svg).startswith(("<svg", "<?xml"))
    assert len(str(svg)) > 10000
    assert svg.summary


def test_plot(tmp_path):
    svg = ulps.plot(bm.digamma, -5, 5, samples=2000)
    check(svg)
    assert "samples=2000" in svg.summary
    assert svg._repr_svg_() == str(svg)
    path = tmp_path / "digamma.svg"
    svg.save(path)
    assert path.read_text() == str(svg)


def test_heatmap():
    u, e = ulps.heatmap("lgamma_p", (1, 1e3), (0.01, 1e3), width=40, height=30, log_x=True, log_y=True)
    check(u)
    check(e)
    assert str(u) != str(e)
    assert "1200 pixels" in u.summary


def test_python_reference():
    svg = ulps.plot(bm.erf, -3, 3, samples=200, reference=mp.erf)
    check(svg)
    assert "worst_ulps=" in svg.summary
    # A reference that raises leaves the points without a value:
    svg = ulps.plot(bm.erf, -3, 3, samples=50, reference=lambda x: 1 / 0)
    assert "undefined=50" in svg.summary


def test_reference_interrupt():
    def ref(x):
        raise KeyboardInterrupt

    with pytest.raises(KeyboardInterrupt):
        ulps.plot(bm.erf, -3, 3, samples=50, reference=ref)


def test_only_bound_functions():
    def digamma(x):
        return bm.digamma(x)

    with pytest.raises(TypeError, match="arbitrary Python callables"):
        ulps.plot(digamma, 1, 2)


def test_unknown_name():
    with pytest.raises(ValueError, match="no ulps reference for 'nonesuch'"):
        ulps.plot("nonesuch", 0, 1)
    with pytest.raises(ValueError, match="no ulps reference for 'digamma'"):
        ulps.heatmap(bm.digamma, (0, 1), (0, 1))


@pytest.mark.skipif(ulps._chrome() is None and importlib.util.find_spec("cairosvg") is None,
                    reason="needs Chrome or cairosvg")
def test_render_to_png(tmp_path):
    path = ulps.render_to_png(ulps.plot(bm.digamma, 1, 2, samples=100), tmp_path / "digamma.png")
    with open(path, "rb") as f:
        assert f.read(8) == b"\x89PNG\r\n\x1a\n"
