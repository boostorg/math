#  Copyright 2026 Nick Thompson
#  Distributed under the Boost Software License, Version 1.0.
#  https://www.boost.org/LICENSE_1_0.txt
import math

import mpmath as mp

import boost.math as bm
from boost.math import ulps


def test_error_classification():
    assert ulps.error(bm.tgamma, mp.gamma, 4.5) < 4
    assert ulps.error(lambda x: 7.0, lambda x: mp.mpf(7) + mp.mpf(2) ** -50, 0.0) == 1  # ulp(7) = 2^-50
    assert ulps.error(lambda x: math.inf, lambda x: mp.mpf(1), 0.0) == ulps.OVERFLOW
    assert ulps.error(lambda x: math.nan, lambda x: mp.mpf(1), 0.0) == ulps.NAN
    assert ulps.error(bm.tgamma, mp.gamma, 200.0) == ulps.HANDLED
    assert math.isnan(ulps.error(bm.tgamma, lambda x: 1 / 0, 2.0))
