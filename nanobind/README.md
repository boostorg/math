<!-- Copyright 2026 Nick Thompson
     Distributed under the Boost Software License, Version 1.0.
     https://www.boost.org/LICENSE_1_0.txt -->
# boost-math

Python bindings for [Boost.Math](https://www.boost.org/libs/math), built with nanobind against Boost.Math's headers.
This is an early development release; names may change.

- `boost.math`: the special functions (gamma, beta, erf, Bessel, elliptic, hypergeometric, polynomials, ...), flat
  and elementwise over NumPy arrays; `threads=` splits large arrays across threads (0: all cores).
- Submodules: `distributions`, `quadrature`, `optimize`, `fft`, `roots`, `statistics`, `interpolators`, `constants`,
  `differentiation`, `expansions`; `boost.math.ulps` plots the error of a function in ulps.
- `errstate`, `seterr`, `geterr`: what domain, overflow and evaluation errors do (raise, warn or ignore).
- Docstrings link to the Boost.Math documentation; the wheel ships type stubs.

```
pip install boost-math
```

```python
import math
import numpy as np
import boost.math as bm
from boost.math.distributions import normal_distribution
from boost.math.quadrature import tanh_sinh

bm.tgamma(4.5), bm.cyl_bessel_j(0, np.linspace(0, 10, 10**6), threads=0)
d = normal_distribution(1.0, 2.0)
d.cdf(0.5), d.quantile(0.975)
integral, error, l1 = tanh_sinh(lambda x: 1 / math.sqrt(x), 0, 1)
with bm.errstate(domain="ignore"):
    bm.lgamma([-1.0, 0.5])  # nan, 0.572...
```

Wheels are built for CPython 3.12 and later (one abi3 wheel per platform) on Linux, macOS and Windows.
Distributed under the [Boost Software License, Version 1.0](https://www.boost.org/LICENSE_1_0.txt).
