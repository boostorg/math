#  Copyright 2026 Nick Thompson
#  Distributed under the Boost Software License, Version 1.0.
#  https://www.boost.org/LICENSE_1_0.txt
"""Dev-only Boost.Math bindings; build instructions are in nanobind/CMakeLists.txt."""

import contextlib

from ._ext import *  # noqa: F401,F403
from . import _ext

__all__ = [n for n in dir(_ext) if not n.startswith("_")] + ["errstate", "geterr", "seterr"]

_CATEGORIES = ("domain", "overflow", "evaluation")
_ACTIONS = ("default", "ignore", "warn", "raise")


def geterr():
    """How this thread handles each category of error; see seterr.

    https://www.boost.org/doc/libs/release/libs/math/doc/html/math_toolkit/error_handling.html
    """
    return dict(zip(_CATEGORIES, (_ACTIONS[a] for a in _ext._geterr())))


def seterr(*, all=None, domain=None, overflow=None, evaluation=None):
    """Sets, for this thread, what errors do: 'raise', 'warn', 'ignore' or 'default', and returns the old settings.

    domain covers pole, indeterminate-result and rounding errors too. By default scalars raise (ValueError, or
    OverflowError and RuntimeError for overflow and evaluation errors), and arrays warn once per call and set the
    failed elements to NaN, except that overflows silently give an infinity.

    https://www.boost.org/doc/libs/release/libs/math/doc/html/math_toolkit/error_handling.html
    """
    old = geterr()
    new = dict(old)
    for category, action in zip(_CATEGORIES, (domain, overflow, evaluation)):
        action = action if action is not None else all
        if action is not None:
            if action not in _ACTIONS:
                raise ValueError("%s=%r: expected one of %s" % (category, action, ", ".join(map(repr, _ACTIONS))))
            new[category] = action
    _ext._seterr(*(_ACTIONS.index(new[c]) for c in _CATEGORIES))
    return old


@contextlib.contextmanager
def errstate(**kwargs):
    """Context manager for seterr: `with errstate(domain='ignore'): ...`.

    https://www.boost.org/doc/libs/release/libs/math/doc/html/math_toolkit/error_handling.html
    """
    old = seterr(**kwargs)
    try:
        yield
    finally:
        seterr(**old)
