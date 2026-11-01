// Copyright 2026 Nick Thompson
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/math/special_functions/jacobi_theta.hpp>

#include "bindings.hpp"
#include "macros.hpp"

void bmpy::bind_jacobi_theta(nb::module_& m) {
   BMPY_R2(jacobi_theta1, "z", "q",
      "The Jacobi theta function θ₁(z, q), nome q in (0, 1); as mpmath.jtheta(1, z, q)."
      BMPY_URL("jacobi_theta/jacobi_theta1"));
   BMPY_R2(jacobi_theta2, "z", "q",
      "The Jacobi theta function θ₂(z, q), nome q in (0, 1)."
      BMPY_URL("jacobi_theta/jacobi_theta2"));
   BMPY_R2(jacobi_theta3, "z", "q", "The Jacobi theta function θ₃(z, q), nome q in (0, 1)."
      BMPY_URL("jacobi_theta/jacobi_theta3"));
   BMPY_R2(jacobi_theta4, "z", "q",
      "The Jacobi theta function θ₄(z, q), nome q in (0, 1)."
      BMPY_URL("jacobi_theta/jacobi_theta4"));
   BMPY_R2(jacobi_theta3m1, "z", "q", "θ₃(z, q) − 1, accurate for small q." BMPY_URL("jacobi_theta/jacobi_theta3"));
   BMPY_R2(jacobi_theta4m1, "z", "q", "θ₄(z, q) − 1, accurate for small q." BMPY_URL("jacobi_theta/jacobi_theta4"));
   BMPY_R2(jacobi_theta1tau, "z", "tau",
      "θ₁(z, q) with q = e^(−π tau), tau real (not e^(iπτ)), accurate as q → 1."
      BMPY_URL("jacobi_theta/jacobi_theta1"));
   BMPY_R2(jacobi_theta2tau, "z", "tau", "θ₂(z, q) with q = e^(−π tau), accurate as q → 1; tau real."
      BMPY_URL("jacobi_theta/jacobi_theta2"));
   BMPY_R2(jacobi_theta3tau, "z", "tau", "θ₃(z, q) with q = e^(−π tau), accurate as q → 1; tau real."
      BMPY_URL("jacobi_theta/jacobi_theta3"));
   BMPY_R2(jacobi_theta4tau, "z", "tau", "θ₄(z, q) with q = e^(−π tau), accurate as q → 1; tau real."
      BMPY_URL("jacobi_theta/jacobi_theta4"));
   BMPY_R2(jacobi_theta3m1tau, "z", "tau", "θ₃(z, e^(−π tau)) − 1." BMPY_URL("jacobi_theta/jacobi_theta3"));
   BMPY_R2(jacobi_theta4m1tau, "z", "tau", "θ₄(z, e^(−π tau)) − 1." BMPY_URL("jacobi_theta/jacobi_theta4"));
}
