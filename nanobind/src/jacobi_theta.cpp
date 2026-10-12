// Copyright 2026 Nick Thompson
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/math/special_functions/jacobi_theta.hpp>

#include "bindings.hpp"
#include "macros.hpp"

void bmpy::bind_jacobi_theta(nb::module_& m) {
   BMPY_R2(jacobi_theta1, "z", "q");
   BMPY_R2(jacobi_theta2, "z", "q");
   BMPY_R2(jacobi_theta3, "z", "q");
   BMPY_R2(jacobi_theta4, "z", "q");
   BMPY_R2(jacobi_theta3m1, "z", "q");
   BMPY_R2(jacobi_theta4m1, "z", "q");
   BMPY_R2(jacobi_theta1tau, "z", "tau");
   BMPY_R2(jacobi_theta2tau, "z", "tau");
   BMPY_R2(jacobi_theta3tau, "z", "tau");
   BMPY_R2(jacobi_theta4tau, "z", "tau");
   BMPY_R2(jacobi_theta3m1tau, "z", "tau");
   BMPY_R2(jacobi_theta4m1tau, "z", "tau");
}
