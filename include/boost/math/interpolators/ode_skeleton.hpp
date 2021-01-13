// Copyright Nick Thompson, 2021
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

// This is a sensible default interpolator for ODE skeletons.
// An ODE skeleton is a list of samples {t_i, v_, dvdt_i} produced by an ODE stepper,
// such as RK45. There are many very sophisticated ODE skeleton interpolators,
// but many require the RHS of v' = f to be available during interpolation,
// or require the stages in a given step (which are discarded).
// These interpolators are generally used for *event* localization, for example,
// when a ball traveling in a fluid strikes the ground.
// The stepper has to be paranoid and try to identify the event using it's interpolator as the solution proceeds.

// This particular interpolator is designed for a different use case I personally encountered:
// Where v' = f is solved on an HPC system, and f is unavailable during visualization and analysis.
// I hope it is useful for more scenarios, but obviously if you have f available or have knowledge of how
// the solution is produced you can exploit it.

// One might wonder what this interpolator actually is.
// It is PQHIP, which is the quintic version of PCHIP.
// It is discussed in Corless, A Graduate Introduction to Numerical Methods, Exercise 8.40.
// In addition, an implementation in matlab is found here:
// http://www.nfillion.com/coderepository/Graduate_Introduction_to_Numerical_Methods/pqhip.m

// This ODE interpolator is thrice continuously differentiable, as is required to analyze the differential geometry of the curve,
// e.g., extracting the curvature and the torsion of the curve, as well as the Frenet frame.
// The computation of the curvature, torsion, and Frenet frame is not difficult once three derivatives are provided,
// but nonetheless to avoid tedium I have exposed those quantities.

#ifndef BOOST_MATH_INTERPOLATORS_ODE_SKELETON_INTERPOLATOR_HPP
#define BOOST_MATH_INTERPOLATORS_ODE_SKELETON_INTERPOLATOR_HPP
#include <memory>
#include <cmath>
#include <boost/math/interpolators/detail/ode_skeleton_detail.hpp>

namespace boost::math::interpolators {

template<class TimeContainer, class SpaceContainer, class TangentContainer>
class ode_skeleton_interpolator
{
public:
    using Real = typename TimeContainer::value_type;
    using SpacePoint = typename SpaceContainer::value_type;
    using TangentPoint = typename TangentContainer::value_type;
    ode_skeleton_interpolator(TimeContainer&& times, SpaceContainer&& v, TangentContainer&& dvdt);

    void operator()(Point& x, Real t) const;

    SpacePoint operator()(Real t) const {
        Point p;
        this->operator()(p, t);
        return p;
    }

    void prime(TangentPoint& dvdt, Real t) const {
        Point v;
        m_imp->eval_with_prime(v, dvdt, t);
    }

    TanentPoint prime(Real t) const {
        TangentPoint p;
        this->prime(p, t);
        return p;
    }

    void eval_with_prime(Point& x, Point& dxdt, Real t) const {
        m_imp->eval_with_prime(x, dxdt, t);
        return;
    }

    std::pair<SpacePoint, TangentPoint> eval_with_prime(Real t) const {
        SpacePoint x;
        TangentPoint dxdt;
        m_imp->eval_with_prime(x, dxdt, t);
        return {x, dxdt};
    }

private:
    std::shared_ptr<detail::ode_skeleton_imp<TimeContainer, SpaceContainer, TangentContainer>> m_imp;
};


}
#endif
