// Copyright Jacob Hass, 2026
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_QUADRATURE_SYMPLECTIC_HPP
#define BOOST_MATH_QUADRATURE_SYMPLECTIC_HPP

#include <vector>
#include <cmath>
#include <map>
#include <string>
#include <stdexcept>
#include <boost/math/policies/policy.hpp>
#include <boost/math/policies/error_handling.hpp>
#include <boost/math/special_functions/fpclassify.hpp>

namespace boost{ namespace math { namespace quadrature { namespace detail {

template <class RealType, class ReturnType, class Func>
std::pair<ReturnType, ReturnType> second_order_yoshida(const ReturnType p0, 
                                                       const ReturnType q0, 
                                                       RealType dt, 
                                                       Func dHdp, 
                                                       Func dHdq)
{
    ReturnType p = p0;
    ReturnType q = q0;

    // Half step in q
    auto dHdp_val = dHdp(p);
    for (unsigned i=0; i < q.size(); i++){
        q[i] = q[i] + dt / 2 * dHdp_val[i];
    }

    // Full step in p
    auto dHdq_val = dHdq(q);
    for (unsigned i=0; i < p.size(); i++){
        p[i] = p[i] - dt * dHdq_val[i];
    }

    // Half step in q
    dHdp_val = dHdp(p);
    for (unsigned i=0; i < q.size(); i++){
        q[i] = q[i] + dt / 2 * dHdp_val[i];
    }

    return std::make_pair(p, q);
}

template <class RealType, class ReturnType, class Func>
std::pair<ReturnType, ReturnType> fourth_order_yoshida(const ReturnType p0, 
                                                       const ReturnType q0, 
                                                       const RealType dt, 
                                                       Func dHdp, 
                                                       Func dHdq)
{
    BOOST_MATH_STD_USING
    
    ReturnType p = p0;
    ReturnType q = q0;

    // RealType x0 = -(std::pow(2, 1/3) / (2 - std::pow(2, 1/3)));
    RealType x1 = 1 / (2 - std::pow(2, 1/3));
    RealType x0 = 1 - 2 * x1; 

    std::vector<RealType> weights = { x1, x0, x1 };

    for (unsigned i=0; i < weights.size(); i++)
    {
        std::tie(p, q) = second_order_yoshida(p, q, weights[i] * dt, dHdp, dHdq);
    }
    
    return std::make_pair(p, q);
}

template <class RealType, class ReturnType, class Func>
std::pair<ReturnType, ReturnType> sixth_order_yoshida(const ReturnType p0, 
                                                      const ReturnType q0, 
                                                      RealType dt, 
                                                      Func dHdp, 
                                                      Func dHdq)
{
    ReturnType p = p0;
    ReturnType q = q0;
    
    // Choosing "System A" solution
    // The following Mathematica command can calculate these coefficients to arbitrary precision
    // FindRoot[{w0+2(w1+w2+w3)==1,
    //           w0^3 + 2(w1^3 + w2^3+w3^3)==0, 
    //           w0^5 + 2(w1^5 + w2^5+w3^5)==0, 
    //           1/6(w0^2w1^3-w0^4*w1) - 1/6(w0^3w1^2-w0*w1^4)+1/6((w0+2w1)^2w2^3-(w0+2w1)(w0^3+2w1^3)w2) - 1/6((w0^3+2w1^3)w2^2-(w0+2w1)w2^4) +1/6((w0+2w1+2w2)^2w3^3-(w0+2w1+2w2)(w0^3+2w1^3+2w2^3)w3)-1/6((w0^3+2w1^3+2w2^3)w3^2-(w0+2w1+2w2)w3^4)==0}, {{w0,1.3151863206839063}, {w1,-1.17767998417887}, {w2,0.235573213359357}, {w3,0.784513610477560}}, WorkingPrecision->64]
    RealType w1 = static_cast<RealType>(-1.17767998417887100694641568096431573463926925263459848447536851379674155618156L);
    RealType w2 = static_cast<RealType>(0.23557321335935813368479318297853460168646808210340111900349313095621471215223L);
    RealType w3 = static_cast<RealType>(0.78451361047755726381949763386634987577682441745149338456794779895125997479548L);
    // w0 = 1.31518632068391121888424972823886251435195350615940796180785516777853373846773
    RealType w0 = 1 - 2 * (w1 + w2 + w3);
    std::vector<RealType> weights = { w3, w2, w1, w0, w1, w2, w3};

    for (unsigned i=0; i < weights.size(); i++)
    {
        std::tie(p, q) = second_order_yoshida(p, q, weights[i] * dt, dHdp, dHdq);
    }
    
    return std::make_pair(p, q);
}

template <class ReturnType, class RealType, class Func, class Policy>
std::pair<std::vector<ReturnType>, std::vector<ReturnType> > integrate_hamiltonian_imp(const ReturnType p0,
                                                                                       const ReturnType q0,
                                                                                       const RealType dt,
                                                                                       const unsigned steps,
                                                                                       Func dHdp,
                                                                                       Func dHdq,
                                                                                       std::string method,
                                                                                       const Policy& pol)
{
    // Not sure how to make this function string nicer
    static const char* function = "boost::math::quadrature::integrate_hamiltonian(p0, q0, %1%, steps, dHdp, dHdq)";

    if ((dt <= 0) || !(boost::math::isfinite)(dt))
    {
        RealType val = (boost::math::policies::raise_domain_error(function, "Time step must be positive and finite but got: dt = %1%.\n", dt, pol));
        std::vector<ReturnType> nan_vec = {ReturnType(val)};
        return std::make_pair(nan_vec, nan_vec);
    }

    // Check if method is available
    std::vector<std::string> available_methods = {"Y6", "Y4", "Y2"};

    typedef std::pair<ReturnType, ReturnType> (*stepperType)(ReturnType, ReturnType, RealType, Func, Func);

    std::map<std::string, stepperType> m{{"Y6", sixth_order_yoshida}, 
                                         {"Y4", fourth_order_yoshida}, 
                                         {"Y2", second_order_yoshida}};
    stepperType stepper = m.at(method);

    std::vector<ReturnType> p(steps);
    std::vector<ReturnType> q(steps);
    p[0] = p0;
    q[0] = q0;

    ReturnType p_current = p0;
    ReturnType q_current = q0;
    for (unsigned i=1; i < steps; i++)
    {
        std::tie(p_current, q_current) = stepper(p_current, q_current, dt, dHdp, dHdq);
        p[i] = p_current;
        q[i] = q_current;
    }
    return std::make_pair(p, q);
}
} // namespace detail

template <class ReturnType, class RealType, class Func, class Policy>
std::pair<std::vector<ReturnType>, std::vector<ReturnType> > integrate_hamiltonian(const ReturnType p0,
                                                                                   const ReturnType q0,
                                                                                   const RealType dt,
                                                                                   const unsigned steps,
                                                                                   Func dHdp,
                                                                                   Func dHdq,
                                                                                   std::string method,
                                                                                   const Policy& pol)
{
    return detail::integrate_hamiltonian_imp(p0, q0, dt, steps, dHdp, dHdq, method, pol); 
}

template <class ReturnType, class RealType, class Func>
std::pair<std::vector<ReturnType>, std::vector<ReturnType> > integrate_hamiltonian(const ReturnType p0,
                                                                                   const ReturnType q0,
                                                                                   const RealType dt,
                                                                                   const unsigned steps,
                                                                                   Func dHdp,
                                                                                   Func dHdq,
                                                                                   std::string method)
{
    return integrate_hamiltonian(p0, q0, dt, steps, dHdp, dHdq, method, boost::math::policies::policy<>()); 
}

template <class ReturnType, class RealType, class Func>
std::pair<std::vector<ReturnType>, std::vector<ReturnType> > integrate_hamiltonian(const ReturnType p0,
                                                                                   const ReturnType q0,
                                                                                   const RealType dt,
                                                                                   const unsigned steps,
                                                                                   Func dHdp,
                                                                                   Func dHdq)
{
    return integrate_hamiltonian(p0, q0, dt, steps, dHdp, dHdq, "Y6", boost::math::policies::policy<>()); 
}
}}}
#endif
