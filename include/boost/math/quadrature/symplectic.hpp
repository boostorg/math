/* Based on "Construction of higher order symplectic integrators" by Haruo Yoshida*/

#ifndef BOOST_MATH_QUADRATURE_SYMPLECTIC_HPP
#define BOOST_MATH_QUADRATURE_SYMPLECTIC_HPP

#include <vector>
#include <cmath>
#include <string>
#include <stdexcept>
#include <boost/math/policies/policy.hpp>
#include <boost/math/policies/error_handling.hpp>
#include <boost/math/special_functions/fpclassify.hpp>

namespace boost{ namespace math {namespace quadrature { 

template <class RealType, class ReturnType, class Func>
std::pair<ReturnType, ReturnType> second_order_yoshida(const ReturnType p0, 
                                                       const ReturnType q0, 
                                                       RealType dt, 
                                                       Func dHdp, 
                                                       Func dHdq)
{
    ReturnType p = p0;
    ReturnType q = q0;

    q = q + dt / 2 * dHdp(p);
    p = p - dt * dHdq(q);
    q = q + dt / 2 * dHdp(p);
    
    return std::make_pair(p, q);
}

template <class RealType, class ReturnType, class Func>
std::pair<ReturnType, ReturnType> fourth_order_yoshida(const ReturnType p0, 
                                                       const ReturnType q0, 
                                                       const RealType dt, 
                                                       Func dHdp, 
                                                       Func dHdq)
{
    ReturnType p = p0;
    ReturnType q = q0;

    RealType x0 = -(std::pow(2, 1/3) / (2 - std::pow(2, 1/3)));
    RealType x1 = 1 / (2 - std::pow(2, 1/3));

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
    RealType w3 = 0.784513610477560;
    RealType w2 = 0.235573213359357;
    RealType w1 = -1.17767998417887;
    RealType w0 = 1 - 2 * (w1 + w2 + w3);
    std::vector<RealType> weights = { w3, w2, w1, w0, w1, w2, w3};

    for (unsigned i=0; i < weights.size(); i++)
    {
        std::tie(p, q) = second_order_yoshida(p, q, weights[i] * dt, dHdp, dHdq);
    }
    
    return std::make_pair(p, q);
}

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

    std::pair<ReturnType, ReturnType> (*stepper)(ReturnType, ReturnType, RealType, Func, Func);

    if (method == "Y6"){
        stepper = sixth_order_yoshida;
    } 
    else if (method == "Y4"){
        stepper = fourth_order_yoshida;
    }
    else if (method == "Y2")
    {
        stepper = second_order_yoshida;
    }
    else
    {
        throw std::domain_error("Method must be in {'Y6', 'Y4', 'Y2'} but got: method = " + method);
    }

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
