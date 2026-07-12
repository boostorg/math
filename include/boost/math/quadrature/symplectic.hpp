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

template<typename...>
using void_t = void;

template<typename T, typename U, typename = void>
struct has_plus : std::false_type {};

template<typename T, typename U>
struct has_plus<T, U, std::void_t<decltype(std::declval<T>() + std::declval<U>())> > : std::true_type {};

template <typename T, typename U>
typename std::enable_if<has_plus<T, U>::value, T>::type
add(T x, U y) 
{
    return x + y;
}

template <typename T, typename U>
typename std::enable_if<!has_plus<T, U>::value, T>::type
add(T vec1, U vec2) 
{
    for (unsigned i=0; i < vec1.size(); i++)
    {
        vec1[i] = vec1[i] + vec2[i];
    }
    return vec1;
}

template<typename...>
using void_t = void;

template<typename T, typename U, typename = void>
struct has_mult : std::false_type {};

template<typename T, typename U>
struct has_mult<T, U, std::void_t<decltype(std::declval<T>() * std::declval<U>())> > : std::true_type {};

template <typename T, typename U>
typename std::enable_if<has_mult<T, U>::value, T>::type
mult_prefactor(T x, U prefactor) 
{
    return x * prefactor;
}

template <typename T, typename U>
typename std::enable_if<!has_plus<T, U>::value, T>::type
mult_prefactor(T vec1, U prefactor) 
{
    for (unsigned i=0; i < vec1.size(); i++)
    {
        vec1[i] = vec1[i] * prefactor;
    }
    return vec1;
}

template <typename RandomAccessContainer, typename RealType, class Func>
std::pair<RandomAccessContainer, RandomAccessContainer> second_order_yoshida(const RandomAccessContainer p0, 
                                                                             const RandomAccessContainer q0, 
                                                                             const RealType dt, 
                                                                             Func dHdp, 
                                                                             Func dHdq)
{
    RandomAccessContainer p = p0;
    RandomAccessContainer q = q0;

    // Half step in q
    RandomAccessContainer dHdp_val = dHdp(p);
    RandomAccessContainer dq = mult_prefactor(dHdp_val, dt / 2);
    q = add(q, dq);

    // Full step in p
    RandomAccessContainer dHdq_val = dHdq(q);
    RandomAccessContainer dp = mult_prefactor(dHdq_val, -dt);
    p = add(p, dp);

    // Half step in q
    dHdp_val = dHdp(p);
    dq = mult_prefactor(dHdp_val, dt / 2);
    q = add(q, dq);

    return std::make_pair(p, q);
}

template <typename RandomAccessContainer, typename RealType, class Func>
std::pair<RandomAccessContainer, RandomAccessContainer> fourth_order_yoshida(const RandomAccessContainer p0, 
                                                                             const RandomAccessContainer q0, 
                                                                             const RealType dt, 
                                                                             Func dHdp, 
                                                                             Func dHdq)
{
    BOOST_MATH_STD_USING

    RandomAccessContainer p = p0;
    RandomAccessContainer q = q0;

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

template <typename RandomAccessContainer, typename RealType, class Func>
std::pair<RandomAccessContainer, RandomAccessContainer> sixth_order_yoshida(const RandomAccessContainer p0, 
                                                                            const RandomAccessContainer q0, 
                                                                            const RealType dt, 
                                                                            Func dHdp, 
                                                                            Func dHdq)
{
    RandomAccessContainer p = p0;
    RandomAccessContainer q = q0;
    
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
    RealType w0 = 1. - 2. * (w1 + w2 + w3);
    std::vector<RealType> weights = { w3, w2, w1, w0, w1, w2, w3};

    for (unsigned i=0; i < weights.size(); i++)
    {
        std::tie(p, q) = second_order_yoshida(p, q, weights[i] * dt, dHdp, dHdq);
    }
    
    return std::make_pair(p, q);
}

template <typename RandomAccessContainer, typename RealType, class Func>
std::pair<RandomAccessContainer, RandomAccessContainer> SRKN_b_order_6(const RandomAccessContainer p0,
                                                                       const RandomAccessContainer q0,
                                                                       const RealType dt, 
                                                                       Func dHdp,
                                                                       Func dHdq)
{ // This method implements SRKN_b^6 in Table 3 here 
  // https://www.sciencedirect.com/science/article/pii/S0377042701004927

    RandomAccessContainer p = p0;
    RandomAccessContainer q = q0;
    
    RealType b1 = static_cast<RealType>(0.0829844064174052);
    RealType b2 = static_cast<RealType>(0.396309801498368);
    RealType b3 = static_cast<RealType>(-0.0390563049223486);
    RealType b4 = 1. - 2. * (b1 + b2 + b3);

    RealType a1 = static_cast<RealType>(0.245298957184271);
    RealType a2 = static_cast<RealType>(0.604872665711080);
    RealType a3 = 0.5 - (a1 + a2);

    std::vector<RealType> b_weights = {b1, b2, b3, b4, b3, b2};
    std::vector<RealType> a_weights = {a1, a2, a3, a3, a2, a1};
    
    RealType a, b;
    for (unsigned int i=0; i < b_weights.size(); i++)
    {
        b = b_weights[i];
        a = a_weights[i];
        
        RandomAccessContainer dHdp_val = dHdp(p);
        RandomAccessContainer dq = mult_prefactor(dHdp_val, dt * b);
        q = add(q, dq);

        RandomAccessContainer dHdq_val = dHdq(q);
        RandomAccessContainer dp = mult_prefactor(dHdq_val, -a * dt);
        p = add(p, dp);
    }
    // Need to do one more step in q
    RandomAccessContainer dHdp_val = dHdp(p);
    RandomAccessContainer dq = mult_prefactor(dHdp_val, dt * b1);
    q = add(q, dq);

    return std::make_pair(p, q);
}

template <typename RandomAccessContainer, typename RealType, class Func, class Policy>
std::pair<std::vector<RandomAccessContainer>, std::vector<RandomAccessContainer> > integrate_hamiltonian_imp(const RandomAccessContainer p0,
                                                                                                             const RandomAccessContainer q0,
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
        std::vector<RandomAccessContainer> nan_vec = {RandomAccessContainer(val)};
        return std::make_pair(nan_vec, nan_vec);
    }

    // Check if method is available
    std::vector<std::string> available_methods = {"Y6", "Y4", "Y2"};

    typedef std::pair<RandomAccessContainer, RandomAccessContainer> (*stepperType)(RandomAccessContainer, RandomAccessContainer, RealType, Func, Func);

    std::map<std::string, stepperType> m{{"Y6", sixth_order_yoshida}, 
                                         {"Y4", fourth_order_yoshida}, 
                                         {"Y2", second_order_yoshida},
                                         {"SRKNB6", SRKN_b_order_6}};
    stepperType stepper = m.at(method);

    std::vector<RandomAccessContainer> p(steps);
    std::vector<RandomAccessContainer> q(steps);
    p[0] = p0;
    q[0] = q0;

    RandomAccessContainer p_current = p0;
    RandomAccessContainer q_current = q0;
    for (unsigned i=1; i < steps; i++)
    {
        std::tie(p_current, q_current) = stepper(p_current, q_current, dt, dHdp, dHdq);
        p[i] = p_current;
        q[i] = q_current;
    }
    return std::make_pair(p, q);
}
} // namespace detail

template <typename RandomAccessContainer, typename RealType, class Func, class Policy>
std::pair<std::vector<RandomAccessContainer>, std::vector<RandomAccessContainer> > integrate_hamiltonian(const RandomAccessContainer p0,
                                                                                                         const RandomAccessContainer q0,
                                                                                                         const RealType dt,
                                                                                                         const unsigned steps,
                                                                                                         Func dHdp,
                                                                                                         Func dHdq,
                                                                                                         std::string method,
                                                                                                         const Policy& pol)
{
    return detail::integrate_hamiltonian_imp(p0, q0, dt, steps, dHdp, dHdq, method, pol); 
}

template <typename RandomAccessContainer, typename RealType, class Func>
std::pair<std::vector<RandomAccessContainer>, std::vector<RandomAccessContainer> > integrate_hamiltonian(const RandomAccessContainer p0,
                                                                                                         const RandomAccessContainer q0,
                                                                                                         const RealType dt,
                                                                                                         const unsigned steps,
                                                                                                         Func dHdp,
                                                                                                         Func dHdq,
                                                                                                         std::string method)
{
    return integrate_hamiltonian(p0, q0, dt, steps, dHdp, dHdq, method, boost::math::policies::policy<>()); 
}

template <typename RandomAccessContainer, typename RealType, class Func>
std::pair<std::vector<RandomAccessContainer>, std::vector<RandomAccessContainer> > integrate_hamiltonian(const RandomAccessContainer p0,
                                                                                                         const RandomAccessContainer q0,
                                                                                                         const RealType dt,
                                                                                                         const unsigned steps,
                                                                                                         Func dHdp,
                                                                                                         Func dHdq)
{
    return integrate_hamiltonian(p0, q0, dt, steps, dHdp, dHdq, "Y6", boost::math::policies::policy<>()); 
}
}}}
#endif
