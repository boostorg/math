// Copyright Jacob Hass, 2026
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

// Disable warnings when defining size_check for numeric types
#ifdef _MSC_VER
#  pragma warning(push)
#  pragma warning (disable : 4100) // 'function': unreferenced parameter
#endif

#ifndef BOOST_MATH_QUADRATURE_SYMPLECTIC_HPP
#define BOOST_MATH_QUADRATURE_SYMPLECTIC_HPP

#include <utility>
#include <vector>
#include <cmath>
#include <map>
#include <string>
#include <stdexcept>
#include <type_traits>
#include <boost/math/policies/policy.hpp>
#include <boost/math/policies/error_handling.hpp>
#include <boost/math/special_functions/fpclassify.hpp>

namespace boost{ namespace math { namespace quadrature { 

BOOST_MATH_EXPORT enum class available_methods 
{
    Y6,
    Y4,
    Y2,
    SRKNB6,
    SRKNB11
};
    
namespace detail {

template<typename...>
using void_t = void;

template<typename T, typename U, typename = void>
struct has_plus : std::false_type {};

template<typename T, typename U>
struct has_plus<T, U, void_t<decltype(std::declval<T>() + std::declval<U>())> > : std::true_type {};

template <typename T, typename U>
BOOST_MATH_INLINE_CONSTEXPR bool has_plus_v = has_plus<T, U>::value;

template <typename T, typename U>
typename std::enable_if_t<has_plus_v<T, U>, void>
add(T& x, const U& y) 
{
    x = x + y;
}

template <typename T, typename U>
typename std::enable_if_t<!has_plus_v<T, U>, void>
add(T& vec1, const U& vec2)
{
    for (std::size_t i=0; i < vec1.size(); i++)
    {
        vec1[i] = vec1[i] + vec2[i];
    }
}

template<typename T, typename U, typename = void>
struct has_mult : std::false_type {};

template<typename T, typename U>
struct has_mult<T, U, void_t<decltype(std::declval<T>() * std::declval<U>())> > : std::true_type {};

template <typename T, typename U>
BOOST_MATH_INLINE_CONSTEXPR bool has_mult_v = has_mult<T, U>::value;

template <typename T, typename U>
typename std::enable_if_t<has_mult_v<T, U>, void>
mult_prefactor(T& x, U prefactor) 
{
    x = x * prefactor;
}

template <typename T, typename U>
typename std::enable_if_t<!has_plus_v<T, U>, void>
mult_prefactor(T& vec1, U prefactor) 
{
    for (std::size_t i=0; i < vec1.size(); i++)
    {
        vec1[i] = vec1[i] * prefactor;
    }
}

template <typename T, typename U, class Policy>
typename std::enable_if_t<!has_plus_v<T, U>, void>
size_check(const T& vec1, const U& vec2, const char* function, const Policy& pol)
{
    if (vec1.size() != vec2.size())
    {
        std::string msg = "Starting vectors (p0, q0) are not the same size but got (" + std::to_string(vec1.size()) + "!=" + std::to_string(vec2.size()) + ")";
        boost::math::policies::raise_evaluation_error(function, msg.c_str(), 0, pol);
    }
}

template <typename T, typename U, class Policy>
typename std::enable_if_t<has_plus_v<T, U>, void>
size_check(const T& vec1, const U& vec2, const char* function, const Policy& pol){    return;    }

template <typename RandomAccessContainer, typename RealType, class Func>
void second_order_yoshida(RandomAccessContainer& p0, RandomAccessContainer& q0, const RealType dt, 
                          const Func dHdp, const Func dHdq)
{

    // Half step in q
    RandomAccessContainer dq = dHdp(p0);
    mult_prefactor(dq, dt / RealType(2));
    add(q0, dq);

    // Full step in p
    RandomAccessContainer dp = dHdq(q0);
    mult_prefactor(dp, -dt);
    add(p0, dp);

    // Half step in q
    dq = dHdp(p0);
    mult_prefactor(dq, dt / RealType(2));
    add(q0, dq);
}

template <typename RandomAccessContainer, typename RealType, class Func>
void fourth_order_yoshida(RandomAccessContainer& p0, RandomAccessContainer& q0, const RealType dt, 
                          const Func dHdp, const Func dHdq)
{
    BOOST_MATH_STD_USING

    RealType x1 = RealType(1) / (RealType(2) - cbrt(RealType(2)));
    RealType x0 = RealType(1) - RealType(2) * x1; 

    std::array<RealType, 3> weights = { x1, x0, x1 };

    for (std::size_t i=0; i < weights.size(); i++)
    {
        second_order_yoshida(p0, q0, weights[i] * dt, dHdp, dHdq);
    }
}

template <typename RandomAccessContainer, typename RealType, class Func>
void sixth_order_yoshida(RandomAccessContainer& p0, RandomAccessContainer& q0, const RealType dt, 
                         const Func dHdp, const Func dHdq)
{
    
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
    RealType w0 = RealType(1) - RealType(2) * (w1 + w2 + w3);
    std::array<RealType, 7> weights = { w3, w2, w1, w0, w1, w2, w3};

    for (std::size_t i=0; i < weights.size(); i++)
    {
        second_order_yoshida(p0, q0, weights[i] * dt, dHdp, dHdq);
    }
}

template <typename RandomAccessContainer, typename RealType, class Func>
void SRKN_b_order_6(RandomAccessContainer& p0, RandomAccessContainer& q0, const RealType dt, 
                    const Func dHdp, const Func dHdq)
{ // This method implements SRKN_b^6 in Table 3 here 
  // https://www.sciencedirect.com/science/article/pii/S0377042701004927
    
    RealType b1 = static_cast<RealType>(0.0829844064174052);
    RealType b2 = static_cast<RealType>(0.396309801498368);
    RealType b3 = static_cast<RealType>(-0.0390563049223486);
    RealType b4 = 1.0 - 2.0 * (b1 + b2 + b3);

    RealType a1 = static_cast<RealType>(0.245298957184271);
    RealType a2 = static_cast<RealType>(0.604872665711080);
    RealType a3 = 0.5 - (a1 + a2);

    std::array<RealType, 6> b_weights = {b1, b2, b3, b4, b3, b2};
    std::array<RealType, 6> a_weights = {a1, a2, a3, a3, a2, a1};
    
    RealType a, b;
    for (std::size_t i=0; i < b_weights.size(); i++)
    {
        b = b_weights[i];
        a = a_weights[i];
        
        RandomAccessContainer dq = dHdp(p0);
        mult_prefactor(dq, dt * b);
        add(q0, dq);

        RandomAccessContainer dp = dHdq(q0);
        mult_prefactor(dp, -a * dt);
        add(p0, dp);
    }
    // Need to do one more step in q
    RandomAccessContainer dq = dHdp(p0);
    mult_prefactor(dq, dt * b1);
    add(q0, dq);
}

template <typename RandomAccessContainer, typename RealType, class Func>
void SRKN_b_order_11(RandomAccessContainer& p0, RandomAccessContainer& q0, const RealType dt, 
                     const Func dHdp, const Func dHdq)
{ // This method implements SRKN_b^11 in Table 3 here 
  // https://www.sciencedirect.com/science/article/pii/S0377042701004927
    
    RealType b1 = static_cast<RealType>(0.0414649985182624);
    RealType b2 = static_cast<RealType>(0.198128671918067);
    RealType b3 = static_cast<RealType>(-0.0400061921041533);
    RealType b4 = static_cast<RealType>(0.0752539843015807);
    RealType b5 = static_cast<RealType>(-0.0115113874206879);
    RealType b6 = 0.5 - (b1 + b2 + b3 + b4 + b5);

    RealType a1 = static_cast<RealType>(0.123229775946271);
    RealType a2 = static_cast<RealType>(0.290553797799558);
    RealType a3 = static_cast<RealType>(-0.127049212625417);
    RealType a4 = static_cast<RealType>(-0.246331761062075);
    RealType a5 = static_cast<RealType>(0.357208872795928);
    RealType a6 = 1.0 - 2.0 * (a1 + a2 + a3 + a4 + a5);

    std::array<RealType, 11> b_weights = {b1, b2, b3, b4, b5, b6, b6, b5, b4, b3, b2};
    std::array<RealType, 11> a_weights = {a1, a2, a3, a4, a5, a6, a5, a4, a3, a2, a1};
    
    RealType a, b;
    for (std::size_t i=0; i < b_weights.size(); i++)
    {
        b = b_weights[i];
        a = a_weights[i];
        
        RandomAccessContainer dq = dHdp(p0);
        mult_prefactor(dq, dt * b);
        add(q0, dq);

        RandomAccessContainer dp = dHdq(q0);
        mult_prefactor(dp, -a * dt);
        add(p0, dp);
    }
    // Need to do one more step in q
    RandomAccessContainer dq = dHdp(p0);
    mult_prefactor(dq, dt * b1);
    add(q0, dq);
}

template <typename RandomAccessContainer, typename RealType, class Func, class Policy>
std::pair<std::vector<RandomAccessContainer>, std::vector<RandomAccessContainer> > integrate_hamiltonian_imp(RandomAccessContainer& p0,
                                                                                                             RandomAccessContainer& q0,
                                                                                                             const RealType dt,
                                                                                                             const unsigned steps,
                                                                                                             const Func& dHdp,
                                                                                                             const Func& dHdq,
                                                                                                             const available_methods& method,
                                                                                                             const Policy& pol)
{
    // Not sure how to make this function string nicer
    static const char* function = "boost::math::quadrature::integrate_hamiltonian(p0, q0, %1%, steps, dHdp, dHdq)";

    if ((dt <= 0) || !(boost::math::isfinite)(dt))
    {
        boost::math::policies::raise_domain_error(function, "Time step must be positive and finite but got: dt = %1%.\n", dt, pol);
    }

    // Check that p0 and q0 have the same size
    size_check(p0, q0, function, pol);
    #ifdef _MSC_VER
        #  pragma warning(pop)
    #endif 

    typedef void (*stepperType)(RandomAccessContainer&, RandomAccessContainer&, RealType, Func, Func);

    stepperType stepper; 
    switch (method) {
        case available_methods::Y6:       stepper = sixth_order_yoshida; break;
        case available_methods::Y4:       stepper = fourth_order_yoshida; break;
        case available_methods::Y2:       stepper = second_order_yoshida; break;
        case available_methods::SRKNB6:   stepper = SRKN_b_order_6; break;
        case available_methods::SRKNB11:  stepper = SRKN_b_order_11; break;
        default: boost::math::policies::raise_domain_error(function, "Incorrect method recieved. Must be in `available_methods` enum class.", 0, pol);
    }

    std::vector<RandomAccessContainer> p(steps);
    std::vector<RandomAccessContainer> q(steps);
    p[0] = p0;
    q[0] = q0;

    for (unsigned i=1; i < steps; i++)
    {
        stepper(p0, q0, dt, dHdp, dHdq);
        p[i] = p0;
        q[i] = q0;
    }
    return std::make_pair(p, q);
}
} // namespace detail

template <typename RandomAccessContainer, typename RealType, class Func, class Policy>
BOOST_MATH_EXPORT std::pair<std::vector<RandomAccessContainer>, std::vector<RandomAccessContainer> > integrate_hamiltonian(RandomAccessContainer p0,
                                                                                                                           RandomAccessContainer q0,
                                                                                                                           const RealType dt,
                                                                                                                           const unsigned steps,
                                                                                                                           Func dHdp,
                                                                                                                           Func dHdq,
                                                                                                                           available_methods method,
                                                                                                                           const Policy& pol)
{
    return detail::integrate_hamiltonian_imp(p0, q0, dt, steps, dHdp, dHdq, method, pol); 
}

template <typename RandomAccessContainer, typename RealType, class Func>
BOOST_MATH_EXPORT std::pair<std::vector<RandomAccessContainer>, std::vector<RandomAccessContainer> > integrate_hamiltonian(RandomAccessContainer p0,
                                                                                                                           RandomAccessContainer q0,
                                                                                                                           const RealType dt,
                                                                                                                           const unsigned steps,
                                                                                                                           Func dHdp,
                                                                                                                           Func dHdq,
                                                                                                                           available_methods method)
{
    return integrate_hamiltonian(p0, q0, dt, steps, dHdp, dHdq, method, boost::math::policies::policy<>()); 
}

template <typename RandomAccessContainer, typename RealType, class Func>
BOOST_MATH_EXPORT std::pair<std::vector<RandomAccessContainer>, std::vector<RandomAccessContainer> > integrate_hamiltonian(RandomAccessContainer p0,
                                                                                                                           RandomAccessContainer q0,
                                                                                                                           const RealType dt,
                                                                                                                           const unsigned steps,
                                                                                                                           Func dHdp,
                                                                                                                           Func dHdq)
{
    return integrate_hamiltonian(p0, q0, dt, steps, dHdp, dHdq, available_methods::Y6, boost::math::policies::policy<>()); 
}
}}}

#endif
