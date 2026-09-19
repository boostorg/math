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
#include <boost/math/special_functions/cbrt.hpp>
#include <boost/math/interpolators/cubic_hermite.hpp>

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

template<typename T, typename = void>
struct has_size : std::false_type {};

template<typename T>
struct has_size<T, void_t<decltype(std::declval<T>().size())> > : std::true_type {};

template <typename T>
BOOST_MATH_INLINE_CONSTEXPR bool has_size_v = has_size<T>::value;

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

template <typename T, typename U>
typename std::enable_if_t<has_plus_v<T, U>, void>
subtract(T& x, const U& y)
{
    x = x - y;
}

template <typename T, typename U>
typename std::enable_if_t<!has_plus_v<T, U>, void>
subtract(T& vec1, const U& vec2)
{
    for (std::size_t i=0; i < vec1.size(); i++)
    {
        vec1[i] = vec1[i] - vec2[i];
    }
}

template <typename T>
typename std::enable_if_t<has_size_v<T>, typename T::value_type>
l2_norm(const T& vec)
{
    typename T::value_type sum = 0;
    for (std::size_t i=0; i < vec.size(); i++)
    {
        sum += vec[i] * vec[i];
    }
    return std::sqrt(sum);
}

template <typename T>
typename std::enable_if_t<!has_size_v<T>, T>
l2_norm(const T& val)
{
    BOOST_MATH_STD_USING;
    return abs(val);
}

template <typename T>
typename std::enable_if_t<has_size_v<T>, typename T::value_type>
avg_l2_norm(const T& vec)
{
    typename T::value_type sum = 0;
    for (std::size_t i=0; i < vec.size(); i++)
    {
        sum += vec[i] * vec[i];
    }
    typename T::value_type size = vec.size();
    return std::sqrt(sum / size);
}

template <typename T>
typename std::enable_if_t<!has_size_v<T>, T>
avg_l2_norm(const T& val)
{
    BOOST_MATH_STD_USING;
    return abs(val);
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
typename std::enable_if_t<!has_mult_v<T, U>, void>
mult_prefactor(T& vec1, U prefactor)
{
    for (std::size_t i=0; i < vec1.size(); i++)
    {
        vec1[i] = vec1[i] * prefactor;
    }
}

template <typename T, typename U, class Policy>
typename std::enable_if_t<(has_size_v<T> && has_size_v<U>), void>
size_check(const T& vec1, const U& vec2, const char* function, const Policy& pol)
{
    if (vec1.size() != vec2.size())
    {
        std::string msg = "Starting vectors (p0, q0) are not the same size but got (" + std::to_string(vec1.size()) + "!=" + std::to_string(vec2.size()) + ")";
        boost::math::policies::raise_evaluation_error(function, msg.c_str(), 0, pol);
    }
}

template <typename T, typename U, class Policy>
typename std::enable_if_t<!(has_size_v<T> && has_size_v<U>), void>
size_check(const T& vec1, const U& vec2, const char* function, const Policy& pol){    return;    }

// Function to initialize array of integer types of size N
template<typename T>
typename std::enable_if_t<!has_size_v<T>, std::vector<T> >
initialize_array(const T& p0, const size_t N)
{
    std::vector<T> initializedArr(N);
    return initializedArr;
}

// Function to initialize vector of vectors of size (N, p0.size())
template<typename T>
typename std::enable_if_t<has_size_v<T>, std::vector<T> >
initialize_array(const T& p0, const size_t N)
{
    std::vector<T> initializedArr(N, T(p0.size(), 0));
    return initializedArr;
}

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
void SRKN_b_6(RandomAccessContainer& p0, RandomAccessContainer& q0, const RealType dt,
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
void SRKN_b_11(RandomAccessContainer& p0, RandomAccessContainer& q0, const RealType dt,
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
        case available_methods::SRKNB6:   stepper = SRKN_b_6; break;
        case available_methods::SRKNB11:  stepper = SRKN_b_11; break;
        default: boost::math::policies::raise_domain_error(function, "Incorrect method recieved. Must be in `available_methods` enum class.", 0, pol);
    }

    std::vector<RandomAccessContainer> p = initialize_array(p0, steps);
    std::vector<RandomAccessContainer> q = initialize_array(q0, steps);
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

template <typename RandomAccessContainer, typename RealType, class stepperType, class Func>
struct error_info
{
    error_info(stepperType stepper_, Func& dHdp_, Func& dHdq_, unsigned order_, RealType atol_, RealType rtol_)
    : stepper(stepper_), dHdp(dHdp_), dHdq(dHdq_), order(order_), atol(atol_), rtol(rtol_)
    {
        safe = RealType(0.9);
        maxScale = RealType(2);
        minScale = RealType(0.2);
        rejected = true;
    }

    // Here, current_p/q represents (p, q) before a single step of size dt is made
    // whereas next_p/q represents (p, q) after a step of size dt
    bool success(RandomAccessContainer& current_p,
                 RandomAccessContainer& current_q,
                 RandomAccessContainer& next_p,
                 RandomAccessContainer& next_q,
                 RealType& dt,
                 RealType& currentError)
    {
        RandomAccessContainer error_p = current_p;
        RandomAccessContainer error_q = current_q;

        // For the error, reference method is using two half steps
        stepper(error_p, error_q, dt / 2, dHdp, dHdq);
        stepper(error_p, error_q, dt / 2, dHdp, dHdq);

        // Get difference in solution between one/two timesteps
        subtract(error_p, next_p);
        RealType error = avg_l2_norm(error_p);

        // Get scale in which to reject error
        RealType currentDistance = l2_norm(current_p);
        RealType newDistance = l2_norm(next_p);
        RealType tol = atol + rtol * std::max(currentDistance, newDistance);
        error /= tol;
        currentError = error;

        RealType scale;
        if (error <= 1.0)
        {
            if (error == 0.0)
            {
                scale = maxScale;
            }
            else
            {
                scale = safe * std::pow(static_cast<RealType>(1 / error), RealType(1) / (order + 1));
                if (scale < minScale) scale = minScale;
                if (scale > maxScale) scale = maxScale;
            }

            if (rejected)
            {
                dt *= std::min(scale, 1.0);
            }
            else
            {
                dt *= scale;
            }
            errorOld = std::max(error, 1e-4);
            rejected = false;
            return true;
        }
        else
        {
            scale = std::max(safe * std::pow(static_cast<RealType>(1 / error), RealType(1) / (order + 1)), minScale);
            dt *= scale;
            rejected = true;
            return false;
        }
    }

    RealType errorOld;
    RealType safe;
    RealType maxScale;
    RealType minScale;
    bool rejected;
    RealType atol;
    RealType rtol;
    unsigned order;
    stepperType stepper;
    Func& dHdp;
    Func& dHdq;
};

template <typename RandomAccessContainer, typename RealType, class Func, class Policy>
std::tuple<std::vector<RealType>, std::vector<RealType>, std::vector<RealType> > integrate_hamiltonian_adaptive(RandomAccessContainer& p0,
                                                                                                               RandomAccessContainer& q0,
                                                                                                               const std::pair<RealType, RealType>& timeInterval,
                                                                                                               const Func& dHdp,
                                                                                                               const Func& dHdq,
                                                                                                               const RealType& atol,
                                                                                                               const RealType& rtol,
                                                                                                               const available_methods& method,
                                                                                                               const Policy& pol)
{
    BOOST_MATH_STD_USING
    // Not sure how to make this function string nicer
    static const char* function = "boost::math::quadrature::integrate_hamiltonian(p0, q0, %1%, steps, dHdp, dHdq)";

    if (!(boost::math::isfinite)(timeInterval.second))
    {
        boost::math::policies::raise_domain_error(function, "Maximum time  must be positive and finite but got: tMax = %1%.\n", timeInterval.second, pol);
    }

    if ((timeInterval.first <= 0))
    {
        boost::math::policies::raise_domain_error(function, "Minimum time  must be positive and finite but got: tMin = %1%.\n", timeInterval.first, pol);
    }

    // Check that p0 and q0 have the same size
    size_check(p0, q0, function, pol);
    #ifdef _MSC_VER
        #  pragma warning(pop)
    #endif

    typedef void (*stepperType)(RandomAccessContainer&, RandomAccessContainer&, RealType, Func, Func);

    unsigned order;
    stepperType stepper;
    switch (method) {
        case available_methods::Y2:       stepper = second_order_yoshida; order = 2; break;
        case available_methods::Y6:       stepper = sixth_order_yoshida; order = 6; break;
        case available_methods::Y4:       stepper = fourth_order_yoshida; order = 4; break;
        case available_methods::SRKNB6:   stepper = SRKN_b_6; order = 4; break;
        case available_methods::SRKNB11:  stepper = SRKN_b_11; order = 6; break;
        default: boost::math::policies::raise_domain_error(function, "Incorrect method recieved. Must be in `available_methods` enum class.", 0, pol);
    }

    std::vector<RandomAccessContainer> p = { p0 };
    std::vector<RandomAccessContainer> q = { q0 };

    std::vector<RealType> time = { timeInterval.first };
    RealType dt = 0.01;

    // To calculate the error, we need to store the current values of p and q seperately
    // from the values updated in the loop. This is because the steppers modify p and q in place
    RandomAccessContainer current_p;
    RandomAccessContainer current_q;

    error_info<RandomAccessContainer, RealType, stepperType, Func> info(
        stepper, dHdp, dHdq, order, atol, rtol);

    RealType currentError = 0;
    std::vector<RealType> error = { currentError };

    while (time.back() < timeInterval.second)
    {
        // Save p0 and q0 for error estimation step
        current_p = p0;
        current_q = q0;

        // Step p0 and q0 forward one step
        stepper(p0, q0, dt, dHdp, dHdq);

        bool stepSuccessful = info.success(current_p, current_q, p0, q0, dt, currentError);
        if (!stepSuccessful)
        {
            // Reject the step and reset p0/q0 to current_p/q
            p0 = current_p;
            q0 = current_q;
        }
        else
        {
            // Accept the step, push back p/q and step t forward
            p.push_back(p0);
            q.push_back(q0);
            time.push_back(time.back() + dt);
            error.push_back(currentError);
        }
    }
    return std::make_tuple(time, p, error);
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

template <typename RealType, class Func, class Policy>
BOOST_MATH_EXPORT boost::math::interpolators::cubic_hermite<std::vector<RealType> > integrate_hamiltonian_cubic(RealType p0,
                                                                                                                RealType q0,
                                                                                                                const RealType dt,
                                                                                                                const unsigned steps,
                                                                                                                Func dHdp,
                                                                                                                Func dHdq,
                                                                                                                available_methods method,
                                                                                                                const Policy& pol)
{
    std::vector<RealType> p;
    std::vector<RealType> q;
    std::tie(p, q) = detail::integrate_hamiltonian_imp(p0, q0, dt, steps, dHdp, dHdq, method, pol);

    std::vector<RealType> dpdt(p.size());
    std::vector<RealType> time(p.size(), 0);
    for (std::size_t i=0; i < p.size(); i++)
    {
        dpdt[i] = -dHdq(p[i]);
        time[i] = i * dt;
    }
    return boost::math::interpolators::cubic_hermite<std::vector<RealType> >(std::move(time), std::move(p), std::move(dpdt));
}

template <typename RealType, class Func>
BOOST_MATH_EXPORT boost::math::interpolators::cubic_hermite<std::vector<RealType> > integrate_hamiltonian_cubic(RealType p0,
                                                                                                                RealType q0,
                                                                                                                const RealType dt,
                                                                                                                const unsigned steps,
                                                                                                                Func dHdp,
                                                                                                                Func dHdq,
                                                                                                                available_methods method)
{
    return integrate_hamiltonian_cubic(p0, q0, dt, steps, dHdp, dHdq, method, boost::math::policies::policy<>());
}

template <typename RealType, class Func>
BOOST_MATH_EXPORT boost::math::interpolators::cubic_hermite<std::vector<RealType> > integrate_hamiltonian_cubic(RealType p0,
                                                                                                                RealType q0,
                                                                                                                const RealType dt,
                                                                                                                const unsigned steps,
                                                                                                                Func dHdp,
                                                                                                                Func dHdq)
{
    return integrate_hamiltonian_cubic(p0, q0, dt, steps, dHdp, dHdq, available_methods::Y6, boost::math::policies::policy<>());
}

}}}

#endif
