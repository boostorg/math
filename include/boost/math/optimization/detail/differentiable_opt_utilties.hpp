#ifndef DIFFERENTIABLE_OPT_UTILITIES_HPP
#define DIFFERENTIABLE_OPT_UTILITIES_HPP
#include <cmath>
#include <random>
#include <type_traits>
#include <vector>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <vector>

namespace boost {
namespace math {
namespace optimization {
template<typename UpdPol>
struct update_policy_real_type;

template<typename UpdPol>
struct update_policy_real_type;

template<template<typename> class UpdPol, typename RealType>
struct update_policy_real_type<UpdPol<RealType>>
{
    using type = RealType;
};

template<typename UpdPol>
using update_policy_real_type_t =
    typename update_policy_real_type<typename std::decay<UpdPol>::type>::type;

template<class RealType>
RealType gradient_norm2(const std::vector<std::reference_wrapper<RealType>>& g)
{
    /* @brief computes 2 norm of a vector of reference wrapped RealTypes
	*/
    RealType sum = RealType(0);
    for (auto& gi : g) {
        RealType val = gi.get();
        sum += val * val;
    }
    return sqrt(sum);
}

template<class RealType>
RealType gradient_norm1(const std::vector<std::reference_wrapper<RealType>>& g)
{
    /* @brief computes 2 norm of a vector of reference wrapped RealTypes
	*/
    RealType sum = RealType(0);
    for (auto& gi : g) {
        RealType val = gi.get();
        sum += abs(val);
    }
    return sqrt(sum);
}

template<typename RealType>
std::vector<RealType> random_vector(size_t n)
{
    /*@brief> generates a random std::vector<RealType> of size n
     * using mt19937 algorithm
     */

    /* TODO: these may need to be marked thread local 
     * in the future
	* 
	* TODO: benchmark.
	*/
    static boost::random::mt19937                             rng{std::random_device{}()};
    static boost::random::uniform_real_distribution<RealType> dist(0.0, 1.0);

    std::vector<RealType> result(n);
    std::generate(result.begin(), result.end(), [&] { return dist(rng); });
    return result;
}

} // namespace optimization
} // namespace math
} // namespace boost
#endif
