#ifndef EIGEN_SUPPORT_HPP
#define EIGEN_SUPPORT_HPP

#include "Eigen/Eigen"

#include "autodiff_reverse.hpp"

namespace boost {
namespace math {
namespace differentiation {
namespace reverse_mode {

template<typename T, size_t order>
struct rvar;
}
} // namespace differentiation
} // namespace math
} // namespace boost
namespace Eigen {
template<typename T, size_t order>
struct NumTraits<boost::math::differentiation::reverse_mode::rvar<T, order>> : NumTraits<T>
{
    typedef boost::math::differentiation::reverse_mode::rvar<T, order> Real;
    typedef boost::math::differentiation::reverse_mode::rvar<T, order> NonInteger;
    typedef boost::math::differentiation::reverse_mode::rvar<T, order> Nested;
    enum {
        IsComplex             = 0,
        IsInteger             = 0,
        IsSigned              = 1,
        RequireInitialization = 1,
        ReadCost              = 1,
        AddCost               = 3,
        MulCost               = 3
    };
};
template<typename T, size_t order, typename BinaryOp>
struct ScalarBinaryOpTraits<boost::math::differentiation::reverse_mode::rvar<T, order>, T, BinaryOp>
{
    typedef boost::math::differentiation::reverse_mode::rvar<T, order> ReturnType;
};

template<typename T, size_t order, typename BinaryOp>
struct ScalarBinaryOpTraits<T, boost::math::differentiation::reverse_mode::rvar<T, order>, BinaryOp>
{
    typedef boost::math::differentiation::reverse_mode::rvar<T, order> ReturnType;
};

}; // namespace Eigen
#endif // EIGENSUPPORT_HPP
