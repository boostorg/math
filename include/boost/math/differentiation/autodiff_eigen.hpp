#ifndef BOOST_MATH_DIFFERENTIATION_AUTODIFF_EIGEN_HPP
#define BOOST_MATH_DIFFERENTIATION_AUTODIFF_EIGEN_HPP

#include <Eigen/Core>
#include <boost/math/differentiation/autodiff.hpp>

namespace Eigen {
template <typename RealType, size_t Order>
struct NumTraits<boost::math::differentiation::autodiff_v1::detail::
                     template fvar<RealType, Order>> : NumTraits<RealType> {
  using fvar =
      boost::math::differentiation::autodiff_v1::detail::template fvar<RealType,
                                                                       Order>;

  enum {
    RequireInitialization = 1,
    ReadCost = 1,
    AddCost = 16,
    MulCost = 16,
  };
};

template <typename RealType, size_t Order, typename BinaryOp, typename A>
struct ScalarBinaryOpTraits<boost::math::differentiation::autodiff_v1::detail::
                                template fvar<RealType, Order>,
                            A, BinaryOp> {
  typedef boost::math::differentiation::autodiff_v1::detail::template fvar<
      RealType, Order>
      ReturnType;
};

template <typename RealType, size_t Order, typename BinaryOp, typename A>
struct ScalarBinaryOpTraits<A,
                            boost::math::differentiation::autodiff_v1::detail::
                                template fvar<RealType, Order>,
                            BinaryOp> {
  typedef boost::math::differentiation::autodiff_v1::detail::template fvar<
      RealType, Order>
      ReturnType;
};

template <typename RealType, size_t Order, typename RealType2, size_t Order2,
          typename BinaryOp>
struct ScalarBinaryOpTraits<
    boost::math::differentiation::autodiff_v1::detail::template fvar<RealType,
                                                                     Order>,
    boost::math::differentiation::autodiff_v1::detail::template fvar<RealType2,
                                                                     Order2>,
    BinaryOp> {
  typedef ScalarBinaryOpTraits<RealType, RealType2, BinaryOp>::ReturnType
      RealReturn;
  const size_t ReturnOrder = (Order > Order2) ? Order : Order2;
  typedef boost::math::differentiation::autodiff_v1::detail::template fvar<
      RealReturn, ReturnOrder>
      ReturnType;
}
}  // namespace Eigen

#endif  // BOOST_MATH_DIFFERENTIATION_AUTODIFF_EIGEN_HPP
