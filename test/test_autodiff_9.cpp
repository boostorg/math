#include <Eigen/Core>
#include <Eigen/Dense>

#include "test_autodiff.hpp"

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

#define BOOST_AUTODIFF_EIGEN_SCALAR_TRAITS(A)                                 \
  template <class RealType, size_t Order, typename BinaryOp>                  \
  struct ScalarBinaryOpTraits<boost::math::differentiation::autodiff_v1::     \
                                  detail::template fvar<RealType, Order>,     \
                              A, BinaryOp> {                                  \
    typedef boost::math::differentiation::autodiff_v1::detail::template fvar< \
        RealType, Order>                                                      \
        ReturnType;                                                           \
  };                                                                          \
  template <class RealType, size_t Order, typename BinaryOp>                  \
  struct ScalarBinaryOpTraits<A,                                              \
                              boost::math::differentiation::autodiff_v1::     \
                                  detail::template fvar<RealType, Order>,     \
                              BinaryOp> {                                     \
    typedef boost::math::differentiation::autodiff_v1::detail::template fvar< \
        RealType, Order>                                                      \
        ReturnType;                                                           \
  };

BOOST_AUTODIFF_EIGEN_SCALAR_TRAITS(float);
BOOST_AUTODIFF_EIGEN_SCALAR_TRAITS(double);
BOOST_AUTODIFF_EIGEN_SCALAR_TRAITS(long double);
BOOST_AUTODIFF_EIGEN_SCALAR_TRAITS(short);
BOOST_AUTODIFF_EIGEN_SCALAR_TRAITS(unsigned short);
BOOST_AUTODIFF_EIGEN_SCALAR_TRAITS(int);
BOOST_AUTODIFF_EIGEN_SCALAR_TRAITS(unsigned int);
BOOST_AUTODIFF_EIGEN_SCALAR_TRAITS(long);
BOOST_AUTODIFF_EIGEN_SCALAR_TRAITS(unsigned long);

}  // namespace Eigen

BOOST_AUTO_TEST_SUITE(test_autodiff_9)

using boost::math::differentiation::autodiff_v1::detail::fvar;

BOOST_AUTO_TEST_CASE_TEMPLATE(eigen_init, T, all_float_types) {
  constexpr int size = 5;

  constexpr std::size_t n = 5;
  typedef fvar<T, n> fTn;
  Eigen::Matrix<fTn, size, 1> x;
  x[0] = fTn(1.5);
  x[1] = make_fvar<T, n - 1>(2.5);
  x[2] = fvar<T, n - 2>(3.5);
  x[3] = 4;
  x[4] = 5.5;

  constexpr std::size_t m = 2;
  typedef fvar<T, m> fTm;
  Eigen::Matrix<fTm, size, 1> y;
  y = x.template cast<fTm>();
  BOOST_CHECK_EQUAL(x[0].derivative(0), y[0].derivative(0));
  BOOST_CHECK_EQUAL(x[1].derivative(0), y[1].derivative(0));
  BOOST_CHECK_EQUAL(x[2].derivative(0), y[2].derivative(0));
  BOOST_CHECK_EQUAL(x[3].derivative(0), y[3].derivative(0));
  BOOST_CHECK_EQUAL(x[4].derivative(0), y[4].derivative(0));

  // constexpr std::size_t p = 3;
  // typedef fvar<T, p> fTp;
  // Eigen::Matrix<fTp, 1, size> z =
  //     Eigen::Matrix<fTn, size, 1>::Random().transpose().template cast<fTp>();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(eigen_general, T, all_float_types) {
  using std::cos;
  using std::exp;
  using std::log;
  using std::pow;
  using std::sin;

  constexpr int dim = 4;
  constexpr std::size_t n = 2;
  constexpr double p = 3.456;

  typedef fvar<T, n> fTn;
  Eigen::Matrix<fTn, dim, 1> x;
  x[0] = make_fvar<T, n>(-1);
  x[1] = make_fvar<T, n>(0);
  x[2] = make_fvar<T, n>(1);
  x[3] = make_fvar<T, n>(5);

  Eigen::Matrix<fTn, dim, 1> y1, y2, y3, y4, y5;
  y1 = x.array().sin();
  y2 = x.array().cos();
  y3 = x.array().exp();
  y4 = x.array().pow(p);
  y5 = x.array().log();

  // Check sin
  BOOST_CHECK_EQUAL(y1[0].derivative(0), sin(x[0].derivative(0)));
  BOOST_CHECK_EQUAL(y1[1].derivative(0), sin(x[1].derivative(0)));
  BOOST_CHECK_EQUAL(y1[2].derivative(0), sin(x[2].derivative(0)));
  BOOST_CHECK_EQUAL(y1[3].derivative(0), sin(x[3].derivative(0)));
  BOOST_CHECK_EQUAL(y1[0].derivative(1), cos(x[0].derivative(0)));
  BOOST_CHECK_EQUAL(y1[1].derivative(1), cos(x[1].derivative(0)));
  BOOST_CHECK_EQUAL(y1[2].derivative(1), cos(x[2].derivative(0)));
  BOOST_CHECK_EQUAL(y1[3].derivative(1), cos(x[3].derivative(0)));
  BOOST_CHECK_EQUAL(y1[0].derivative(2), -sin(x[0].derivative(0)));
  BOOST_CHECK_EQUAL(y1[1].derivative(2), -sin(x[1].derivative(0)));
  BOOST_CHECK_EQUAL(y1[2].derivative(2), -sin(x[2].derivative(0)));
  BOOST_CHECK_EQUAL(y1[3].derivative(2), -sin(x[3].derivative(0)));

  // Check cos
  BOOST_CHECK_EQUAL(y2[0].derivative(0), cos(x[0].derivative(0)));
  BOOST_CHECK_EQUAL(y2[1].derivative(0), cos(x[1].derivative(0)));
  BOOST_CHECK_EQUAL(y2[2].derivative(0), cos(x[2].derivative(0)));
  BOOST_CHECK_EQUAL(y2[3].derivative(0), cos(x[3].derivative(0)));
  BOOST_CHECK_EQUAL(y2[0].derivative(1), -sin(x[0].derivative(0)));
  BOOST_CHECK_EQUAL(y2[1].derivative(1), -sin(x[1].derivative(0)));
  BOOST_CHECK_EQUAL(y2[2].derivative(1), -sin(x[2].derivative(0)));
  BOOST_CHECK_EQUAL(y2[3].derivative(1), -sin(x[3].derivative(0)));
  BOOST_CHECK_EQUAL(y2[0].derivative(2), -cos(x[0].derivative(0)));
  BOOST_CHECK_EQUAL(y2[1].derivative(2), -cos(x[1].derivative(0)));
  BOOST_CHECK_EQUAL(y2[2].derivative(2), -cos(x[2].derivative(0)));
  BOOST_CHECK_EQUAL(y2[3].derivative(2), -cos(x[3].derivative(0)));

  // Check exp
  BOOST_CHECK_EQUAL(y3[0].derivative(0), exp(x[0].derivative(0)));
  BOOST_CHECK_EQUAL(y3[1].derivative(0), exp(x[1].derivative(0)));
  BOOST_CHECK_EQUAL(y3[2].derivative(0), exp(x[2].derivative(0)));
  BOOST_CHECK_EQUAL(y3[3].derivative(0), exp(x[3].derivative(0)));
  BOOST_CHECK_EQUAL(y3[0].derivative(1), exp(x[0].derivative(0)));
  BOOST_CHECK_EQUAL(y3[1].derivative(1), exp(x[1].derivative(0)));
  BOOST_CHECK_EQUAL(y3[2].derivative(1), exp(x[2].derivative(0)));
  BOOST_CHECK_EQUAL(y3[3].derivative(1), exp(x[3].derivative(0)));
  BOOST_CHECK_EQUAL(y3[0].derivative(2), exp(x[0].derivative(0)));
  BOOST_CHECK_EQUAL(y3[1].derivative(2), exp(x[1].derivative(0)));
  BOOST_CHECK_EQUAL(y3[2].derivative(2), exp(x[2].derivative(0)));
  BOOST_CHECK_EQUAL(y3[3].derivative(2), exp(x[3].derivative(0)));

  // Check pow (without negative or zero)
  T powTol = 1e-4;
  BOOST_CHECK_CLOSE(y4[2].derivative(0), pow(x[2].derivative(0), p), powTol);
  BOOST_CHECK_CLOSE(y4[3].derivative(0), pow(x[3].derivative(0), p), powTol);
  BOOST_CHECK_CLOSE(y4[2].derivative(1), p * pow(x[2].derivative(0), p - 1),
                    powTol);
  BOOST_CHECK_CLOSE(y4[3].derivative(1), p * pow(x[3].derivative(0), p - 1),
                    powTol);
  BOOST_CHECK_CLOSE(y4[2].derivative(2),
                    (p - 1) * p * pow(x[2].derivative(0), p - 2), powTol);
  BOOST_CHECK_CLOSE(y4[3].derivative(2),
                    (p - 1) * p * pow(x[3].derivative(0), p - 2), powTol);

  // Check log (without negative or zero)
  T logTol = 1e-5;
  BOOST_CHECK_CLOSE(y5[2].derivative(0), log(x[2].derivative(0)), logTol);
  BOOST_CHECK_CLOSE(y5[3].derivative(0), log(x[3].derivative(0)), logTol);
  BOOST_CHECK_CLOSE(y5[2].derivative(1), 1 / x[2].derivative(0), logTol);
  BOOST_CHECK_CLOSE(y5[3].derivative(1), 1 / x[3].derivative(0), logTol);
  BOOST_CHECK_CLOSE(y5[2].derivative(2), -1 / pow(x[2].derivative(0), 2),
                    logTol);
  BOOST_CHECK_CLOSE(y5[3].derivative(2), -1 / pow(x[3].derivative(0), 2),
                    logTol);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(eigen_scalar, T, all_float_types) {
  constexpr int dim = 4;
  constexpr size_t n = 4;

  typedef fvar<T, n> fTn;
  fTn x = make_fvar<T, n>(4);
  Eigen::Matrix<fTn, dim, 1> X;
  Eigen::Matrix<fTn, dim, 1> Z;
  Eigen::Matrix<fTn, dim, dim> I = Eigen::Matrix<fTn, dim, dim>::Identity();
  X[0] = x;
  Z[0] = x;
  X[1] = 2;
  Z[1] = x;
  X[2] = x;
  Z[2] = -3;
  X[3] = 4;
  Z[3] = 5;

  // y = x*x + 2*x - 3*x + 20
  fTn y1 = X.transpose() * Z;
  fTn y2 = Z.transpose() * I * X;

  BOOST_CHECK_EQUAL(y1.derivative(0), x * x + 2 * x - 3 * x + 20);
  BOOST_CHECK_EQUAL(y1.derivative(1), 2 * x - 1);
  BOOST_CHECK_EQUAL(y1.derivative(2), 2);
  BOOST_CHECK_EQUAL(y1.derivative(3), 0);
  BOOST_CHECK_EQUAL(y1.derivative(4), 0);
  BOOST_CHECK_EQUAL(y2.derivative(0), x * x + 2 * x - 3 * x + 20);
  BOOST_CHECK_EQUAL(y2.derivative(1), 2 * x - 1);
  BOOST_CHECK_EQUAL(y2.derivative(2), 2);
  BOOST_CHECK_EQUAL(y2.derivative(3), 0);
  BOOST_CHECK_EQUAL(y2.derivative(4), 0);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(eigen_vector, T, all_float_types) {
  // Note: Can't handle taking gradient all at once. That would require the same
  // approach as `make_ftuple`, which has upper limits.
  using std::cos;
  using std::sin;

  constexpr int dim = 3;
  constexpr size_t n = 4;

  typedef fvar<T, n> fTn;
  fTn x = make_fvar<T, n>(5);
  T xD0 = x.derivative(0);
  Eigen::Matrix<fTn, dim, 1> X;
  X[0] = 1;
  X[1] = x;
  X[2] = x * x;
  Eigen::Matrix<fTn, dim, dim> M;
  M(0, 0) = 1;
  M(0, 1) = 2;
  M(0, 2) = x;
  M(1, 0) = 1 / x;
  M(1, 1) = 4;
  M(1, 2) = 0;
  M(2, 0) = 5;
  M(2, 1) = sin(x);
  M(2, 2) = cos(x * x);

  Eigen::Matrix<fTn, dim, 1> Y = M * X;

  T powTol = 1e-4;

  // Y[0] = 1 + 2*x + x*x*x
  BOOST_CHECK_CLOSE(Y[0].derivative(0), 1 + 2 * xD0 + pow(xD0, 3), powTol);
  BOOST_CHECK_EQUAL(Y[0].derivative(1), 2 + 3 * xD0 * xD0);
  BOOST_CHECK_EQUAL(Y[0].derivative(2), 6 * xD0);
  BOOST_CHECK_EQUAL(Y[0].derivative(3), 6);
  BOOST_CHECK_EQUAL(Y[0].derivative(4), 0);

  // Y[1] = 1/x + 4*x + 0
  BOOST_CHECK_EQUAL(Y[1].derivative(0), 1 / xD0 + 4 * xD0);
  BOOST_CHECK_EQUAL(Y[1].derivative(1), -1 / (xD0 * xD0) + 4);
  BOOST_CHECK_CLOSE(Y[1].derivative(2), 2 / (pow(xD0, 3)), powTol);
  BOOST_CHECK_CLOSE(Y[1].derivative(3), -6 / (pow(xD0, 3) * xD0), powTol);
  BOOST_CHECK_CLOSE(Y[1].derivative(4), 24 / (pow(xD0, 3) * xD0 * xD0), powTol);

  // Y[2] = 5 + x*sin(x) + x*x*cos(x*x)
  BOOST_CHECK_EQUAL(Y[2].derivative(0),
                    5 + xD0 * sin(xD0) + xD0 * xD0 * cos(xD0 * xD0));
  BOOST_CHECK_CLOSE(Y[2].derivative(1),
                    2 * xD0 * cos(xD0 * xD0) -
                        2 * pow(xD0, 3) * sin(xD0 * xD0) + sin(xD0) +
                        xD0 * cos(xD0),
                    powTol);
  BOOST_CHECK_CLOSE(Y[2].derivative(2),
                    -xD0 * (10 * xD0 * sin(xD0 * xD0) + sin(xD0)) +
                        (2 - 4 * pow(xD0, 4)) * cos(xD0 * xD0) + 2 * cos(xD0),
                    powTol);
  BOOST_CHECK_CLOSE(
      Y[2].derivative(3),
      -24 * xD0 * sin(xD0 * xD0) + 8 * pow(xD0, 5) * sin(xD0 * xD0) -
          36 * pow(xD0, 3) * cos(xD0 * xD0) - 3 * sin(xD0) - xD0 * cos(xD0),
      powTol);
  BOOST_CHECK_CLOSE(
      Y[2].derivative(4),
      -24 * sin(xD0 * xD0) + 112 * pow(xD0, 4) * sin(xD0 * xD0) +
          4 * (4 * pow(xD0, 4) - 39) * xD0 * xD0 * cos(xD0 * xD0) +
          xD0 * sin(xD0) - 4 * cos(xD0),
      powTol);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(eigen_determinant, T, all_float_types) {
  using std::cos;
  using std::sin;
  constexpr int dim = 4;
  constexpr size_t n = 1;

  typedef fvar<T, n> fTn;
  fTn x = make_fvar<T, n>(3);
  T xD0 = x.derivative(0);
  Eigen::Matrix<fTn, dim, dim> M = 10 * Eigen::Matrix<fTn, dim, dim>::Random();
  M(0, 3) = x;
  M(1, 0) = 3 * x;
  M(1, 2) = 1 / x;
  M(2, 2) = x * x;
  M(3, 1) = sin(x);
  fTn y = M.determinant();

  Eigen::Matrix<T, dim, dim> dMdx = Eigen::Matrix<T, dim, dim>::Zero();
  dMdx(0, 3) = 1;
  dMdx(1, 0) = 3;
  dMdx(1, 2) = -1 / (xD0 * xD0);
  dMdx(2, 2) = 2 * xD0;
  dMdx(3, 1) = cos(xD0);

  T detTol = 1e-3;

  T ans = y.derivative(0) *
          (M.inverse() * dMdx.template cast<fTn>()).trace().derivative(0);
  BOOST_CHECK_CLOSE(y.derivative(1), ans, detTol);
}

BOOST_AUTO_TEST_SUITE_END()
