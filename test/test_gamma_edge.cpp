// Copyright 2024 Christopher Kormanyos
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/math/special_functions/gamma.hpp>
#include <boost/core/lightweight_test.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>

namespace local
{
  template<typename NumericType>
  auto is_close_fraction(const NumericType& a,
                         const NumericType& b,
                         const NumericType& tol) noexcept -> bool
  {
    using std::fabs;

    auto result_is_ok = bool { };

    if(b == static_cast<NumericType>(0))
    {
      result_is_ok = (fabs(a - b) < tol);
    }
    else
    {
      const auto delta = fabs(1 - (a / b));

      result_is_ok = (delta < tol);
    }

    return result_is_ok;
  }

  auto tgamma_under_cbrt_epsilon() -> void
  {
    // This test is intended to hit the lines:

    // template <class T, class Policy>
    // T gamma_imp(T z, const Policy& pol, const lanczos::undefined_lanczos&)
    // ...
    // ...near the comment:
    //      Special case for ultra-small z:

    using local_float_type = boost::multiprecision::number<boost::multiprecision::cpp_dec_float<250>, boost::multiprecision::et_off>;

    static_assert(   (std::numeric_limits<local_float_type>::digits10 >= 248)
                  && (std::numeric_limits<local_float_type>::digits10 <= 252), "Error: Multiprecision wrong number of digits");

    // Table[N[Gamma[n (10^-84)], 260], {n, 1, 10, 1}]
    using local_data_array_type = std::array<local_float_type, static_cast<std::size_t>(UINT8_C(10))>;

    const local_data_array_type ctrl_data =
    {{
      static_cast<local_float_type>("9.9999999999999999999999999999999999999999999999999999999999999999999999999999999999942278433509846713939348790991759756895784066406007640119423276511513227322233532906404199270358122304076394840169355222697867950624167125421902178604121309579597843800360126781E+83"),
      static_cast<local_float_type>("4.9999999999999999999999999999999999999999999999999999999999999999999999999999999999942278433509846713939348790991759756895784066406007640119423276511513227322233532906503104869890919559615934405319418693491786302696988534466221756972734972784545720974832491521E+83"),
      static_cast<local_float_type>("3.3333333333333333333333333333333333333333333333333333333333333333333333333333333333275611766843180047272682124325093090229117399739340973452756609844846560655566866239935343802757050148488807303802815497619037988103143276843874668674681969322826931482456693779E+83"),
      static_cast<local_float_type>("2.4999999999999999999999999999999999999999999999999999999999999999999999999999999999942278433509846713939348790991759756895784066406007640119423276511513227322233532906700916068956514070695013535619545635079623006842631352554860913709962299194441475323232733555E+83"),
      static_cast<local_float_type>("1.9999999999999999999999999999999999999999999999999999999999999999999999999999999999942278433509846713939348790991759756895784066406007640119423276511513227322233532906799821668489311326234553100769609105873541358915452761599180492078575962399389352497160610849E+83"),
      static_cast<local_float_type>("1.6666666666666666666666666666666666666666666666666666666666666666666666666666666666608945100176513380606015457658426423562450733072674306786089943178179893988900199573565393934688775248440759332586339243334126377654940837310166737113856292271003896337573658994E+83"),
      static_cast<local_float_type>("1.4285714285714285714285714285714285714285714285714285714285714285714285714285714285656564147795560999653634505277474042610069780691721925833708990797227513036519247192711918581840620123027917945355450333175663777346809865402105363101517574523570821130186163706E+83"),
      static_cast<local_float_type>("1.2499999999999999999999999999999999999999999999999999999999999999999999999999999999942278433509846713939348790991759756895784066406007640119423276511513227322233532907096538467087703092853171796219799518255296415133916988732139227184416952014232984017855267840E+83"),
      static_cast<local_float_type>("1.1111111111111111111111111111111111111111111111111111111111111111111111111111111111053389544620957825050459902102870868006895177517118751230534387622624338433344644018306555177731611459503822472480974100160325878317849508887569916664141726330291972302168272984E+83"),
      static_cast<local_float_type>("9.9999999999999999999999999999999999999999999999999999999999999999999999999999999999422784335098467139393487909917597568957840664060076401194232765115132273222335329072943496661532976039322509265199264598431331192795598068207783839216442784241287383640775600912E+82"),
    }};

    unsigned index = 1U;

    const local_float_type little { "1E-84" };
    const local_float_type my_tol { std::numeric_limits<local_float_type>::epsilon() * 256 };

    for(const auto& ctrl : ctrl_data)
    {
      const auto x_small = static_cast<local_float_type>(static_cast<local_float_type>(index) * little);

      ++index;

      const auto g_val   = boost::math::tgamma(x_small);

      const auto result_tgamma_x_small_is_ok = is_close_fraction(g_val, ctrl, my_tol);

      BOOST_TEST(result_tgamma_x_small_is_ok);

      if(!result_tgamma_x_small_is_ok)
      {
        break;
      }
    }
  }

  auto tgamma_undefined_lanczos_known_error() -> void
  {
    // This test is intended to hit the lines:

    // template <class T, class Policy>
    // T gamma_imp(T z, const Policy& pol, const lanczos::undefined_lanczos&)
    // ...
    // ...for edge cases that raise errors such as domain error.

    using local_float_type = boost::multiprecision::number<boost::multiprecision::cpp_dec_float<250>, boost::multiprecision::et_off>;

    {
      local_float_type zero { 0 };

      bool domain_error_is_ok { false };

      try
      {
        boost::math::tgamma(zero);
      }
      catch(std::domain_error& err)
      {
        static_cast<void>(err.what());

        domain_error_is_ok = true;
      }

      BOOST_TEST(domain_error_is_ok);
    }

    {
      local_float_type my_nan = std::numeric_limits<local_float_type>::quiet_NaN();

      bool domain_error_is_ok { false };

      try
      {
        boost::math::tgamma(my_nan);
      }
      catch(std::domain_error& err)
      {
        static_cast<void>(err.what());

        domain_error_is_ok = true;
      }

      BOOST_TEST(domain_error_is_ok);
    }

    {
      local_float_type my_inf = -std::numeric_limits<local_float_type>::infinity();

      bool domain_error_is_ok { false };

      try
      {
        boost::math::tgamma(my_inf);
      }
      catch(std::domain_error& err)
      {
        static_cast<void>(err.what());

        domain_error_is_ok = true;
      }

      BOOST_TEST(domain_error_is_ok);
    }
  }
}

auto main() -> int
{
  local::tgamma_under_cbrt_epsilon();
  local::tgamma_undefined_lanczos_known_error();

  const auto result_is_ok = (boost::report_errors() == 0);

  return (result_is_ok ? 0 : -1);
}
