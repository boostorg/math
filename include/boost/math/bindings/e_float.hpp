//  Copyright John Maddock 2012 - 2021.
//  Copyright Christopher Kormanyos 2016 - 2021.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

// Wrapper that works with e_float-2021, see also
// https://github.com/ckormanyos/e_float-2021

#ifndef E_FLOAT_2017_08_18_HPP_
  #define E_FLOAT_2017_08_18_HPP_

  #include <cstdint>
  #include <sstream>
  #include <tuple>
  #include <type_traits>

  #include <e_float/e_float.h>
  #include <e_float/e_float_functions.h>

  #include <boost/config.hpp>
  #include <boost/multiprecision/number.hpp>

  namespace boost { namespace math { namespace ef {

  // Forward declaration of the e_float multiple precision class.
  // This class binds native ::e_float to boost::multiprecsion::e_float.
  class e_float;

  } } }

  // Define the number category as a floating-point kind
  // for the e_float. This is needed for properly
  // interacting as a backend with boost::muliprecision.
  template<>
  struct boost::multiprecision::number_category<boost::math::ef::e_float>
    : public std::integral_constant<int, boost::multiprecision::number_kind_floating_point> { };

  namespace boost { namespace math { namespace ef {

  // This is the e_float multiple precision class.
  class e_float
  {
  public:
    using   signed_types = std::tuple<  signed char,   signed short,   signed int,   signed long,   signed long long, std::intmax_t>;
    using unsigned_types = std::tuple<unsigned char, unsigned short, unsigned int, unsigned long, unsigned long long, std::uintmax_t>;
    using float_types    = std::tuple<float, double, long double>;
    using exponent_type  = std::int64_t;

    e_float() : m_value() { }

    explicit e_float(const ::e_float& rep) : m_value(rep) { }

    e_float(const e_float& other) : m_value(other.m_value) { }
    e_float(e_float&& other) : m_value(static_cast<::e_float&&>(other.m_value)) { }

    template<typename UnsignedIntegralType,
             typename std::enable_if<(   (std::is_integral<UnsignedIntegralType>::value == true)
                                      && (std::is_unsigned<UnsignedIntegralType>::value == true))>::type const* = nullptr>
    e_float(UnsignedIntegralType u) : m_value(::e_float(std::uint64_t(u))) { }

    template<typename SignedIntegralType,
             typename std::enable_if<(   (std::is_integral<SignedIntegralType>::value == true)
                                      && (std::is_signed  <SignedIntegralType>::value == true))>::type const* = nullptr>
    e_float(SignedIntegralType n) : m_value(::e_float(std::int64_t(n))) { }

    template<typename FloatingPointType,
             typename std::enable_if<std::is_floating_point<FloatingPointType>::value == true>::type const* = nullptr>
    e_float(FloatingPointType f) : m_value(::e_float(static_cast<long double>(f))) { }

    e_float(const char* c) : m_value(c) { }

    e_float(const std::string& str) : m_value(str) { }

    ~e_float() { }

    e_float& operator=(const e_float& other)
    {
      if(this != &other)
      {
        m_value = other.m_value;
      }

      return *this;
    }

    e_float& operator=(e_float&& other)
    {
      m_value.operator=(static_cast<::e_float&&>(other.m_value));

      return *this;
    }

    template<typename ArithmeticType,
             typename std::enable_if<std::is_arithmetic<ArithmeticType>::value == true>::type const* = nullptr>
    e_float& operator=(const ArithmeticType& x)
    {
      m_value = ::e_float(x);

      return *this;
    }

    e_float& operator=(const std::string& str_rep)  { m_value = ::e_float(str_rep);  return *this; }
    e_float& operator=(const char*        char_ptr) { m_value = ::e_float(char_ptr); return *this; }

    void swap(e_float& other_mp_cpp_backend)
    {
      m_value.swap(other_mp_cpp_backend.m_value);
    }

          ::e_float&  representation()       { return m_value; }
    const ::e_float&  representation() const { return m_value; }
    const ::e_float& crepresentation() const { return m_value; }

    std::string str(std::streamsize number_of_digits, const std::ios::fmtflags format_flags) const
    {
      std::string        my_result_str;
      std::stringstream  my_stream_str;

      my_stream_str.flags(format_flags);

      static_cast<void>(my_stream_str.precision(number_of_digits));

      m_value.wr_string(my_result_str, my_stream_str);

      return my_result_str;
    }

    void negate()
    {
      m_value.negate();
    }

    int compare(const e_float& other_mp_cpp_backend) const
    {
      return static_cast<int>(m_value.compare(other_mp_cpp_backend.m_value));
    }

    template<typename ArithmeticType,
             typename std::enable_if<std::is_arithmetic<ArithmeticType>::value == true>::type const* = nullptr>
    int compare(ArithmeticType x)
    {
      return static_cast<int>(m_value.compare(::e_float(x)));
    }

  private:
    ::e_float m_value;

    e_float& operator=(const ::e_float&) = delete;
  };

  inline void eval_add(e_float& result, const e_float& x)
  {
    result.representation() += x.crepresentation();
  }

  inline void eval_subtract(e_float& result, const e_float& x)
  {
    result.representation() -= x.crepresentation();
  }

  inline void eval_multiply(e_float& result, const e_float& x)
  {
    result.representation() *= x.crepresentation();
  }

  inline void eval_multiply(e_float& result, const e_float& x, const e_float& y)
  {
    result.representation() = x.crepresentation() * y.crepresentation();
  }

  template<typename SignedIntegralType,
           typename std::enable_if<(   (std::is_integral<SignedIntegralType>::value == true)
                                    && (std::is_signed  <SignedIntegralType>::value == true))>::type const* = nullptr>
  void eval_multiply(e_float& result, const SignedIntegralType& n)
  {
    result.representation().mul_signed_long_long(static_cast<std::int64_t>(n));
  }

  inline void eval_divide(e_float& result, const e_float& x)
  {
    result.representation() /= x.crepresentation();
  }

  template<typename SignedIntegralType,
           typename std::enable_if<(   (std::is_integral<SignedIntegralType>::value == true)
                                    && (std::is_signed  <SignedIntegralType>::value == true))>::type const* = nullptr>
  void eval_divide(e_float& result, const SignedIntegralType& n)
  {
    result.representation().div_signed_long_long(static_cast<std::int64_t>(n));
  }

  inline bool eval_eq(const e_float& a, const e_float& b)
  {
    return (a.compare(b) == 0);
  }

  template<typename ArithmeticType,
           typename std::enable_if<std::is_arithmetic<ArithmeticType>::value == true>::type const* = nullptr>
  bool eval_eq(const e_float& a, const ArithmeticType& b)
  {
    return (a.compare(b) == 0);
  }

  template<typename ArithmeticType,
           typename std::enable_if<std::is_arithmetic<ArithmeticType>::value == true>::type const* = nullptr>
  bool eval_eq(const ArithmeticType& a, const e_float& b)
  {
    return (e_float(a).compare(b) == 0);
  }

  inline bool eval_gt(const e_float& a, const e_float& b)
  {
    return (a.compare(b) == 1);
  }

  template<typename ArithmeticType,
           typename std::enable_if<std::is_arithmetic<ArithmeticType>::value == true>::type const* = nullptr>
  bool eval_gt(const e_float& a, const ArithmeticType& b)
  {
    return (a.compare(b) == 1);
  }

  template<typename ArithmeticType,
           typename std::enable_if<std::is_arithmetic<ArithmeticType>::value>::type const* = nullptr>
  bool eval_gt(const ArithmeticType& a, const e_float& b)
  {
    return (e_float(a).compare(b) == 1);
  }

  inline bool eval_lt(const e_float& a, const e_float& b)
  {
    return (a.compare(b) == -1);
  }

  template<typename ArithmeticType,
           typename std::enable_if<std::is_arithmetic<ArithmeticType>::value == true>::type const* = nullptr>
  bool eval_lt(const e_float& a, const ArithmeticType& b)
  {
    return (a.compare(b) == -1);
  }

  template<typename ArithmeticType,
           typename std::enable_if<std::is_arithmetic<ArithmeticType>::value == true>::type const* = nullptr>
  bool eval_lt(const ArithmeticType& a, const e_float& b)
  {
    return (e_float(a).compare(b) == -1);
  }

  inline bool eval_is_zero(const e_float& x)
  {
    return x.crepresentation().iszero();
  }

  inline int eval_get_sign(const e_float& x)
  {
    if     (x.crepresentation().iszero()) { return  0; }
    else if(x.crepresentation().isneg ()) { return -1; }
    else                                  { return  1; }
  }

  inline void eval_convert_to(unsigned long long* result,
                              const e_float& val)
  {
    *result = (val.crepresentation()).extract_unsigned_long_long();
  }

  inline void eval_convert_to(signed long long* result,
                              const e_float& val)
  {
    *result = (val.crepresentation()).extract_signed_long_long();
  }

  inline void eval_convert_to(long double* result,
                              const e_float& val)
  {
    *result = (val.crepresentation()).extract_long_double();
  }

  inline void eval_frexp(      e_float&                         result,
                         const e_float&                         x,
                               typename e_float::exponent_type* expptr)
  {
    using local_exponent_type = int;

    local_exponent_type exp2;

    result.representation() = ::ef::frexp(x.crepresentation(), &exp2);

    *expptr = static_cast<typename e_float::exponent_type>(exp2);
  }

  inline void eval_frexp(      e_float& result,
                         const e_float& x,
                         int*           expptr,
                         typename std::enable_if<std::is_same<typename e_float::exponent_type, int>::value == false>::type const* = nullptr)
  {
    result.representation() = ::ef::frexp(x.crepresentation(), expptr);
  }

  inline void eval_ldexp(      e_float& result,
                         const e_float& x,
                         const typename e_float::exponent_type exp_value)
  {
    using local_exponent_type = int;

    result.representation() = ::ef::ldexp(x.crepresentation(), local_exponent_type(exp_value));
  }

  inline void eval_ldexp(      e_float& result,
                         const e_float& x,
                         int            exp_value,
                         typename std::enable_if<std::is_same<typename e_float::exponent_type, int>::value == false>::type const* = nullptr)
  {
    using local_exponent_type = int;

    result.representation() = ::ef::ldexp(x.crepresentation(), local_exponent_type(exp_value));
  }

  inline e_float::exponent_type eval_ilogb(const e_float& val)
  {
    if(val.crepresentation().iszero())
    {
      return (std::numeric_limits<e_float::exponent_type>::min)();
    }

    if ((val.crepresentation().isinf)())
    {
      return (std::numeric_limits<e_float::exponent_type>::max)();
    }

    if ((val.crepresentation().isnan)())
    {
      #if defined(FP_ILOGBNAN)
      return FP_ILOGBNAN;
      #else
      return (std::numeric_limits<e_float::exponent_type>::max)();
      #endif
    }

   // Set the result to base-10 exponent (order) of the value.
   return val.crepresentation().order();
  }

  template<typename IntegralType,
           typename std::enable_if<(   (std::is_fundamental<IntegralType>::value == true)
                                    && (std::is_integral   <IntegralType>::value == true))>::type const* = nullptr>
  void eval_scalbn(e_float& result, const e_float& val, IntegralType e)
  {
    const std::int64_t my_e = static_cast<std::int64_t>(e);

    const e_float t(::e_float(1.0, my_e));

    eval_multiply(result, val, t);
  }

  inline void eval_floor(      e_float& result,
                         const e_float& x)
  {
    result.representation() = ::ef::floor(x.crepresentation());
  }

  inline void eval_ceil(      e_float& result,
                        const e_float& x)
  {
    result.representation() = ::ef::ceil(x.crepresentation());
  }

  inline int eval_fpclassify(const e_float& x)
  {
    if     ((x.crepresentation().isinf)()) { return FP_INFINITE; }
    else if((x.crepresentation().isnan)()) { return FP_NAN; }
    else if( x.crepresentation().iszero()) { return FP_ZERO; }
    else                                   { return FP_NORMAL; }
  }

  inline void eval_trunc(      e_float& result,
                         const e_float& x)
  {
    result.representation() = ::ef::integer_part(x.crepresentation());
  }

  inline void eval_abs(      e_float& result,
                       const e_float& x)
  {
    result.representation() = x.crepresentation();

    if(result.crepresentation().isneg())
    {
      result.representation().negate();
    }
  }

  inline void eval_fabs(      e_float& result,
                        const e_float& x)
  {
    result.representation() = x.crepresentation();

    if(result.crepresentation().isneg())
    {
      result.representation().negate();
    }
  }

  inline void eval_sqrt(      e_float& result,
                        const e_float& x)
  {
    result.representation() = ::ef::sqrt(x.crepresentation());
  }

  inline void eval_sin(      e_float& result,
                       const e_float& x)
  {
    result.representation() = ::ef::sin(x.crepresentation());
  }

  inline void eval_cos(      e_float& result,
                       const e_float& x)
  {
    result.representation() = ::ef::cos(x.crepresentation());
  }

  inline void eval_tan(      e_float& result,
                       const e_float& x)
  {
    result.representation() = ::ef::tan(x.crepresentation());
  }

  inline void eval_asin(      e_float& result,
                        const e_float& x)
  {
    result.representation() = ::ef::asin(x.crepresentation());
  }

  inline void eval_acos(      e_float& result,
                        const e_float& x)
  {
    result.representation() = ::ef::acos(x.crepresentation());
  }

  inline void eval_atan(      e_float& result,
                        const e_float& x)
  {
    result.representation() = ::ef::atan(x.crepresentation());
  }

  inline void eval_atan2(      e_float& result,
                         const e_float& y,
                         const e_float& x)
  {
    result.representation() = ::ef::atan2(y.crepresentation(), x.crepresentation());
  }

  inline void eval_log(      e_float& result,
                       const e_float& x)
  {
    result.representation() = ::ef::log(x.crepresentation());
  }

  inline void eval_log10(      e_float& result,
                         const e_float& x)
  {
    result.representation() = ::ef::log10(x.crepresentation());
  }

  inline void eval_exp(      e_float& result,
                       const e_float& x)
  {
    result.representation() = ::ef::exp(x.crepresentation());
  }

  inline void eval_sinh(      e_float& result,
                        const e_float& x)
  {
    result.representation() = ::ef::sinh(x.crepresentation());
  }

  inline void eval_cosh(      e_float& result,
                        const e_float& x)
  {
    result.representation() = ::ef::cosh(x.crepresentation());
  }

  inline void eval_tanh(      e_float& result,
                        const e_float& x)
  {
    result.representation() = ::ef::tanh(x.crepresentation());
  }

  inline void eval_fmod(      e_float& result,
                        const e_float& x,
                        const e_float& y)
  {
    if(y.crepresentation().iszero())
    {
      result.representation() = ::ef::zero();
    }
    else
    {
      // Calculate the fractional part of x such that:
      //   x = (integer_part * y) + fractional_part,
      // where
      //   |fractional_part| < |y|,
      // and fractional_part has the same sign as x.

      const ::e_float integer_part = ::ef::floor(x.crepresentation() / y.crepresentation());

      result.representation() =
        x.crepresentation() - (integer_part * y.crepresentation());

      if(x.crepresentation().isneg() != y.crepresentation().isneg())
      {
        result.representation() -= y.crepresentation();
      }
    }
  }

  inline void eval_pow(      e_float& result,
                       const e_float& x,
                       const e_float& a)
  {
    result.representation() = ::ef::pow(x.crepresentation(), a.crepresentation());
  }

  } } } // namespace boost::math::ef

  namespace boost { namespace math { namespace policies {

  // Specialization of the precision structure.
  template<typename ThisPolicy,
           const boost::multiprecision::expression_template_option ExpressionTemplates>
  struct precision<boost::multiprecision::number<boost::math::ef::e_float,
                                                 ExpressionTemplates>,
                   ThisPolicy>
  {
    using precision_type = typename ThisPolicy::precision_type;

    using local_digits_2 = digits2<((::e_float::ef_digits10 + 1LL) * 1000LL) / 301LL>;

    using type = typename std::conditional<((local_digits_2::value <= precision_type::value) || (precision_type::value <= 0)),
                                           local_digits_2,
                                           precision_type>::type;
  };

  } } } // namespaces boost::math::policies

  namespace std
  {
    template<const boost::multiprecision::expression_template_option ExpressionTemplates>
    class numeric_limits<boost::multiprecision::number<boost::math::ef::e_float,
                                                       ExpressionTemplates>>
    {
    public:
      static constexpr bool is_specialized = true;
      static constexpr bool is_signed      = true;
      static constexpr bool is_integer     = false;
      static constexpr bool is_exact       = false;
      static constexpr bool is_bounded     = true;
      static constexpr bool is_modulo      = false;
      static constexpr bool is_iec559      = false;
      static constexpr int  digits         = ::e_float::ef_digits10;
      static constexpr int  digits10       = ::e_float::ef_digits10;
      static constexpr int  max_digits10   = ::e_float::ef_digits10 + 1;

      static constexpr std::int64_t max_exponent   = ::e_float::ef_max_exp;
      static constexpr std::int64_t max_exponent10 = ::e_float::ef_max_exp10;
      static constexpr std::int64_t min_exponent   = ::e_float::ef_min_exp;
      static constexpr std::int64_t min_exponent10 = ::e_float::ef_min_exp10;

      static constexpr int                     radix             = 10;
      static constexpr std::float_round_style  round_style       = std::round_indeterminate;
      static constexpr bool                    has_infinity      = true;
      static constexpr bool                    has_quiet_NaN     = true;
      static constexpr bool                    has_signaling_NaN = false;
      static constexpr std::float_denorm_style has_denorm        = std::denorm_absent;
      static constexpr bool                    has_denorm_loss   = false;
      static constexpr bool                    traps             = false;
      static constexpr bool                    tinyness_before   = false;

      static constexpr boost::multiprecision::number<boost::math::ef::e_float, ExpressionTemplates> (min)        () { return boost::math::ef::e_float(ef::value_min()); }
      static constexpr boost::multiprecision::number<boost::math::ef::e_float, ExpressionTemplates> (max)        () { return boost::math::ef::e_float(ef::value_max()); }
      static constexpr boost::multiprecision::number<boost::math::ef::e_float, ExpressionTemplates> lowest       () { return boost::math::ef::e_float(ef::zero()); }
      static constexpr boost::multiprecision::number<boost::math::ef::e_float, ExpressionTemplates> epsilon      () { return boost::math::ef::e_float(ef::value_eps()); }
      static constexpr boost::multiprecision::number<boost::math::ef::e_float, ExpressionTemplates> round_error  () { return boost::math::ef::e_float(ef::half()); }
      static constexpr boost::multiprecision::number<boost::math::ef::e_float, ExpressionTemplates> infinity     () { return boost::math::ef::e_float(ef::value_inf()); }
      static constexpr boost::multiprecision::number<boost::math::ef::e_float, ExpressionTemplates> quiet_NaN    () { return boost::math::ef::e_float(ef::value_nan()); }
      static constexpr boost::multiprecision::number<boost::math::ef::e_float, ExpressionTemplates> signaling_NaN() { return boost::math::ef::e_float(ef::zero()); }
      static constexpr boost::multiprecision::number<boost::math::ef::e_float, ExpressionTemplates> denorm_min   () { return boost::math::ef::e_float(ef::zero()); }
    };
  } // namespace std

#endif // E_FLOAT_2017_08_18_HPP_
