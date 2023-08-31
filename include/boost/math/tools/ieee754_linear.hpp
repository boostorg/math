//  (C) Copyright Ryan Elandt 2023.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_TOOLS_IEEE754_LINEAR_HPP
#define BOOST_MATH_TOOLS_IEEE754_LINEAR_HPP

#ifdef _MSC_VER
#pragma once
#endif

#include <type_traits>

#include <boost/math/tools/config.hpp>  // For BOOST_MATH_FLOAT128_TYPE
#include <boost/cstdfloat.hpp>  // For numeric_limits support for 128 bit floats

namespace boost {
namespace math {
namespace tools {
namespace detail {
namespace ieee754_linear {

   // The `IsFloat32` class contains a static constexpr method `value()` that
   // returns true if the type T is a 32 bit floating point type. This duck type
   // test for 32 bit float types returns true for the following types:
   //   - `float`
   //   - `_Float32` (NOTE: not a type alias for `float`)
   //   - `std_real_concept` when emulating a 32 bit float with EMULATE32
   //   - other types that seem to be 32 bit floats
   class IsFloat32 {
   public:
      template <typename T>
      static constexpr bool value() {
         return std::numeric_limits<T>::is_iec559 &&
                std::numeric_limits<T>::radix == 2 &&
                std::numeric_limits<T>::digits == 24 &&  // Mantissa has 23 bits + 1 implicit bit
                sizeof(T) == 4;
      }
   };

   // The `IsFloat64` class contains a static constexpr method `value()` that
   // returns true if the type T is a 64 bit floating point type. This duck type
   // test for 64 bit float types returns true for the following types:
   //   - `double`
   //   - `_Float64` (NOTE: not a type alias for `double`)
   //   - `std_real_concept` when emulating a 64 bit float with EMULATE64
   //   - other types that seem to be 64 bit floats
   class IsFloat64 {
   public:
      template <typename T>
      static constexpr bool value() {
         return std::numeric_limits<T>::is_iec559 &&
                std::numeric_limits<T>::radix == 2 &&
                std::numeric_limits<T>::digits == 53 &&  // Mantissa has 52 bits + 1 implicit bit
                sizeof(T) == 8;
      }
   };

   // The `IsFloat128` class contains a static constexpr method `value()` that
   // returns true if the type T is a 128 bit floating point type. This duck type
   // test for 128 bit float types returns true for the following types:
   //   - `__float128`
   //   - `_Float128` (NOTE: not a type alias for `__float128`)
   //   - `std_real_concept` when emulating a 128 bit float with EMULATE128
   //   - other types that seem to be 128 bit floats
   class IsFloat128 {
   public:
      template <typename T>
      static constexpr bool value() {
         return std::numeric_limits<T>::is_iec559 &&
                std::numeric_limits<T>::radix == 2 &&
                std::numeric_limits<T>::digits == 113 &&  // Mantissa has 112 bits + 1 implicit bit
                sizeof(T) == 16;
      }
   };

   // The `Layout_IEEE754Linear` class represents IEEE 754 floating point types
   // for which increasing float values (with the same sign) have increasing 
   // values in bit space. And for which there is a corresponding unsigned integer
   // type U that has the same number of bits as T. Types that satisfy these
   //  conditions include 32, 64, and 128 bit floats. The table below shows select
   // numbers for `float` (single precision).
   //
   //   0            |  0 00000000 00000000000000000000000  |  positive zero
   //   1.4013e-45   |  0 00000000 00000000000000000000001  |  std::numeric_limits<float>::denorm_min()
   //   1.17549e-38  |  0 00000001 00000000000000000000000  |  std::numeric_limits<float>::min()
   //   1.19209e-07  |  0 01101000 00000000000000000000000  |  std::numeric_limits<float>::epsilon()
   //   1            |  0 01111111 00000000000000000000000  |  positive one
   //   3.40282e+38  |  0 11111110 11111111111111111111111  |  std::numeric_limits<float>::max()
   //   inf          |  0 11111111 00000000000000000000000  |  std::numeric_limits<float>::infinity()
   //   nan          |  0 11111111 10000000000000000000000  |  std::numeric_limits<float>::quiet_NaN() 
   //
   // NOTE: the 80 bit `long double` is not "linear" due to the "integer part". See:
   //       https://en.wikipedia.org/wiki/Extended_precision#x86_extended_precision_format
   //
   template <typename T, typename U>
   class Layout_IEEE754Linear {
      static_assert(std::numeric_limits<T>::is_iec559, "Type must be IEEE 754 floating point.");
      static_assert(std::is_unsigned<U>::value, "U must be an unsigned integer type.");
      static_assert(sizeof(T) == sizeof(U), "Type and uint size must be the same.");

   public:
      using type_float = T;  // The linear floating point type
   };

   // The `Layout_Unspecified` class represents floating point types for which
   // are not supported by the `Layout_IEEE754Linear` class.
   template <typename T>
   class Layout_Unspecified {
   public:
      using type_float = T;  // The linear floating point type
   };

   // The `LayoutIdentifier` class identifies the layout type for a floating
   // point type T. The layout type is one of the following:
   //   - `Layout_IEEE754Linear<T, U>` for 32, 64, and 128 bit floats
   //   - `Layout_Unspecified<T>` for other types
   class LayoutIdentifier {
   public:
      // Layout: 32 bit linear
      template <typename T>
      static typename std::enable_if<IsFloat32::value<T>(), Layout_IEEE754Linear<T, std::uint32_t>>::type
      get_layout() { return Layout_IEEE754Linear<T, std::uint32_t>(); }

      // Layout: 64 bit linear
      template <typename T>
      static typename std::enable_if<IsFloat64::value<T>(), Layout_IEEE754Linear<T, std::uint64_t>>::type
      get_layout() { return Layout_IEEE754Linear<T, std::uint64_t>(); }

      // Layout: 128 bit linear
#if defined(BOOST_HAS_INT128) && defined(BOOST_HAS_FLOAT128)
      // NOTE: returns Layout_IEEE754Linear<BOOST_MATH_FLOAT128_TYPE, ...> instead of Layout_IEEE754Linear<T, ...>
      //       in order to work with non-trivial types that wrap trivial 128 bit floating point types.
      template <typename T>
      static typename std::enable_if<IsFloat128::value<T>(), Layout_IEEE754Linear<BOOST_MATH_FLOAT128_TYPE, boost::uint128_type>>::type
      get_layout() { return Layout_IEEE754Linear<BOOST_MATH_FLOAT128_TYPE, boost::uint128_type>(); }

      template <typename T>
      static constexpr bool is_layout_nonlinear() {
         return !IsFloat32::value<T>() &&
                !IsFloat64::value<T>() &&
                !IsFloat128::value<T>();
      }
#else
      template <typename T>
      static constexpr bool is_layout_nonlinear() {
         return !IsFloat32::value<T>() &&
                !IsFloat64::value<T>();
      }
#endif

      // Layout: unspecified
      template <typename T>
      static typename std::enable_if<is_layout_nonlinear<T>(), Layout_Unspecified<T>>::type
      get_layout() { return Layout_Unspecified<T>(); }
   };

   // Base class for the `BitsInflated` and `BitsDeflated` classes.
   //
   // This class stores the bit values of a floating point type `T` as a
   // sign bit as a unsigned integer magnitude. For example, a floating
   // point type with 50 discrete values (zero is counted twice) between
   // -Inf and +Inf is shown below.
   //
   //    -24 -20 -16 -12 -8  -4   0   4   8   12  16  20  24
   //     .   .   .   .   .   .   .   .   .   .   .   .   .
   //     |-----------------------|-----------------------|
   //   -Inf                      0                     +Inf
   //
   template <typename T, typename U>
   class BitsFromZero: public Layout_IEEE754Linear<T, U> {
   public:
      bool sign_bit() const { return sign_bit_; }
      const U& mag() const { return mag_; }
      U& mag() { return mag_; }

   protected:
      BitsFromZero(const bool sign, const U mag) : sign_bit_(sign), mag_(mag) {}

      BitsFromZero(const T x) { 
         // The sign bit is 1, all other bits are 0
         constexpr U bits_sign_mask = U(1) << (sizeof(U) * 8 - 1);

         U bits_;
         std::memcpy(&bits_, &x, sizeof(U));
         sign_bit_ = bits_ & bits_sign_mask;
         mag_ = bits_ & ~bits_sign_mask;
      }

      void flip_sign_bit() { sign_bit_ = !sign_bit_; }

   private:
      bool sign_bit_;
      U mag_;
   };

   // The inflated bitspace representation of a floating point type `T`.
   //
   // This representation always includes denormal numbers, even if denorm
   // support is turned off. For example, a floating point type with 50
   // discrete values (zero is counted twice) between -Inf and +Inf is shown
   // below with *** marking denormal numbers.
   //
   //    -24 -20 -16 -12 -8  -4   0   4   8   12  16  20  24
   //     .   .   .   .   .   .   .   .   .   .   .   .   .
   //     |-------------------|***|***|-------------------|
   //   -Inf                      0                     +Inf
   //
   template <typename T, typename U>
   class BitsInflated : public BitsFromZero<T, U> {
   public:
      BitsInflated(const T x) : BitsFromZero<T, U>(x) {}
      BitsInflated(const bool sign, const U mag) : BitsFromZero<T, U>(sign, mag) {}
      
      T reinterpret_as_float() const {
         T f_out;
         std::memcpy(&f_out, &this->mag(), sizeof(T));
         return this->sign_bit() ? -f_out : f_out;
      }

      void divide_by_2() { this->mag() >>= 1; }
   };

   // The deflated bitspace representation of a floating point type `T`. 
   //
   // Denorm Support ON:
   //   This representation is identical to the inflated representation.
   //
   // Denorm Support OFF:
   //    This representation shifts the bitspace representation toward `0`
   //    to remove gaps in the bitspace caused by denormal numbers. For 
   //    example, consider a floating point type with 50 discrete values
   //    (zero is counted twice) between -Inf and +Inf is shown below with
   //    *** marking denormal numbers shown below.
   // 
   //                  Inflated representation:
   //    -24 -20 -16 -12 -8  -4   0   4   8   12  16  20  24
   //     .   .   .   .   .   .   .   .   .   .   .   .   .
   //     |-------------------|***|***|-------------------|
   //   -Inf                      0                     +Inf
   // _________________________________________________________________
   //
   //                  Deflated representation:
   //    -24 -20 -16 -12 -8  -4   0   4   8   12  16  20  24
   //     .   .   .   .   .   .   .   .   .   .   .   .   .
   //        |-------------------|||-------------------|
   //      -Inf                   0                  +Inf
   //
   template <typename T, typename U>
   class BitsDeflated : public BitsFromZero<T, U> {
   public:
      BitsDeflated(const bool sign, const U mag) : BitsFromZero<T, U>(sign, mag) {}

      T static_cast_int_value_to_float() const {
         T mag_float = static_cast<T>(this->mag());
         return this->sign_bit() ? -mag_float : mag_float;
      }

      // Move over `n` in bitspace
      BitsDeflated<T,U> operator+(int n) const {
         return BitsDeflated<T,U>(this->sign_bit(), this->mag() + n);
      }

      BitsDeflated<T,U> operator-(const BitsDeflated<T,U>& y) const {
         auto y_copy = y;
         y_copy.flip_sign_bit();
         return *this + y_copy;
      }

      BitsDeflated<T,U> operator+(const BitsDeflated<T,U>& y) const {
         // Gaurantee that y has the larger magnitude
         if (y.mag() < this->mag()) { return y + *this; }

         // Call *this x
         const BitsDeflated<T, U>& x = *this;

         const U mag_x = x.mag();
         const U mag_y = y.mag();

         // Calculate the deflated sum in bits (always positive)
         U bits_sum = (x.sign_bit() == y.sign_bit()) ? (mag_y + mag_x) 
                                                     : (mag_y - mag_x);

         // Sum always has the sign of the bigger magnitude (y)
         return BitsDeflated<T,U>(y.sign_bit(), bits_sum);
      }

      BitsDeflated<T,U>& operator>>(const int n) {
         this->mag() >>= n;
         return *this;
      }
   };

   // The `Denormflator` converts between inflated and deflated bitspace representations
   // to compensate for gaps in the bitspace caused by denormal numbers. The `deflate`
   // operation shifts the bitspace representation toward `0`. The `inflate` operation
   // shifts the bitspace representation away from `0`. These operations only have 
   // an effect if denorm support is turned off. If denorm support is turned on, then
   // the shift is zero. The effect of these operations is illustrated in various
   // cartoons below.
   //
   template <typename T, typename U>
   class Denormflator : public Layout_IEEE754Linear<T, U> {
   public:
      Denormflator() : has_denorm_(calc_has_denorm()) {}

      BitsDeflated<T, U> deflate(const BitsInflated<T, U>& bit_view) const {
         const U penalty = bits_denorm_shift();
         const U mag_un = bit_view.mag();
         const U mag_sh = (has_denorm_ | (mag_un < penalty)) ? mag_un : mag_un - penalty;
         return BitsDeflated<T, U>(bit_view.sign_bit(), mag_sh);
      }

      BitsInflated<T,U> inflate(const BitsDeflated<T,U>& b) const {
         const U mag = b.mag();
         const U mag_inflated = (has_denorm_ | (mag == 0)) ? mag : mag + bits_denorm_shift();
         return BitsInflated<T, U>(b.sign_bit(), mag_inflated);
      }

   private:
      // Denorm penalty
      static constexpr U bits_denorm_shift() { return (U(1) << (std::numeric_limits<T>::digits - 1)) - 1; }
      
      static bool calc_has_denorm() {
         return boost::math::detail::get_smallest_value<T>() != (std::numeric_limits<T>::min)();
      }

      bool has_denorm_;
   };

   // Check if T has a member function named "value".
   template <typename T>
   class has_value_member
   {
      template <typename U>
      static auto test(int) -> decltype(std::declval<U>().value(), std::true_type());

      template <typename U>
      static std::false_type test(...);

   public:
      static constexpr bool value = decltype(test<T>(0))::value;
   };

   // Allows one to static_cast from `boost::math::concepts::std_real_concept` to
   // type `__float`
   class StaticCast {
   public:
      template <typename T, typename V>
      static typename std::enable_if<has_value_member<V>::value, T>::type value(V x) {
         return static_cast<T>(x.value());
      }

      template <typename T, typename V>
      static typename std::enable_if<!has_value_member<V>::value, T>::type value(V x) {
         return static_cast<T>(x);
      }
   };

}  // namespace ieee754_linear
}  // namespace detail
}  // namespace tools
}  // namespace math
}  // namespace boost   

#endif  // BOOST_MATH_TOOLS_IEEE754_LINEAR_HPP
