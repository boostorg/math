//  Copyright (c) 2024 Matt Borland
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
//  Regular use of <type_traits> is not compatible with CUDA
//  Adds aliases to unify the support
//  Also adds convience overloads like is_same_v so we don't have to wait for C++17

#ifndef BOOST_MATH_TOOLS_TYPE_TRAITS
#define BOOST_MATH_TOOLS_TYPE_TRAITS

#include <boost/math/tools/config.hpp>

#ifdef BOOST_MATH_ENABLE_CUDA

#include <cuda/std/type_traits>

namespace boost {
namespace math {

// Helper classes
BOOST_MATH_EXPORT using cuda::std::integral_constant;
BOOST_MATH_EXPORT using cuda::std::true_type;
BOOST_MATH_EXPORT using cuda::std::false_type;

// Primary type categories
BOOST_MATH_EXPORT using cuda::std::is_void;
BOOST_MATH_EXPORT using cuda::std::is_null_pointer;
BOOST_MATH_EXPORT using cuda::std::is_integral;
BOOST_MATH_EXPORT using cuda::std::is_floating_point;
BOOST_MATH_EXPORT using cuda::std::is_array;
BOOST_MATH_EXPORT using cuda::std::is_enum;
BOOST_MATH_EXPORT using cuda::std::is_union;
BOOST_MATH_EXPORT using cuda::std::is_class;
BOOST_MATH_EXPORT using cuda::std::is_function;
BOOST_MATH_EXPORT using cuda::std::is_pointer;
BOOST_MATH_EXPORT using cuda::std::is_lvalue_reference;
BOOST_MATH_EXPORT using cuda::std::is_rvalue_reference;
BOOST_MATH_EXPORT using cuda::std::is_member_object_pointer;
BOOST_MATH_EXPORT using cuda::std::is_member_function_pointer;

// Composite Type Categories
BOOST_MATH_EXPORT using cuda::std::is_fundamental;
BOOST_MATH_EXPORT using cuda::std::is_arithmetic;
BOOST_MATH_EXPORT using cuda::std::is_scalar;
BOOST_MATH_EXPORT using cuda::std::is_object;
BOOST_MATH_EXPORT using cuda::std::is_compound;
BOOST_MATH_EXPORT using cuda::std::is_reference;
BOOST_MATH_EXPORT using cuda::std::is_member_pointer;

// Type properties
BOOST_MATH_EXPORT using cuda::std::is_const;
BOOST_MATH_EXPORT using cuda::std::is_volatile;
BOOST_MATH_EXPORT using cuda::std::is_trivially_copyable;
BOOST_MATH_EXPORT using cuda::std::is_standard_layout;
BOOST_MATH_EXPORT using cuda::std::is_empty;
BOOST_MATH_EXPORT using cuda::std::is_polymorphic;
BOOST_MATH_EXPORT using cuda::std::is_abstract;
BOOST_MATH_EXPORT using cuda::std::is_final;
BOOST_MATH_EXPORT using cuda::std::is_signed;
BOOST_MATH_EXPORT using cuda::std::is_unsigned; 

// Supported Operations
BOOST_MATH_EXPORT using cuda::std::is_constructible;
BOOST_MATH_EXPORT using cuda::std::is_trivially_constructible;
BOOST_MATH_EXPORT using cuda::std::is_nothrow_constructible;

BOOST_MATH_EXPORT using cuda::std::is_default_constructible;
BOOST_MATH_EXPORT using cuda::std::is_trivially_default_constructible;
BOOST_MATH_EXPORT using cuda::std::is_nothrow_default_constructible;

BOOST_MATH_EXPORT using cuda::std::is_copy_constructible;
BOOST_MATH_EXPORT using cuda::std::is_trivially_copy_constructible;
BOOST_MATH_EXPORT using cuda::std::is_nothrow_copy_constructible;

BOOST_MATH_EXPORT using cuda::std::is_move_constructible;
BOOST_MATH_EXPORT using cuda::std::is_trivially_move_constructible;
BOOST_MATH_EXPORT using cuda::std::is_nothrow_move_constructible;

BOOST_MATH_EXPORT using cuda::std::is_assignable;
BOOST_MATH_EXPORT using cuda::std::is_trivially_assignable;
BOOST_MATH_EXPORT using cuda::std::is_nothrow_assignable;

BOOST_MATH_EXPORT using cuda::std::is_copy_assignable;
BOOST_MATH_EXPORT using cuda::std::is_trivially_copy_assignable;
BOOST_MATH_EXPORT using cuda::std::is_nothrow_copy_assignable;

BOOST_MATH_EXPORT using cuda::std::is_move_assignable;
BOOST_MATH_EXPORT using cuda::std::is_trivially_move_assignable;
BOOST_MATH_EXPORT using cuda::std::is_nothrow_move_assignable;

BOOST_MATH_EXPORT using cuda::std::is_destructible;
BOOST_MATH_EXPORT using cuda::std::is_trivially_destructible;
BOOST_MATH_EXPORT using cuda::std::is_nothrow_destructible;

BOOST_MATH_EXPORT using cuda::std::has_virtual_destructor;

// Property Queries
BOOST_MATH_EXPORT using cuda::std::alignment_of;
BOOST_MATH_EXPORT using cuda::std::rank;
BOOST_MATH_EXPORT using cuda::std::extent;

// Type Relationships
BOOST_MATH_EXPORT using cuda::std::is_same;
BOOST_MATH_EXPORT using cuda::std::is_base_of;
BOOST_MATH_EXPORT using cuda::std::is_convertible;

// Const-volatility specifiers
BOOST_MATH_EXPORT using cuda::std::remove_cv;
BOOST_MATH_EXPORT using cuda::std::remove_cv_t;
BOOST_MATH_EXPORT using cuda::std::remove_const;
BOOST_MATH_EXPORT using cuda::std::remove_const_t;
BOOST_MATH_EXPORT using cuda::std::remove_volatile;
BOOST_MATH_EXPORT using cuda::std::remove_volatile_t;
BOOST_MATH_EXPORT using cuda::std::add_cv;
BOOST_MATH_EXPORT using cuda::std::add_cv_t;
BOOST_MATH_EXPORT using cuda::std::add_const;
BOOST_MATH_EXPORT using cuda::std::add_const_t;
BOOST_MATH_EXPORT using cuda::std::add_volatile;
BOOST_MATH_EXPORT using cuda::std::add_volatile_t;

// References
BOOST_MATH_EXPORT using cuda::std::remove_reference;
BOOST_MATH_EXPORT using cuda::std::remove_reference_t;
BOOST_MATH_EXPORT using cuda::std::add_lvalue_reference;
BOOST_MATH_EXPORT using cuda::std::add_lvalue_reference_t;
BOOST_MATH_EXPORT using cuda::std::add_rvalue_reference;
BOOST_MATH_EXPORT using cuda::std::add_rvalue_reference_t;

// Pointers
BOOST_MATH_EXPORT using cuda::std::remove_pointer;
BOOST_MATH_EXPORT using cuda::std::remove_pointer_t;
BOOST_MATH_EXPORT using cuda::std::add_pointer;
BOOST_MATH_EXPORT using cuda::std::add_pointer_t;

// Sign Modifiers
BOOST_MATH_EXPORT using cuda::std::make_signed;
BOOST_MATH_EXPORT using cuda::std::make_signed_t;
BOOST_MATH_EXPORT using cuda::std::make_unsigned;
BOOST_MATH_EXPORT using cuda::std::make_unsigned_t;

// Arrays
BOOST_MATH_EXPORT using cuda::std::remove_extent;
BOOST_MATH_EXPORT using cuda::std::remove_extent_t;
BOOST_MATH_EXPORT using cuda::std::remove_all_extents;
BOOST_MATH_EXPORT using cuda::std::remove_all_extents_t;

// Misc transformations
BOOST_MATH_EXPORT using cuda::std::decay;
BOOST_MATH_EXPORT using cuda::std::decay_t;
BOOST_MATH_EXPORT using cuda::std::enable_if;
BOOST_MATH_EXPORT using cuda::std::enable_if_t;
BOOST_MATH_EXPORT using cuda::std::conditional;
BOOST_MATH_EXPORT using cuda::std::conditional_t;
BOOST_MATH_EXPORT using cuda::std::common_type;
BOOST_MATH_EXPORT using cuda::std::common_type_t;
BOOST_MATH_EXPORT using cuda::std::underlying_type;
BOOST_MATH_EXPORT using cuda::std::underlying_type_t;

#else // STD versions

#ifndef BOOST_MATH_BUILD_MODULE
#include <type_traits>
#endif

namespace boost {
namespace math {

// Helper classes
BOOST_MATH_EXPORT using std::integral_constant;
BOOST_MATH_EXPORT using std::true_type;
BOOST_MATH_EXPORT using std::false_type;

// Primary type categories
BOOST_MATH_EXPORT using std::is_void;
BOOST_MATH_EXPORT using std::is_null_pointer;
BOOST_MATH_EXPORT using std::is_integral;
BOOST_MATH_EXPORT using std::is_floating_point;
BOOST_MATH_EXPORT using std::is_array;
BOOST_MATH_EXPORT using std::is_enum;
BOOST_MATH_EXPORT using std::is_union;
BOOST_MATH_EXPORT using std::is_class;
BOOST_MATH_EXPORT using std::is_function;
BOOST_MATH_EXPORT using std::is_pointer;
BOOST_MATH_EXPORT using std::is_lvalue_reference;
BOOST_MATH_EXPORT using std::is_rvalue_reference;
BOOST_MATH_EXPORT using std::is_member_object_pointer;
BOOST_MATH_EXPORT using std::is_member_function_pointer;

// Composite Type Categories
BOOST_MATH_EXPORT using std::is_fundamental;
BOOST_MATH_EXPORT using std::is_arithmetic;
BOOST_MATH_EXPORT using std::is_scalar;
BOOST_MATH_EXPORT using std::is_object;
BOOST_MATH_EXPORT using std::is_compound;
BOOST_MATH_EXPORT using std::is_reference;
BOOST_MATH_EXPORT using std::is_member_pointer;

// Type properties
BOOST_MATH_EXPORT using std::is_const;
BOOST_MATH_EXPORT using std::is_volatile;
BOOST_MATH_EXPORT using std::is_trivially_copyable;
BOOST_MATH_EXPORT using std::is_standard_layout;
BOOST_MATH_EXPORT using std::is_empty;
BOOST_MATH_EXPORT using std::is_polymorphic;
BOOST_MATH_EXPORT using std::is_abstract;
BOOST_MATH_EXPORT using std::is_final;
BOOST_MATH_EXPORT using std::is_signed;
BOOST_MATH_EXPORT using std::is_unsigned; 

// Supported Operations
BOOST_MATH_EXPORT using std::is_constructible;
BOOST_MATH_EXPORT using std::is_trivially_constructible;
BOOST_MATH_EXPORT using std::is_nothrow_constructible;

BOOST_MATH_EXPORT using std::is_default_constructible;
BOOST_MATH_EXPORT using std::is_trivially_default_constructible;
BOOST_MATH_EXPORT using std::is_nothrow_default_constructible;

BOOST_MATH_EXPORT using std::is_copy_constructible;
BOOST_MATH_EXPORT using std::is_trivially_copy_constructible;
BOOST_MATH_EXPORT using std::is_nothrow_copy_constructible;

BOOST_MATH_EXPORT using std::is_move_constructible;
BOOST_MATH_EXPORT using std::is_trivially_move_constructible;
BOOST_MATH_EXPORT using std::is_nothrow_move_constructible;

BOOST_MATH_EXPORT using std::is_assignable;
BOOST_MATH_EXPORT using std::is_trivially_assignable;
BOOST_MATH_EXPORT using std::is_nothrow_assignable;

BOOST_MATH_EXPORT using std::is_copy_assignable;
BOOST_MATH_EXPORT using std::is_trivially_copy_assignable;
BOOST_MATH_EXPORT using std::is_nothrow_copy_assignable;

BOOST_MATH_EXPORT using std::is_move_assignable;
BOOST_MATH_EXPORT using std::is_trivially_move_assignable;
BOOST_MATH_EXPORT using std::is_nothrow_move_assignable;

BOOST_MATH_EXPORT using std::is_destructible;
BOOST_MATH_EXPORT using std::is_trivially_destructible;
BOOST_MATH_EXPORT using std::is_nothrow_destructible;

BOOST_MATH_EXPORT using std::has_virtual_destructor;

// Property Queries
BOOST_MATH_EXPORT using std::alignment_of;
BOOST_MATH_EXPORT using std::rank;
BOOST_MATH_EXPORT using std::extent;

// Type Relationships
BOOST_MATH_EXPORT using std::is_same;
BOOST_MATH_EXPORT using std::is_base_of;
BOOST_MATH_EXPORT using std::is_convertible;

// Const-volatility specifiers
BOOST_MATH_EXPORT using std::remove_cv;
BOOST_MATH_EXPORT using std::remove_cv_t;
BOOST_MATH_EXPORT using std::remove_const;
BOOST_MATH_EXPORT using std::remove_const_t;
BOOST_MATH_EXPORT using std::remove_volatile;
BOOST_MATH_EXPORT using std::remove_volatile_t;
BOOST_MATH_EXPORT using std::add_cv;
BOOST_MATH_EXPORT using std::add_cv_t;
BOOST_MATH_EXPORT using std::add_const;
BOOST_MATH_EXPORT using std::add_const_t;
BOOST_MATH_EXPORT using std::add_volatile;
BOOST_MATH_EXPORT using std::add_volatile_t;

// References
BOOST_MATH_EXPORT using std::remove_reference;
BOOST_MATH_EXPORT using std::remove_reference_t;
BOOST_MATH_EXPORT using std::add_lvalue_reference;
BOOST_MATH_EXPORT using std::add_lvalue_reference_t;
BOOST_MATH_EXPORT using std::add_rvalue_reference;
BOOST_MATH_EXPORT using std::add_rvalue_reference_t;

// Pointers
BOOST_MATH_EXPORT using std::remove_pointer;
BOOST_MATH_EXPORT using std::remove_pointer_t;
BOOST_MATH_EXPORT using std::add_pointer;
BOOST_MATH_EXPORT using std::add_pointer_t;

// Sign Modifiers
BOOST_MATH_EXPORT using std::make_signed;
BOOST_MATH_EXPORT using std::make_signed_t;
BOOST_MATH_EXPORT using std::make_unsigned;
BOOST_MATH_EXPORT using std::make_unsigned_t;

// Arrays
BOOST_MATH_EXPORT using std::remove_extent;
BOOST_MATH_EXPORT using std::remove_extent_t;
BOOST_MATH_EXPORT using std::remove_all_extents;
BOOST_MATH_EXPORT using std::remove_all_extents_t;

// Misc transformations
BOOST_MATH_EXPORT using std::decay;
BOOST_MATH_EXPORT using std::decay_t;
BOOST_MATH_EXPORT using std::enable_if;
BOOST_MATH_EXPORT using std::enable_if_t;
BOOST_MATH_EXPORT using std::conditional;
BOOST_MATH_EXPORT using std::conditional_t;
BOOST_MATH_EXPORT using std::common_type;
BOOST_MATH_EXPORT using std::common_type_t;
BOOST_MATH_EXPORT using std::underlying_type;
BOOST_MATH_EXPORT using std::underlying_type_t;

#endif 

BOOST_MATH_EXPORT template <bool B>
using bool_constant = boost::math::integral_constant<bool, B>;

BOOST_MATH_EXPORT template <typename T>
BOOST_MATH_INLINE_CONSTEXPR bool is_void_v = boost::math::is_void<T>::value;

BOOST_MATH_EXPORT template <typename T>
BOOST_MATH_INLINE_CONSTEXPR bool is_null_pointer_v = boost::math::is_null_pointer<T>::value;

BOOST_MATH_EXPORT template <typename T>
BOOST_MATH_INLINE_CONSTEXPR bool is_integral_v = boost::math::is_integral<T>::value;

BOOST_MATH_EXPORT template <typename T>
BOOST_MATH_INLINE_CONSTEXPR bool is_floating_point_v = boost::math::is_floating_point<T>::value;

BOOST_MATH_EXPORT template <typename T>
BOOST_MATH_INLINE_CONSTEXPR bool is_array_v = boost::math::is_array<T>::value;

BOOST_MATH_EXPORT template <typename T>
BOOST_MATH_INLINE_CONSTEXPR bool is_enum_v = boost::math::is_enum<T>::value;

BOOST_MATH_EXPORT template <typename T>
BOOST_MATH_INLINE_CONSTEXPR bool is_union_v = boost::math::is_union<T>::value;

BOOST_MATH_EXPORT template <typename T>
BOOST_MATH_INLINE_CONSTEXPR bool is_class_v = boost::math::is_class<T>::value;

BOOST_MATH_EXPORT template <typename T>
BOOST_MATH_INLINE_CONSTEXPR bool is_function_v = boost::math::is_function<T>::value;

BOOST_MATH_EXPORT template <typename T>
BOOST_MATH_INLINE_CONSTEXPR bool is_pointer_v = boost::math::is_pointer<T>::value;

BOOST_MATH_EXPORT template <typename T>
BOOST_MATH_INLINE_CONSTEXPR bool is_lvalue_reference_v = boost::math::is_lvalue_reference<T>::value;

BOOST_MATH_EXPORT template <typename T>
BOOST_MATH_INLINE_CONSTEXPR bool is_rvalue_reference_v = boost::math::is_rvalue_reference<T>::value;

BOOST_MATH_EXPORT template <typename T>
BOOST_MATH_INLINE_CONSTEXPR bool is_member_object_pointer_v = boost::math::is_member_object_pointer<T>::value;

BOOST_MATH_EXPORT template <typename T>
BOOST_MATH_INLINE_CONSTEXPR bool is_member_function_pointer_v = boost::math::is_member_function_pointer<T>::value;

BOOST_MATH_EXPORT template <typename T>
BOOST_MATH_INLINE_CONSTEXPR bool is_fundamental_v = boost::math::is_fundamental<T>::value;

BOOST_MATH_EXPORT template <typename T>
BOOST_MATH_INLINE_CONSTEXPR bool is_arithmetic_v = boost::math::is_arithmetic<T>::value;

BOOST_MATH_EXPORT template <typename T>
BOOST_MATH_INLINE_CONSTEXPR bool is_scalar_v = boost::math::is_scalar<T>::value;

BOOST_MATH_EXPORT template <typename T>
BOOST_MATH_INLINE_CONSTEXPR bool is_object_v = boost::math::is_object<T>::value;

BOOST_MATH_EXPORT template <typename T>
BOOST_MATH_INLINE_CONSTEXPR bool is_compound_v = boost::math::is_compound<T>::value;

BOOST_MATH_EXPORT template <typename T>
BOOST_MATH_INLINE_CONSTEXPR bool is_reference_v = boost::math::is_reference<T>::value;

BOOST_MATH_EXPORT template <typename T>
BOOST_MATH_INLINE_CONSTEXPR bool is_member_pointer_v = boost::math::is_member_pointer<T>::value;

BOOST_MATH_EXPORT template <typename T>
BOOST_MATH_INLINE_CONSTEXPR bool is_const_v = boost::math::is_const<T>::value;

BOOST_MATH_EXPORT template <typename T>
BOOST_MATH_INLINE_CONSTEXPR bool is_volatile_v = boost::math::is_volatile<T>::value;

BOOST_MATH_EXPORT template <typename T>
BOOST_MATH_INLINE_CONSTEXPR bool is_trivially_copyable_v = boost::math::is_trivially_copyable<T>::value;

BOOST_MATH_EXPORT template <typename T>
BOOST_MATH_INLINE_CONSTEXPR bool is_standard_layout_v = boost::math::is_standard_layout<T>::value;

BOOST_MATH_EXPORT template <typename T>
BOOST_MATH_INLINE_CONSTEXPR bool is_empty_v = boost::math::is_empty<T>::value;

BOOST_MATH_EXPORT template <typename T>
BOOST_MATH_INLINE_CONSTEXPR bool is_polymorphic_v = boost::math::is_polymorphic<T>::value;

BOOST_MATH_EXPORT template <typename T>
BOOST_MATH_INLINE_CONSTEXPR bool is_abstract_v = boost::math::is_abstract<T>::value;

BOOST_MATH_EXPORT template <typename T>
BOOST_MATH_INLINE_CONSTEXPR bool is_final_v = boost::math::is_final<T>::value;

BOOST_MATH_EXPORT template <typename T>
BOOST_MATH_INLINE_CONSTEXPR bool is_signed_v = boost::math::is_signed<T>::value;

BOOST_MATH_EXPORT template <typename T>
BOOST_MATH_INLINE_CONSTEXPR bool is_unsigned_v = boost::math::is_unsigned<T>::value;

BOOST_MATH_EXPORT template <typename T, typename... Args>
BOOST_MATH_INLINE_CONSTEXPR bool is_constructible_v = boost::math::is_constructible<T, Args...>::value;

BOOST_MATH_EXPORT template <typename T, typename... Args>
BOOST_MATH_INLINE_CONSTEXPR bool is_trivially_constructible_v = boost::math::is_trivially_constructible<T, Args...>::value;

BOOST_MATH_EXPORT template <typename T, typename... Args>
BOOST_MATH_INLINE_CONSTEXPR bool is_nothrow_constructible_v = boost::math::is_nothrow_constructible<T, Args...>::value;

BOOST_MATH_EXPORT template <typename T>
BOOST_MATH_INLINE_CONSTEXPR bool is_default_constructible_v = boost::math::is_default_constructible<T>::value;

BOOST_MATH_EXPORT template <typename T>
BOOST_MATH_INLINE_CONSTEXPR bool is_trivially_default_constructible_v = boost::math::is_trivially_default_constructible<T>::value;

BOOST_MATH_EXPORT template <typename T>
BOOST_MATH_INLINE_CONSTEXPR bool is_nothrow_default_constructible_v = boost::math::is_nothrow_default_constructible<T>::value;

BOOST_MATH_EXPORT template <typename T>
BOOST_MATH_INLINE_CONSTEXPR bool is_copy_constructible_v = boost::math::is_copy_constructible<T>::value;

BOOST_MATH_EXPORT template <typename T>
BOOST_MATH_INLINE_CONSTEXPR bool is_trivially_copy_constructible_v = boost::math::is_trivially_copy_constructible<T>::value;

BOOST_MATH_EXPORT template <typename T>
BOOST_MATH_INLINE_CONSTEXPR bool is_nothrow_copy_constructible_v = boost::math::is_nothrow_copy_constructible<T>::value;

BOOST_MATH_EXPORT template <typename T>
BOOST_MATH_INLINE_CONSTEXPR bool is_move_constructible_v = boost::math::is_move_constructible<T>::value;

BOOST_MATH_EXPORT template <typename T>
BOOST_MATH_INLINE_CONSTEXPR bool is_trivially_move_constructible_v = boost::math::is_trivially_move_constructible<T>::value;

BOOST_MATH_EXPORT template <typename T>
BOOST_MATH_INLINE_CONSTEXPR bool is_nothrow_move_constructible_v = boost::math::is_nothrow_move_constructible<T>::value;

BOOST_MATH_EXPORT template <typename T, typename U>
BOOST_MATH_INLINE_CONSTEXPR bool is_assignable_v = boost::math::is_assignable<T, U>::value;

BOOST_MATH_EXPORT template <typename T, typename U>
BOOST_MATH_INLINE_CONSTEXPR bool is_trivially_assignable_v = boost::math::is_trivially_assignable<T, U>::value;

BOOST_MATH_EXPORT template <typename T, typename U>
BOOST_MATH_INLINE_CONSTEXPR bool is_nothrow_assignable_v = boost::math::is_nothrow_assignable<T, U>::value;

BOOST_MATH_EXPORT template <typename T>
BOOST_MATH_INLINE_CONSTEXPR bool is_copy_assignable_v = boost::math::is_copy_assignable<T>::value;

BOOST_MATH_EXPORT template <typename T>
BOOST_MATH_INLINE_CONSTEXPR bool is_trivially_copy_assignable_v = boost::math::is_trivially_copy_assignable<T>::value;

BOOST_MATH_EXPORT template <typename T>
BOOST_MATH_INLINE_CONSTEXPR bool is_nothrow_copy_assignable_v = boost::math::is_nothrow_copy_assignable<T>::value;

BOOST_MATH_EXPORT template <typename T>
BOOST_MATH_INLINE_CONSTEXPR bool is_move_assignable_v = boost::math::is_move_assignable<T>::value;

BOOST_MATH_EXPORT template <typename T>
BOOST_MATH_INLINE_CONSTEXPR bool is_trivially_move_assignable_v = boost::math::is_trivially_move_assignable<T>::value;

BOOST_MATH_EXPORT template <typename T>
BOOST_MATH_INLINE_CONSTEXPR bool is_nothrow_move_assignable_v = boost::math::is_nothrow_move_assignable<T>::value;

BOOST_MATH_EXPORT template <typename T>
BOOST_MATH_INLINE_CONSTEXPR bool is_destructible_v = boost::math::is_destructible<T>::value;

BOOST_MATH_EXPORT template <typename T>
BOOST_MATH_INLINE_CONSTEXPR bool is_trivially_destructible_v = boost::math::is_trivially_destructible<T>::value;

BOOST_MATH_EXPORT template <typename T>
BOOST_MATH_INLINE_CONSTEXPR bool is_nothrow_destructible_v = boost::math::is_nothrow_destructible<T>::value;

BOOST_MATH_EXPORT template <typename T>
BOOST_MATH_INLINE_CONSTEXPR bool has_virtual_destructor_v = boost::math::has_virtual_destructor<T>::value;

BOOST_MATH_EXPORT template <typename T, typename U>
BOOST_MATH_INLINE_CONSTEXPR bool is_same_v = boost::math::is_same<T, U>::value;

BOOST_MATH_EXPORT template <typename T, typename U>
BOOST_MATH_INLINE_CONSTEXPR bool is_base_of_v = boost::math::is_base_of<T, U>::value;

} // namespace math
} // namespace boost

#endif // BOOST_MATH_TOOLS_TYPE_TRAITS
