//  Boost common_factor.hpp header file  -------------------------------------//

//  (C) Copyright Daryle Walker, Stephen Cleary, Paul Moore 2001.  Permission
//  to copy, use, modify, sell and distribute this software is granted provided
//  this copyright notice appears in all copies.  This software is provided "as
//  is" without express or implied warranty, and with no claim as to its
//  suitability for any purpose. 

//  See http://www.boost.org for updates, documentation, and revision history. 

#ifndef BOOST_MATH_COMMON_FACTOR_HPP
#define BOOST_MATH_COMMON_FACTOR_HPP

#include <boost/math_fwd.hpp>  // self include

#include <boost/config.hpp>  // for BOOST_STATIC_CONSTANT, etc.
#include <boost/limits.hpp>  // for std::numeric_limits


namespace boost
{
namespace math
{


//  Forward declarations for function templates  -----------------------------//

template < typename IntegerType >
    IntegerType  gcd( IntegerType const &a, IntegerType const &b );

template < typename IntegerType >
    IntegerType  lcm( IntegerType const &a, IntegerType const &b );


//  Greatest common divisor evaluator class declaration  ---------------------//

template < typename IntegerType >
class gcd_evaluator
{
public:
    // Types
    typedef IntegerType  result_type, first_argument_type, second_argument_type;

    // Function object interface
    result_type  operator ()( first_argument_type const &a,
     second_argument_type const &b ) const;

};  // boost::math::gcd_evaluator


//  Least common multiple evaluator class declaration  -----------------------//

template < typename IntegerType >
class lcm_evaluator
{
public:
    // Types
    typedef IntegerType  result_type, first_argument_type, second_argument_type;

    // Function object interface
    result_type  operator ()( first_argument_type const &a,
     second_argument_type const &b ) const;

};  // boost::math::lcm_evaluator


//  Implementation details  --------------------------------------------------//

namespace detail
{
#ifndef BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION
    // Build GCD with Euclid's recursive algorithm
    template < unsigned long Value1, unsigned long Value2 >
    struct static_gcd_helper_t
    {
    private:
        BOOST_STATIC_CONSTANT( unsigned long, new_value1 = Value2 );
        BOOST_STATIC_CONSTANT( unsigned long, new_value2 = Value1 % Value2 );

        #ifndef __BORLANDC__
        #define BOOST_DETAIL_GCD_HELPER_VAL(Value)  Value
        #else
        typedef static_gcd_helper_t  self_type;
        #define BOOST_DETAIL_GCD_HELPER_VAL(Value)  (self_type:: Value )
        #endif

        typedef static_gcd_helper_t< BOOST_DETAIL_GCD_HELPER_VAL(new_value1),
         BOOST_DETAIL_GCD_HELPER_VAL(new_value2) >  next_step_type;

        #undef BOOST_DETAIL_GCD_HELPER_VAL

    public:
        BOOST_STATIC_CONSTANT( unsigned long, value = next_step_type::value );
    };

    // Non-recursive case
    template < unsigned long Value1 >
    struct static_gcd_helper_t< Value1, 0UL >
    {
        BOOST_STATIC_CONSTANT( unsigned long, value = Value1 );
    };
#else
    // Use inner class template workaround from Peter Dimov
    template < unsigned long Value1 >
    struct static_gcd_helper2_t
    {
        template < unsigned long Value2 >
        struct helper
        {
            BOOST_STATIC_CONSTANT( unsigned long, value
             = static_gcd_helper2_t<Value2>::helper<Value1 % Value2>::value );
        };

        template <  >
        struct helper< 0UL >
        {
            BOOST_STATIC_CONSTANT( unsigned long, value = Value1 );
        };
    };

    // Special case
    template <  >
    struct static_gcd_helper2_t< 0UL >
    {
        template < unsigned long Value2 >
        struct helper
        {
            BOOST_STATIC_CONSTANT( unsigned long, value = Value2 );
        };
    };

    // Build the GCD from the above template(s)
    template < unsigned long Value1, unsigned long Value2 >
    struct static_gcd_helper_t
    {
        BOOST_STATIC_CONSTANT( unsigned long, value
         = static_gcd_helper2_t<Value1>::BOOST_NESTED_TEMPLATE
         helper<Value2>::value );
    };
#endif

#ifndef BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION
    // Build the LCM from the GCD
    template < unsigned long Value1, unsigned long Value2 >
    struct static_lcm_helper_t
    {
        typedef static_gcd_helper_t<Value1, Value2>  gcd_type;

        BOOST_STATIC_CONSTANT( unsigned long, value = Value1 / gcd_type::value
         * Value2 );
    };

    // Special case for zero-GCD values
    template < >
    struct static_lcm_helper_t< 0UL, 0UL >
    {
        BOOST_STATIC_CONSTANT( unsigned long, value = 0UL );
    };
#else
    // Adapt GCD's inner class template workaround for LCM
    template < unsigned long Value1 >
    struct static_lcm_helper2_t
    {
        template < unsigned long Value2 >
        struct helper
        {
            typedef static_gcd_helper_t<Value1, Value2>  gcd_type;

            BOOST_STATIC_CONSTANT( unsigned long, value = Value1
             / gcd_type::value * Value2 );
        };

        template <  >
        struct helper< 0UL >
        {
            BOOST_STATIC_CONSTANT( unsigned long, value = 0UL );
        };
    };

    // Special case
    template <  >
    struct static_lcm_helper2_t< 0UL >
    {
        template < unsigned long Value2 >
        struct helper
        {
            BOOST_STATIC_CONSTANT( unsigned long, value = 0UL );
        };
    };

    // Build the LCM from the above template(s)
    template < unsigned long Value1, unsigned long Value2 >
    struct static_lcm_helper_t
    {
        BOOST_STATIC_CONSTANT( unsigned long, value
         = static_lcm_helper2_t<Value1>::BOOST_NESTED_TEMPLATE
         helper<Value2>::value );
    };
#endif

    // Greatest common divisor for rings (including unsigned integers)
    template < typename RingType >
    RingType
    gcd_euclidean
    (
        RingType  a,
        RingType  b
    )
    {
        // Avoid repeated construction
        #ifndef __BORLANDC__
        RingType const  zero = static_cast<RingType>( 0 );
        #else
        RingType  zero = static_cast<RingType>( 0 );
        #endif

        // Reduce by GCD-remainder property [GCD(a,b) == GCD(b,a MOD b)]
        while ( true )
        {
            if ( a == zero )
                return b;
            b %= a;

            if ( b == zero )
                return a;
            a %= b;
        }
    }

    // Greatest common divisor for (signed) integers
    template < typename IntegerType >
    inline
    IntegerType
    gcd_integer
    (
        IntegerType const &  a,
        IntegerType const &  b
    )
    {
        // Avoid repeated construction
        IntegerType const  zero = static_cast<IntegerType>( 0 );
        IntegerType const  result = gcd_euclidean( a, b );

        return ( result < zero ) ? -result : result;
    }

    // Least common multiple for rings (including unsigned integers)
    template < typename RingType >
    inline
    RingType
    lcm_euclidean
    (
        RingType const &  a,
        RingType const &  b
    )
    {
        RingType const  zero = static_cast<RingType>( 0 );
        RingType const  temp = gcd_euclidean( a, b );

        return ( temp != zero ) ? ( a / temp * b ) : zero;
    }

    // Least common multiple for (signed) integers
    template < typename IntegerType >
    inline
    IntegerType
    lcm_integer
    (
        IntegerType const &  a,
        IntegerType const &  b
    )
    {
        // Avoid repeated construction
        IntegerType const  zero = static_cast<IntegerType>( 0 );
        IntegerType const  result = lcm_euclidean( a, b );

        return ( result < zero ) ? -result : result;
    }

    // Function objects to find the best way of computing GCD or LCM
#ifndef BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION
    template < typename T, bool IsSpecialized, bool IsSigned >
    struct gcd_optimal_evaluator_helper_t
    {
        T  operator ()( T const &a, T const &b )
        {
            return gcd_euclidean( a, b );
        }
    };

    template < typename T >
    struct gcd_optimal_evaluator_helper_t< T, true, true >
    {
        T  operator ()( T const &a, T const &b )
        {
            return gcd_integer( a, b );
        }
    };
#else
    template < bool IsSpecialized, bool IsSigned >
    struct gcd_optimal_evaluator_helper2_t
    {
        template < typename T >
        struct helper
        {
            T  operator ()( T const &a, T const &b )
            {
                return gcd_euclidean( a, b );
            }
        };
    };

    template < >
    struct gcd_optimal_evaluator_helper2_t< true, true >
    {
        template < typename T >
        struct helper
        {
            T  operator ()( T const &a, T const &b )
            {
                return gcd_integer( a, b );
            }
        };
    };

    template < typename T, bool IsSpecialized, bool IsSigned >
    struct gcd_optimal_evaluator_helper_t
        : gcd_optimal_evaluator_helper2_t<IsSpecialized, IsSigned>
           ::BOOST_NESTED_TEMPLATE helper<T>
    {
    };
#endif

    template < typename T >
    struct gcd_optimal_evaluator
    {
        T  operator ()( T const &a, T const &b )
        {
            typedef ::std::numeric_limits<T>  limits_type;

            typedef gcd_optimal_evaluator_helper_t<T,
             limits_type::is_specialized, limits_type::is_signed>  helper_type;

            helper_type  solver;

            return solver( a, b );
        }
    };

#ifndef BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION
    template < typename T, bool IsSpecialized, bool IsSigned >
    struct lcm_optimal_evaluator_helper_t
    {
        T  operator ()( T const &a, T const &b )
        {
            return lcm_euclidean( a, b );
        }
    };

    template < typename T >
    struct lcm_optimal_evaluator_helper_t< T, true, true >
    {
        T  operator ()( T const &a, T const &b )
        {
            return lcm_integer( a, b );
        }
    };
#else
    template < bool IsSpecialized, bool IsSigned >
    struct lcm_optimal_evaluator_helper2_t
    {
        template < typename T >
        struct helper
        {
            T  operator ()( T const &a, T const &b )
            {
                return lcm_euclidean( a, b );
            }
        };
    };

    template < >
    struct lcm_optimal_evaluator_helper2_t< true, true >
    {
        template < typename T >
        struct helper
        {
            T  operator ()( T const &a, T const &b )
            {
                return lcm_integer( a, b );
            }
        };
    };

    template < typename T, bool IsSpecialized, bool IsSigned >
    struct lcm_optimal_evaluator_helper_t
        : lcm_optimal_evaluator_helper2_t<IsSpecialized, IsSigned>
           ::BOOST_NESTED_TEMPLATE helper<T>
    {
    };
#endif

    template < typename T >
    struct lcm_optimal_evaluator
    {
        T  operator ()( T const &a, T const &b )
        {
            typedef ::std::numeric_limits<T>  limits_type;

            typedef lcm_optimal_evaluator_helper_t<T,
             limits_type::is_specialized, limits_type::is_signed>  helper_type;

            helper_type  solver;

            return solver( a, b );
        }
    };

    // Functions to find the GCD or LCM in the best way
    template < typename T >
    inline
    T
    gcd_optimal
    (
        T const &  a,
        T const &  b
    )
    {
        gcd_optimal_evaluator<T>  solver;

        return solver( a, b );
    }

    template < typename T >
    inline
    T
    lcm_optimal
    (
        T const &  a,
        T const &  b
    )
    {
        lcm_optimal_evaluator<T>  solver;

        return solver( a, b );
    }

}  // namespace detail


//  Compile-time greatest common divisor evaluator class declaration  --------//

template < unsigned long Value1, unsigned long Value2 >
struct static_gcd
{
    BOOST_STATIC_CONSTANT( unsigned long, value
     = (detail::static_gcd_helper_t<Value1, Value2>::value) );

};  // boost::math::static_gcd


//  Compile-time least common multiple evaluator class declaration  ----------//

template < unsigned long Value1, unsigned long Value2 >
struct static_lcm
{
    BOOST_STATIC_CONSTANT( unsigned long, value
     = (detail::static_lcm_helper_t<Value1, Value2>::value) );

};  // boost::math::static_lcm


//  Greatest common divisor evaluator member function definition  ------------//

template < typename IntegerType >
inline
typename gcd_evaluator<IntegerType>::result_type
gcd_evaluator<IntegerType>::operator ()
(
    first_argument_type const &   a,
    second_argument_type const &  b
) const
{
    return detail::gcd_optimal( a, b );
}


//  Least common multiple evaluator member function definition  --------------//

template < typename IntegerType >
inline
typename lcm_evaluator<IntegerType>::result_type
lcm_evaluator<IntegerType>::operator ()
(
    first_argument_type const &   a,
    second_argument_type const &  b
) const
{
    return detail::lcm_optimal( a, b );
}


//  Greatest common divisor and least common multiple function definitions  --//

template < typename IntegerType >
inline
IntegerType
gcd
(
    IntegerType const &  a,
    IntegerType const &  b
)
{
    gcd_evaluator<IntegerType>  solver;

    return solver( a, b );
}

template < typename IntegerType >
inline
IntegerType
lcm
(
    IntegerType const &  a,
    IntegerType const &  b
)
{
    lcm_evaluator<IntegerType>  solver;

    return solver( a, b );
}


}  // namespace math
}  // namespace boost


#endif  // BOOST_MATH_COMMON_FACTOR_HPP
