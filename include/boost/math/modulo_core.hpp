//  Boost math/modulo_core.hpp header file  ----------------------------------//

//  Copyright 2002 Daryle Walker.  Use, modification, and distribution are
//  subject to the Boost Software License, Version 1.0.  (See accompanying file
//  LICENSE_1_0.txt or a copy at <http://www.boost.org/LICENSE_1_0.txt>.)

//  See <http://www.boost.org/libs/math/> for the library's home page.

#ifndef BOOST_MATH_MODULO_CORE_HPP
#define BOOST_MATH_MODULO_CORE_HPP

#include <algorithm>   // for std::swap, std::min, std::max
#include <climits>     // for INT_MAX
#include <stdexcept>   // for std::domain_error
#include <utility>     // for std::pair

#include <boost/config.hpp>         // for BOOST_STATIC_CONSTANT, etc.
#include <boost/operators.hpp>      // for boost::equality_comparable, etc.
#include <boost/static_assert.hpp>  // for BOOST_STATIC_ASSERT

#include <boost/math/common_factor.hpp>  // for boost::math::static_lcm, etc.


namespace boost
{
namespace math
{


//  Forward declarations  ----------------------------------------------------//

template < int Modulus >
    class modulo;

template < int Modulus >
    typename modulo<Modulus>::value_type  residue( modulo<Modulus> const &x );

template < int Modulus >
    modulo<Modulus>  pow( modulo<Modulus> x, unsigned e );
template < int Modulus >
    modulo<Modulus>  pow( modulo<Modulus> x, int e );

template < int Modulus >
    void  swap( modulo<Modulus> &a, modulo<Modulus> &b );

template < int Modulus1, int Modulus2 >
    modulo<Modulus1 * Modulus2>
      chinese_remainder( modulo<Modulus1> x, modulo<Modulus2> y );


//  Modular arithmetic class declaration  ------------------------------------//

template < int Modulus >
class modulo
    : equality_comparable2< modulo<Modulus>, int
    , bitwise1< modulo<Modulus>
    , totally_ordered1< modulo<Modulus>
    , unit_steppable< modulo<Modulus>
    , arithmetic1< modulo<Modulus> > > > > >
{
    BOOST_STATIC_ASSERT( Modulus >= 2 );

    typedef modulo<Modulus>  self_type;

public:
    // Types
    typedef int                      value_type;
    typedef value_type self_type::*  safe_bool;

    // Template parameters
    BOOST_STATIC_CONSTANT( value_type, modulus = Modulus );

    // Constructors (use automatic destructor and copy constructor)
    explicit  modulo( value_type v = 0 );

    #ifndef BOOST_NO_MEMBER_TEMPLATES
    template < value_type Modulus2 >
    modulo( modulo<Modulus2> const &o )
        : r_( residue(o) % self_type::modulus )
    {
        // The destination modulus has to be a factor of the source modulus.
        // This gives a many-to-one mapping.  The regular copy constructor
        // (identical source and destination moduli) is an one-to-one mapping.
        // The one-to-many case (destination modulus a multiple of the source
        // modulus) and the many-to-many case (any other two moduli) are
        // ambiguous and don't have constructors.
        BOOST_STATIC_ASSERT( (Modulus2 % self_type::modulus) == 0 );
    }
    #endif

    // Assigning operations
    void  swap( self_type &o );
    void  assign( value_type v );

    // Status operations
    bool  is_invertible() const;

    // Self-reversal operations
    void  negate_self();
    void  reciprocate_self();
    void  complement_self();
    void  not_self();

    // Equivalence comparisons with real numbers
    bool  operator ==( value_type v ) const;

    // Boolean logic operators (use automatic && and ||)
    bool  operator !() const;
          operator safe_bool() const;

    // Multi-valued logic operators
    self_type  operator ~() const;

    self_type &  operator &=( self_type const &o );
    self_type &  operator |=( self_type const &o );
    self_type &  operator ^=( self_type const &o );

    bool  operator <( self_type const &o ) const;

    // Operators for multi-valued logic and modulo arithmetic
    bool  operator ==( self_type const &o ) const;

    self_type &  operator ++();
    self_type &  operator --();

    self_type  operator +() const;

    // Modulo arithmetic operators
    self_type  operator -() const;

    self_type &  operator +=( self_type const &o );
    self_type &  operator -=( self_type const &o );
    self_type &  operator *=( self_type const &o );
    self_type &  operator /=( self_type const &o );

private:
    // Friends
    friend  value_type residue<Modulus>( self_type const & );

    // Data members
    value_type  r_;

};  // boost::math::modulo


//  Modular arithmetic class function definitions  ---------------------------//

#ifndef BOOST_NO_INCLASS_MEMBER_INITIALIZATION
template < int Modulus >
typename modulo<Modulus>::value_type const
  modulo<Modulus>::modulus;
#endif

template < int Modulus >
inline
modulo<Modulus>::modulo
(
    value_type  v  // = 0
)
    : r_( (v < 0) ? (self_type::modulus
     - static_cast<value_type>( static_cast<unsigned>(-v)
     % static_cast<unsigned>(self_type::modulus) )) : (v % self_type::modulus) )
{
}

template < int Modulus >
inline
void
modulo<Modulus>::swap
(
    self_type &  o
)
{
    using std::swap;

    swap( this->r_, o.r_ );
}

template < int Modulus >
void
modulo<Modulus>::assign
(
    value_type  v
)
{
    if ( v >= 0 )
    {
        this->r_ = v % self_type::modulus;
    }
    else
    {
        unsigned  vau = static_cast<unsigned>( -v );

        vau %= static_cast<unsigned>( self_type::modulus );
        this->r_ = self_type::modulus - static_cast<value_type>( vau );
    }
}

template < int Modulus >
inline
bool
modulo<Modulus>::is_invertible
(
) const
{
    return gcd( this->r_, self_type::modulus ) == 1;
}

template < int Modulus >
inline
void
modulo<Modulus>::negate_self
(
)
{
    if ( this->r_ )
    {
        this->r_ = self_type::modulus - this->r_;
    }
}

template < int Modulus >
void
modulo<Modulus>::reciprocate_self
(
)
{
    // We'll do the extended GCD algorithm, keeping only the parts we need,
    // and using this modulo type to substitute for negative value handling.
    typedef std::pair<value_type, self_type>  mgcd_t;

    mgcd_t  r2( self_type::modulus, self_type(0) );
    mgcd_t  r1( this->r_, self_type(1) );

    while ( r1.first )
    {
        self_type const  q( r2.first / r1.first );
        mgcd_t const     r0( r2.first % r1.first, r2.second - q * r1.second );

        r2 = r1;
        r1 = r0;
    }

    if ( r2.first != 1 )
    {
        throw std::domain_error( "noninvertible argument" );
    }

    this->r_ = r2.second.r_;
}

template < int Modulus >
inline
void
modulo<Modulus>::complement_self
(
)
{
    // Use NOT( this ) == MAX_VAL - |this| == MODULUS - 1 - |this|
    this->negate_self();
    this->operator --();
}

template < int Modulus >
inline
void
modulo<Modulus>::not_self
(
)
{
    // C-like Boolean NOT: nonzero -> 0 and 0 -> 1
    this->r_ = !this->r_;
}

template < int Modulus >
inline
bool
modulo<Modulus>::operator ==
(
    value_type  v
) const
{
    self_type const  vm( v );

    return this->r_ == vm.r_;
}

template < int Modulus >
inline
bool
modulo<Modulus>::operator !
(
) const
{
    return !this->operator safe_bool();
}

// I couldn't get this to compile with the "safe_bool" typedef.  Help!
template < int Modulus >
inline
modulo<Modulus>::operator int modulo<Modulus>::*
(
) const
{
    return this->r_ ? &self_type::r_ : 0;
}

template < int Modulus >
inline
modulo<Modulus>
modulo<Modulus>::operator ~
(
) const
{
    self_type  temp = *this;

    temp.complement_self();
    return temp;
}

template < int Modulus >
inline
modulo<Modulus> &
modulo<Modulus>::operator &=
(
    modulo<Modulus> const &  o
)
{
    using std::min;

    // Carry out: this AND other == MIN( |this|, |other| )
    this->r_ = min( this->r_, o.r_ );
    return *this;
}

template < int Modulus >
inline
modulo<Modulus> &
modulo<Modulus>::operator |=
(
    modulo<Modulus> const &  o
)
{
    using std::max;

    // Carry out: this OR other == MAX( |this|, |other| )
    this->r_ = max( this->r_, o.r_ );
    return *this;
}

template < int Modulus >
inline
modulo<Modulus> &
modulo<Modulus>::operator ^=
(
    modulo<Modulus> const &  o
)
{
    // Carry out: this XOR other
    // == ( this AND NOT(other) ) OR ( NOT(this) AND other )
    self_type  not_this = this->operator ~();

    this->operator &=( ~o );
    not_this.operator &=( o );
    this->operator |=( not_this );
    return *this;
}

template < int Modulus >
inline
bool
modulo<Modulus>::operator <
(
    modulo<Modulus> const &  o
) const
{
    return this->r_ < o.r_;
}

template < int Modulus >
inline
bool
modulo<Modulus>::operator ==
(
    self_type const &  o
) const
{
    return this->r_ == o.r_;
}

template < int Modulus >
modulo<Modulus> &
modulo<Modulus>::operator ++
(
)
{
    if ( ++this->r_ >= self_type::modulus )
    {
        // Wrap around if residue was (Modulus - 1)
        // i.e. r_ jumps from Modulus to 0
        this->r_ = 0;
    }
    return *this;
}

template < int Modulus >
modulo<Modulus> &
modulo<Modulus>::operator --
(
)
{
    if ( !this->r_-- )
    {
        // Wrap around if residue was zero
        // i.e. r_ jumps from -1 to (Modulus - 1)
        this->r_ += self_type::modulus;
    }
    return *this;
}

template < int Modulus >
inline
modulo<Modulus>
modulo<Modulus>::operator +
(
) const
{
    return *this;
}

template < int Modulus >
inline
modulo<Modulus>
modulo<Modulus>::operator -
(
) const
{
    self_type  temp = *this;

    temp.negate_self();
    return temp;
}

template < int Modulus >
modulo<Modulus> &
modulo<Modulus>::operator +=
(
    modulo<Modulus> const &  o
)
{
    value_type const  limit = self_type::modulus - this->r_;

    if ( o.r_ < limit )
    {
        this->r_ += o.r_;
    }
    else
    {
        // Simulate the wrap-around
        this->r_ = o.r_ - limit;
    }
    return *this;
}

template < int Modulus >
inline
modulo<Modulus> &
modulo<Modulus>::operator -=
(
    modulo<Modulus> const &  o
)
{
    self_type  temp = o;

    temp.negate_self();
    return this->operator +=( temp );
}

template < int Modulus >
modulo<Modulus> &
modulo<Modulus>::operator *=
(
    modulo<Modulus> const &  o
)
{
    // Only nonzero values are changed by multiplication
    if ( this->r_ )
    {
        value_type  trial_q = o.r_, trial_r = 0;

        // We could just multiply the residues together then reduce
        // that product by the modulus, put this piecewise method
        // works even if that product is too large to be represented!
        while ( trial_q )
        {
            value_type const  diff = self_type::modulus - trial_r;
            value_type const  diff_q = diff / this->r_;
            value_type const  diff_r = diff % this->r_;

            if ( trial_q <= diff_q )
            {
                break;
            }
            else
            {
                trial_q -= diff_q;
                if ( diff_r )
                {
                    --trial_q;
                }

                trial_r = this->r_ - diff_r;
            }
        }

        this->assign( trial_r + trial_q * this->r_ );
    }

    return *this;
}

template < int Modulus >
inline
modulo<Modulus> &
modulo<Modulus>::operator /=
(
    modulo<Modulus> const &  o
)
{
    self_type  temp = o;

    temp.reciprocate_self();
    return this->operator *=( temp );
}


//  Modular arithmetic operating function definitions  -----------------------//

template < int Modulus >
inline
typename modulo<Modulus>::value_type
residue
(
    modulo<Modulus> const &  x
)
{
    return x.r_;
}

template < int Modulus >
modulo<Modulus>
pow
(
    modulo<Modulus>  x,
    unsigned         e
)
{
    if ( !x && !e )
    {
        throw std::domain_error( "zeroth power of the zero congruence class"
         " is undefined" );
    }

    modulo<Modulus>  temp( 1 );

    for ( ; e ; e /= 2u, x *= x )
    {
        if ( e % 2u != 0 )
        {
            temp *= x;
        }
    }

    return temp;
}

template < int Modulus >
modulo<Modulus>
pow
(
    modulo<Modulus>  x,
    int              e
)
{
    unsigned  en;

    if ( e < 0 )
    {
        x.reciprocate_self();
        en = static_cast<unsigned>( -e );
    }
    else
    {
        en = static_cast<unsigned>( e );
    }

    return pow( x, en );
}

template < int Modulus >
inline
void
swap
(
    modulo<Modulus> &  a,
    modulo<Modulus> &  b
)
{
    a.swap( b );
}


//  Chinese remainder theorem solver function  -------------------------------//

template < int Modulus1, int Modulus2 >
modulo<Modulus1 * Modulus2>
chinese_remainder
(
    modulo<Modulus1>  x,
    modulo<Modulus2>  y
)
{
    // Avoid overflow
    BOOST_STATIC_ASSERT( (INT_MAX / Modulus1) >= Modulus2 );

    // Make sure the moduli are coprime
    BOOST_STATIC_ASSERT( (Modulus1 * Modulus2)
     == (static_lcm<Modulus1, Modulus2>::value) );

    // Types
    typedef modulo<Modulus1 * Modulus2>  result_type;

    // The algorithm
    result_type  xx( residue(x) ), yy( residue(y) );

    xx *= result_type( Modulus2 );
    yy *= result_type( Modulus1 );

    modulo<Modulus1>  m2_as_m1( Modulus2 );
    modulo<Modulus2>  m1_as_m2( Modulus1 );

    m2_as_m1.reciprocate_self();
    m1_as_m2.reciprocate_self();

    xx *= result_type( residue(m2_as_m1) );
    yy *= result_type( residue(m1_as_m2) );

    xx += yy;
    return xx;
}


}  // namespace math
}  // namespace boost


#endif  // BOOST_MATH_MODULO_CORE_HPP
