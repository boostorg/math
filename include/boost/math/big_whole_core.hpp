//  Boost math/big_whole_core.hpp header file  -------------------------------//

//  Copyright 2004 Daryle Walker.  Use, modification, and distribution are
//  subject to the Boost Software License, Version 1.0.  (See accompanying file
//  LICENSE_1_0.txt or a copy at <http://www.boost.org/LICENSE_1_0.txt>.)

//  See <http://www.boost.org/libs/math/> for the library's home page.

#ifndef BOOST_MATH_BIG_WHOLE_CORE_HPP
#define BOOST_MATH_BIG_WHOLE_CORE_HPP

#include <boost/math_fwd.hpp>  // self include

#include <boost/cstdint.hpp>  // for boost::uintmax_t

#include <algorithm>  // for std::swap, std::min, std::max
#include <cstddef>    // for NULL, std::size_t
#include <limits>     // for std::numeric_limits
#include <memory>     // for std::auto_ptr
#include <valarray>   // for std::valarray, std::slice


namespace boost
{
namespace math
{


//  Class/template/function/operator forward declarations  -------------------//

class big_whole;

void  swap( big_whole &a, big_whole &b );

big_whole  operator !( big_whole const &w );

bool  operator ==( big_whole const &lhs, big_whole const &rhs );
bool  operator !=( big_whole const &lhs, big_whole const &rhs );
bool  operator < ( big_whole const &lhs, big_whole const &rhs );
bool  operator > ( big_whole const &lhs, big_whole const &rhs );
bool  operator <=( big_whole const &lhs, big_whole const &rhs );
bool  operator >=( big_whole const &lhs, big_whole const &rhs );

big_whole  operator &( big_whole const &lhs, big_whole const &rhs );
big_whole  operator |( big_whole const &lhs, big_whole const &rhs );
big_whole  operator ^( big_whole const &lhs, big_whole const &rhs );

big_whole &  operator &=( big_whole &lhs, big_whole const &rhs );
big_whole &  operator |=( big_whole &lhs, big_whole const &rhs );
big_whole &  operator ^=( big_whole &lhs, big_whole const &rhs );


//  Arbitrary-length whole-number object class declaration  ------------------//

class big_whole
{
    struct dummy { dummy *d; };

    typedef big_whole                  self_type;
    typedef dummy * dummy::*           bool_type;
    typedef std::size_t                size_type;
    typedef std::valarray<bool>        mask_type;
    typedef std::valarray<size_type>  index_type;

public:
    // Lifetime management
    big_whole();
    big_whole( self_type const &other );

    big_whole( uintmax_t v );

    explicit  big_whole( mask_type const &b );
    explicit  big_whole( index_type const &i );

    // Object-mutating operations
    void  swap( self_type &other );

    void  assign( self_type const &other );
    void  assign( uintmax_t v );

    void  reconfigure( mask_type const &b );
    void  reconfigure( index_type const &i );

    // Value-accessing operations
    uintmax_t   to_uintmax() const;
    mask_type   to_bit_vector() const;
    index_type  to_bit_indices() const;

    bool  is_even() const;

    // Bit-twiddling operations
    void  reset();
    void  reset( size_type from, size_type to );
    void  reset( size_type i );

    void  set( size_type from, size_type to );
    void  set( size_type i );

    void  flip( size_type from, size_type to );
    void  flip( size_type i );

    void  bit_assign( size_type from, size_type to, bool value );
    void  bit_assign( size_type i, bool value );

    void  bits_assign( size_type from, size_type to, self_type const &values );

    // Bit-inspecting operations
    size_type  length() const;

    size_type  count() const;
    bool       any() const;
    bool       none() const;

    bool       test( size_type i ) const;
    self_type  tests( size_type from, size_type to ) const;

    self_type  reverse( size_type cap ) const;
    self_type  reverse() const;

    // Self-operator mutators
    void  not_self();

    // Object-accessing operations
    int  compare( self_type const &other ) const;

    // Operators
    self_type &  operator =( self_type const &rhs );

    operator bool_type() const;

protected:
    // Member types
    typedef unsigned int              word_type;
    typedef std::valarray<word_type>    va_type;
    typedef std::auto_ptr<va_type>      ap_type;

    typedef std::numeric_limits<word_type>  limits_type;

    // Helper functions
    static  size_type  words_for_bits( size_type bits );

    static  va_type  uintmax_to_va( uintmax_t v );
    static  va_type  bit_vector_to_va( mask_type const &b );
    static  va_type  bit_indices_to_va( index_type const &i );

    static  size_type  wlength( va_type const &v );
    static  size_type  blength( word_type w );

    static  word_type  count_set_bits_for_word( word_type w );

    static  int  do_compare( va_type const &lhs, va_type const &rhs );

private:
    // The non-assignment operators are the primary routines
    friend  self_type operator &( self_type const &lhs, self_type const &rhs );
    friend  self_type operator |( self_type const &lhs, self_type const &rhs );
    friend  self_type operator ^( self_type const &lhs, self_type const &rhs );

    // More member types
    enum bit_operation { reset_bit, set_bit, flip_bit };

    // Object-specific helper functions
    void  bit_change( size_type from, size_type to, bit_operation op );
    void  bit_change( size_type i, bit_operation op );

    // Member data
    ap_type  x_;

};  // boost::math::big_whole


//  Arbitrary-length whole-number constructor definitions  -------------------//

inline
big_whole::big_whole
(
)
    : x_()
{
}

inline
big_whole::big_whole
(
    big_whole const &  other
)
    : x_( other.x_.get() ? new va_type(*other.x_) : NULL )
{
}

inline
big_whole::big_whole
(
    uintmax_t  v
)
    : x_( new va_type(self_type::uintmax_to_va( v )) )
{
}

inline
big_whole::big_whole
(
    std::valarray<bool> const &  b
)
    : x_( new va_type(self_type::bit_vector_to_va( b )) )
{
}

inline
big_whole::big_whole
(
    std::valarray<std::size_t> const &  i
)
    : x_( new va_type(self_type::bit_indices_to_va( i )) )
{
}


//  Arbitrary-length whole-number object-mutating member definitions  --------//

void
big_whole::swap
(
    big_whole &  other
)
{
    ap_type  temp( other.x_.release() );

    other.x_.reset( this->x_.release() );
    this->x_.reset( temp.release() );
}

inline
void
big_whole::assign
(
    big_whole const &  other
)
{
    this->x_.reset( other.x_.get() ? new va_type(*other.x_) : NULL );
}

inline
void
big_whole::assign
(
    uintmax_t  v
)
{
    self_type  temp( v );

    this->swap( temp );
}

inline
void
big_whole::reconfigure
(
    std::valarray<bool> const &  b
)
{
    self_type  temp( b );

    this->swap( temp );
}

inline
void
big_whole::reconfigure
(
    std::valarray<std::size_t> const &  i
)
{
    self_type  temp( i );

    this->swap( temp );
}


//  Arbitrary-length whole-number value-accessing member definitions  --------//

uintmax_t
big_whole::to_uintmax
(
) const
{
    // NOTE: this test should be optimized during compile-time

    if ( std::numeric_limits<uintmax_t>::digits > limits_type::digits )
    {
        uintmax_t  temp = 0;

        if ( va_type const * const  v = this->x_.get() )
        {
            for ( size_type  i = v->size() ; i > 0 ; --i )
            {
                temp <<= limits_type::digits;
                temp |= static_cast<uintmax_t>( (*v)[i - 1] );
            }
        }

        return temp;
    }
    else
    {
        return ( this->x_.get() && this->x_->size() ) ? (*this->x_)[ 0 ] : 0;
    }
}

std::valarray<bool>
big_whole::to_bit_vector
(
) const
{
    index_type const  indices = this->to_bit_indices();
    mask_type         temp;

    if ( indices.size() )
    {
        temp.resize( indices.max() + 1, false );
        temp[ indices ] = true;
    }

    return temp;
}

std::valarray<std::size_t>
big_whole::to_bit_indices
(
) const
{
    index_type  temp;

    if ( va_type const * const  v = this->x_.get() )
    {
        if ( size_type const  s = v->size() )
        {
            int const        d = limits_type::digits;
            index_type       temp2( s * d );
            word_type const  mask = 1;
            size_type        bits_used = 0;

            for ( size_type  i = 0, ii = 0 ; i < s ; ++i, ii += d )
            {
                for ( int  j = 0, jj = ii ; j < d ; ++j, ++jj )
                {
                    if ( (mask << j) & (*v)[i] )
                    {
                        temp2[ bits_used++ ] = jj;
                    }
                }
            }

            if ( bits_used )
            {
                temp.resize( bits_used );
                temp = temp2[ std::slice(0, bits_used, 1) ];
            }
        }
    }

    return temp;
}

bool
big_whole::is_even
(
) const
{
    return !this->x_.get() || !this->x_->size() || ( 0u == (1u
     & ( *this->x_ )[ 0 ]) );
}


//  Arbitrary-length whole-number bit-twiddling member definitions  ----------//

inline
void
big_whole::reset
(
)
{
    this->x_.reset();
}

inline
void
big_whole::reset
(
    std::size_t  from,
    std::size_t  to
)
{
    this->bit_change( from, to, reset_bit );
}

inline
void
big_whole::reset
(
    std::size_t  i
)
{
    this->bit_change( i, reset_bit );
}

inline
void
big_whole::set
(
    std::size_t  from,
    std::size_t  to
)
{
    this->bit_change( from, to, set_bit );
}

inline
void
big_whole::set
(
    std::size_t  i
)
{
    this->bit_change( i, set_bit );
}

inline
void
big_whole::flip
(
    std::size_t  from,
    std::size_t  to
)
{
    this->bit_change( from, to, flip_bit );
}

inline
void
big_whole::flip
(
    std::size_t  i
)
{
    this->bit_change( i, flip_bit );
}

inline
void
big_whole::bit_assign
(
    std::size_t  from,
    std::size_t  to,
    bool         value
)
{
    this->bit_change( from, to, value ? set_bit : reset_bit );
}

inline
void
big_whole::bit_assign
(
    std::size_t  i,
    bool         value
)
{
    this->bit_change( i, value ? set_bit : reset_bit );
}

void
big_whole::bits_assign
(
    std::size_t        from,
    std::size_t        to,
    big_whole const &  values
)
{
    if ( from > to )
    {
        std::swap( from, to );
    }

    index_type const  t_ids = this->to_bit_indices();
    index_type const  v_ids = values.to_bit_indices();

    index_type const   pre_t_ids = t_ids[ t_ids < from ];
    index_type const   new_v_ids = v_ids[ v_ids <= (to - from) ] + from;
    index_type const  post_t_ids = t_ids[ t_ids > to ];

    size_type const   pre_s =  pre_t_ids.size();
    size_type const   new_s =  new_v_ids.size();
    size_type const  post_s = post_t_ids.size();

    index_type   new_t_ids( pre_s + new_s + post_s );
    size_type *  new_t_ptr = &new_t_ids[ 0 ];

    for ( size_type  i = 0 ; i < pre_s ; ++i, ++new_t_ptr )
    {
        *new_t_ptr = pre_t_ids[ i ];
    }

    for ( size_type  j = 0 ; j < new_s ; ++j, ++new_t_ptr )
    {
        *new_t_ptr = new_v_ids[ j ];
    }

    for ( size_type  k = 0 ; k < post_s ; ++k, ++new_t_ptr )
    {
        *new_t_ptr = post_t_ids[ k ];
    }

    this->reconfigure( new_t_ids );
}


//  Arbitrary-length whole-number bit-inspecting member definitions  ---------//

std::size_t
big_whole::length
(
) const
{
    if ( va_type const * const  v = this->x_.get() )
    {
        if ( size_type const  o = self_type::wlength(*v) )
        {
            size_type const  oo = o - 1;

            return self_type::blength( (*v)[oo] ) + oo * limits_type::digits;
        }
    }

    return 0;
}

std::size_t
big_whole::count
(
) const
{
    if ( va_type const * const  v = this->x_.get() )
    {
        return v->size() ? v->apply( self_type::count_set_bits_for_word ).sum()
         : 0u;
    }
    else
    {
        return 0;
    }
}

inline
bool
big_whole::any
(
) const
{
    return this->x_.get() && this->x_->size() && this->x_->max();
}

inline
bool
big_whole::none
(
) const
{
    return !this->any();
}

bool
big_whole::test
(
    std::size_t  i
) const
{
    if ( va_type const * const  v = this->x_.get() )
    {
        size_type const  wi = i / limits_type::digits;
        size_type const  wj = i % limits_type::digits;

        return ( v->size() > wi ) ? ( (*v)[wi] & (1u << wj) ) : false;
    }
    else
    {
        return false;
    }
}

big_whole
big_whole::tests
(
    std::size_t  from,
    std::size_t  to
) const
{
    if ( from > to )
    {
        std::swap( from, to );
    }

    index_type const      ids = this->to_bit_indices();
    index_type        new_ids = ids[ (ids >= from) && (ids <= to) ];

    if ( new_ids.size() )
    {
        new_ids -= from;
    }

    return self_type( new_ids );
}

big_whole
big_whole::reverse
(
    std::size_t  cap
) const
{
    index_type const      ids = this->to_bit_indices();
    index_type const  new_ids = ids[ ids <= cap ];

    return self_type( cap - new_ids );
}

inline
big_whole
big_whole::reverse
(
) const
{
    return this->any() ? this->reverse( this->length() - 1 ) : *this;
}


//  Arbitrary-length whole-number self-operator member definitions  ----------//

inline
void
big_whole::not_self
(
)
{
    this->x_.reset( this->any() ? NULL : new va_type(1u, 1) );
}


//  Arbitrary-length whole-number object-accessing member definitions  -------//

int
big_whole::compare
(
    big_whole const &  other
) const
{
    if ( va_type const * const  t = this->x_.get() )
    {
        if ( va_type const * const  o = other.x_.get() )
        {
            return self_type::do_compare( *t, *o );
        }
        else
        {
            return this->any() ? +1 : 0;
        }
    }
    else
    {
        return other.any() ? -1 : 0;
    }
}


//  Arbitrary-length whole-number member operator definitions  ---------------//

inline
big_whole &
big_whole::operator =
(
    big_whole const &  rhs
)
{
    return this->assign( rhs ), *this;
}

inline
big_whole::operator big_whole::bool_type
(
) const
{
    return this->any() ? &dummy::d : NULL;
}


//  Arbitrary-length whole-number helper static member function definitions  -//

inline
std::size_t
big_whole::words_for_bits
(
    std::size_t  bits
)
{
    return ( bits / limits_type::digits ) + static_cast< std::size_t >( 0
     != (bits % limits_type::digits) );
}

big_whole::va_type
big_whole::uintmax_to_va
(
    uintmax_t  v
)
{
    int const  d = std::numeric_limits<uintmax_t>::digits;

    // NOTE: this test should be optimized during compile-time

    if ( d > limits_type::digits )
    {
        size_type const  s = self_type::words_for_bits( d );
        va_type          temp( s );
        size_type        words_used = 0;

        for ( ; v && (words_used < s) ; ++words_used )
        {
            temp[ words_used ] = static_cast< word_type >( v );
            v >>= limits_type::digits;
        }

        return words_used ? temp[ std::slice(0, words_used, 1) ] : va_type();
    }
    else
    {
        return va_type( static_cast<word_type>(v), 1 );
    }
}

big_whole::va_type
big_whole::bit_vector_to_va
(
    std::valarray<bool> const &  b
)
{
    va_type  temp;

    if ( size_type const  bs = b.size() )
    {
        word_type const  m = 1;
        int const        d = limits_type::digits;

        temp.resize( self_type::words_for_bits(bs), 0u );

        for ( size_type  i = 0 ; i < bs ; ++i )
        {
            if ( b[i] )
            {
                size_type const  wi = i / d;
                size_type const  wj = i % d;

                temp[ wi ] |= ( m << wj );
            }
        }
    }

    return temp;
}

big_whole::va_type
big_whole::bit_indices_to_va
(
    std::valarray<std::size_t> const &  i
)
{
    mask_type  temp;

    if ( i.size() )
    {
        temp.resize( i.max() + 1, false );
        temp[ i ] = true;
    }

    return self_type::bit_vector_to_va( temp );
}

std::size_t
big_whole::wlength
(
    big_whole::va_type const &  v
)
{
    size_type  i = v.size();

    while ( i && !v[i - 1] )
    {
        --i;
    }

    return i;
}

std::size_t
big_whole::blength
(
    big_whole::word_type  w
)
{
    int  i = limits_type::digits;

    while ( i && !(w & ( 1u << (i - 1) )) )
    {
        --i;
    }

    return i;
}

big_whole::word_type
big_whole::count_set_bits_for_word
(
    big_whole::word_type  w
)
{
    word_type  c = 0;

    for ( int i = 0 ; w && (i < limits_type::digits) ; ++i, w >>= 1 )
    {
        c += ( 0u != (w & 1u) );
    }

    return c;
}

int
big_whole::do_compare
(
    big_whole::va_type const &  lhs,
    big_whole::va_type const &  rhs
)
{
    size_type const  l_size = lhs.size();
    size_type const  r_size = rhs.size();

    for ( size_type  i = l_size ; i > r_size ; --i )
    {
        if ( lhs[i - 1] )
        {
            return +1;
        }
    }

    for ( size_type  j = r_size ; j > l_size ; --j )
    {
        if ( rhs[j - 1] )
        {
            return -1;
        }
    }

    for ( size_type  k = std::min(l_size, r_size) ; k > 0 ; --k )
    {
        size_type const  kk = k - 1;
        word_type const  lw = lhs[ kk ];
        word_type const  rw = rhs[ kk ];

        if ( lw < rw )
        {
            return -1;
        }
        else if ( lw > rw )
        {
            return +1;
        }
    }

    return 0;
}


//  Arbitrary-length whole-number helper member function definitions  --------//

void
big_whole::bit_change
(
    std::size_t               from,
    std::size_t               to,
    big_whole::bit_operation  op
)
{
    using std::slice;

    if ( from > to )
    {
        std::swap( from, to );
    }

    size_type const  fi = from / limits_type::digits;
    size_type const  fj = from % limits_type::digits;
    size_type const  ti = to / limits_type::digits;
    size_type const  tj = to % limits_type::digits;

    va_type *  v = this->x_.get();
    size_type  v_size = v ? v->size() : 0;

    word_type const  fm = ~( (1u << fj) - 1u );
    word_type const  tm = ( 1u << tj ) | ( (1u << tj) - 1u );

    if ( (ti >= v_size) && (op != reset_bit) )
    {
        ap_type  p( new va_type(ti + 1) );

        if ( v && v_size )
        {
            ( *p )[ slice(0, v_size, 1) ] = *v;
        }

        this->x_.reset( p.release() );
        v = this->x_.get();
        v_size = v->size();
    }

    if ( fi == ti )
    {
        word_type const  ftm = fm & tm;

        switch ( op )
        {
        case flip_bit:
            ( *v )[ fi ] ^= ftm;
            break;

        case set_bit:
            ( *v )[ fi ] |= ftm;
            break;

        case reset_bit:
        default:
            if ( ti < v_size )
            {
                ( *v )[ fi ] &= ~ftm;
            }
            break;
        }
    }
    else
    {
        // lowest affected word
        switch ( op )
        {
        case flip_bit:
            ( *v )[ fi ] ^= fm;
            break;

        case set_bit:
            ( *v )[ fi ] |= fm;
            break;

        case reset_bit:
        default:
            if ( fi < v_size )
            {
                ( *v )[ fi ] &= ~fm;
            }
            break;
        }

        // middle affected word(s), if any
        size_type const  start = fi + 1;
        size_type const  stop = std::min( v_size, ti );

        if ( stop > start )
        {
            slice const  ids( start, stop - start, 1 );

            switch ( op )
            {
            case flip_bit:
                ( *v )[ ids ] = static_cast<va_type>( (*v)[ids] )
                 ^ limits_type::max();
                break;

            case set_bit:
                ( *v )[ ids ] = limits_type::max();
                break;

            case reset_bit:
            default:
                ( *v )[ ids ] = limits_type::min();
                break;
            }
        }

        // highest affected word
        switch ( op )
        {
        case flip_bit:
            ( *v )[ ti ] ^= tm;
            break;

        case set_bit:
            ( *v )[ ti ] |= tm;
            break;

        case reset_bit:
        default:
            if ( ti < v_size )
            {
                ( *v )[ ti ] &= ~tm;
            }
            break;
        }
    }
}

void
big_whole::bit_change
(
    std::size_t               i,
    big_whole::bit_operation  op
)
{
    size_type const      ii = i / limits_type::digits;
    va_type *             v = this->x_.get();
    size_type        v_size = v ? v->size() : 0;

    if ( (ii >= v_size) && (op != reset_bit) )
    {
        ap_type  p( new va_type(ii + 1) );

        if ( v && v_size )
        {
            ( *p )[ std::slice(0, v_size, 1) ] = *v;
        }

        this->x_.reset( p.release() );
        v = this->x_.get();
        v_size = v->size();
    }

    size_type const  ij = i % limits_type::digits;
    word_type const   m = 1u << ij;

    switch ( op )
    {
    case flip_bit:
        ( *v )[ ii ] ^= m;
        break;

    case set_bit:
        ( *v )[ ii ] |= m;
        break;

    case reset_bit:
    default:
        if ( ii < v_size )
        {
            ( *v )[ ii ] &= ~m;
        }
        break;
    }
}


//  Arbitrary-length whole-number non-member function definitions  -----------//

inline
void
swap
(
    big_whole &  a,
    big_whole &  b
)
{
    a.swap( b );
}


//  Arbitrary-length whole-number non-member operator definitions  -----------//

inline
big_whole
operator !
(
    big_whole const &  w
)
{
    big_whole  temp( w );

    return temp.not_self(), temp;
}

#define BOOST_PRIVATE_COMPARE_OP( Op )   \
    inline                               \
    bool                                 \
    operator Op                          \
    (                                    \
        big_whole const &  lhs,          \
        big_whole const &  rhs           \
    )                                    \
    {                                    \
        return lhs.compare( rhs ) Op 0;  \
    }

BOOST_PRIVATE_COMPARE_OP( == )
BOOST_PRIVATE_COMPARE_OP( != )
BOOST_PRIVATE_COMPARE_OP( < )
BOOST_PRIVATE_COMPARE_OP( > )
BOOST_PRIVATE_COMPARE_OP( <= )
BOOST_PRIVATE_COMPARE_OP( >= )

#undef BOOST_PRIVATE_COMPARE_OP

big_whole
operator &
(
    big_whole const &  lhs,
    big_whole const &  rhs
)
{
    typedef big_whole::va_type      va_type;
    typedef big_whole::size_type  size_type;

    if ( va_type const * const  l = lhs.x_.get() )
    {
        if ( va_type const * const  r = rhs.x_.get() )
        {
            big_whole  temp;

            temp.x_.reset( new va_type(std::min( l->size(), r->size() )) );

            va_type &  t = *temp.x_;

            if ( size_type const  s = t.size() )
            {
                std::slice const  sl( 0, s, 1 );

                t  = ( *l )[ sl ];
                t &= ( *r )[ sl ];
            }

            return temp;
        }
        else
        {
            return rhs;
        }
    }
    else
    {
        return lhs;
    }
}

big_whole
operator |
(
    big_whole const &  lhs,
    big_whole const &  rhs
)
{
    typedef big_whole::va_type      va_type;
    typedef big_whole::size_type  size_type;

    if ( va_type const * const  l = lhs.x_.get() )
    {
        if ( va_type const * const  r = rhs.x_.get() )
        {
            using std::slice;

            big_whole        temp;
            size_type const  ls = l->size();
            size_type const  rs = r->size();

            temp.x_.reset( new va_type(std::max( ls, rs )) );

            va_type &  t = *temp.x_;

            t[ slice(0, ls, 1) ]  = *l;
            t[ slice(0, rs, 1) ] |= *r;

            return temp;
        }
        else
        {
            return lhs;
        }
    }
    else
    {
        return rhs;
    }
}

big_whole
operator ^
(
    big_whole const &  lhs,
    big_whole const &  rhs
)
{
    typedef big_whole::va_type      va_type;
    typedef big_whole::size_type  size_type;

    if ( va_type const * const  l = lhs.x_.get() )
    {
        if ( va_type const * const  r = rhs.x_.get() )
        {
            using std::slice;

            big_whole        temp;
            size_type const  ls = l->size();
            size_type const  rs = r->size();

            temp.x_.reset( new va_type(std::max( ls, rs )) );

            va_type &  t = *temp.x_;

            t[ slice(0, ls, 1) ]  = *l;
            t[ slice(0, rs, 1) ] ^= *r;

            return temp;
        }
        else
        {
            return lhs;
        }
    }
    else
    {
        return rhs;
    }
}

#define BOOST_PRIVATE_MIXED_ASSIGN_OP( Op )  \
    inline                                   \
    big_whole &                              \
    operator Op##=                           \
    (                                        \
        big_whole &        lhs,              \
        big_whole const &  rhs               \
    )                                        \
    {                                        \
        return lhs = lhs Op rhs;             \
    }

BOOST_PRIVATE_MIXED_ASSIGN_OP( & )
BOOST_PRIVATE_MIXED_ASSIGN_OP( | )
BOOST_PRIVATE_MIXED_ASSIGN_OP( ^ )

#undef BOOST_PRIVATE_MIXED_ASSIGN_OP


}  // namespace math
}  // namespace boost


namespace std
{


//  Specialization of numeric-limits declaration  ----------------------------//

template < >
class numeric_limits< ::boost::math::big_whole >
{
    typedef ::boost::math::big_whole  type;

    // NOTE: functions whose definitions aren't given here are never defined

public:
    static  bool const  is_specialized = true;
    static  type  min() throw()  { return type(); }
    static  type  max() throw();

    static  int const  digits = 0;
    static  int const  digits10 = 0;
    static  bool const  is_signed = false;
    static  bool const  is_integer = true;
    static  bool const  is_exact = true;
    static  int const  radix = 2;
    static  type  epsilon() throw();
    static  type  round_error() throw();

    static  int const  min_exponent = 0;
    static  int const  min_exponent10 = 0;
    static  int const  max_exponent = 0;
    static  int const  max_exponent10 = 0;

    static  bool const  has_infinity = false;
    static  bool const  has_quiet_NaN = false;
    static  bool const  has_signaling_NaN = false;
    static  float_denorm_style const  has_denorm = denorm_absent;
    static  bool const  has_denorm_loss = false;
    static  type  infinity() throw();
    static  type  quiet_NaN() throw();
    static  type  signaling_NaN() throw();
    static  type  denorm_min() throw();

    static  bool const  is_iec559 = false;
    static  bool const  is_bounded = false;
    static  bool const  is_modulo = false;

    static  bool const  traps = false;
    static  bool const  tinyness_before = false;
    static  float_round_style const  round_style = round_toward_zero;

};  // std::numeric_limits<boost::math::big_whole>


}  // namespace std


#endif  // BOOST_MATH_BIG_WHOLE_CORE_HPP
