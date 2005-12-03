//  Boost math/cayley_element.hpp header file  -------------------------------//

//  Copyright 2005 Daryle Walker.  Use, modification, and distribution are
//  subject to the Boost Software License, Version 1.0.  (See accompanying file
//  LICENSE_1_0.txt or a copy at <http://www.boost.org/LICENSE_1_0.txt>.)

//  See <http://www.boost.org/libs/math/> for the library's home page.

/** \file
    \brief Cayley-Dickson hypercomplex element classes

    The main class, boost::math::cayley_element, represents a unit for one
    hypercomplex basis.  The class boost::math::negated_cayley_element and
    class template boost::math::scaled_cayley_element represent the results
    of most types of manipulations with cayley_element, like negation or
    multiplication.  Operations like addition, subtraction, and the
    transcedentals are not included because they need a general hypercomplex
    type to represent the results.  Such hypercomplex types should include:
    conversion constructors for the element types, mixed arithmetic and equality
    operators, transcedental functions, and element-on-element and
    element-on-scalar addition and subtraction operators.  The last category
    of operators could be in a sub-namespace where it has to be explicitly
    injected into the computing workspaces, so those operators won't clash with
    another hypercomplex type that does the same service.  With a single such
    hypercomplex type, users should be working only with cayley_element and use
    various arithmetic operations to implicitly utilize the other types:

    \code
typedef boost::math::cayley_element  e;
using namespace  hypercomplex;
using namespace  hypercomplex::mixed_ops;
my_hypercomplex_t const  x = ( 2.0 + e(1) - 3.0 * e(4) ) / -e(13);
    \endcode
 */

#ifndef BOOST_MATH_CAYLEY_ELEMENT_HPP
#define BOOST_MATH_CAYLEY_ELEMENT_HPP

#include <algorithm>  // for std::max, std::swap
#include <cmath>      // for std::abs, std::atan2, std::pow
#include <cstddef>    // for std::size_t, NULL
#include <cstdlib>    // for std::abs


namespace boost
{
namespace math
{


//  Forward declarations  ----------------------------------------------------//

// Cayley-Dickson hypercomplex element types
class cayley_element;

class negated_cayley_element;

template < typename RealType >
    class scaled_cayley_element;

// Equality operators
//! Compares two scaled elements for equality
template < typename T >
  bool  operator ==( scaled_cayley_element<T> const &l,
   scaled_cayley_element<T> const &r );
//! Compares two signed elements for equality
bool    operator ==( negated_cayley_element const &l,
 negated_cayley_element const &r );
//! Compares two basic elements for equality
bool    operator ==( cayley_element const &l, cayley_element const &r );

//! Compares two scaled elements for anti-equality
template < typename T >
  bool  operator !=( scaled_cayley_element<T> const &l,
   scaled_cayley_element<T> const &r );
//! Compares two signed elements for anti-equality
bool    operator !=( negated_cayley_element const &l,
 negated_cayley_element const &r );
//! Compares two basic elements for anti-equality
bool    operator !=( cayley_element const &l, cayley_element const &r );

// Unary operators
//! Generates the logical negation of a scaled element
template < typename T >
  bool  operator !( scaled_cayley_element<T> const &x );
//! Generates the logical negation of a signed element
bool    operator !( negated_cayley_element const &x );
//! Generates the logical negation of a basic element
bool    operator !( cayley_element const &x );

//! Generates a copy of a scaled element
template < typename T >
  scaled_cayley_element<T>  operator +( scaled_cayley_element<T> const &x );
//! Generates a copy of a signed element
negated_cayley_element      operator +( negated_cayley_element const &x );
//! Generates a copy of a basic element
cayley_element              operator +( cayley_element const &x );

//! Generates the additive inverse of a scaled element
template < typename T >
  scaled_cayley_element<T>  operator -( scaled_cayley_element<T> const &x );
//! Generates the additive inverse of a signed element
negated_cayley_element      operator -( negated_cayley_element const &x );
//! Generates the additive inverse of a basic element
negated_cayley_element      operator -( cayley_element const &x );

// Condition functions
//! Finds the dot product of two basic elements
int  dot_product( cayley_element const &l, cayley_element const &r );
//! Finds the dot product of two signed elements
int  dot_product( negated_cayley_element const &l,
      negated_cayley_element const &r );
//! Finds the dot product of two scaled elements
template < typename T >
  T  dot_product( scaled_cayley_element<T> const &l,
      scaled_cayley_element<T> const &r );

//! Finds the hypercomplex conjugate of a basic element
negated_cayley_element      conj( cayley_element const &x );
//! Finds the hypercomplex conjugate of a signed element
negated_cayley_element      conj( negated_cayley_element const &x );
//! Finds the hypercomplex conjugate of a scaled element
template < typename T >
  scaled_cayley_element<T>  conj( scaled_cayley_element<T> const &x );

//! Finds the magnitude of a basic element
int  abs( cayley_element const &x );
//! Finds the magnitude of a signed element
int  abs( negated_cayley_element const &x );
//! Finds the magnitude of a scaled element
template < typename T >
  T  abs( scaled_cayley_element<T> const &x );

//! Finds the vector angle of a basic element
double  arg( cayley_element const &x );
//! Finds the vector angle of a signed element
double  arg( negated_cayley_element const &x );
//! Finds the vector angle of a scaled element
template < typename T >
  T     arg( scaled_cayley_element<T> const &x );

//! Finds the Cayley-norm of a basic element
int  norm( cayley_element const &x );
//! Finds the Cayley-norm of a signed element
int  norm( negated_cayley_element const &x );
//! Finds the Cayley-norm of a scaled element
template < typename T >
  T  norm( scaled_cayley_element<T> const &x );

//! Finds the sign of a basic element
cayley_element              sgn( cayley_element const &x );
//! Finds the sign of a signed element
negated_cayley_element      sgn( negated_cayley_element const &x );
//! Finds the sign of a scaled element
template < typename T >
  scaled_cayley_element<T>  sgn( scaled_cayley_element<T> const &x );

//! Finds the reciprocal of a basic element
negated_cayley_element      reciprocal( cayley_element const &x );
//! Finds the reciprocal of a signed element
negated_cayley_element      reciprocal( negated_cayley_element const &x );
//! Finds the reciprocal of a scaled element
template < typename T >
  scaled_cayley_element<T>  reciprocal( scaled_cayley_element<T> const &x );

//! Finds the infinity-norm of a basic element
int  sup( cayley_element const &x );
//! Finds the infinity-norm of a signed element
int  sup( negated_cayley_element const &x );
//! Finds the infinity-norm of a scaled element
template < typename T >
  T  sup( scaled_cayley_element<T> const &x );

//! Finds the 1-norm of a basic element
int  l1( cayley_element const &x );
//! Finds the 1-norm of a signed element
int  l1( negated_cayley_element const &x );
//! Finds the 1-norm of a scaled element
template < typename T >
  T  l1( scaled_cayley_element<T> const &x );

// Scalar-multiplicative operators
//! Generates the product of a scalar and a scaled element
template < typename T >
  scaled_cayley_element<T>  operator *( T const &l,
   scaled_cayley_element<T> const &r );
//! Generates an amplified copy of a basic element
template < typename T >
  scaled_cayley_element<T>  operator *( scaled_cayley_element<T> const &l,
   T const &r );
//! Generates the quotient of a scalar divided by a scaled element
template < typename T >
  scaled_cayley_element<T>  operator /( T const &l,
   scaled_cayley_element<T> const &r );
//! Generates an attenuated copy of a scaled element
template < typename T >
  scaled_cayley_element<T>  operator /( scaled_cayley_element<T> const &l,
   T const &r );

//! Generates the product of a scalar and a signed element
template < typename T >
  scaled_cayley_element<T>  operator *( T const &l,
   negated_cayley_element const &r );
//! Generates an up-scaled copy of a signed element
template < typename T >
  scaled_cayley_element<T>  operator *( negated_cayley_element const &l,
   T const &r );
//! Generates the quotient of a scalar divided by a signed element
template < typename T >
  scaled_cayley_element<T>  operator /( T const &l,
   negated_cayley_element const &r );
//! Generates a down-scaled copy of a signed element
template < typename T >
  scaled_cayley_element<T>  operator /( negated_cayley_element const &l,
   T const &r );

//! Generates the product of a scalar and a basic element
template < typename T >
  scaled_cayley_element<T>  operator *( T const &l, cayley_element const &r );
//! Generates an up-scaled copy of a basic element
template < typename T >
  scaled_cayley_element<T>  operator *( cayley_element const &l, T const &r );
//! Generates the quotient of a scalar divided by a basic element
template < typename T >
  scaled_cayley_element<T>  operator /( T const &l, cayley_element const &r );
//! Generates a down-scaled copy of a basic element
template < typename T >
  scaled_cayley_element<T>  operator /( cayley_element const &l, T const &r );

// Element-multiplicative operators
//! Generates the product of two signed elements
negated_cayley_element  operator *( negated_cayley_element const &l,
 negated_cayley_element const &r );
//! Generates the product of a signed element and a basic element
negated_cayley_element  operator *( negated_cayley_element const &l,
 cayley_element const &r );
//! Generates the product of a basic element and a signed element
negated_cayley_element  operator *( cayley_element const &l,
 negated_cayley_element const &r );
//! Generates the product of two basic elements
negated_cayley_element  operator *( cayley_element const &l,
 cayley_element const &r );

//! Generates the product of two scaled elements
template < typename T >
  scaled_cayley_element<T>  operator *( scaled_cayley_element<T> const &l,
   scaled_cayley_element<T> const &r );
//! Generates the product of a scaled element and a basic element
template < typename T >
  scaled_cayley_element<T>  operator *( scaled_cayley_element<T> const &l,
   cayley_element const &r );
//! Generates the product of a basic element and a scaled element
template < typename T >
  scaled_cayley_element<T>  operator *( cayley_element const &l,
   scaled_cayley_element<T> const &r );
//! Generates the product of a scaled element and a signed element
template < typename T >
  scaled_cayley_element<T>  operator *( scaled_cayley_element<T> const &l,
   negated_cayley_element const &r );
//! Generates the product of a signed element and a scaled element
template < typename T >
  scaled_cayley_element<T>  operator *( negated_cayley_element const &l,
   scaled_cayley_element<T> const &r );

//! Generates the quotient of a basic element divided by another
negated_cayley_element  operator /( cayley_element const &l,
 cayley_element const &r );
//! Generates the quotient of a basic element divided by a signed element
negated_cayley_element  operator /( cayley_element const &l,
 negated_cayley_element const &r );
//! Generates the quotient of a signed element divided by a basic element
negated_cayley_element  operator /( negated_cayley_element const &l,
 cayley_element const &r );
//! Generates the quotient of a signed element divided by another
negated_cayley_element  operator /( negated_cayley_element const &l,
 negated_cayley_element const &r );

//! Generates the quotient of a scaled element divided by another
template < typename T >
  scaled_cayley_element<T>  operator /( scaled_cayley_element<T> const &l,
   scaled_cayley_element<T> const &r );
//! Generates the quotient of a scaled element divided by a basic element
template < typename T >
  scaled_cayley_element<T>  operator /( scaled_cayley_element<T> const &l,
   cayley_element const &r );
//! Generates the quotient of a basic element divided by a scaled element
template < typename T >
  scaled_cayley_element<T>  operator /( cayley_element const &l,
   scaled_cayley_element<T> const &r );
//! Generates the quotient of a scaled element divided by a signed element
template < typename T >
  scaled_cayley_element<T>  operator /( scaled_cayley_element<T> const &l,
   negated_cayley_element const &r );
//! Generates the quotient of a signed element divided by a scaled element
template < typename T >
  scaled_cayley_element<T>  operator /( negated_cayley_element const &l,
   scaled_cayley_element<T> const &r );

// Shift operators
//! Generates a basic element copy with a higher basis
cayley_element  operator <<( cayley_element const &l, std::size_t r );
//! Generates a basic element copy with a lower basis
cayley_element  operator >>( cayley_element const &l, std::size_t r );

//! Generates a signed element copy with a higher basis
negated_cayley_element  operator <<( negated_cayley_element const &l,
 std::size_t r );
//! Generates a signed element copy with a lower basis
negated_cayley_element  operator >>( negated_cayley_element const &l,
 std::size_t r );

//! Generates a scaled element copy with a higher basis
template < typename T >
  scaled_cayley_element<T>  operator <<( scaled_cayley_element<T> const &l,
   std::size_t r );
//! Generates a scaled element copy with a lower basis
template < typename T >
  scaled_cayley_element<T>  operator >>( scaled_cayley_element<T> const &l,
   std::size_t r );

// Component functions
//! Extracts the real component of a basic element
int  real( cayley_element const &x );
//! Extracts the real component of a signed element
int  real( negated_cayley_element const &x );
//! Extracts the real component of a scaled element
template < typename T >
  T  real( scaled_cayley_element<T> const &x );

//! Extracts the imaginary (complex) component of a basic element
int  imag( cayley_element const &x );
//! Extracts the imaginary (complex) component of a signed element
int  imag( negated_cayley_element const &x );
//! Extracts the imaginary (complex) component of a scaled element
template < typename T >
  T  imag( scaled_cayley_element<T> const &x );

//! Extracts the unreal components of a scaled element
template < typename T >
  scaled_cayley_element<T>  unreal( scaled_cayley_element<T> const &x );
//! Extracts the unreal components of a signed element
scaled_cayley_element<int>  unreal( negated_cayley_element const &x );
//! Extracts the unreal components of a basic element
scaled_cayley_element<int>  unreal( cayley_element const &x );

// Integer-power functions
//! Raises a signed element to an integer power
negated_cayley_element      pow( negated_cayley_element const &b, int e );
//! Raises a basic element to an integer power
negated_cayley_element      pow( cayley_element const &b, int e );
//! Raises a scaled element to an integer power
template < typename T >
  scaled_cayley_element<T>  pow( scaled_cayley_element<T> const &b, int e );


//  Cayley algebra element classes declarations  -----------------------------//

//! Represents a Cayley-Dickson hypercomplex element unit
/** Represents a unit (\e i.e. length and direction +1) for a Cayley-Dickson
    hypercomplex element.  The desired basis is determined at run-time.
 */
class cayley_element
{
    struct dummy { dummy *d; };
    typedef dummy * dummy::*  bool_type;
    typedef cayley_element    self_type;

public:
    // Types
    //! Type for enumerating hypercomplex bases
    /** Represents an index for an individual basis.  0 is for the real numbers,
        1 for imaginary units, 2 and 3 for the quaternion units, 4 through 7 for
        the octonion units, \e etc.
     */
    typedef std::size_t  index_type;
    //! Type for enumerating Cayley-Dickson construction levels
    /** Represents a level that a basis can belong to.  0 is for real numbers,
        1 for complex numbers, 2 for quaternions, 3 for octonions, \e etc.
     */
    typedef std::size_t  rung_type;

    // Index <-> range conversion functions
    //! Compute the minimum rung a basis can belong to
    static  rung_type  minimum_rung_for_index( index_type i );
    //! Compute the maximum basis index for a given rung
    static  index_type  maximum_index_for_rung( rung_type r );

    // Lifetime management (use automatic destructor and copy constructor)
    //! Create an element with the given basis
    explicit  cayley_element( index_type b );

    // Accessors
    //! Returns this element's basis
    index_type  basis() const;

    //! Returns the (minimum) rung for this element's basis
    rung_type  minimum_rung() const;

    // Mutators
    //! Changes the element's basis
    void  basis( index_type b );

    // Operators (use automatic assignment)
    //! Converts this element to a Boolean value
    operator bool_type() const;

    //! Changes the basis to a higher index
    self_type &  operator <<=( std::size_t r );
    //! Changes the basis to a lower index
    self_type &  operator >>=( std::size_t r );

    // Self-serving mutators for unary operators & condition operations
    //! Applies *this = +(*this)
    self_type &  same_self();

    //! Applies *this = sgn( *this );
    self_type &  sign_self();

private:
    // Member data
    index_type  b_;

};  // boost::math::cayley_element

//! Represents hypercomplex units with a sign (positive or negative)
/** Represents an extension of boost::math::cayley_element by giving run-time
    control over a unit's direction, either positive or negative.  (The length
    stays at 1.)
 */
class negated_cayley_element
{
    struct dummy { dummy *d; };
    typedef dummy * dummy::*        bool_type;
    typedef negated_cayley_element  self_type;

public:
    // Types
    //! Type for enumerating hypercomplex bases
    /** \see boost::math::cayley_element::index_type */
    typedef cayley_element::index_type  index_type;
    //! Type for enumerating Cayley-Dickson construction levels
    /** \see boost::math::cayley_element::rung_type */
    typedef cayley_element::rung_type  rung_type;

    // Lifetime management (use automatic destructor and copy constructor)
    //! Create an element with the given basis and sign
    explicit  negated_cayley_element( index_type b, bool is_negative = false );
    //! Convert a basic element
    negated_cayley_element( cayley_element const &c );

    // Accessors
    //! Returns this element's basis
    index_type  basis() const;
    //! Returns whether or not this element is negative
    bool        negative() const;

    //! Returns the (minimum) rung for this element's basis
    rung_type  minimum_rung() const;

    // Mutators
    //! Changes the element's basis
    void  basis( index_type b );
    //! Changes the element's sign
    void  negative( bool is_negative );

    // Operators (use automatic assignment)
    //! Converts this element to a Boolean value
    operator bool_type() const;

    //! Multiplies this element by another
    self_type &  operator *=( self_type const &r );
    //! Divides this element by another
    self_type &  operator /=( self_type const &r );

    //! Changes the basis to a higher index
    self_type &  operator <<=( std::size_t r );
    //! Changes the basis to a lower index
    self_type &  operator >>=( std::size_t r );

    // Self-serving mutators for unary operators & condition operations
    //! Applies *this = +(*this);
    self_type &  same_self();
    //! Applies *this = -(*this);
    self_type &  negate_self();

    //! Applies *this = conj( *this );
    self_type &  conjugate_self();
    //! Applies *this = sgn( *this );
    self_type &  sign_self();
    //! Applies *this = reciprocal( *this );
    self_type &  reciprocate_self();

private:
    // Member data
    cayley_element  e_;
    bool            n_;

};  // boost::math::negated_cayley_element

//! Represents hypercomplex elements as a unit and a real scale
/** Represents an extension of boost::math::cayley_element by giving run-time
    control over an element's length \e via a scale factor.  This may include
    direction if the scale's type allows negative values.

    \param RealType  The type used for representing scale factors.  It should
                     represent (a subset of) real numbers.  It must support
                     default construction, where said default value must
                     represent zero.  It must support Boolean conversion, where
                     zero represents \c False and non-zero values represent
                     \c True.  If the type only represents nonnegative values,
                     then it cannot be used for element-multiplication,
                     element-division, reciprocation, negation, conjugation, or
                     power.  If the type uses the "quotient & remainder" style
                     of division (instead of "'exact' quotient" style), then it
                     cannot be used for division, reciprocation, or power (with
                     negative exponents).  Otherwise, \c RealType should support
                     all the arithmetic and equality operators.  Operations with
                     \c scaled_cayley_element are \em not guaranteed if a
                     computation based on \c RealType generates an overflow,
                     underflow, not-a-number, infinity, or other irregular
                     result.
 */
template < typename RealType >
class scaled_cayley_element
{
    struct dummy { dummy *d; };
    typedef dummy * dummy::*                 bool_type;
    typedef scaled_cayley_element<RealType>  self_type;

public:
    // Template parameters
    //! Type for the scalar coefficient
    /** Represents the scale factor applied to a hypercomplex element.  This
        alias exposes the type for later meta-programming.
        \see std::complex::value_type
     */
    typedef RealType  value_type;

    // Other types
    //! Type for enumerating hypercomplex bases
    /** \see boost::math::cayley_element::index_type */
    typedef cayley_element::index_type  index_type;
    //! Type for enumerating Cayley-Dickson construction levels
    /** \see boost::math::cayley_element::rung_type */
    typedef cayley_element::rung_type  rung_type;

    // Lifetime management (use automatic destructor and copy constructor)
    //! Create an element with the given basis and scale
    scaled_cayley_element( index_type b, value_type scale );
    //! Cross-version conversion
    template < typename RealType2 >
        scaled_cayley_element( scaled_cayley_element<RealType2> const &c );
    //! Convert a signed element
    scaled_cayley_element( negated_cayley_element const &c );
    //! Convert a basic element
    scaled_cayley_element( cayley_element const &c );
    //! Default constructor; convert a real scalar
    scaled_cayley_element( value_type r = value_type() );

    // Accessors
    //! Returns this element's basis
    index_type  basis() const;
    //! Returns this element's scale factor
    value_type  scale() const;

    //! Returns the (minimum) rung for this element's basis
    rung_type  minimum_rung() const;

    // Mutators
    //! Changes the element's basis
    void  basis( index_type b );
    //! Changes the element's scale (length and/or orientation)
    void  scale( value_type s );

    // Operators (use automatic assignment)
    //! Converts this element to a Boolean value
    operator bool_type() const;

    //! Multiplies this element by a scalar
    self_type &  operator *=( value_type const &r );
    //! Divides this element by a scalar
    self_type &  operator /=( value_type const &r );

    //! Multiplies this element by another
    self_type &  operator *=( self_type const &r );
    //! Divides this element by another
    self_type &  operator /=( self_type const &r );

    //! Changes the basis to a higher index
    self_type &  operator <<=( std::size_t r );
    //! Changes the basis to a lower index
    self_type &  operator >>=( std::size_t r );

    // Self-serving mutators for unary operators & condition operations
    //! Applies *this = !(*this);
    self_type &  not_self();
    //! Applies *this = +(*this);
    self_type &  same_self();
    //! Applies *this = -(*this);
    self_type &  negate_self();

    //! Applies *this = conj( *this );
    self_type &  conjugate_self();
    //! Applies *this = sgn( *this );
    self_type &  sign_self();
    //! Applies *this = reciprocal( *this );
    self_type &  reciprocate_self();

private:
    // Member data
    cayley_element  e_;
    value_type      s_;

};  // boost::math::scaled_cayley_element


//  Cayley algebra elements component member function definitions  -----------//

/** \return  The index for this element's basis */
inline  cayley_element::index_type
cayley_element::basis() const
{ return this->b_; }

/** \param b  The index for the new basis
    \post     <tt>this->basis() == \a b</tt>
    \see      index_type
 */
inline  void
cayley_element::basis( index_type b )
{ this->b_ = b; }

/** \see  cayley_element::basis() */
inline  negated_cayley_element::index_type
negated_cayley_element::basis() const
{ return this->e_.basis(); }

/** \see  cayley_element::basis(index_type) */
inline  void
negated_cayley_element::basis( index_type b )
{ this->e_.basis( b ); }

/** \return  \c True if this element is negative, \c False if positive */
inline  bool
negated_cayley_element::negative() const
{ return this->n_; }

/** \param is_negative  The new sign, using \c False for positive (original
                        orientation) and \c True for negative (reverse
                        orientation)
    \post               <tt>this->negative() == \a is_negative</tt>
 */
inline  void
negated_cayley_element::negative( bool is_negative )
{ this->n_ = is_negative; }

/** \see  cayley_element::basis() */
template < typename T >
inline  typename scaled_cayley_element<T>::index_type
scaled_cayley_element<T>::basis() const
{ return this->e_.basis(); }

/** \see  cayley_element::basis(index_type) */
template < typename T >
inline  void
scaled_cayley_element<T>::basis( index_type b )
{ this->e_.basis( b ); }

/** \return  The coefficient for this element */
template < typename T >
inline  typename scaled_cayley_element<T>::value_type
scaled_cayley_element<T>::scale() const
{ return this->s_; }

/** \param s  The new coefficient
    \pre      \a s should not be irregular (\e e.g. not-a-number)
    \post     <tt>this->scale() == \a s</tt>
    \see      value_type
 */
template < typename T >
inline  void
scaled_cayley_element<T>::scale( value_type s )
{ this->s_ = s; }


//  Cayley algebra elements rung/index member function definitions  ----------//

/** Computes the minimum rung for a given basis.  Groupings are based on powers
    of 2, so the result is the exponent for the smallest power of 2 greater than
    the given basis index.

    \param i  The basis index to be analyzed.
    \return   The minimum rung \a i can belong to.
 */
cayley_element::rung_type
cayley_element::minimum_rung_for_index( index_type i )
{ rung_type r = 0; for ( ; i ; i >>=1 ) ++r; return r; }

/** Computes the maximum basis index for a given rung.  Since groupings are
    based on powers of 2 (zero-based), the result is one less than 2 taken to
    the power of the rung value.

    \param r  The rung level to be analyzed.
    \return   The maximum basis index \a r can support.
 */
inline  cayley_element::index_type
cayley_element::maximum_index_for_rung( rung_type r )
{ return ( static_cast<index_type>(1) << r ) - 1u; }

/** \return  The smallest rung level that can support this element's basis */
inline  cayley_element::rung_type
cayley_element::minimum_rung() const
{ return self_type::minimum_rung_for_index( this->b_ ); }

/** \see  cayley_element::minimum_rung() */
inline  negated_cayley_element::rung_type
negated_cayley_element::minimum_rung() const
{ return this->e_.minimum_rung(); }

/** \see  cayley_element::minimum_rung() */
template < typename T >
inline  typename scaled_cayley_element<T>::rung_type
scaled_cayley_element<T>::minimum_rung() const
{ return this->e_.minimum_rung(); }


//  Cayley algebra elements Boolean conversion member operator definitions  --//

/** Performs Boolean conversion.  Since all objects of this type have an
    implicit coefficient of 1, all objects convert to \c True.  It enables an
    element to be used in a Boolean context for certain C++ expressions or
    statements.

    \return  \c True if the coefficient is non-zero, \c False if it is.  The
             actual return type is a unspecified built-in numeric or pointer
             type using the built-in Boolean rules (non-zero/null to \c True,
             zero/null to \c False).
 */
inline
cayley_element::operator bool_type() const 
{ return &self_type::dummy::d; }

/** Performs Boolean conversion.  Since all objects of this type have an
    implicit coefficient of either -1 or +1, all objects convert to \c True.

    \see cayley_element::operator bool_type()
 */
inline
negated_cayley_element::operator bool_type() const 
{ return &self_type::dummy::d; }

/** Performs Boolean conversion.

    \return  <tt>this->scale() != 0</tt>

    \see cayley_element::operator bool_type()
 */
template < typename T >
inline
scaled_cayley_element<T>::operator bool_type() const
{ return /*this->s_*/ ( this->s_ != T() ) ? &self_type::dummy::d : NULL; }


//  Cayley algebra elements constructor definitions  -------------------------//

/** Constructs an element with the given basis.  It is a unit for that basis, so
    it has a coefficent of 1 (by definition).

    \param b  The index of the basis to use

    \post  <tt>this->basis() == \a b</tt>

    \see index_type
 */
inline
cayley_element::cayley_element
(
    index_type  b
)
    : b_( b )
{}

/** Constructs an element with the given basis and optional sign.  The element
    is a unit, but in a selectable orientation, giving a coefficent of either +1
    or -1.  Not using the \a is_negative parameter has an effect similar to the
    ::boost::math::cayley_element::cayley_element(index_type) constructor.

    \param b            The index of the basis to use
    \param is_negative  Whether or not the element should be negated unit.  If
                        it is omitted, it defaults to \c False.  A \c False
                        value results in a positive unit, use \c True for a
                        negative unit (\e i.e. reversed orientation).

    \post  <tt>this->basis() == \a b && this->negative() == \a is_negative</tt>

    \see index_type
 */
inline
negated_cayley_element::negated_cayley_element
(
    index_type  b,
    bool        is_negative  // = false
)
    : e_( b ), n_( is_negative )
{}

/** Constructs an element by converting a basic unit element.  This unit keeps
    the same basis and coefficent (of +1).

    \param c  The basic unit element to convert

    \post  <tt>this->basis() == c.basis() && not this->negative()</tt>

    \see boost::math::cayley_element
 */
inline
negated_cayley_element::negated_cayley_element
(
    cayley_element const &  c
)
    : e_( c ), n_( false )
{}

/** Constructs an element with the given basis and scale coefficent.

    \param b      The index of the basis to use
    \param scale  The scalar length of the element.  It may be negative, if the
                  type allows it, to reverse the orientation.

    \post  <tt>this->basis() == \a b && this->scale() == \a scale</tt>

    \see index_type
    \see value_type
 */
template < typename T >
inline
scaled_cayley_element<T>::scaled_cayley_element
(
    index_type  b,
    value_type  scale
)
    : e_( b ), s_( scale )
{}

/** Constructs an element by converting from a scaled element with another
    scalar type.  This element copies the basis index and coefficent.  By the
    rules of C++, this constructor template is \e not used for the copy
    constructor.  The automatically-defined copy constructor is used when the
    source and new elements have the same \c value_type.

    \param c  The scaled element to convert

    \pre  The \c value_type for \a c must be convertible to this element's
          \c value_type, and the conversion should be as value-preserving
          as possible.

    \post  <tt>this->basis() == c.basis() && this->scale() == c.scale()</tt>
 */
template < typename T >
template < typename U >
inline
scaled_cayley_element<T>::scaled_cayley_element
(
    scaled_cayley_element<U> const &  c
)
    : e_( c.basis() ), s_( c.scale() )
{}

/** Constructs an element by converting a signed unit element.  This element
    keeps the same basis and coefficent (of either -1 or +1).

    \param c  The signed unit element to convert

    \pre  if \c value_type does not support negative values, then c.negative()
          must not be \c True

    \post  <tt>this->basis() == c.basis() && this->scale() ==
           value_type(c.negative() ? -1 : +1)</tt>

    \see boost::math::negated_cayley_element
 */
template < typename T >
inline
scaled_cayley_element<T>::scaled_cayley_element
(
    negated_cayley_element const &  c
)
    : e_( c.basis() ), s_( c.negative() ? -1 : +1 )
{}

/** Constructs an element by converting a basic unit element.  This element
    keeps the same basis and coefficent (of +1).

    \param c  The basic unit element to convert

    \post  <tt>this->basis() == c.basis() && this->scale()
            == value_type( +1 )</tt>

    \see boost::math::cayley_element
 */
template < typename T >
inline
scaled_cayley_element<T>::scaled_cayley_element
(
    cayley_element const &  c
)
    : e_( c ), s_( 1 )
{}

/** Constructs an element by converting a real scalar.  It also acts as the
    default constructor, with a basis index and coefficent of zero in that case.

    \param r  The scalar to convert.  If it is omitted, it defaults to zero.

    \post  <tt>this->basis() == 0 && this->scale() == \a r</tt>

    \see value_type
 */
template < typename T >
inline
scaled_cayley_element<T>::scaled_cayley_element
(
    value_type  r  // = value_type()
)
    : e_( 0 ), s_( r )
{}


//  Cayley algebra elements equality operator function definitions  ----------//

/** Checks if two elements are equal, in basis and coefficent.  Two zero
    elements compare as equal, even if their bases are different.

    \param l  The first element to be compared.
    \param r  The second element to be compared.

    \return  \c True if both coefficents are zero (no matter the bases) or if
             <tt>l.scale() == r.scale() && l.basis() == r.basis()</tt>,
             otherwise \c False.

    \relates  scaled_cayley_element
 */
template < typename T >
inline  bool
operator ==
(
    scaled_cayley_element<T> const &  l,
    scaled_cayley_element<T> const &  r
)
{
    // if both scales are zero, then they're equal, no matter the bases
    return ( l.scale() == r.scale() ) && ( /*!l.scale()*/(l.scale() == T()) || (l.basis()
     == r.basis()) );
}

/** \overload

    \relates  negated_cayley_element
    \see  bool operator ==(scaled_cayley_element<T> const &,
           scaled_cayley_element<T> const &)
 */
inline  bool
operator ==( negated_cayley_element const &l, negated_cayley_element const &r )
{ return ( l.basis() == r.basis() ) && ( l.negative() == r.negative() ); }

/** \overload

    \relates  cayley_element
    \see  bool operator ==(scaled_cayley_element<T> const &,
           scaled_cayley_element<T> const &)
 */
inline  bool
operator ==( cayley_element const &l, cayley_element const &r )
{ return l.basis() == r.basis(); }

/** Checks if two elements are not equal, in basis and/or coefficent.  However,
    two zero elements compare as equal, even if their bases are different.

    \param l  The first element to be compared.
    \param r  The second element to be compared.

    \return  <tt>!( \a l == \a r )</tt>

    \relates  scaled_cayley_element
    \see  bool operator ==(scaled_cayley_element<T> const &,
           scaled_cayley_element<T> const &)
 */
template < typename T >
inline  bool
operator !=
(
    scaled_cayley_element<T> const &  l,
    scaled_cayley_element<T> const &  r
)
{ return !( l == r ); }

/** \overload

    \relates  negated_cayley_element
    \see  bool operator !=(scaled_cayley_element<T> const &,
           scaled_cayley_element<T> const &)
 */
inline  bool
operator !=( negated_cayley_element const &l, negated_cayley_element const &r )
{ return !( l == r ); }

/** \overload

    \relates  cayley_element
    \see  bool operator !=(scaled_cayley_element<T> const &,
           scaled_cayley_element<T> const &)
 */
inline  bool
operator !=( cayley_element const &l, cayley_element const &r )
{ return !( l == r ); }


//  Cayley algebra elements unary operator function definitions  -------------//

/** Gives the opposite of an element's non-zero state, \e i.e. the opposite of
    its Boolean conversion.  Note that the return type is \c bool, so
    <tt>!!x</tt> cannot be used to return the original value.

    \param x  The element to be analyzed.

    \return  \c False if the coefficent of \a x is non-zero, \c True otherwise.

    \relates  scaled_cayley_element
 */
template < typename T >
inline  bool
operator !( scaled_cayley_element<T> const &x )
{ return x ? false : true; }

/** \overload

    \relates  negated_cayley_element
    \see  bool operator !(scaled_cayley_element<T> const &)
 */
inline  bool
operator !( negated_cayley_element const &x )
{ return false; }

/** \overload

    \relates  cayley_element
    \see  bool operator !(scaled_cayley_element<T> const &)
 */
inline  bool
operator !( cayley_element const &x )
{ return false; }

/** Applies the identity function to an element.

    \param x  The element to be analyzed.

    \return  \a x

    \relates  scaled_cayley_element
 */
template < typename T >
inline  scaled_cayley_element<T>
operator +( scaled_cayley_element<T> const &x )
{ return x; }

/** \overload

    \relates  negated_cayley_element
    \see  scaled_cayley_element<T> operator +(scaled_cayley_element<T> const &)
 */
inline  negated_cayley_element
operator +( negated_cayley_element const &x )
{ return x; }

/** \overload

    \relates  cayley_element
    \see  scaled_cayley_element<T> operator +(scaled_cayley_element<T> const &)
 */
inline  cayley_element
operator +( cayley_element const &x )
{ return x; }

/** Negates an element.  The sum of an element and its negative should be zero.
    Note that this library does \em not provide additive binary operators, since
    the result will be (in general) outside the range of a single elment.

    \pre  The \c value_type for \a x must support negation through an unary
          <tt>operator -</tt>, and the actual coefficent of \a x must have an
          additive inverse that is representable in \c value_type.

    \param x  The element to be analyzed.

    \return  An element \a y such that \a x and \a y have the same basis index,
             but their coefficents are additive inverses of each other.

    \relates  scaled_cayley_element
 */
template < typename T >
inline  scaled_cayley_element<T>
operator -( scaled_cayley_element<T> const &x )
{ return scaled_cayley_element<T>( x.basis(), -x.scale() ); }

/** \overload

    \relates  negated_cayley_element
    \see  scaled_cayley_element<T> operator -(scaled_cayley_element<T> const &)
 */
inline  negated_cayley_element
operator -( negated_cayley_element const &x )
{ return negated_cayley_element( x.basis(), !x.negative() ); }

/** \overload

    \relates  cayley_element
    \see  scaled_cayley_element<T> operator -(scaled_cayley_element<T> const &)
 */
inline  negated_cayley_element
operator -( cayley_element const &x )
{ return negated_cayley_element( x.basis(), true ); }


//  Cayley algebra elements unary-action member function definitions  --------//

/** Transforms this element to itself (\e i.e. no change)

    \return  A reference to this element

    \post  <tt>*this == +old_this</tt>

    \see  cayley_element operator +(cayley_element const &)
 */
inline  cayley_element &
cayley_element::same_self()
{ this->b_ = +this->b_; return *this; }

/** \see  cayley_element::same_self() */
inline  negated_cayley_element &
negated_cayley_element::same_self()
{ this->e_.same_self(); return *this; }

/** Transforms this element to its negative (\e i.e. additive inverse)

    \return  A reference to this element

    \post  <tt>*this == -old_this</tt>

    \see  negated_cayley_element operator -(negated_cayley_element const &)
 */
inline  negated_cayley_element &
negated_cayley_element::negate_self()
{ this->n_ = !this->n_; return *this; }

/** Transforms this element to its logical negation.  The canoical \c False
    value for the coefficent is zero, for \c True it's +1.

    \return  A reference to this element

    \post  <tt>static_cast\<bool\>(*this) == !static_cast\<bool\>(old_this)</tt>

    \see  bool operator !(scaled_cayley_element<T> const &)
 */
template < typename T >
inline  scaled_cayley_element<T> &
scaled_cayley_element<T>::not_self()
{ this->s_ = !this->s_; return *this; }

/** \see  cayley_element::same_self() */
template < typename T >
inline  scaled_cayley_element<T> &
scaled_cayley_element<T>::same_self()
{ this->e_.same_self(); this->s_ = +this->s_; return *this; }

/** \see  negated_cayley_element::negate_self() */
template < typename T >
inline  scaled_cayley_element<T> &
scaled_cayley_element<T>::negate_self()
{ this->s_ = -this->s_; return *this; }


//  Cayley algebra elements condition function definitions  ------------------//

/** Computes the Euclidean inner-product of two elements.  It works
    component-wise like the vector dot-product.  It is commutative and the
    result is a scalar.

    \param l  The left-side factor
    \param r  The right-side factor

    \return  \f$l \cdot r = \sum_i {l_i r_i} = \frac{ \bar{l}r
             + \bar{r}l }{ 2 }\f$

    \relates  cayley_element
 */
inline  int
dot_product( cayley_element const &l, cayley_element const &r )
{ return l.basis() == r.basis(); }

/** \overload

    \relates  negated_cayley_element
    \see  int dot_product(cayley_element const &, cayley_element const &)
 */
inline  int
dot_product( negated_cayley_element const &l, negated_cayley_element const &r )
{
    return ( l.basis() == r.basis() ) ? ( (l.negative() == r.negative()) ? +1
     : -1 ) : 0;
}

/** \overload

    \relates  scaled_cayley_element
    \see  int dot_product(cayley_element const &, cayley_element const &)
 */
template < typename T >
inline  T
dot_product
(
    scaled_cayley_element<T> const &  l,
    scaled_cayley_element<T> const &  r
)
{ return ( l.basis() == r.basis() ) ? ( l.scale() * r.scale() ) : T(); }

/** Computes the conjugate of an element.  The conjugate of a hypercomplex
    number negates all the unreal components.  It is a (linear) involution,
    \e i.e. it's its own inverse (\f$\bar{\bar{x}} = x\f$), that is a core
    element of Cayley-Dickson construction.

    \param x  The element to be analyzed

    \return  \f$\bar{x}\f$, which is the identity function for a real \a x
             (\e i.e. \f$x \textrm{ if } x \in \mathbb{R}\f$), otherwise it
             expands as \f$\overline{(x_l \oplus x_h)} = (\overline{x_l} \oplus
             {-x_h})\f$

    \relates  cayley_element
 */
inline  negated_cayley_element
conj( cayley_element const &x )
{ return x.basis() ? -x : +x; }

/** \overload

    \relates  negated_cayley_element
    \see  negated_cayley_element conj(cayley_element const &)
 */
inline  negated_cayley_element
conj( negated_cayley_element const &x )
{ return x.basis() ? -x : +x; }

/** \overload

    \relates  scaled_cayley_element
    \see  negated_cayley_element conj(cayley_element const &)
 */
template < typename T >
inline  scaled_cayley_element<T>
conj( scaled_cayley_element<T> const &x )
{ return x.basis() ? -x : +x; }

/** Computes the absolute value of an element.  If the components were
    coordinates of a point in Euclidean hyperspace, then the absolute value
    would be the distance from the origin, also the length of a vector from the
    origin to the element's point.  It is also the 2-norm.  By definition, a
    unit element will always have a length of one.

    \param x  The element to be analyzed

    \return  \f$|x| = \|x\|_2 = \sqrt{\sum_i {|x_i|^2}}\f$

    \relates  cayley_element
 */
inline  int
abs( cayley_element const & )
{ return 1; }

/** \overload
    Even if the sign is negative (\e i.e. the orientation is reversed), the
    magnitude of a unit element remains 1.

    \relates  negated_cayley_element
    \see  int abs(cayley_element const &)
 */
inline  int
abs( negated_cayley_element const & )
{ return 1; }

/** \overload
    The magnitude is equivalent to the absolute value of the scale.

    \pre  The \c value_type for \a x must support a function named "abs" that
          takes one argument, which needs to convert from \c value_type, and
          returns a \c value_type.  The function should return the absolute
          value of its argument.  The function needs to be found in the "std"
          name-space or discovered by argument-dependent name lookup on
          \c value_type.

    \relates  scaled_cayley_element
    \see  int abs(cayley_element const &)
 */
template < typename T >
inline  T
abs( scaled_cayley_element<T> const &x )
{ using std::abs; return abs( x.scale() ); }

/** Computes the angle of an element.  Treating the element like a point/vector
    (see the  int abs(cayley_element const &)  notes), the angle is taken from
    the real unit vector towards the element vector, measured on the plane the
    two vectors share.  Therefore, the angle ranges from zero to \f$\pi\f$.  A
    purely real element has an angle of zero if positive, or \f$\pi\f$ if
    negative.  A purely unreal element has an angle of \f$\pi / 2\f$.  A zero
    element is considered to have an angle of zero.

    \param x  The element to be analyzed

    \return  \f$\arg x = \arccos {\frac{\Re(x)}{|x|}}
                       = \arctan {\frac{|\mathfrak{Ur}(x)|}{\Re(x)}}\f$

    \relates  cayley_element
 */
inline  double
arg( cayley_element const &x )
{ return x.basis() ? std::atan2( 1.0, 0.0 ) : std::atan2( 0.0, 1.0 ); }

/** \overload

    \relates  negated_cayley_element
    \see  double arg(cayley_element const &)
 */
inline  double
arg( negated_cayley_element const &x )
{
    return x.basis() ? std::atan2( 1.0, 0.0 ) : std::atan2( 0.0, x.negative()
     ? -1.0 : +1.0 );
}

/** \overload

    \pre  The \c value_type for \a x must support a function named "atan2" that
          takes two arguments, which need to convert from \c value_type, and
          returns a \c value_type.  The function should return the arc-tangent
          value of its arguments, where the first argument is the y-coordinate
          and the second argument is the x-coordinate.  (Unlike the
          single-argument "atan" function with y/x, this function works when the
          x-coordinate is zero.)  The function needs to be found in the "std"
          name-space or discovered by argument-dependent name lookup on
          \c value_type.

    \relates  scaled_cayley_element
    \see  double arg(cayley_element const &)
 */
template < typename T >
inline  T
arg( scaled_cayley_element<T> const &x )
{
    using std::atan2;

    return /*x.scale()*/ ( x.scale() != T() ) ? ( x.basis() ? atan2(abs( x ), T()) : atan2(T(),
     x.scale()) ) : T();
}

/** Computes the Cayley norm of an element.  This norm does \em not follow the
    rules of conventional norms.  It happens to equal the square of the
    element's magnitude.  It is the product of an element and its conjugate (in
    either order, so it equals the Cayley-norm of the conjugate).

    \param x  The element to be analyzed

    \return  \f$|x|^2 = x \bar{x} = \sum_i {|x_i|^2} = x \cdot x\f$

    \relates  cayley_element
 */
inline  int
norm( cayley_element const & )
{ return 1; }

/** \overload

    \relates  negated_cayley_element
    \see  int norm(cayley_element const &)
 */
inline  int
norm( negated_cayley_element const & )
{ return 1; }

/** \overload

    \relates  scaled_cayley_element
    \see  int norm(cayley_element const &)
 */
template < typename T >
inline  T
norm( scaled_cayley_element<T> const &x )
{ return x.scale() * x.scale(); }

/** Computes the sign of an element.  The sign of a hypercomplex value retains
    the original's direction, but has a length of one.  A zero element has a
    sign of zero.  Note that unit and zero elements are unchanged via this
    function.

    \param x  The element to be analyzed

    \return  \f$\mathrm{sgn}(x) = \frac{x}{|x|}\f$

    \relates  cayley_element
 */
inline  cayley_element
sgn( cayley_element const &x )
{ return x; }

/** \overload

    \relates  negated_cayley_element
    \see  cayley_element sgn(cayley_element const &)
 */
inline  negated_cayley_element
sgn( negated_cayley_element const &x )
{ return x; }

/** \overload

    \relates  scaled_cayley_element
    \see  cayley_element sgn(cayley_element const &)
 */
template < typename T >
inline  scaled_cayley_element<T>
sgn( scaled_cayley_element<T> const &x )
{ return /*x.scale()*/ ( x.scale() != T() ) ? negated_cayley_element( x.basis(), x.scale() < T() ) : x; }

/** Computes the multiplicative inverse of an element.

    \param x  The element to be analyzed

    \return  \f$x^{-1} = \frac{1}{x} = \frac{\bar{x}}{|x|^2}\f$

    \relates  cayley_element
 */
inline  negated_cayley_element
reciprocal( cayley_element const &x )
{ return conj( x ); }

/** \overload

    \relates  negated_cayley_element
    \see  negated_cayley_element reciprocal(cayley_element const &)
 */
inline  negated_cayley_element
reciprocal( negated_cayley_element const &x )
{ return conj( x ); }

/** \overload

    \pre  The coefficient for \a x must be invertible.  For instance, it cannot
          be zero.

    \relates  scaled_cayley_element
    \see  negated_cayley_element reciprocal(cayley_element const &)
 */
template < typename T >
inline  scaled_cayley_element<T>
reciprocal( scaled_cayley_element<T> const &x )
{
    return scaled_cayley_element<T>( x.basis(), static_cast<T>(x.basis() ? -1
     : +1) / x.scale() );
}

/** Computes the supremum-, or infinity-, norm of an element.  It is equal to
    the largest absolute value of a component.

    \param x  The element to be analyzed

    \return  \f$\|x\|_\infty = \lim_{n \rightarrow \infty} \sqrt[n]{\sum_i
             {|x_i|^n}} = \max_i |x_i|\f$

    \relates  cayley_element
 */
inline  int
sup( cayley_element const &x )
{ return abs( x ); }

/** \overload

    \relates  negated_cayley_element
    \see  int sup(cayley_element const &)
 */
inline  int
sup( negated_cayley_element const &x )
{ return abs( x ); }

/** \overload

    \relates  scaled_cayley_element
    \see  int sup(cayley_element const &)
 */
template < typename T >
inline  T
sup( scaled_cayley_element<T> const &x )
{ return abs( x ); }

/** Computes the 1-norm of an element.  It is equal to the sum of the absolute
    values of each component.

    \param x  The element to be analyzed

    \return  \f$\|x\|_1 = \sum_i |x_i|\f$

    \relates  cayley_element
 */
inline  int
l1( cayley_element const &x )
{ return abs( x ); }

/** \overload

    \relates  negated_cayley_element
    \see  int l1(cayley_element const &)
 */
inline  int
l1( negated_cayley_element const &x )
{ return abs( x ); }

/** \overload

    \relates  scaled_cayley_element
    \see  int l1(cayley_element const &)
 */
template < typename T >
inline  T
l1( scaled_cayley_element<T> const &x )
{ return abs( x ); }


//  Cayley algebra elements condition member function definitions  -----------//

/** Transforms this element to its unit element

    \return  A reference to this element

    \post  <tt>*this == sgn( old_this )</tt>

    \see  cayley_element sgn(cayley_element const &)
 */
inline  cayley_element &
cayley_element::sign_self()
{ return *this; }

/** Transforms this element to its conjugate

    \return  A reference to this element

    \post  <tt>*this == conj( old_this )</tt>

    \see  negated_cayley_element conj(negated_cayley_element const &)
 */
inline  negated_cayley_element &
negated_cayley_element::conjugate_self()
{ return this->basis() ? this->negate_self() : this->same_self(); }

/** \see  cayley_element::sign_self() */
inline  negated_cayley_element &
negated_cayley_element::sign_self()
{ return *this; }

/** Transforms this element to its reciprocal (\e i.e. multiplicative inverse)

    \return  A reference to this element

    \post  <tt>*this == reciprocal( old_this )</tt>

    \see  negated_cayley_element reciprocal(negated_cayley_element const &)
 */
inline  negated_cayley_element &
negated_cayley_element::reciprocate_self()
{ return this->conjugate_self(); }

/** \see  negated_cayley_element::conjugate_self() */
template < typename T >
inline  scaled_cayley_element<T> &
scaled_cayley_element<T>::conjugate_self()
{ return this->basis() ? this->negate_self() : this->same_self(); }

/** \see  cayley_element::sign_self() */
template < typename T >
scaled_cayley_element<T> &
scaled_cayley_element<T>::sign_self()
{
    if ( /*this->s_*/ this->s_ != T() )
    {
        this->s_ = static_cast<T>( (this->s_ < T()) ? -1 : +1 );
    }

    return *this;
}

/** \pre  <tt>*this != 0</tt>

    \see  negated_cayley_element::reciprocate_self()
 */
template < typename T >
scaled_cayley_element<T> &
scaled_cayley_element<T>::reciprocate_self()
{
    T const  ss = this->s_;

    this->s_ /= ss;  // separate in case "ss * ss" overflows
    this->s_ /= ss;

    return this->conjugate_self();
}


//  Cayley algebra elements scalar-multiplicative operator definitions  ------//

/** Multiplies a scalar by a scaled element.  The product has the same basis as
    the element, but its scale is the scalar multiplied by the element's scale.

    \param l  The multiplicand
    \param r  The multiplier

    \return  The product of \a l and \a r

    \relates  scaled_cayley_element
 */
template < typename T >
inline  scaled_cayley_element<T>
operator *( T const &l, scaled_cayley_element<T> const &r )
{ return scaled_cayley_element<T>( r.basis(), l * r.scale() ); }

/** Multiplies a scaled element by a scalar.  The product has the same basis as
    the element, but its scale is the element's scaled multiplied by the scalar.

    \param l  The multiplicand
    \param r  The multiplier

    \return  The product of \a l and \a r

    \relates  scaled_cayley_element
 */
template < typename T >
inline  scaled_cayley_element<T>
operator *( scaled_cayley_element<T> const &l, T const &r )
{ return scaled_cayley_element<T>( l.basis(), l.scale() * r ); }

/** Divides a scalar by a scaled element.  The quotient is the same as
    multiplying the scalar by the reciprocal of the element.

    \param l  The dividend
    \param r  The divisor

    \pre  <tt>\a r != 0</tt>

    \return  The quotient of \a l and \a r

    \relates  scaled_cayley_element
 */
template < typename T >
inline  scaled_cayley_element<T>
operator /( T const &l, scaled_cayley_element<T> const &r )
{ return l * reciprocal( r ); }

/** Divides a scaled element by a scalar.  The quotient has the same basis as
    the element, but its scale is the element's scaled divided by the scalar.

    \param l  The dividend
    \param r  The divisor

    \pre  <tt>\a r != 0</tt>

    \return  The quotient of \a l and \a r

    \relates  scaled_cayley_element
 */
template < typename T >
inline  scaled_cayley_element<T>
operator /( scaled_cayley_element<T> const &l, T const &r )
{ return scaled_cayley_element<T>( l.basis(), l.scale() / r ); }

/** \overload

    \relates  negated_cayley_element
    \see  scaled_cayley_element<T> operator *(T const &,
           scaled_cayley_element<T> const &)
 */
template < typename T >
inline  scaled_cayley_element<T>
operator *( T const &l, negated_cayley_element const &r )
{ return scaled_cayley_element<T>( r.basis(), r.negative() ? -l : +l ); }

/** \overload

    \relates  negated_cayley_element
    \see  scaled_cayley_element<T> operator *(scaled_cayley_element<T> const &,
           T const &)
 */
template < typename T >
inline  scaled_cayley_element<T>
operator *( negated_cayley_element const &l, T const &r )
{ return scaled_cayley_element<T>( l.basis(), l.negative() ? -r : +r ); }

/** \overload

    \relates  negated_cayley_element
    \see  scaled_cayley_element<T> operator /(T const &,
           scaled_cayley_element<T> const &)
 */
template < typename T >
inline  scaled_cayley_element<T>
operator /( T const &l, negated_cayley_element const &r )
{ return l * reciprocal( r ); }

/** \overload

    \pre  <tt>\a r != 0</tt>

    \relates  negated_cayley_element
    \see  scaled_cayley_element<T> operator /(scaled_cayley_element<T> const &,
           T const &)
 */
template < typename T >
inline  scaled_cayley_element<T>
operator /( negated_cayley_element const &l, T const &r )
{
    return scaled_cayley_element<T>( l.basis(), static_cast<T>(l.negative()
     ? -1 : +1) / r );
}

/** \overload

    \relates  cayley_element
    \see  scaled_cayley_element<T> operator *(T const &,
           scaled_cayley_element<T> const &)
 */
template < typename T >
inline  scaled_cayley_element<T>
operator *( T const &l, cayley_element const &r )
{ return scaled_cayley_element<T>( r.basis(), l ); }

/** \overload

    \relates  cayley_element
    \see  scaled_cayley_element<T> operator *(scaled_cayley_element<T> const &,
           T const &)
 */
template < typename T >
inline  scaled_cayley_element<T>
operator *( cayley_element const &l, T const &r )
{ return scaled_cayley_element<T>( l.basis(), r ); }

/** \overload

    \relates  cayley_element
    \see  scaled_cayley_element<T> operator /(T const &,
           scaled_cayley_element<T> const &)
 */
template < typename T >
inline  scaled_cayley_element<T>
operator /( T const &l, cayley_element const &r )
{ return l * reciprocal( r ); }

/** \overload

    \pre  <tt>\a r != 0</tt>

    \relates  cayley_element
    \see  scaled_cayley_element<T> operator /(scaled_cayley_element<T> const &,
           T const &)
 */
template < typename T >
inline  scaled_cayley_element<T>
operator /( cayley_element const &l, T const &r )
{ return scaled_cayley_element<T>( l.basis(), static_cast<T>(1) / r ); }


//  Cayley algebra elements multiplicative operator definitions  -------------//

/** Computes the product of two signed elements.  (This is the core routine
    instead of the basic element version because the product can be a negative
    in general, so incorporating the signs in the routine is faster.)  It uses
    an optimized breakdown of the Cayley-Dickson expansion rule for
    multiplication.

    \param l  The multiplicand
    \param r  The multiplier

    \return  The product

    \relates  negated_cayley_element
 */
negated_cayley_element
operator *( negated_cayley_element const &l, negated_cayley_element const &r )
{
    cayley_element::index_type  offset = 0u, li = l.basis(), ri = r.basis();
    bool                        is_neg = ( l.negative() != r.negative() );

    // The base case is simple real-scalar multiplication, but we
    // don't explicitly do it since every coefficient is 1, which
    // give coefficient products of 1.  Otherwise, we got rung-n
    // Cayley numbers.  Each one can be divided into 2 rung-(n-1)
    // halves.  The product is, with l=(ll;lh), r=(rl;rh), built
    // l * r = ( ll * rl - rh' * lh ; lh * rl'+ rh * ll )
    // where x' is the conjugate of x.  Since each factor has only
    // one non-zero component, exactly one of the 4 sub-products
    // is non-zero.  We proceed only with the non-zero sub-product
    // until we have two reals left.

    for ( cayley_element::rung_type n
     = cayley_element::minimum_rung_for_index(std::max( li, ri )) ; n ; --n )
    {
        cayley_element::index_type const  half_size = 1u << ( n - 1u );
        bool const                        lower_l = ( li < half_size );
        bool const                        lower_r = ( ri < half_size );

        // shift from upper half to lower half
        if ( !lower_r )  ri -= half_size;
        if ( !lower_l )  li -= half_size;

        // upper v. lower of result depends on halves used
        if ( lower_l != lower_r )  offset += half_size;

        // using upper-left half always conjugates the right-half
        if ( !lower_l && ri )  is_neg = !is_neg;

        // using upper-right half always reverses multiply order
        if ( !lower_r )  std::swap( li, ri );

        // one of the sub-products is negated
        if ( !lower_l && !lower_r )  is_neg = !is_neg;
    }

    return negated_cayley_element( offset, is_neg );
}

/** \overload

    \relatesalso  negated_cayley_element
    \relatesalso  cayley_element
    \see  negated_cayley_element operator *(negated_cayley_element const &,
           negated_cayley_element const &)
 */
inline  negated_cayley_element
operator *( negated_cayley_element const &l, cayley_element const &r )
{ return l * negated_cayley_element( r ); }

/** \overload

    \relatesalso  cayley_element
    \relatesalso  negated_cayley_element
    \see  negated_cayley_element operator *(negated_cayley_element const &,
           negated_cayley_element const &)
 */
inline  negated_cayley_element
operator *( cayley_element const &l, negated_cayley_element const &r )
{ return negated_cayley_element( l ) * r; }

/** \overload

    \relates  cayley_element
    \see  negated_cayley_element operator *(negated_cayley_element const &,
           negated_cayley_element const &)
 */
inline  negated_cayley_element
operator *( cayley_element const &l, cayley_element const &r )
{ return negated_cayley_element( l ) * negated_cayley_element( r ); }

/** \overload

    The product is the combination of the product of the scales and the
    products of the unscaled elements.  The scale of the final product is
    negated if the element multiplication indicates a negated result.

    \relates  scaled_cayley_element
    \see  negated_cayley_element operator *(negated_cayley_element const &,
           negated_cayley_element const &)
 */
template < typename T >
scaled_cayley_element<T>
operator *
(
    scaled_cayley_element<T> const &  l,
    scaled_cayley_element<T> const &  r
)
{
    T const                       s = l.scale() * r.scale();
    cayley_element const          ll( l.basis() ), rr( r.basis() );
    negated_cayley_element const  e = ll * rr;

    return scaled_cayley_element<T>( e.basis(), e.negative() ? -s : +s );
}

/** \overload

    \relatesalso  scaled_cayley_element
    \relatesalso  cayley_element
    \see  negated_cayley_element operator *(negated_cayley_element const &,
           negated_cayley_element const &)
 */
template < typename T >
inline  scaled_cayley_element<T>
operator *( scaled_cayley_element<T> const &l, cayley_element const &r )
{ return l * scaled_cayley_element<T>( r ); }

/** \overload

    \relatesalso  cayley_element
    \relatesalso  scaled_cayley_element
    \see  negated_cayley_element operator *(negated_cayley_element const &,
           negated_cayley_element const &)
 */
template < typename T >
inline  scaled_cayley_element<T>
operator *( cayley_element const &l, scaled_cayley_element<T> const &r )
{ return scaled_cayley_element<T>( l ) * r; }

/** \overload

    \relatesalso  scaled_cayley_element
    \relatesalso  negated_cayley_element
    \see  negated_cayley_element operator *(negated_cayley_element const &,
           negated_cayley_element const &)
 */
template < typename T >
inline  scaled_cayley_element<T>
operator *( scaled_cayley_element<T> const &l, negated_cayley_element const &r )
{ return l * scaled_cayley_element<T>( r ); }

/** \overload

    \relatesalso  negated_cayley_element
    \relatesalso  scaled_cayley_element
    \see  negated_cayley_element operator *(negated_cayley_element const &,
           negated_cayley_element const &)
 */
template < typename T >
inline  scaled_cayley_element<T>
operator *( negated_cayley_element const &l, scaled_cayley_element<T> const &r )
{ return scaled_cayley_element<T>( l ) * r; }

/** Computes the quotient of two basic elements.  It is the same as multiplying
    by the reciprocal of the divisor.

    \param l  The dividend
    \param r  The divisor

    \return  The quotient

    \relates  cayley_element
 */
inline  negated_cayley_element
operator /( cayley_element const &l, cayley_element const &r )
{ return l * reciprocal( r ); }

/** \overload

    \relatesalso  cayley_element
    \relatesalso  negated_cayley_element
    \see  cayley_element operator /(cayley_element const &,
           cayley_element const &)
 */
inline  negated_cayley_element
operator /( cayley_element const &l, negated_cayley_element const &r )
{ return l * reciprocal( r ); }

/** \overload

    \relatesalso  negated_cayley_element
    \relatesalso  cayley_element
    \see  cayley_element operator /(cayley_element const &,
           cayley_element const &)
 */
inline  negated_cayley_element
operator /( negated_cayley_element const &l, cayley_element const &r )
{ return l * reciprocal( r ); }

/** \overload

    \relates  negated_cayley_element
    \see  cayley_element operator /(cayley_element const &,
           cayley_element const &)
 */
inline  negated_cayley_element
operator /( negated_cayley_element const &l, negated_cayley_element const &r )
{ return l * reciprocal( r ); }

/** \overload

    \relates  scaled_cayley_element
    \see  cayley_element operator /(cayley_element const &,
           cayley_element const &)
 */
template < typename T >
inline  scaled_cayley_element<T>
operator /
(
    scaled_cayley_element<T> const &  l,
    scaled_cayley_element<T> const &  r
)
{ return l * reciprocal( r ); }

/** \overload

    \relatesalso  scaled_cayley_element
    \relatesalso  cayley_element
    \see  cayley_element operator /(cayley_element const &,
           cayley_element const &)
 */
template < typename T >
inline  scaled_cayley_element<T>
operator /( scaled_cayley_element<T> const &l, cayley_element const &r )
{ return l / scaled_cayley_element<T>( r ); }

/** \overload

    \relatesalso  cayley_element
    \relatesalso  scaled_cayley_element
    \see  cayley_element operator /(cayley_element const &,
           cayley_element const &)
 */
template < typename T >
inline  scaled_cayley_element<T>
operator /( cayley_element const &l, scaled_cayley_element<T> const &r )
{ return scaled_cayley_element<T>( l ) / r; }

/** \overload

    \relatesalso  scaled_cayley_element
    \relatesalso  negated_cayley_element
    \see  cayley_element operator /(cayley_element const &,
           cayley_element const &)
 */
template < typename T >
inline  scaled_cayley_element<T>
operator /( scaled_cayley_element<T> const &l, negated_cayley_element const &r )
{ return l / scaled_cayley_element<T>( r ); }

/** \overload

    \relatesalso  negated_cayley_element
    \relatesalso  scaled_cayley_element
    \see  cayley_element operator /(cayley_element const &,
           cayley_element const &)
 */
template < typename T >
inline  scaled_cayley_element<T>
operator /( negated_cayley_element const &l, scaled_cayley_element<T> const &r )
{ return scaled_cayley_element<T>( l ) / r; }


//  Cayley algebra elements multiplicative member operator definitions  ------//

/** Multiplies this element by another.  This element is the multiplicand and
    later stores the product.  (Order matters for element multiplication.)

    \param r  The multiplier

    \post  <tt>*this == old_this * \a r</tt>

    \return  A reference to this element
 */
inline  negated_cayley_element &
negated_cayley_element::operator *=( negated_cayley_element const &r )
{ return *this = *this * r; }

/** Divides this element by another.  This element is the dividend and later
    stores the quotient.

    \param r  The divisor

    \post  <tt>*this == old_this / \a r</tt>

    \return  A reference to this element
 */
inline  negated_cayley_element &
negated_cayley_element::operator /=( negated_cayley_element const &r )
{ return *this = *this / r; }

/** Multiplies this element by a real scalar.  This element is the multiplicand
    and later stores the product.

    \param r  The multiplier

    \post  <tt>*this == old_this * \a r</tt>

    \return  A reference to this element

    \see  value_type
 */
template < typename T >
inline  scaled_cayley_element<T> &
scaled_cayley_element<T>::operator *=( value_type const &r )
{ this->s_ *= r; return *this; }

/** Divides this element by a real scalar.  This element is the dividend and
    later stores the quotient.

    \param r  The divisor

    \pre  <tt>\a r != 0</tt>

    \post  <tt>*this == old_this / \a r</tt>

    \return  A reference to this element

    \see  value_type
 */
template < typename T >
inline  scaled_cayley_element<T> &
scaled_cayley_element<T>::operator /=( value_type const &r )
{ this->s_ /= r; return *this; }

/** \see  negated_cayley_element::operator *=(negated_cayley_element const &) */
template < typename T >
inline  scaled_cayley_element<T> &
scaled_cayley_element<T>::operator *=( self_type const &r )
{ return *this = *this * r; }

/** \see  negated_cayley_element::operator /=(negated_cayley_element const &)

    \pre  <tt>\a r != 0</tt>
 */
template < typename T >
inline  scaled_cayley_element<T> &
scaled_cayley_element<T>::operator /=( self_type const &r )
{ return *this = *this / r; }


//  Cayley algebra elements shift operator definitions  ----------------------//

/** Copies \a l such that the new basis is up-shifted by \a r.

    \param l  The element to be shifted
    \param r  The shift amount (must be nonnegative)

    \return  A new element \a x such that <tt>x.basis() == l.basis() + r</tt>

    \relates  cayley_element
 */
inline  cayley_element
operator <<( cayley_element const &l, std::size_t r )
{ return cayley_element( l.basis() + r ); }

/** Copies \a l such that the new basis is down-shifted by \a r.

    \pre  <tt>r \<= l.basis()</tt>

    \param l  The element to be shifted
    \param r  The shift amount (must be nonnegative)

    \return  A new element \a x such that <tt>x.basis() == l.basis() - r</tt>

    \relates  cayley_element
 */
inline  cayley_element
operator >>( cayley_element const &l, std::size_t r )
{ return cayley_element( l.basis() - r ); }

/** \overload

    \relates  negated_cayley_element
    \see  cayley_element operator <<(cayley_element const &, std::size_t)
 */
inline  negated_cayley_element
operator <<( negated_cayley_element const &l, std::size_t r )
{ return negated_cayley_element( l.basis() + r, l.negative() ); }

/** \overload

    \relates  negated_cayley_element
    \see  cayley_element operator >>(cayley_element const &, std::size_t)
 */
inline  negated_cayley_element
operator >>( negated_cayley_element const &l, std::size_t r )
{ return negated_cayley_element( l.basis() - r, l.negative() ); }

/** \overload

    \relates  scaled_cayley_element
    \see  cayley_element operator <<(cayley_element const &, std::size_t)
 */
template < typename T >
inline  scaled_cayley_element<T>
operator <<( scaled_cayley_element<T> const &l, std::size_t r )
{ return scaled_cayley_element<T>( l.basis() + r, l.scale() ); }

/** \overload

    \relates  scaled_cayley_element
    \see  cayley_element operator >>(cayley_element const &, std::size_t)
 */
template < typename T >
inline  scaled_cayley_element<T>
operator >>( scaled_cayley_element<T> const &l, std::size_t r )
{ return scaled_cayley_element<T>( l.basis() - r, l.scale() ); }


//  Cayley algebra elements shift member operator definitions  ---------------//

/** Changes the element's basis by increasing the basis index.

    \param r  The desired index increase

    \pre  <tt>this->basis() + \a r
           \<= std::numeric_limits\<std::size_t\>::max()</tt>

    \post  <tt>this->basis() == old_this.basis() + \a r</tt>

    \return  A reference to this element
 */
inline  cayley_element &
cayley_element::operator <<=( std::size_t r )
{ this->b_ += r; return *this; }

/** Changes the element's basis by decreasing the basis index.

    \param r  The desired index decrease

    \pre  <tt>this->basis() >= \a r</tt>

    \post  <tt>this->basis() == old_this.basis() - \a r</tt>

    \return  A reference to this element
 */
inline  cayley_element &
cayley_element::operator >>=( std::size_t r )
{ this->b_ -= r; return *this; }

/** \see  cayley_element::operator<<=(std::size_t) */
inline  negated_cayley_element &
negated_cayley_element::operator <<=( std::size_t r )
{ this->e_ <<= r; return *this; }

/** \see  cayley_element::operator>>=(std::size_t) */
inline  negated_cayley_element &
negated_cayley_element::operator >>=( std::size_t r )
{ this->e_ >>= r; return *this; }

/** \see  cayley_element::operator<<=(std::size_t) */
template < typename T >
inline  scaled_cayley_element<T> &
scaled_cayley_element<T>::operator <<=( std::size_t r )
{ this->e_ <<= r; return *this; }

/** \see  cayley_element::operator>>=(std::size_t) */
template < typename T >
inline  scaled_cayley_element<T> &
scaled_cayley_element<T>::operator >>=( std::size_t r )
{ this->e_ >>= r; return *this; }


//  Cayley algebra elements component function definitions  ------------------//

/** Since a unit element is purely real or purely unreal, this function returns
    either 1 or 0.

    \param x  The element to be analyzed

    \return  \f$x_0 = \Re(x) = x \cdot 1 = \frac{x + \bar{x}}{2}\f$

    \relates  cayley_element
 */
inline  int
real( cayley_element const &x )
{ return !x.basis(); }

/** \overload
    Now signs are considered, the result may be -1, 0, or +1.

    \relates  negated_cayley_element
    \see  int real(cayley_element const &)
 */
inline  int
real( negated_cayley_element const &x )
{ return x.basis() ? 0 : ( x.negative() ? -1 : +1 ); }

/** \overload
    Now the scale is the coefficent, so the result is the scale or zero.

    \relates  scaled_cayley_element
    \see  int real(cayley_element const &)
 */
template < typename T >
inline  T
real( scaled_cayley_element<T> const &x )
{ return x.basis() ? T() : x.scale(); }

/** Similar to real(\a x), this function returns either 1 or 0.

    \param x  The element to be analyzed

    \return  \f$x_1 = \Im(x) = x \cdot \imath\f$

    \relates  cayley_element
    \see  int real(cayley_element const &)
 */
inline  int
imag( cayley_element const &x )
{ return x.basis() == 1u; }

/** \overload

    \relates  negated_cayley_element
    \see  int imag(cayley_element const &)
 */
inline  int
imag( negated_cayley_element const &x )
{ return ( x.basis() == 1u ) ? ( x.negative() ? -1 : +1 ) : 0; }

/** \overload

    \relates  scaled_cayley_element
    \see  int imag(cayley_element const &)
 */
template < typename T >
inline  T
imag( scaled_cayley_element<T> const &x )
{ return ( x.basis() == 1u ) ? x.scale() : T(); }

/** Since an element is purely real or purely unreal, this function returns
    either zero or \a x.

    \param x  The element to be analyzed

    \return  \f$\mathfrak{Ur}(x) = \frac{x - \bar{x}}{2}\f$

    \relates  scaled_cayley_element
 */
template < typename T >
inline  scaled_cayley_element<T>
unreal( scaled_cayley_element<T> const &x )
{ return x.basis() ? x : scaled_cayley_element<T>(); }

/** \overload

    \relates  negated_cayley_element
    \see  scaled_cayley_element<T> unreal(scaled_cayley_element<T> const &)
 */
inline  scaled_cayley_element<int>
unreal( negated_cayley_element const &x )
{ return unreal( static_cast< scaled_cayley_element<int> >(x) ); }

/** \overload

    \relates  cayley_element
    \see  scaled_cayley_element<T> unreal(scaled_cayley_element<T> const &)
 */
inline  scaled_cayley_element<int>
unreal( cayley_element const &x )
{ return unreal( static_cast< scaled_cayley_element<int> >(x) ); }


//  Cayley algebra elements integer-power function definitions  --------------//

/** Raises an element to an integer power, only supporting "int" as the
    exponent's type.  Use this function in piecemeal for integer exponents
    that require a type with an expanded range than "int".  Powers of 1 stay 1;
    powers of -1 alternate between +1 and -1; and powers of unreal unit elements
    alternate sign and alternate basis with a cycle of 4.

    \param b  The element to be exponentiated
    \param e  The exponent

    \return  \f$b^e\f$

    \relates  negated_cayley_element
 */
negated_cayley_element
pow( negated_cayley_element const &b, int e )
{
    bool const  final_negate = ( b.negative() && (e % 2) );

    if ( cayley_element::index_type const  n = b.basis() )
    {
        bool const              s = ( e < 0 );
        unsigned const          a = static_cast<unsigned>( s ? -e : +e );
        negated_cayley_element  r( (a & 1u) ? n : 0u, a & 2u );

        if ( s )
        {
            r.reciprocate_self();
        }

        return final_negate ? -r : r;
    }
    else
    {
        // have +1 or -1
        return negated_cayley_element( 0u, final_negate );
    }
}

/** \overload
    Without signs, there are fewer cycles to consider.

    \relates  cayley_element
    \see  negated_cayley_element pow(negated_cayley_element const &, int)
 */
inline  negated_cayley_element
pow( cayley_element const &b, int e )
{
    return pow( negated_cayley_element(b), e );
}

/** \overload
    The basis and sign cycle as before, but a non-unit or non-zero scale
    strictly changes as \f$|b|^e\f$ (without cycles).

    \relates  scaled_cayley_element
    \see  negated_cayley_element pow(negated_cayley_element const &, int)
 */
template < typename T >
inline  scaled_cayley_element<T>
pow( scaled_cayley_element<T> const &b, int e )
{
    using std::pow;

    return pow( b.scale(), e ) * pow( cayley_element(b.basis()), e );
}


}  // namespace math
}  // namespace boost


#endif  // BOOST_MATH_CAYLEY_ELEMENT_HPP
