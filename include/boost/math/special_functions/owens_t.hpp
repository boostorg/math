// (C) Benjamin Sobotta 2012

//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_OWENS_T_HPP
#define BOOST_OWENS_T_HPP

// Reference:
// Mike Patefield, David Tandy
// FAST AND ACCURATE CALCULATION OF OWEN'S T-FUNCTION
// Journal of Statistical Software, 5 (5), 1-25

#ifdef _MSC_VER
#  pragma once
#endif

#include <boost/config/no_tr1/cmath.hpp>
#include <boost/math/special_functions/erf.hpp>
#include <boost/math/special_functions/expm1.hpp>
#include <boost/throw_exception.hpp>
#include <boost/assert.hpp>
#include <boost/math/constants/constants.hpp>

#include <stdexcept>

namespace boost
{
   namespace math
   {
      namespace detail
      {
         // owens_t_znorm1(x) = P(-oo<Z<=x)-0.5 with Z being normally distributed.
         template<typename RealType>
         inline RealType owens_t_znorm1(const RealType x)
         {
            using namespace boost::math::constants;
            return erf(x/root_two<RealType>())*half<RealType>();
         } // RealType owens_t_znorm1(const RealType x)

         // owens_t_znorm2(x) = P(x<=Z<oo) with Z being normally distributed.
         template<typename RealType>
         inline RealType owens_t_znorm2(const RealType x)
         {
            using namespace boost::math::constants;
            return erfc(x/root_two<RealType>())*half<RealType>();
         } // RealType owens_t_znorm2(const RealType x)

         // Auxiliary function, it computes an array key that is used to determine
         // the specific computation method for Owen's T and the order thereof
         // used in owens_t_dispatch.
         template<typename RealType>
         inline unsigned short owens_t_compute_code(const RealType h, const RealType a)
         {
            static const RealType hrange[] =
            {0.02, 0.06, 0.09, 0.125, 0.26, 0.4,  0.6,  1.6,  1.7,  2.33,  2.4,  3.36, 3.4,  4.8};

            static const RealType arange[] = {0.025, 0.09, 0.15, 0.36, 0.5, 0.9, 0.99999};
            /*
            original select array from paper:
            1, 1, 2,13,13,13,13,13,13,13,13,16,16,16, 9
            1, 2, 2, 3, 3, 5, 5,14,14,15,15,16,16,16, 9
            2, 2, 3, 3, 3, 5, 5,15,15,15,15,16,16,16,10
            2, 2, 3, 5, 5, 5, 5, 7, 7,16,16,16,16,16,10
            2, 3, 3, 5, 5, 6, 6, 8, 8,17,17,17,12,12,11
            2, 3, 5, 5, 5, 6, 6, 8, 8,17,17,17,12,12,12
            2, 3, 4, 4, 6, 6, 8, 8,17,17,17,17,17,12,12
            2, 3, 4, 4, 6, 6,18,18,18,18,17,17,17,12,12
            */                  
            // subtract one because the array is written in FORTRAN in mind - in C arrays start @ zero
            static const unsigned short select[] =
            {
               0,    0 ,   1  , 12   ,12 ,  12  , 12  , 12 ,  12  , 12  , 12  , 15  , 15 ,  15  ,  8,
               0  ,  1  ,  1   , 2 ,   2   , 4  ,  4  , 13 ,  13  , 14  , 14 ,  15  , 15  , 15  ,  8,
               1  ,  1   , 2 ,   2  ,  2  ,  4   , 4  , 14  , 14 ,  14  , 14 ,  15  , 15 ,  15  ,  9,
               1  ,  1   , 2 ,   4  ,  4  ,  4   , 4  ,  6  ,  6 ,  15  , 15 ,  15 ,  15 ,  15  ,  9,
               1  ,  2   , 2  ,  4  ,  4  ,  5   , 5  ,  7  ,  7  , 16   ,16 ,  16 ,  11 ,  11 ,  10,
               1  ,  2   , 4  ,  4   , 4  ,  5   , 5  ,  7  ,  7  , 16  , 16 ,  16 ,  11  , 11 ,  11,
               1  ,  2   , 3  ,  3  ,  5  ,  5   , 7  ,  7  , 16 ,  16  , 16 ,  16 ,  16  , 11 ,  11,
               1  ,  2   , 3   , 3   , 5  ,  5 ,  17  , 17  , 17 ,  17  , 16 ,  16 ,  16 ,  11 ,  11
            };

            unsigned short ihint = 14, iaint = 7;
            for(unsigned short i = 0; i != 14; i++)
            {
               if( h <= hrange[i] )
               {
                  ihint = i;
                  break;
               }
            } // for(unsigned short i = 0; i != 14; i++)

            for(unsigned short i = 0; i != 7; i++)
            {
               if( a <= arange[i] )
               {
                  iaint = i;
                  break;
               }
            } // for(unsigned short i = 0; i != 7; i++)

            // interprete select array as 8x15 matrix
            return select[iaint*15 + ihint];

         } // unsigned short owens_t_compute_code(const RealType h, const RealType a)

         // compute the value of Owen's T function with method T1 from the reference paper
         template<typename RealType>
         inline RealType owens_t_T1(const RealType h, const RealType a, const unsigned short m)
         {
            BOOST_MATH_STD_USING
            using namespace boost::math::constants;

            const RealType hs = -h*h*half<RealType>();
            const RealType dhs = exp( hs );
            const RealType as = a*a;

            unsigned short j=1;
            RealType jj = 1;
            RealType aj = a / (static_cast<RealType>(2)*pi<RealType>());
            RealType dj = expm1( hs );
            RealType gj = hs*dhs;

            RealType val = atan( a ) / (static_cast<RealType>(2)*pi<RealType>());

            while( true )
            {
               val += dj*aj/jj;

               if( m <= j )
                  break;

               j++;
               jj += static_cast<RealType>(2);
               aj *= as;
               dj = gj - dj;
               gj *= hs / static_cast<RealType>(j);
            } // while( true )

            return val;
         } // RealType owens_t_T1(const RealType h, const RealType a, const unsigned short m)

         // compute the value of Owen's T function with method T2 from the reference paper
         template<typename RealType>
         inline RealType owens_t_T2(const RealType h, const RealType a, const unsigned short m, const RealType ah)
         {
            BOOST_MATH_STD_USING
            using namespace boost::math::constants;

            const unsigned short maxii = m+m+1;
            const RealType hs = h*h;
            const RealType as = -a*a;
            const RealType y = static_cast<RealType>(1) / hs;

            unsigned short ii = 1;
            RealType val = 0;
            RealType vi = a * exp( -ah*ah*half<RealType>() ) / root_two_pi<RealType>();
            RealType z = owens_t_znorm1(ah)/h;

            while( true )
            {
               val += z;
               if( maxii <= ii )
               {
                  val *= exp( -hs*half<RealType>() ) / root_two_pi<RealType>();
                  break;
               } // if( maxii <= ii )
               z = y * ( vi - static_cast<RealType>(ii) * z );
               vi *= as;
               ii += 2;
            } // while( true )

            return val;
         } // RealType owens_t_T2(const RealType h, const RealType a, const unsigned short m, const RealType ah)

         // compute the value of Owen's T function with method T3 from the reference paper
         template<typename RealType>
         inline RealType owens_t_T3(const RealType h, const RealType a, const unsigned short m, const RealType ah)
         {
            BOOST_MATH_STD_USING
            using namespace boost::math::constants;

            static const RealType c2[] =
            {
               0.99999999999999987510,
               -0.99999999999988796462,      0.99999999998290743652,
               -0.99999999896282500134,      0.99999996660459362918,
               -0.99999933986272476760,      0.99999125611136965852,
               -0.99991777624463387686,      0.99942835555870132569,
               -0.99697311720723000295,      0.98751448037275303682,
               -0.95915857980572882813,      0.89246305511006708555,
               -0.76893425990463999675,      0.58893528468484693250,
               -0.38380345160440256652,      0.20317601701045299653,
               -0.82813631607004984866E-01,  0.24167984735759576523E-01,
               -0.44676566663971825242E-02,  0.39141169402373836468E-03
            };

            const RealType as = a*a;
            const RealType hs = h*h;
            const RealType y = static_cast<RealType>(1)/hs;

            RealType ii = 1;
            unsigned short i = 0;
            RealType vi = a * exp( -ah*ah*half<RealType>() ) / root_two_pi<RealType>();
            RealType zi = owens_t_znorm1(ah)/h;
            RealType val = 0;

            while( true )
            {
               BOOST_ASSERT(i < 21);
               val += zi*c2[i];
               if( m <= i ) // if( m < i+1 )
               {
                  val *= exp( -hs*half<RealType>() ) / root_two_pi<RealType>();
                  break;
               } // if( m < i )
               zi = y * (ii*zi - vi);
               vi *= as;
               ii += 2;
               i++;
            } // while( true )

            return val;
         } // RealType owens_t_T3(const RealType h, const RealType a, const unsigned short m, const RealType ah)

         // compute the value of Owen's T function with method T4 from the reference paper
         template<typename RealType>
         inline RealType owens_t_T4(const RealType h, const RealType a, const unsigned short m)
         {
            BOOST_MATH_STD_USING
            using namespace boost::math::constants;

            const unsigned short maxii = m+m+1;
            const RealType hs = h*h;
            const RealType as = -a*a;

            unsigned short ii = 1;
            RealType ai = a * exp( -hs*(static_cast<RealType>(1)-as)*half<RealType>() ) / (static_cast<RealType>(2)*pi<RealType>());
            RealType yi = 1;
            RealType val = 0;

            while( true )
            {
               val += ai*yi;
               if( maxii <= ii )
                  break;
               ii += 2;
               yi = (static_cast<RealType>(1)-hs*yi) / static_cast<RealType>(ii);
               ai *= as;
            } // while( true )

            return val;
         } // RealType owens_t_T4(const RealType h, const RealType a, const unsigned short m)

         // compute the value of Owen's T function with method T5 from the reference paper
         template<typename RealType>
         inline RealType owens_t_T5(const RealType h, const RealType a, const unsigned short m)
         {
            BOOST_MATH_STD_USING
            static const RealType pts[] = {0.35082039676451715489E-02,
               0.31279042338030753740E-01,  0.85266826283219451090E-01,
               0.16245071730812277011,      0.25851196049125434828,
               0.36807553840697533536,      0.48501092905604697475,
               0.60277514152618576821,      0.71477884217753226516,
               0.81475510988760098605,      0.89711029755948965867,
               0.95723808085944261843,      0.99178832974629703586};
            static const RealType wts[] = { 0.18831438115323502887E-01,
               0.18567086243977649478E-01,  0.18042093461223385584E-01,
               0.17263829606398753364E-01,  0.16243219975989856730E-01,
               0.14994592034116704829E-01,  0.13535474469662088392E-01,
               0.11886351605820165233E-01,  0.10070377242777431897E-01,
               0.81130545742299586629E-02,  0.60419009528470238773E-02,
               0.38862217010742057883E-02,  0.16793031084546090448E-02};

            const RealType as = a*a;
            const RealType hs = -h*h*boost::math::constants::half<RealType>();

            RealType val = 0;
            for(unsigned short i = 0; i < m; ++i)
            {
               BOOST_ASSERT(i < 13);
               const RealType r = static_cast<RealType>(1) + as*pts[i];
               val += wts[i] * exp( hs*r ) / r;
            } // for(unsigned short i = 0; i < m; ++i)

            return val*a;
         } // RealType owens_t_T5(const RealType h, const RealType a, const unsigned short m)

         // compute the value of Owen's T function with method T6 from the reference paper
         template<typename RealType>
         inline RealType owens_t_T6(const RealType h, const RealType a)
         {
            BOOST_MATH_STD_USING
            using namespace boost::math::constants;

            const RealType normh = owens_t_znorm2( h );
            const RealType y = static_cast<RealType>(1) - a;
            const RealType r = atan2(y, ( static_cast<RealType>(1) + a ) );

            RealType val = normh * ( static_cast<RealType>(1) - normh ) * half<RealType>();

            if( r != 0 )
               val -= r * exp( -y*h*h*half<RealType>()/r ) / (static_cast<RealType>(2)*pi<RealType>());

            return val;
         } // RealType owens_t_T6(const RealType h, const RealType a, const unsigned short m)

         // This routine dispatches the call to one of six subroutines, depending on the values
         // of h and a.
         // preconditions: h >= 0, 0<=a<=1, ah=a*h
         template<typename RealType>
         inline RealType owens_t_dispatch(const RealType h, const RealType a, const RealType ah)
         {
            static const unsigned short meth[] = {1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 3, 4, 4, 4, 4, 5, 6}; // 18 entries
            static const unsigned short ord[] = {2, 3, 4, 5, 7, 10, 12, 18, 10, 20, 30, 20, 4, 7, 8, 20, 13, 0}; // 18 entries

            RealType val = 0; // avoid compiler warnings, 0 will be overwritten in any case

            const unsigned short icode = owens_t_compute_code(h,a);
            const unsigned short m = ord[icode];

            // determine the appropriate method, T1 ... T6
            switch( meth[icode] )
            {
            case 1: // T1
               val = owens_t_T1(h,a,m);
               break;
            case 2: // T2
               val = owens_t_T2(h,a,m,ah);
               break;
            case 3: // T3
               val = owens_t_T3(h,a,m,ah);
               break;
            case 4: // T4
               val = owens_t_T4(h,a,m);
               break;
            case 5: // T5
               val = owens_t_T5(h,a,m);
               break;
            case 6: // T6
               val = owens_t_T6(h,a);
               break;
            default:
               BOOST_THROW_EXCEPTION(std::logic_error("selection routine in Owen's T function failed"));
            }
            return val;
         } // RealType owens_t_dispatch(RealType h, RealType a, RealType ah)

         // compute Owen's T function, T(h,a), for arbitrary values of h and a
         template<typename RealType, class Policy>
         inline RealType owens_t(RealType h, RealType a, const Policy&)
         {
            BOOST_MATH_STD_USING
            // exploit that T(-h,a) == T(h,a)
            h = fabs(h);

            // Use equation (2) in the paper to remap the arguments
            // such that h>=0 and 0<=a<=1 for the call of the actual
            // computation routine.

            const RealType fabs_a = fabs(a);
            const RealType fabs_ah = fabs_a*h;

            RealType val = 0.0; // avoid compiler warnings, 0.0 will be overwritten in any case

            if(fabs_a <= 1)
            {
               val = detail::owens_t_dispatch(h, fabs_a, fabs_ah);
            } // if(fabs_a <= 1.0)
            else 
            {
               if( h <= 0.67 )
               {
                  const RealType normh = detail::owens_t_znorm1(h);
                  const RealType normah = detail::owens_t_znorm1(fabs_ah);
                  val = static_cast<RealType>(1)/static_cast<RealType>(4) - normh*normah -
                     detail::owens_t_dispatch(fabs_ah, static_cast<RealType>(1) / fabs_a, h);
               } // if( h <= 0.67 )
               else
               {
                  const RealType normh = detail::owens_t_znorm2(h);
                  const RealType normah = detail::owens_t_znorm2(fabs_ah);
                  val = constants::half<RealType>()*(normh+normah) - normh*normah -
                     detail::owens_t_dispatch(fabs_ah, static_cast<RealType>(1) / fabs_a, h);
               } // else [if( h <= 0.67 )]
            } // else [if(fabs_a <= 1)]

            // exploit that T(h,-a) == -T(h,a)
            if(a < 0)
            {
               return -val;
            } // if(a < 0)

            return val;
         } // RealType owens_t(RealType h, RealType a)

      } // namespace detail

      template <class T1, class T2, class Policy>
      inline typename tools::promote_args<T1, T2>::type owens_t(T1 h, T2 a, const Policy& pol)
      {
         typedef typename tools::promote_args<T1, T2>::type result_type;
         typedef typename policies::evaluation<result_type, Policy>::type value_type;
         return policies::checked_narrowing_cast<result_type, Policy>(detail::owens_t(static_cast<value_type>(h), static_cast<value_type>(a), pol), "boost::math::ellint_2<%1%>(%1%,%1%)");
      }

      template <class T1, class T2>
      inline typename tools::promote_args<T1, T2>::type owens_t(T1 h, T2 a)
      {
         return owens_t(h, a, policies::policy<>());
      }


   } // namespace math
} // namespace boost

#endif
// EOF
