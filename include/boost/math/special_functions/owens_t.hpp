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
            return erf(x*one_div_root_two<RealType>())*half<RealType>();
         } // RealType owens_t_znorm1(const RealType x)

         // owens_t_znorm2(x) = P(x<=Z<oo) with Z being normally distributed.
         template<typename RealType>
         inline RealType owens_t_znorm2(const RealType x)
         {
            using namespace boost::math::constants;
            return erfc(x*one_div_root_two<RealType>())*half<RealType>();
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

         template<typename RealType>
         inline unsigned short owens_t_get_order(const unsigned short icode, RealType)
         {
            static const unsigned short ord[] = {2, 3, 4, 5, 7, 10, 12, 18, 10, 20, 30, 0, 4, 7, 8, 20, 0, 0}; // 18 entries

            BOOST_ASSERT(icode<18);

            return ord[icode];
         } // unsigned short owens_t_get_order(const unsigned short icode, RealType)

#ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS

        template<>
        inline unsigned short owens_t_get_order(const unsigned short icode, long double)
        {
           // method ================>>>       {1, 1, 1, 1, 1,  1,  1,  1,  2,  2,  2,  3, 4,  4,  4,  4,  5, 6}
           static const unsigned short ord[] = {3, 4, 5, 6, 8, 11, 13, 19, 10, 20, 30,  0, 7, 10, 11, 23,  0, 0}; // 18 entries

          BOOST_ASSERT(icode<18);

          return ord[icode];
        } // unsigned short owens_t_get_order(const unsigned short icode, long double)

#endif // BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS

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
            RealType aj = a * one_div_two_pi<RealType>();
            RealType dj = expm1( hs );
            RealType gj = hs*dhs;

            RealType val = atan( a ) * one_div_two_pi<RealType>();

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
            RealType vi = a * exp( -ah*ah*half<RealType>() ) * one_div_root_two_pi<RealType>();
            RealType z = owens_t_znorm1(ah)/h;

            while( true )
            {
               val += z;
               if( maxii <= ii )
               {
                  val *= exp( -hs*half<RealType>() ) * one_div_root_two_pi<RealType>();
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
         inline RealType owens_t_T3(const RealType h, const RealType a, const RealType ah)
         {
            BOOST_MATH_STD_USING
            using namespace boost::math::constants;

	    const unsigned short m = 20;

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
            RealType vi = a * exp( -ah*ah*half<RealType>() ) * one_div_root_two_pi<RealType>();
            RealType zi = owens_t_znorm1(ah)/h;
            RealType val = 0;

            while( true )
            {
               BOOST_ASSERT(i < 21);
               val += zi*c2[i];
               if( m <= i ) // if( m < i+1 )
               {
                  val *= exp( -hs*half<RealType>() ) * one_div_root_two_pi<RealType>();
                  break;
               } // if( m < i )
               zi = y * (ii*zi - vi);
               vi *= as;
               ii += 2;
               i++;
            } // while( true )

            return val;
         } // RealType owens_t_T3(const RealType h, const RealType a, const RealType ah)

#ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS

        // compute the value of Owen's T function with method T3 from the reference paper
        template<>
        inline long double owens_t_T3(const long double h, const long double a, const long double ah)
        {
          BOOST_MATH_STD_USING
          using namespace boost::math::constants;
          
          const unsigned short m = 30;

          static const long double c2[] =
          {
              0.99999999999999999999999729978162447266851932041876728736094298092917625009873l,
            -0.99999999999999999999467056379678391810626533251885323416799874878563998732905968l,
              0.99999999999999999824849349313270659391127814689133077036298754586814091034842536l,
            -0.9999999999999997703859616213643405880166422891953033591551179153879839440241685l,
              0.99999999999998394883415238173334565554173013941245103172035286759201504179038147l,
            -0.9999999999993063616095509371081203145247992197457263066869044528823599399470977l,
              0.9999999999797336340409464429599229870590160411238245275855903767652432017766116267l,
            -0.999999999574958412069046680119051639753412378037565521359444170241346845522403274l,
              0.9999999933226234193375324943920160947158239076786103108097456617750134812033362048l,
            -0.9999999188923242461073033481053037468263536806742737922476636768006622772762168467l,
              0.9999992195143483674402853783549420883055129680082932629160081128947764415749728967l,
            -0.999993935137206712830997921913316971472227199741857386575097250553105958772041501l,
              0.99996135597690552745362392866517133091672395614263398912807169603795088421057688716l,
            -0.99979556366513946026406788969630293820987757758641211293079784585126692672425362469l,
              0.999092789629617100153486251423850590051366661947344315423226082520411961968929483l,
            -0.996593837411918202119308620432614600338157335862888580671450938858935084316004769854l,
              0.98910017138386127038463510314625339359073956513420458166238478926511821146316469589567l,
            -0.970078558040693314521331982203762771512160168582494513347846407314584943870399016019l,
              0.92911438683263187495758525500033707204091967947532160289872782771388170647150321633673l,
            -0.8542058695956156057286980736842905011429254735181323743367879525470479126968822863l,
              0.73796526033030091233118357742803709382964420335559408722681794195743240930748630755l,
            -0.58523469882837394570128599003785154144164680587615878645171632791404210655891158l,
              0.415997776145676306165661663581868460503874205343014196580122174949645271353372263l,
            -0.2588210875241943574388730510317252236407805082485246378222935376279663808416534365l,
              0.1375535825163892648504646951500265585055789019410617565727090346559210218472356689l,
            -0.0607952766325955730493900985022020434830339794955745989150270485056436844239206648l,
              0.0216337683299871528059836483840390514275488679530797294557060229266785853764115l,
            -0.00593405693455186729876995814181203900550014220428843483927218267309209471516256l,
              0.0011743414818332946510474576182739210553333860106811865963485870668929503649964142l,
            -1.489155613350368934073453260689881330166342484405529981510694514036264969925132e-4l,
              9.072354320794357587710929507988814669454281514268844884841547607134260303118208e-6l
          };

          const long double as = a*a;
          const long double hs = h*h;
          const long double y = 1.l/hs;

          long double ii = 1;
          unsigned short i = 0;
          long double vi = a * exp( -ah*ah*half<long double>() ) * one_div_root_two_pi<long double>();
          long double zi = owens_t_znorm1(ah)/h;
          long double val = 0;

          while( true )
          {
              BOOST_ASSERT(i < 31);
              val += zi*c2[i];
              if( m <= i ) // if( m < i+1 )
              {
                val *= exp( -hs*half<long double>() ) * one_div_root_two_pi<long double>();
                break;
              } // if( m < i )
              zi = y * (ii*zi - vi);
              vi *= as;
              ii += 2;
              i++;
          } // while( true )

          return val;
        } // long double owens_t_T3(const long double h, const long double a, const long double ah)

#endif // BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS

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
            RealType ai = a * exp( -hs*(static_cast<RealType>(1)-as)*half<RealType>() ) * one_div_two_pi<RealType>();
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
         inline RealType owens_t_T5(const RealType h, const RealType a)
         {
            BOOST_MATH_STD_USING
            /*
               NOTICE:
               - The pts[] array contains the squares (!) of the abscissas, i.e. the roots of the Legendre
                 polynomial P_n(x), instead of the plain roots as required in Gauss-Legendre
                 quadrature, because T5(h,a,m) contains only x^2 terms.
               - The wts[] array contains the weights for Gauss-Legendre quadrature scaled with a factor
                 of 1/(2*pi) according to T5(h,a,m).
             */

            const unsigned short m = 13;
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
         } // RealType owens_t_T5(const RealType h, const RealType a)

#ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS

        // compute the value of Owen's T function with method T5 from the reference paper
        template<>
        inline long double owens_t_T5(const long double h, const long double a)
        {
          BOOST_MATH_STD_USING
            /*
              NOTICE:
              - The pts[] array contains the squares (!) of the abscissas, i.e. the roots of the Legendre
              polynomial P_n(x), instead of the plain roots as required in Gauss-Legendre
              quadrature, because T5(h,a,m) contains only x^2 terms.
              - The wts[] array contains the weights for Gauss-Legendre quadrature scaled with a factor
              of 1/(2*pi) according to T5(h,a,m).
            */

          const unsigned short m = 19;
          static const long double pts[] = {0.0016634282895983227941l,
                                            0.014904509242697054183l,
                                            0.04103478879005817919l,
                                            0.079359853513391511008l,
                                            0.1288612130237615133l,
                                            0.18822336642448518856l,
                                            0.25586876186122962384l,
                                            0.32999972011807857222l,
                                            0.40864620815774761438l,
                                            0.48971819306044782365l,
                                            0.57106118513245543894l,
                                            0.6505134942981533829l,
                                            0.72596367859928091618l,
                                            0.79540665919549865924l,
                                            0.85699701386308739244l,
                                            0.90909804422384697594l,
                                            0.95032536436570154409l,
                                            0.97958418733152273717l,
                                            0.99610366384229088321l};
          static const long double wts[] = {0.012975111395684900835l,
                                            0.012888764187499150078l,
                                            0.012716644398857307844l,
                                            0.012459897461364705691l,
                                            0.012120231988292330388l,
                                            0.011699908404856841158l,
                                            0.011201723906897224448l,
                                            0.010628993848522759853l,
                                            0.0099855296835573320047l,
                                            0.0092756136096132857933l,
                                            0.0085039700881139589055l,
                                            0.0076757344408814561254l,
                                            0.0067964187616556459109l,
                                            0.005871875456524750363l,
                                            0.0049082589542498110071l,
                                            0.0039119870792519721409l,
                                            0.0028897090921170700834l,
                                            0.0018483371329504443947l,
                                            0.00079623320100438873578l};

          const long double as = a*a;
          const long double hs = -h*h*boost::math::constants::half<long double>();

          long double val = 0;
          for(unsigned short i = 0; i < m; ++i)
            {
              BOOST_ASSERT(i < 19);
              const long double r = 1.l + as*pts[i];
              val += wts[i] * exp( hs*r ) / r;
            } // for(unsigned short i = 0; i < m; ++i)

          return val*a;
        } // long double owens_t_T5(const long double h, const long double a)

#endif // BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS


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
               val -= r * exp( -y*h*h*half<RealType>()/r ) * one_div_two_pi<RealType>();

            return val;
         } // RealType owens_t_T6(const RealType h, const RealType a, const unsigned short m)

         // This routine dispatches the call to one of six subroutines, depending on the values
         // of h and a.
         // preconditions: h >= 0, 0<=a<=1, ah=a*h
         template<typename RealType>
         inline RealType owens_t_dispatch(const RealType h, const RealType a, const RealType ah)
         {
            static const unsigned short meth[] = {1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 3, 4, 4, 4, 4, 5, 6}; // 18 entries
            //static const unsigned short ord[] = {2, 3, 4, 5, 7, 10, 12, 18, 10, 20, 30, 20, 4, 7, 8, 20, 13, 0}; // 18 entries

            RealType val = 0; // avoid compiler warnings, 0 will be overwritten in any case

            const unsigned short icode = owens_t_compute_code(h,a);
            const unsigned short m = owens_t_get_order(icode, val /* just a dummy for the type */);

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
               val = owens_t_T3(h,a,ah);
               break;
            case 4: // T4
               val = owens_t_T4(h,a,m);
               break;
            case 5: // T5
               val = owens_t_T5(h,a);
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
               val = owens_t_dispatch(h, fabs_a, fabs_ah);
            } // if(fabs_a <= 1.0)
            else 
            {
               if( h <= 0.67 )
               {
                  const RealType normh = owens_t_znorm1(h);
                  const RealType normah = owens_t_znorm1(fabs_ah);
                  val = static_cast<RealType>(1)/static_cast<RealType>(4) - normh*normah -
                     owens_t_dispatch(fabs_ah, static_cast<RealType>(1 / fabs_a), h);
               } // if( h <= 0.67 )
               else
               {
                  const RealType normh = detail::owens_t_znorm2(h);
                  const RealType normah = detail::owens_t_znorm2(fabs_ah);
                  val = constants::half<RealType>()*(normh+normah) - normh*normah -
                     owens_t_dispatch(fabs_ah, static_cast<RealType>(1 / fabs_a), h);
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
         return policies::checked_narrowing_cast<result_type, Policy>(detail::owens_t(static_cast<value_type>(h), static_cast<value_type>(a), pol), "boost::math::owens_t<%1%>(%1%,%1%)");
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
