//  (C) Copyright John Maddock 2007 - 2025.
//  (C) Copyright Christopher Kormanyos 2025.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#define BOOST_MATH_OVERFLOW_ERROR_POLICY ignore_error
#include <boost/math/concepts/real_concept.hpp>
#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/math/special_functions/math_fwd.hpp>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/type_traits/is_floating_point.hpp>
#include <boost/array.hpp>
#include "functor.hpp"

#include "handle_test_result.hpp"
#include "table_type.hpp"

#ifndef SC_
#  define SC_(x) static_cast<typename table_type<T>::type>(BOOST_JOIN(x, L))
#endif

template <class Real, class T>
void do_test_cyl_bessel_j(const T& data, const char* type_name, const char* test_name)
{
#if !(defined(ERROR_REPORTING_MODE) && !defined(BESSEL_J_FUNCTION_TO_TEST))
   typedef Real                   value_type;

   typedef value_type (*pg)(value_type, value_type);
#ifdef BESSEL_J_FUNCTION_TO_TEST
   pg funcp = BESSEL_J_FUNCTION_TO_TEST;
#elif defined(BOOST_MATH_NO_DEDUCED_FUNCTION_POINTERS)
   pg funcp = boost::math::cyl_bessel_j<value_type, value_type>;
#else
   pg funcp = boost::math::cyl_bessel_j;
#endif

   boost::math::tools::test_result<value_type> result;

   std::cout << "Testing " << test_name << " with type " << type_name
      << "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";

   //
   // test cyl_bessel_j against data:
   //
   result = boost::math::tools::test_hetero<Real>(
      data, 
      bind_func<Real>(funcp, 0, 1), 
      extract_result<Real>(2));
   handle_test_result(result, data[result.worst()], result.worst(), type_name, "cyl_bessel_j", test_name);
   std::cout << std::endl;

#endif
}

template <class T>
T cyl_bessel_j_int_wrapper(T v, T x)
{
#ifdef BESSEL_JN_FUNCTION_TO_TEST
   return static_cast<T>(BESSEL_JN_FUNCTION_TO_TEST(boost::math::itrunc(v), x));
#else
   return static_cast<T>(boost::math::cyl_bessel_j(boost::math::itrunc(v), x));
#endif
}


template <class Real, class T>
void do_test_cyl_bessel_j_int(const T& data, const char* type_name, const char* test_name)
{
#if !(defined(ERROR_REPORTING_MODE) && !defined(BESSEL_JN_FUNCTION_TO_TEST))
   typedef Real                   value_type;

   typedef value_type (*pg)(value_type, value_type);
#if defined(BOOST_MATH_NO_DEDUCED_FUNCTION_POINTERS)
   pg funcp = cyl_bessel_j_int_wrapper<value_type>;
#else
   pg funcp = cyl_bessel_j_int_wrapper;
#endif

   boost::math::tools::test_result<value_type> result;

   std::cout << "Testing " << test_name << " with type " << type_name
      << "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";

   //
   // test cyl_bessel_j against data:
   //
   result = boost::math::tools::test_hetero<Real>(
      data, 
      bind_func<Real>(funcp, 0, 1), 
      extract_result<Real>(2));
   handle_test_result(result, data[result.worst()], result.worst(), type_name, "cyl_bessel_j (integer orders)", test_name);
   std::cout << std::endl;
#endif
}

template <class Real, class T>
void do_test_sph_bessel_j(const T& data, const char* type_name, const char* test_name)
{
#if !(defined(ERROR_REPORTING_MODE) && !defined(BESSEL_JS_FUNCTION_TO_TEST))
   typedef Real                   value_type;

   typedef value_type (*pg)(unsigned, value_type);
#ifdef BESSEL_JS_FUNCTION_TO_TEST
   pg funcp = BESSEL_JS_FUNCTION_TO_TEST;
#elif defined(BOOST_MATH_NO_DEDUCED_FUNCTION_POINTERS)
   pg funcp = boost::math::sph_bessel<value_type>;
#else
   pg funcp = boost::math::sph_bessel;
#endif

   boost::math::tools::test_result<value_type> result;

   std::cout << "Testing " << test_name << " with type " << type_name
      << "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";

   //
   // test sph_bessel against data:
   //
   result = boost::math::tools::test_hetero<Real>(
      data, 
      bind_func_int1<Real>(funcp, 0, 1), 
      extract_result<Real>(2));
   handle_test_result(result, data[result.worst()], result.worst(), type_name, "sph_bessel", test_name);
   std::cout << std::endl;
#endif
}

template <class T>
void test_bessel(T, const char* name)
{
   using table_entry_type = typename table_type<T>::type;

   //
   // The actual test data is rather verbose, so it's in a separate file
   //
   // The contents are as follows, each row of data contains
   // three items, input value a, input value b and erf(a, b):
   // 
    // function values calculated on http://functions.wolfram.com/
    static const std::array<std::array<table_entry_type, 3>, 8> j0_data = {{
       { { SC_(0.0), SC_(0.0), SC_(1.0) } },
       { { SC_(0.0), SC_(1.0), SC_(0.7651976865579665514497175261026632209093) } },
       { { SC_(0.0), SC_(-2.0), SC_(0.2238907791412356680518274546499486258252) } },
       { { SC_(0.0), SC_(4.0), SC_(-0.3971498098638473722865907684516980419756) } },
       { { SC_(0.0), SC_(-8.0), SC_(0.1716508071375539060908694078519720010684) } },
        {{ SC_(0.0), SC_(1e-05), SC_(0.999999999975000000000156249999999565972) }},
        {{ SC_(0.0), SC_(1e-10), SC_(0.999999999999999999997500000000000000000) }},
        {{ SC_(0.0), SC_(-1e+01), SC_(-0.2459357644513483351977608624853287538296) }},
    }};
    static const std::array<std::array<table_entry_type, 3>, 6> j0_tricky = {{
        // Big numbers make the accuracy of std::sin the limiting factor:
       { { SC_(0.0), SC_(1e+03), SC_(0.02478668615242017456133073111569370878617) } },
       { { SC_(0.0), SC_(1e+05), SC_(-0.001719201116235972192570601477073201747532) } },
        // test at the roots:
       { { SC_(0.0), SC_(2.4048252105712890625) /*T(2521642.0) / (1024 * 1024)*/, SC_(1.80208819970046790002973759410972422387259992955354630042138e-7) } },
       { { SC_(0.0), SC_(5.52007770538330078125) /*T(5788221.0) / (1024 * 1024)*/, SC_(-1.37774249380686777043369399806210229535671843632174587432454e-7) } },
       { { SC_(0.0), SC_(8.65372753143310546875) /*T(9074091.0) / (1024 * 1024)*/, SC_(1.03553057441100845081018471279571355857520645127532785991335e-7) } },
       { { SC_(0.0), SC_(11.791534423828125) /*T(12364320.0) / (1024 * 1024)*/, SC_(-3.53017140778223781420794006033810387155048392363051866610931e-9) } }
    }};    

    static const std::array<std::array<table_entry_type, 3>, 8> j1_data = {{
       { { SC_(1.0), SC_(0.0), SC_(0.0) } },
       { { SC_(1.0), SC_(1.0), SC_(0.4400505857449335159596822037189149131274) } },
       { { SC_(1.0), SC_(-2.0), SC_(-0.5767248077568733872024482422691370869203) } },
       { { SC_(1.0), SC_(4.0), SC_(-6.604332802354913614318542080327502872742e-02) } },
       { { SC_(1.0), SC_(-8.0), SC_(-0.2346363468539146243812766515904546115488) } },
       { { SC_(1.0), SC_(1e-05), SC_(4.999999999937500000000260416666666124132e-06) } },
       { { SC_(1.0), SC_(1e-10), SC_(4.999999999999999999993750000000000000000e-11) } },
       { { SC_(1.0), SC_(-1e+01), SC_(-4.347274616886143666974876802585928830627e-02) } },
    }};
    static const std::array<std::array<table_entry_type, 3>, 5> j1_tricky = {{
        // Big numbers make the accuracy of std::sin the limiting factor:
       { { SC_(1.0), SC_(1e+03), SC_(4.728311907089523917576071901216916285418e-03) } },
       { { SC_(1.0), SC_(1e+05), SC_(1.846757562882567716362123967114215743694e-03) } },
        // test zeros:
       { { SC_(1.0), SC_(3.8317050933837890625) /*T(4017834) / (1024 * 1024)*/, SC_(3.53149033321258645807835062770856949751958513973522222203044e-7) } },
       { { SC_(1.0), SC_(7.01558589935302734375) /*T(7356375) / (1024 * 1024)*/, SC_(-2.31227973111067286051984021150135526024117175836722748404342e-7) } },
       { { SC_(1.0), SC_(10.1734676361083984375) /*T(10667654) / (1024 * 1024)*/, SC_(1.24591331097191900488116495350277530373473085499043086981229e-7) } },
    }};

    static const std::array<std::array<table_entry_type, 3>, 17> jn_data = {{
        // This first one is a modified test case from https://svn.boost.org/trac/boost/ticket/2733
       { { SC_(-1.0), SC_(1.25), SC_(-0.510623260319880467069474837274910375352924050139633057168856) } },
       { { SC_(2.0), SC_(0.0), SC_(0.0) } },
       { { SC_(-2.0), SC_(0.0), SC_(0.0) } },
       { { SC_(2.0), SC_(1e-02), SC_(1.249989583365885362413250958437642113452e-05) } },
       { { SC_(5.0), SC_(10.0), SC_(-0.2340615281867936404436949416457777864635) } },
       { { SC_(5.0), SC_(-10.0), SC_(0.2340615281867936404436949416457777864635) } },
       { { SC_(-5.0), SC_(1e+06), SC_(7.259643842453285052375779970433848914846e-04) } },
       { { SC_(5.0), SC_(1e+06), SC_(-0.000725964384245328505237577997043384891484649290328285235308619) } },
       { { SC_(-5.0), SC_(-1.0), SC_(2.497577302112344313750655409880451981584e-04) } },
       { { SC_(10.0), SC_(10.0), SC_(0.2074861066333588576972787235187534280327) } },
       { { SC_(10.0), SC_(-10.0), SC_(0.2074861066333588576972787235187534280327) } },
       { { SC_(10.0), SC_(-5.0), SC_(1.467802647310474131107532232606627020895e-03) } },
       { { SC_(-10.0), SC_(1e+06), SC_(-3.310793117604488741264958559035744460210e-04) } },
       { { SC_(10.0), SC_(1e+06), SC_(-0.000331079311760448874126495855903574446020957243277028930713243) } },
        {{ SC_(1e+02), SC_(8e+01), SC_(4.606553064823477354141298259169874909670e-06) }},
        {{ SC_(1e+03), SC_(1e+05), SC_(1.283178112502480365195139312635384057363e-03) }},
        { { SC_(10.0), SC_(1e-100), SC_(2.69114445546737213403880070546737213403880070546737213403880e-1010) } },
    }};
    do_test_cyl_bessel_j<T>(j0_data, name, "Bessel J0: Mathworld Data");
    do_test_cyl_bessel_j<T>(j0_tricky, name, "Bessel J0: Mathworld Data (Tricky cases)");
    do_test_cyl_bessel_j<T>(j1_data, name, "Bessel J1: Mathworld Data");
    do_test_cyl_bessel_j<T>(j1_tricky, name, "Bessel J1: Mathworld Data (tricky cases)");
    do_test_cyl_bessel_j<T>(jn_data, name, "Bessel JN: Mathworld Data");

    do_test_cyl_bessel_j_int<T>(j0_data, name, "Bessel J0: Mathworld Data (Integer Version)");
    do_test_cyl_bessel_j_int<T>(j0_tricky, name, "Bessel J0: Mathworld Data (Tricky cases) (Integer Version)");
    do_test_cyl_bessel_j_int<T>(j1_data, name, "Bessel J1: Mathworld Data (Integer Version)");
    do_test_cyl_bessel_j_int<T>(j1_tricky, name, "Bessel J1: Mathworld Data (tricky cases) (Integer Version)");
    do_test_cyl_bessel_j_int<T>(jn_data, name, "Bessel JN: Mathworld Data (Integer Version)");

    static const std::array<std::array<table_entry_type, 3>, 20> jv_data = {{
        //SC_(-2.4), {{ SC_(0.0), std::numeric_limits<T>::infinity() }},
       { { SC_(22.5), SC_(0.0), SC_(0.0) } },
       { { SC_(2.3994140625) /*2457.0 / 1024*/, SC_(0.0009765625) /* 1 / 1024*/, SC_(3.80739920118603335646474073457326714709615200130620574875292e-9) } },
        {{ SC_(5.5), SC_(3.1416015625) /* 3217/1024*/, SC_(0.0281933076257506091621579544064767140470089107926550720453038) }},
        {{ SC_(-5.5), SC_(3.1416015625) /* 3217/1024*/, SC_(-2.55820064470647911823175836997490971806135336759164272675969) }},
        {{ SC_(-5.5), SC_(1e+04), SC_(2.449843111985605522111159013846599118397e-03) }},
        {{ SC_(5.5), SC_(1e+04), SC_(0.00759343502722670361395585198154817047185480147294665270646578) }},
        {{ SC_(5.5), SC_(1e+06), SC_(-0.000747424248595630177396350688505919533097973148718960064663632) }},
        {{ SC_(5.125), SC_(1e+06), SC_(-0.000776600124835704280633640911329691642748783663198207360238214) }},
        {{ SC_(5.875), SC_(1e+06), SC_(-0.000466322721115193071631008581529503095819705088484386434589780) }},
        { { SC_(0.5), SC_(101.0), SC_(0.0358874487875643822020496677692429287863419555699447066226409) } },
        {{ SC_(-5.5), SC_(1e+04), SC_(0.00244984311198560552211115901384659911839737686676766460822577) }},
        {{ SC_(-5.5), SC_(1e+06), SC_(0.000279243200433579511095229508894156656558211060453622750659554) }},
        { { SC_(-0.5), SC_(101.0), SC_(0.0708184798097594268482290389188138201440114881159344944791454) } },
        {{ SC_(-10.0002994537353515625) /* -10486074 / (1024*1024)*/, SC_(0.0009765625) /* 1/1024*/, SC_(1.41474013160494695750009004222225969090304185981836460288562e35) }},
        { { SC_(-10.0002994537353515625) /* -10486074 / (1024*1024)*/, SC_(15.0), SC_(-0.0902239288885423309568944543848111461724911781719692852541489) } },
        {{ SC_(-10.0002994537353515625) /* -10486074 / (1024*1024)*/, SC_(100.0), SC_(-0.05476136603168065513386371539426045507795139476742228638) }},
        {{ SC_(-10.0002994537353515625) /* -10486074 / (1024*1024)*/, SC_(20000.0), SC_(-0.00556869085445857782456414284057389040183758546505700058) }},
        // Bug report https://svn.boost.org/trac/boost/ticket/4812:
        {{ SC_(1.5), SC_(7.845703125) /* 8034/1024*/, SC_(0.0339477646369710610146236955872928005087352629422508823945264) }},
        {{ SC_(8.5), SC_(12.566370614359172953850573533118011536788677597500423283899778369231265625144835994512139301368468271928592346053) /*Pi * 4*/, SC_(0.0436807946352780974532519564114026730332781693877984686758680) }},
        {{ SC_(-8.5), SC_(12.566370614359172953850573533118011536788677597500423283899778369231265625144835994512139301368468271928592346053) /*Pi * 4*/, SC_(-0.257086543428224355151772807588810984369026142375675714560864) }},
    }};
    do_test_cyl_bessel_j<T>(jv_data, name, "Bessel J: Mathworld Data");
    static const std::array<std::array<table_entry_type, 3>, 4> jv_large_data = {{
        // Bug report https://svn.boost.org/trac/boost/ticket/5560:
        {{ SC_(-0.5), SC_(1.2458993688871959419388378518880931736878259938089494331010226962863582408064841833232475731084062642684629e-206) /*static_cast<T>(std::ldexp(0.5, -683))*/, SC_(7.14823099969225685526188875418476476336424046896822867989728e102) }},
        { { SC_(256.0), SC_(512.0), SC_(0.00671672065717513246956991122723250578101154313313749938944675) } },
        { { SC_(-256.0), SC_(8.0), SC_(1.46866142030022704638298523775638527553596432641223316232692e-353) } },
        { { SC_(-2.5), SC_(4.0), SC_(-0.0145679476685218007666785535204236327832335803441449596297004) } },
    }};
    if(static_cast<T>(jv_large_data[0][1]) != 0)
      do_test_cyl_bessel_j<T>(jv_large_data, name, "Bessel J: Mathworld Data (large values)");

#include "bessel_j_int_data.ipp"
    do_test_cyl_bessel_j<T>(bessel_j_int_data, name, "Bessel JN: Random Data");

#include "bessel_j_data.ipp"
    do_test_cyl_bessel_j<T>(bessel_j_data, name, "Bessel J: Random Data");

#include "bessel_j_large_data.ipp"
    do_test_cyl_bessel_j<T>(bessel_j_large_data, name, "Bessel J: Random Data (Tricky large values)");

#include "sph_bessel_data.ipp"
    do_test_sph_bessel_j<T>(sph_bessel_data, name, "Bessel j: Random Data");

    //
    // Some special cases:
    //

    // Specific tests for Git-issue1292.
    {
        using local_ctrl_array_type = std::array<table_entry_type, std::size_t { UINT8_C(43) }>;

        // Table[N[BesselJ[3, n/10], 40], {n, 9, 51, 1}]
        static const local_ctrl_array_type ctrl_data =
        {{
            SC_(0.01443402847586617545767791623904539755731), SC_(0.01956335398266840591890532162175150825451), SC_(0.02569452861246328174726417617756888741432), SC_(0.03287433692499494270867882730165246683837),
            SC_(0.04113582571991693187673486447516908751463), SC_(0.05049771328895129623567992727476043273558), SC_(0.06096395114113963064394955997646387979571), SC_(0.07252344333261900300034928368068248675877),
            SC_(0.08514992694801526415321095754253909148633), SC_(0.09880201565861918291536618746528733463749), SC_(0.1134234066389601112649841240858617923591), SC_(0.1289432494744020510987933329692398352700),
            SC_(0.1452766740542063665759023355570418120750), SC_(0.1623254728332874543121706910035271736854), SC_(0.1799789312775334540800304157279732327839), SC_(0.1981147987975668248498434552081155790183),
            SC_(0.2166003910391135247666890035159637217168), SC_(0.2352938130489638091015220916013129483423), SC_(0.2540452915872273499615464996563039918262), SC_(0.2726986037216204380267188592437356599939),
            SC_(0.2910925878291867784836313080855848616815), SC_(0.3090627222552516436182601949468331494291), SC_(0.3264427561473409695937042738575781129080), SC_(0.3430663764006682009386373318558777864023),
            SC_(0.3587688942275418259451574456258027163924), SC_(0.3733889346000900583527754127339797472980), SC_(0.3867701117168813668578718121131100327218), SC_(0.3987626737105880326848194417650226836608),
            SC_(0.4092251000454309977422936498249743734653), SC_(0.4180256354477855744864458808409348352597), SC_(0.4250437447674560017637404058105525727991), SC_(0.4301714738756219403581834788533355563393),
            SC_(0.4333147025616927046073022200802734463060), SC_(0.4343942763872007823091130214493427347554), SC_(0.4333470055809823422144251313032973397899), SC_(0.4301265203055088083605755042771532591535),
            SC_(0.4247039729774556002468140098011553543390), SC_(0.4170685797734672711167804755454067582755), SC_(0.4072279949807128989552790124633945783765), SC_(0.3952085134465309348696666123753181072022),
            SC_(0.3810550980268886849843356923521907577982), SC_(0.3648312306136669944635769493587219791343), SC_(0.3466185870197064968846647990300282094299)
        }};

        const T tolerance { 128 * boost::math::tools::epsilon<T>() };

        int n_val { 9 };

        for(const T& ctrl: ctrl_data)
        {
            const T x_val { static_cast<T>(static_cast<T>(n_val) / 10) };

            const T jn_val { boost::math::cyl_bessel_j(3, x_val) };

            ++n_val;

            BOOST_CHECK_CLOSE_FRACTION(jn_val, ctrl, tolerance);
        }
    }

    // More specific tests for Git-issue1292.
    {
        using local_ctrl_array_type = std::array<table_entry_type, std::size_t { UINT8_C(43) }>;

        // Table[N[BesselJ[31 / 10, n/10], 40], {n, 9, 51, 1}]
        static const local_ctrl_array_type ctrl_data =
        {{
            SC_(0.01175139795214295170487105485346781171863), SC_(0.01610092560641321584451701371378836908343), SC_(0.02135659148701787280713314897691402425560), SC_(0.02757316602094671775387912184375438913864),
            SC_(0.03479372470942323424236519547108889796879), SC_(0.04304884612770317349181083608393216357649), SC_(0.05235595486839477302545989170267853515412), SC_(0.06271881444009116864560457340917173135492),
            SC_(0.07412717428914531462554459978389560408816), SC_(0.08655657401374878955779092959288219537482), SC_(0.09996830657138141192006478769855556316995), SC_(0.1143095409011066041623799431143340837007),
            SC_(0.1295136029333754366226206053365805922136), SC_(0.1455004124749064242445094865982308151833), SC_(0.1621770719619318324763655518038365467780), SC_(0.1794386015949089605595848208470049166853),
            SC_(0.1971688139222575484595939469080319406355), SC_(0.2152413195483352761119610246805962611319), SC_(0.2335206543185642795903443565627832597552), SC_(0.2518635170977563283446184976264924077045),
            SC_(0.2701201061202457234361463802484836927464), SC_(0.2881355408650536940851830262098172944660), SC_(0.3057513555072237123913927136839853900197), SC_(0.3228070492275195774383850177821175151772),
            SC_(0.3391416780352496831991017115498371039332), SC_(0.3545954722799571667706920253581728149226), SC_(0.3690114637024459111815176585531943143085), SC_(0.3822371057078773326211699424651752096338),
            SC_(0.3941258705356621024731438508712990073773), SC_(0.4045388071531719615179931506197074798366), SC_(0.4133460440118967760814988295083688541739), SC_(0.4204282212729730452147271372029342292548),
            SC_(0.4256778377298582292448895379334671298297), SC_(0.4290004984236535775073755577697874033289), SC_(0.4303160498540635572666996674883665298558), SC_(0.4295595907277138677145484170304736319185),
            SC_(0.4266823473457215250850626209433240464185), SC_(0.4216524040030023831246842156638018726147), SC_(0.4144552801406973126588418824207056743782), SC_(0.4050943474472003316564108983454900189419),
            SC_(0.3935910816286283019261303507307813773886), SC_(0.3799851451515116989978414766085547417497), SC_(0.3643342988837623802358278078893847672838)
        }};

        const T tolerance { 128 * boost::math::tools::epsilon<T>() };

        const T vu_val { static_cast<T>(static_cast<T>(31) / 10) };

        int n_val { 9 };

        for(const T& ctrl : ctrl_data)
        {
            const T x_val { static_cast<T>(static_cast<T>(n_val) / 10) };

            const T jn_val { boost::math::cyl_bessel_j(vu_val, x_val) };

            ++n_val;

            BOOST_CHECK_CLOSE_FRACTION(jn_val, ctrl, tolerance);
        }
    }

    BOOST_CHECK_EQUAL(boost::math::sph_bessel(0, T(0)), T(1));
    BOOST_CHECK_EQUAL(boost::math::sph_bessel(1, T(0)), T(0));
    BOOST_CHECK_EQUAL(boost::math::sph_bessel(100000, T(0)), T(0));

    //
    // Special cases that are errors:
    //
    BOOST_MATH_CHECK_THROW(boost::math::cyl_bessel_j(T(-2.5), T(0)), std::domain_error);
    BOOST_MATH_CHECK_THROW(boost::math::cyl_bessel_j(T(-2.5), T(-2)), std::domain_error);
    BOOST_MATH_CHECK_THROW(boost::math::cyl_bessel_j(T(2.5), T(-2)), std::domain_error);

    //
    // special cases for code coverage:
    //
    T tolerance = boost::math::tools::epsilon<T>() * 2000;
    BOOST_CHECK_CLOSE_FRACTION(boost::math::sph_bessel(200, T(0.5)), SC_(3.070403008048099934928128420285169174541102108657574230431e-497), tolerance);
    BOOST_MATH_CHECK_THROW(boost::math::sph_bessel(2, T(-2.0)), std::domain_error);
    BOOST_CHECK_EQUAL(boost::math::cyl_bessel_j(T(0), T(2.5)), boost::math::cyl_bessel_j(T(0), T(-2.5)));
    BOOST_CHECK_EQUAL(boost::math::cyl_bessel_j(T(1), T(2.5)), -boost::math::cyl_bessel_j(T(1), T(-2.5)));
    #ifndef SYCL_LANGUAGE_VERSION
    BOOST_CHECK_CLOSE_FRACTION(boost::math::cyl_bessel_j(364, T(38.5)), SC_(1.793940496519190500748409872348034004417458734118663909894e-309), tolerance);
    #endif
    //
    // Special cases at infinity:
    //
    BOOST_IF_CONSTEXPR (std::numeric_limits<T>::has_infinity)
    {
       BOOST_CHECK_EQUAL(boost::math::cyl_bessel_j(T(0), std::numeric_limits<T>::infinity()), T(0));
       BOOST_CHECK_EQUAL(boost::math::cyl_bessel_j(T(2), std::numeric_limits<T>::infinity()), T(0));
       BOOST_CHECK_EQUAL(boost::math::cyl_bessel_j(T(1.25), std::numeric_limits<T>::infinity()), T(0));
       BOOST_CHECK_EQUAL(boost::math::cyl_bessel_j(T(-1.25), std::numeric_limits<T>::infinity()), T(0));
       BOOST_CHECK_EQUAL(boost::math::sph_bessel(0, std::numeric_limits<T>::infinity()), T(0));
       BOOST_CHECK_EQUAL(boost::math::sph_bessel(1, std::numeric_limits<T>::infinity()), T(0));
       BOOST_CHECK_EQUAL(boost::math::sph_bessel(2, std::numeric_limits<T>::infinity()), T(0));
       BOOST_CHECK_EQUAL(boost::math::sph_bessel(3, std::numeric_limits<T>::infinity()), T(0));

       BOOST_CHECK_EQUAL(boost::math::cyl_bessel_j(T(0), -std::numeric_limits<T>::infinity()), T(0));
       BOOST_CHECK_EQUAL(boost::math::cyl_bessel_j(T(2), -std::numeric_limits<T>::infinity()), T(0));
    }
}
