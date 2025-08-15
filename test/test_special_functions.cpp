#include "test_autodiff_reverse.hpp"
#include <boost/math/tools/test_value.hpp>
#include <boost/utility/identity_type.hpp>
BOOST_AUTO_TEST_SUITE(test_special_functions)

using namespace rdiff;
BOOST_AUTO_TEST_CASE_TEMPLATE(test_tgamma, T, all_float_types)
{
    const T              eps = 1000 * boost::math::tools::epsilon<T>(); // percent
    constexpr unsigned   m   = 4; // 10 is too many derivatives for reverse mode
    const T              cx  = 3;
    // Mathematica: N[Table[D[Gamma[x],{x,n}] /. x->3, {n, 0, 10}], 52]
    std::array<T, m + 1> answers{
        {BOOST_MATH_TEST_VALUE(T, 2.0),
         BOOST_MATH_TEST_VALUE(T, 1.845568670196934278786975819835195137915681328120153),
         BOOST_MATH_TEST_VALUE(T, 2.492929991902693057942510065508124245503778067273315),
         BOOST_MATH_TEST_VALUE(T, 3.449965013523673365279327178241708777509009968597547),
         BOOST_MATH_TEST_VALUE(T, 5.521798578098737512443417699412265532987916790978887)}};
    rvar<T, m> x = make_rvar<T, m>(cx);
    rvar<T, m> y = tgamma(x);

    auto g1 = grad_nd<1>(y, &x);
    auto g2 = grad_nd<2>(y, &x);
    auto g3 = grad_nd<3>(y, &x);
    auto g4 = grad_nd<4>(y, &x);

    BOOST_CHECK_CLOSE(y.item(), answers[0], eps);
    BOOST_CHECK_CLOSE(g1[0]->item(), answers[1], eps);
    BOOST_CHECK_CLOSE(g2[0][0]->item(), answers[2], eps);
    BOOST_CHECK_CLOSE(g3[0][0][0]->item(), answers[3], eps);
    BOOST_CHECK_CLOSE(g4[0][0][0][0], answers[4], eps);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(tgamma2_test, T, all_float_types)
{
    //const T eps = 5000 * boost::math::tools::epsilon<T>(); // ok for non-multiprecision
    const T              eps = 500000 * boost::math::tools::epsilon<T>(); // percent
    constexpr unsigned   m   = 4; // 10 is too many derivatives for reverse mode
    const T              cx  = T(-1.5);
    // Mathematica: N[Table[D[Gamma[x],{x,n}] /. x->-3/2, {n, 0, 10}], 52]
    std::array<T, m + 1> answers{{
        BOOST_MATH_TEST_VALUE(T, 2.363271801207354703064223311121526910396732608163183),
        BOOST_MATH_TEST_VALUE(T, 1.661750260668596505586468565464938761014714509096807),
        BOOST_MATH_TEST_VALUE(T, 23.33417984355457252918927856618603412638766668207679),
        BOOST_MATH_TEST_VALUE(T, 47.02130025080143055642555842335081335790754507072526),
        BOOST_MATH_TEST_VALUE(T, 1148.336052788822231948472800239024335856568111484074),
    }};
    rvar<T, m>           x = make_rvar<T, m>(cx);
    rvar<T, m>           y = rdiff::tgamma(1.0 * x);

    auto g1 = grad_nd<1>(y, &x);
    auto g2 = grad_nd<2>(y, &x);
    auto g3 = grad_nd<3>(y, &x);
    auto g4 = grad_nd<4>(y, &x);

    BOOST_CHECK_CLOSE(y.item(), answers[0], eps);
    BOOST_CHECK_CLOSE(g1[0]->item(), answers[1], eps);
    BOOST_CHECK_CLOSE(g2[0][0]->item(), answers[2], eps);
    BOOST_CHECK_CLOSE(g3[0][0][0]->item(), answers[3], eps);
    BOOST_CHECK_CLOSE(g4[0][0][0][0], answers[4], eps);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(lgamma_test, T, all_float_types)
{
    const T              eps = 1000 * boost::math::tools::epsilon<T>(); // percent
    constexpr unsigned   m   = 4;
    const T              cx  = 3;
    // Mathematica: N[Table[D[LogGamma[x],{x,n}] /. x->3, {n, 0, 10}], 52]
    std::array<T, m + 1> answers{{
        BOOST_MATH_TEST_VALUE(T, 0.6931471805599453094172321214581765680755001343602553),
        BOOST_MATH_TEST_VALUE(T, 0.9227843350984671393934879099175975689578406640600764),
        BOOST_MATH_TEST_VALUE(T, 0.3949340668482264364724151666460251892189499012067984),
        BOOST_MATH_TEST_VALUE(T, -0.1541138063191885707994763230228999815299725846809978),
        BOOST_MATH_TEST_VALUE(T, 0.1189394022668291490960221792470074166485057115123614),
    }};

    rvar<T, m> x = make_rvar<T, m>(cx);
    rvar<T, m> y = lgamma(x);

    auto g1 = grad_nd<1>(y, &x);
    auto g2 = grad_nd<2>(y, &x);
    auto g3 = grad_nd<3>(y, &x);
    auto g4 = grad_nd<4>(y, &x);

    BOOST_CHECK_CLOSE(y.item(), answers[0], eps);
    BOOST_CHECK_CLOSE(g1[0]->item(), answers[1], eps);
    BOOST_CHECK_CLOSE(g2[0][0]->item(), answers[2], eps);
    BOOST_CHECK_CLOSE(g3[0][0][0]->item(), answers[3], eps);
    BOOST_CHECK_CLOSE(g4[0][0][0][0], answers[4], eps);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(trigamma_test, T, all_float_types)
{
    const T              eps = 1000 * boost::math::tools::epsilon<T>(); // percent
    constexpr unsigned   m   = 4;
    const T              cx  = 3;
    // Mathematica: N[Table[D[TriGamma[x],{x,n}] /. x->3, {n, 0, 10}], 52]
    std::array<T, m + 1> answers{
        {BOOST_MATH_TEST_VALUE(
             T, 0.3949340668482264364724151666460251892189499012067984377355582293700074704),
         BOOST_MATH_TEST_VALUE(
             T, -0.154113806319188570799476323022899981529972584680997763584543110683676411),
         BOOST_MATH_TEST_VALUE(
             T, 0.1189394022668291490960221792470074166485057115123614460978572926647236971),
         BOOST_MATH_TEST_VALUE(
             T, -0.136266123440878231952771674968820033369942068045907487380624269691286154),
         BOOST_MATH_TEST_VALUE(
             T, 0.2061674381338967657421515749104633482180988039424274210890396805198619482)

        }};

    rvar<T, m> x  = make_rvar<T, m>(cx);
    rvar<T, m> y  = rdiff::trigamma(x);
    auto       g1 = grad_nd<1>(y, &x);
    auto       g2 = grad_nd<2>(y, &x);
    auto       g3 = grad_nd<3>(y, &x);
    auto       g4 = grad_nd<4>(y, &x);

    BOOST_CHECK_CLOSE(y.item(), answers[0], eps);
    BOOST_CHECK_CLOSE(g1[0]->item(), answers[1], eps);
    BOOST_CHECK_CLOSE(g2[0][0]->item(), answers[2], eps);
    BOOST_CHECK_CLOSE(g3[0][0][0]->item(), answers[3], eps);
    BOOST_CHECK_CLOSE(g4[0][0][0][0], answers[4], eps);
}
BOOST_AUTO_TEST_CASE_TEMPLATE(gamma_ratio_test, T, all_float_types)
{
    const T            eps = 10000 * boost::math::tools::epsilon<T>(); // percent
    constexpr unsigned m   = 4;
    const T            cx  = 2;
    const T            cy  = 3;
    // Mathematica: N[Table[D[Gamma[x]/Gamma[y],{x,n}] /. {x->2, y->3}, {n, 0, 4}], 52]

    const T gamma_ratio_2_3 = 0.5;
    std::array<T, m>
        df_x{BOOST_MATH_TEST_VALUE(T, 0.2113921675492335696967439549587987844789203320300382),
             BOOST_MATH_TEST_VALUE(T, 0.4118403304264396947888835614182322768970241847882905),
             BOOST_MATH_TEST_VALUE(T, 0.2447307577412587991365064524330787790317162149669510),
             BOOST_MATH_TEST_VALUE(T, 0.8909881290421667798378415199869088251835467678108199)};
    std::array<T, m>
        df_y{BOOST_MATH_TEST_VALUE(T, -0.4613921675492335696967439549587987844789203320300382),
             BOOST_MATH_TEST_VALUE(T, 0.2282984311274468280132123497310058721569946156115303),
             BOOST_MATH_TEST_VALUE(T, 0.2308256574719042764437644864605804725961017246349102),
             BOOST_MATH_TEST_VALUE(T, -0.7562811950086087527214514759242552399052000145125505)};

    rvar<T, m> x            = make_rvar<T, m>(cx);
    rvar<T, m> y            = make_rvar<T, m>(cy);
    rvar<T, m> f_rvar_rvar  = tgamma_ratio(x, y);
    rvar<T, m> f_rvar_float = tgamma_ratio(x, cy);
    rvar<T, m> f_float_rvar = tgamma_ratio(cx, y);

    BOOST_CHECK_CLOSE(f_rvar_rvar.item(), gamma_ratio_2_3, eps);
    BOOST_CHECK_CLOSE(f_float_rvar.item(), gamma_ratio_2_3, eps);
    BOOST_CHECK_CLOSE(f_rvar_float.item(), gamma_ratio_2_3, eps);

    auto g1x = grad_nd<1>(f_rvar_rvar, &x);
    auto g2x = grad_nd<2>(f_rvar_rvar, &x);
    auto g3x = grad_nd<3>(f_rvar_rvar, &x);
    auto g4x = grad_nd<4>(f_rvar_rvar, &x);

    BOOST_CHECK_CLOSE(g1x[0]->item(), df_x[0], eps);
    BOOST_CHECK_CLOSE(g2x[0][0]->item(), df_x[1], eps);
    BOOST_CHECK_CLOSE(g3x[0][0][0]->item(), df_x[2], eps);
    BOOST_CHECK_CLOSE(g4x[0][0][0][0], df_x[3], eps);

    auto g1y = grad_nd<1>(f_rvar_rvar, &y);
    auto g2y = grad_nd<2>(f_rvar_rvar, &y);
    auto g3y = grad_nd<3>(f_rvar_rvar, &y);
    auto g4y = grad_nd<4>(f_rvar_rvar, &y);

    BOOST_CHECK_CLOSE(g1y[0]->item(), df_y[0], eps);
    BOOST_CHECK_CLOSE(g2y[0][0]->item(), df_y[1], eps);
    BOOST_CHECK_CLOSE(g3y[0][0][0]->item(), df_y[2], eps);
    BOOST_CHECK_CLOSE(g4y[0][0][0][0], df_y[3], eps);

    auto g1y_left_const = grad_nd<1>(f_float_rvar, &y);
    auto g2y_left_const = grad_nd<2>(f_float_rvar, &y);
    auto g3y_left_const = grad_nd<3>(f_float_rvar, &y);
    auto g4y_left_const = grad_nd<4>(f_float_rvar, &y);

    BOOST_CHECK_CLOSE(g1y_left_const[0]->item(), df_y[0], eps);
    BOOST_CHECK_CLOSE(g2y_left_const[0][0]->item(), df_y[1], eps);
    BOOST_CHECK_CLOSE(g3y_left_const[0][0][0]->item(), df_y[2], eps);
    BOOST_CHECK_CLOSE(g4y_left_const[0][0][0][0], df_y[3], eps);

    auto g1x_right_const = grad_nd<1>(f_rvar_float, &x);
    auto g2x_right_const = grad_nd<2>(f_rvar_float, &x);
    auto g3x_right_const = grad_nd<3>(f_rvar_float, &x);
    auto g4x_right_const = grad_nd<4>(f_rvar_float, &x);

    BOOST_CHECK_CLOSE(g1x_right_const[0]->item(), df_x[0], eps);
    BOOST_CHECK_CLOSE(g2x_right_const[0][0]->item(), df_x[1], eps);
    BOOST_CHECK_CLOSE(g3x_right_const[0][0][0]->item(), df_x[2], eps);
    BOOST_CHECK_CLOSE(g4x_right_const[0][0][0][0], df_x[3], eps);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(tgamma_delta_ratio_test, T, all_float_types)
{
    const T            eps             = 1e5 * boost::math::tools::epsilon<T>(); // percent
    constexpr unsigned m               = 4;
    const T            cx              = 2;
    const T            cy              = 0.1;
    // Mathematica: N[Table[D[Gamma[x]/Gamma[x+y],{x,n}] /. {x->2, y->0.1}, {n, 0, 4}], 52]
    const T            gamma_ratio_x_y = 0.95557909646525255223365090810180899602197122;
    std::array<T, m> df_x{BOOST_MATH_TEST_VALUE(T, -0.05977303350010646571844880597482195837913128),
                          BOOST_MATH_TEST_VALUE(T, 0.0401284967562269597550347310948856584574501),
                          BOOST_MATH_TEST_VALUE(T, -0.0503369437770023045591873014080458783558063),
                          BOOST_MATH_TEST_VALUE(T, 0.091481457097213527148292723095043376439440)};

    std::array<T, m> df_y{BOOST_MATH_TEST_VALUE(T, -0.4637769064331622567147529229387456253036133),
                          BOOST_MATH_TEST_VALUE(T, -0.35480830287673262213106039091536927521578898),
                          BOOST_MATH_TEST_VALUE(T, 1.0779782470986784425646874793787772212820541),
                          BOOST_MATH_TEST_VALUE(T, -0.7728508571145974103710087287530096388223087)};

    rvar<T, m> x           = make_rvar<T, m>(cx);
    rvar<T, m> y           = make_rvar<T, m>(cy);
    rvar<T, m> f_rvar_rvar = tgamma_delta_ratio(x, y);

    rvar<T, m> f_rvar_float = tgamma_delta_ratio(x, cy);
    rvar<T, m> f_float_rvar = tgamma_delta_ratio(cx, y);

    BOOST_CHECK_CLOSE(f_rvar_rvar.item(), gamma_ratio_x_y, eps);
    BOOST_CHECK_CLOSE(f_float_rvar.item(), gamma_ratio_x_y, eps);
    BOOST_CHECK_CLOSE(f_rvar_float.item(), gamma_ratio_x_y, eps);

    auto g1x = grad_nd<1>(f_rvar_rvar, &x);
    auto g2x = grad_nd<2>(f_rvar_rvar, &x);
    auto g3x = grad_nd<3>(f_rvar_rvar, &x);
    auto g4x = grad_nd<4>(f_rvar_rvar, &x);

    BOOST_CHECK_CLOSE(g1x[0]->item(), df_x[0], eps);
    BOOST_CHECK_CLOSE(g2x[0][0]->item(), df_x[1], eps);
    BOOST_CHECK_CLOSE(g3x[0][0][0]->item(), df_x[2], eps);
    BOOST_CHECK_CLOSE(g4x[0][0][0][0], df_x[3], eps);

    auto g1y = grad_nd<1>(f_rvar_rvar, &y);
    auto g2y = grad_nd<2>(f_rvar_rvar, &y);
    auto g3y = grad_nd<3>(f_rvar_rvar, &y);
    auto g4y = grad_nd<4>(f_rvar_rvar, &y);

    BOOST_CHECK_CLOSE(g1y[0]->item(), df_y[0], eps);
    BOOST_CHECK_CLOSE(g2y[0][0]->item(), df_y[1], eps);
    BOOST_CHECK_CLOSE(g3y[0][0][0]->item(), df_y[2], eps);
    BOOST_CHECK_CLOSE(g4y[0][0][0][0], df_y[3], eps);

    auto g1y_left_const = grad_nd<1>(f_float_rvar, &y);
    auto g2y_left_const = grad_nd<2>(f_float_rvar, &y);
    auto g3y_left_const = grad_nd<3>(f_float_rvar, &y);
    auto g4y_left_const = grad_nd<4>(f_float_rvar, &y);

    BOOST_CHECK_CLOSE(g1y_left_const[0]->item(), df_y[0], eps);
    BOOST_CHECK_CLOSE(g2y_left_const[0][0]->item(), df_y[1], eps);
    BOOST_CHECK_CLOSE(g3y_left_const[0][0][0]->item(), df_y[2], eps);
    BOOST_CHECK_CLOSE(g4y_left_const[0][0][0][0], df_y[3], eps);

    auto g1x_right_const = grad_nd<1>(f_rvar_float, &x);
    auto g2x_right_const = grad_nd<2>(f_rvar_float, &x);
    auto g3x_right_const = grad_nd<3>(f_rvar_float, &x);
    auto g4x_right_const = grad_nd<4>(f_rvar_float, &x);

    BOOST_CHECK_CLOSE(g1x_right_const[0]->item(), df_x[0], eps);
    BOOST_CHECK_CLOSE(g2x_right_const[0][0]->item(), df_x[1], eps);
    BOOST_CHECK_CLOSE(g3x_right_const[0][0][0]->item(), df_x[2], eps);
    BOOST_CHECK_CLOSE(g4x_right_const[0][0][0][0], df_x[3], eps);
}

BOOST_AUTO_TEST_SUITE_END()
