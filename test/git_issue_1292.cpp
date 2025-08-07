#include <boost/math/special_functions/bessel.hpp>

#include <boost/core/lightweight_test.hpp>

#include <array>
#include <cstddef>
#include <cstdint>

template<typename FloatType>
auto test_n() -> void
{
  using local_float_type = FloatType;

  using local_ctrl_array_type = std::array<local_float_type, std::size_t { UINT8_C(43) }>;

  // Table[N[BesselJ[3, n/10], 40], {n, 9, 51, 1}]
  const local_ctrl_array_type ctrl_data =
  {{
    static_cast<local_float_type>(0.01443402847586617545767791623904539755731L),
    static_cast<local_float_type>(0.01956335398266840591890532162175150825451L),
    static_cast<local_float_type>(0.02569452861246328174726417617756888741432L),
    static_cast<local_float_type>(0.03287433692499494270867882730165246683837L),
    static_cast<local_float_type>(0.04113582571991693187673486447516908751463L),
    static_cast<local_float_type>(0.05049771328895129623567992727476043273558L),
    static_cast<local_float_type>(0.06096395114113963064394955997646387979571L),
    static_cast<local_float_type>(0.07252344333261900300034928368068248675877L),
    static_cast<local_float_type>(0.08514992694801526415321095754253909148633L),
    static_cast<local_float_type>(0.09880201565861918291536618746528733463749L),
    static_cast<local_float_type>(0.1134234066389601112649841240858617923591L),
    static_cast<local_float_type>(0.1289432494744020510987933329692398352700L),
    static_cast<local_float_type>(0.1452766740542063665759023355570418120750L),
    static_cast<local_float_type>(0.1623254728332874543121706910035271736854L),
    static_cast<local_float_type>(0.1799789312775334540800304157279732327839L),
    static_cast<local_float_type>(0.1981147987975668248498434552081155790183L),
    static_cast<local_float_type>(0.2166003910391135247666890035159637217168L),
    static_cast<local_float_type>(0.2352938130489638091015220916013129483423L),
    static_cast<local_float_type>(0.2540452915872273499615464996563039918262L),
    static_cast<local_float_type>(0.2726986037216204380267188592437356599939L),
    static_cast<local_float_type>(0.2910925878291867784836313080855848616815L),
    static_cast<local_float_type>(0.3090627222552516436182601949468331494291L),
    static_cast<local_float_type>(0.3264427561473409695937042738575781129080L),
    static_cast<local_float_type>(0.3430663764006682009386373318558777864023L),
    static_cast<local_float_type>(0.3587688942275418259451574456258027163924L),
    static_cast<local_float_type>(0.3733889346000900583527754127339797472980L),
    static_cast<local_float_type>(0.3867701117168813668578718121131100327218L),
    static_cast<local_float_type>(0.3987626737105880326848194417650226836608L),
    static_cast<local_float_type>(0.4092251000454309977422936498249743734653L),
    static_cast<local_float_type>(0.4180256354477855744864458808409348352597L),
    static_cast<local_float_type>(0.4250437447674560017637404058105525727991L),
    static_cast<local_float_type>(0.4301714738756219403581834788533355563393L),
    static_cast<local_float_type>(0.4333147025616927046073022200802734463060L),
    static_cast<local_float_type>(0.4343942763872007823091130214493427347554L),
    static_cast<local_float_type>(0.4333470055809823422144251313032973397899L),
    static_cast<local_float_type>(0.4301265203055088083605755042771532591535L),
    static_cast<local_float_type>(0.4247039729774556002468140098011553543390L),
    static_cast<local_float_type>(0.4170685797734672711167804755454067582755L),
    static_cast<local_float_type>(0.4072279949807128989552790124633945783765L),
    static_cast<local_float_type>(0.3952085134465309348696666123753181072022L),
    static_cast<local_float_type>(0.3810550980268886849843356923521907577982L),
    static_cast<local_float_type>(0.3648312306136669944635769493587219791343L),
    static_cast<local_float_type>(0.3466185870197064968846647990300282094299L)
  }};

  int n_val { 9 };

  for(std::size_t index { UINT8_C(0) }; index < std::tuple_size<local_ctrl_array_type>::value; ++index)
  {
    const local_float_type x_val { static_cast<local_float_type>(static_cast<local_float_type>(n_val) / 10) };

    const local_float_type jn_val { boost::math::cyl_bessel_j(3, x_val) };

    ++n_val;

    using std::fabs;

    const local_float_type ratio { jn_val / ctrl_data[index] };
    const local_float_type delta { fabs(1 - ratio) };

    BOOST_TEST(delta < 128 * std::numeric_limits<local_float_type>::epsilon());
  }
}

template<typename FloatType>
auto test_vu() -> void
{
  using local_float_type = FloatType;

  using local_ctrl_array_type = std::array<local_float_type, std::size_t { UINT8_C(43) }>;

  // Table[N[BesselJ[31 / 10, n/10], 40], {n, 9, 51, 1}]
  const local_ctrl_array_type ctrl_data =
  {{
    static_cast<local_float_type>(0.01175139795214295170487105485346781171863L),
    static_cast<local_float_type>(0.01610092560641321584451701371378836908343L),
    static_cast<local_float_type>(0.02135659148701787280713314897691402425560L),
    static_cast<local_float_type>(0.02757316602094671775387912184375438913864L),
    static_cast<local_float_type>(0.03479372470942323424236519547108889796879L),
    static_cast<local_float_type>(0.04304884612770317349181083608393216357649L),
    static_cast<local_float_type>(0.05235595486839477302545989170267853515412L),
    static_cast<local_float_type>(0.06271881444009116864560457340917173135492L),
    static_cast<local_float_type>(0.07412717428914531462554459978389560408816L),
    static_cast<local_float_type>(0.08655657401374878955779092959288219537482L),
    static_cast<local_float_type>(0.09996830657138141192006478769855556316995L),
    static_cast<local_float_type>(0.1143095409011066041623799431143340837007L),
    static_cast<local_float_type>(0.1295136029333754366226206053365805922136L),
    static_cast<local_float_type>(0.1455004124749064242445094865982308151833L),
    static_cast<local_float_type>(0.1621770719619318324763655518038365467780L),
    static_cast<local_float_type>(0.1794386015949089605595848208470049166853L),
    static_cast<local_float_type>(0.1971688139222575484595939469080319406355L),
    static_cast<local_float_type>(0.2152413195483352761119610246805962611319L),
    static_cast<local_float_type>(0.2335206543185642795903443565627832597552L),
    static_cast<local_float_type>(0.2518635170977563283446184976264924077045L),
    static_cast<local_float_type>(0.2701201061202457234361463802484836927464L),
    static_cast<local_float_type>(0.2881355408650536940851830262098172944660L),
    static_cast<local_float_type>(0.3057513555072237123913927136839853900197L),
    static_cast<local_float_type>(0.3228070492275195774383850177821175151772L),
    static_cast<local_float_type>(0.3391416780352496831991017115498371039332L),
    static_cast<local_float_type>(0.3545954722799571667706920253581728149226L),
    static_cast<local_float_type>(0.3690114637024459111815176585531943143085L),
    static_cast<local_float_type>(0.3822371057078773326211699424651752096338L),
    static_cast<local_float_type>(0.3941258705356621024731438508712990073773L),
    static_cast<local_float_type>(0.4045388071531719615179931506197074798366L),
    static_cast<local_float_type>(0.4133460440118967760814988295083688541739L),
    static_cast<local_float_type>(0.4204282212729730452147271372029342292548L),
    static_cast<local_float_type>(0.4256778377298582292448895379334671298297L),
    static_cast<local_float_type>(0.4290004984236535775073755577697874033289L),
    static_cast<local_float_type>(0.4303160498540635572666996674883665298558L),
    static_cast<local_float_type>(0.4295595907277138677145484170304736319185L),
    static_cast<local_float_type>(0.4266823473457215250850626209433240464185L),
    static_cast<local_float_type>(0.4216524040030023831246842156638018726147L),
    static_cast<local_float_type>(0.4144552801406973126588418824207056743782L),
    static_cast<local_float_type>(0.4050943474472003316564108983454900189419L),
    static_cast<local_float_type>(0.3935910816286283019261303507307813773886L),
    static_cast<local_float_type>(0.3799851451515116989978414766085547417497L),
    static_cast<local_float_type>(0.3643342988837623802358278078893847672838L)
  }};

  int n_val { 9 };

  const local_float_type vu_val { static_cast<local_float_type>(static_cast<local_float_type>(31) / 10) };

  for(std::size_t index { UINT8_C(0) }; index < std::tuple_size<local_ctrl_array_type>::value; ++index)
  {
    const local_float_type x_val { static_cast<local_float_type>(static_cast<local_float_type>(n_val) / 10) };

    const local_float_type jn_val { boost::math::cyl_bessel_j(vu_val, x_val) };

    ++n_val;

    using std::fabs;

    const local_float_type ratio { jn_val / ctrl_data[index] };
    const local_float_type delta { fabs(1 - ratio) };

    BOOST_TEST(delta < 128 * std::numeric_limits<local_float_type>::epsilon());
  }
}

auto main() -> int
{
  test_n<float>();
  test_n<double>();
  test_n<long double>();

  test_vu<float>();
  test_vu<double>();
  test_vu<long double>();

  return boost::report_errors();
}
