//  (C) Copyright Matt Borland 2022.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <array>
#include <boost/math/tools/color_maps.hpp>

template <typename Real>
void test()
{
    // Test for constexpr-ness
    constexpr boost::math::tools::smooth_cool_warm_color_map<Real> cm;
    static_assert(cm(0.5)[0] > 0);

    // Same rgb values as plasma
    constexpr std::array<Real, 256> r {0.18500126283629117,0.1911189044113562,0.19712864908128114,0.2030416233571155,0.20887268107458867,0.21463072331602523,0.22032157730150803,0.2259575605729931,0.23155039796588678,0.23710641646603325,0.24262940495715968,0.24812510063567056,0.25359583439603345,0.259047386098491,0.26447987886309154,0.2698981057405734,0.2753043830681906,0.28069900692141886,0.28608641029403303,0.2914638822295278,0.29683511768078935,0.3021991776471102,0.307559383374685,0.31291664864816615,0.31826993831935324,0.3236225760753923,0.3289709965431237,0.33431870303843003,0.33966373680896456,0.3450078077262057,0.35034847407826647,0.3556874153529715,0.36102581201779066,0.3663615341003603,0.37169632394652985,0.3770276964703286,0.382357946820387,0.38768438207138256,0.3930062860759083,0.3983273716748813,0.4036425817932003,0.40895431552409317,0.41426001622424247,0.4195623952577179,0.4248563616485503,0.4301441775161412,0.4354268819265687,0.4406999402134,0.44596644643576516,0.45122228710713624,0.4564697629199718,0.4617056654598785,0.4669310023070961,0.47214475738324385,0.4773458407394621,0.48253461632433253,0.4877090328365206,0.49286905392912533,0.49801338714995985,0.5031412603622233,0.5082531828521651,0.5133462124834312,0.5184220260592395,0.5234778559822143,0.5285143492877827,0.5335306956538576,0.5385244789001884,0.5434980954707963,0.5484481027774952,0.5533750875909487,0.5582766534405778,0.5631554062153734,0.5680081900745759,0.5728333975138157,0.5776342744306816,0.5824076737546731,0.5871533984089229,0.5918701205096012,0.5965604367679973,0.6012201835310419,0.6058512183412395,0.6104521263820015,0.6150225220246975,0.6195628616073623,0.6240719801812181,0.6285513580399021,0.6329989221793823,0.6374162963015114,0.6418004809769666,0.6461536800309777,0.650476135769243,0.6547651168217998,0.6590245256639604,0.6632516551223815,0.6674477291126537,0.6716108031251352,0.6757429658756737,0.6798454057970371,0.6839159522860683,0.6879568949603344,0.6919652822236589,0.6959450042288188,0.6998941469583252,0.7038136863269572,0.7077036076142078,0.7115648684573126,0.7153974874768007,0.7192003105763611,0.7229763232162184,0.7267244865470727,0.7304437173416366,0.7341377518333918,0.7378037499188499,0.7414454244046753,0.7450587348350122,0.7486483346580772,0.7522104746641185,0.755748680542265,0.7592629017974487,0.7627520603974456,0.766216825373388,0.7696582925955131,0.773077107895631,0.7764713899964052,0.7798428195281752,0.7831920638854646,0.7865182091233281,0.7898228934175395,0.7931051428523151,0.7963656261858302,0.7996034565968984,0.8028201310987728,0.8060159189746301,0.8091884512009491,0.8123413974246995,0.8154719410234268,0.8185825533178849,0.8216705444481089,0.8247373765508464,0.8277822680188923,0.8308060259837033,0.8338091392025193,0.8367899778785102,0.8397489706421578,0.8426854482941508,0.8456008173688659,0.8484934838710892,0.8513631273192183,0.8542102862752687,0.8570334035212495,0.8598348399896373,0.8626110724885258,0.8653643579173034,0.8680923147108822,0.8707954368901418,0.8734743717670276,0.8761266641066531,0.8787535738894046,0.8813526771710705,0.8839261632148947,0.8864707162026111,0.8889878746971145,0.8914771762860599,0.8939362684731502,0.896366404995977,0.8987661885216094,0.9011348456797987,0.9034710462469783,0.9057762624345291,0.9080481396866358,0.9102863527112363,0.9124909968821041,0.9146608551755517,0.9167960662705978,0.9188933135427245,0.9209542277904884,0.9229781644169556,0.9249630587811338,0.9269098852966564,0.928815526839391,0.9306829366142825,0.9325070229302802,0.9342897265833519,0.9360290542859102,0.9377243972740386,0.9393770933973468,0.9409842661949558,0.9425456696997273,0.9440594310662802,0.9455271344382453,0.9469451155697922,0.9483136189409812,0.9496329509923357,0.9509013387619204,0.9521193680507284,0.9532833326102016,0.9543957893212712,0.9554530885899511,0.9564564274893136,0.9574039803720095,0.9582930778369493,0.9591252451658657,0.9598979045875052,0.9606124721276129,0.9612663555645452,0.9618587615061226,0.9623887555724971,0.9628559105126825,0.9632594414845331,0.9635970174364671,0.9638686962955997,0.964073330664618,0.9642088982457305,0.9642762050566095,0.9642722292530789,0.9641978734579372,0.9640500241140576,0.9638285003224122,0.9635333683894992,0.9631614522931872,0.962712538770477,0.9621866932888067,0.9615816256657105,0.960894345402275,0.9601234461397462,0.959270237581007,0.958333448735586,0.9573126236007495,0.9562062303425974,0.9550135580216242,0.9537297610630108,0.9523531724953679,0.9508874276947369,0.9493337129963176,0.9476874510448859,0.9459406787196255,0.944103300221897,0.9421780061800054,0.9401426587007371,0.9380209885848168,0.9358000817994506,0.9334816723713976,0.9310720730839354,0.9285606186650813,0.9259594465149981,0.92327053647064,0.9204816890380406,0.9176163146663037,0.9146765909646813,0.9116700174610778,0.9086068264551928,0.9055122790291467,0.9024220298268155,0.8994080090565665,0.8965597700566809,0.894058310302958};
    constexpr std::array<Real, 256> g {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,3.874898278890354e-05,0.007471928506941776,0.015271185178250529,0.023435206026773358,0.03197423156761071,0.04087712538541078,0.04941897732161648,0.057376591831635744,0.06489209701289092,0.07205664723234939,0.07893485649247015,0.08557896376179144,0.092021751309049,0.09829643543437619,0.10442077591051005,0.11041573939018576,0.11629982277048792,0.12207972438632506,0.12777061405330994,0.13337789223781238,0.1389126551987294,0.14437902230001867,0.14978346343836335,0.15513277880704496,0.16042886120342012,0.1656775394015491,0.1708811060968512,0.17604431343185956,0.1811675458797579,0.1862549283640475,0.19131016034079534,0.19633283758483017,0.20132751808066945,0.20629365678022926,0.21123542314263571,0.2161517783484883,0.22104657452100646,0.22592119959006707,0.23077446982485766,0.23561210642337616,0.24043047593352906,0.24523511752950322,0.2500223780012234,0.25479746355598837,0.25956029374277906,0.26431022207991034,0.2690513274100934,0.2737818849425285,0.27850483541745347,0.28321834594018785,0.28792637941737553,0.29262694015300356,0.2973217868017699,0.30201349640843383,0.3067011675000855,0.311386513443518,0.3160685293738308,0.3207499291278433,0.3254307091326806,0.33011024027988417,0.3347921147847294,0.3394741977280383,0.3441601465946952,0.3488478095363264,0.3535397891569865,0.3582349147529724,0.36293457844520105,0.3676412168097567,0.3723536223699429,0.3770733241910323,0.38179911504484987,0.3865345827145719,0.3912774119706224,0.3960300112396703,0.40079369284363137,0.4055661645442056,0.4103509864051946,0.4151458393212769,0.41995530212234555,0.4247759985436038,0.429609255031397,0.4344584280288323,0.4393211443376441,0.4441999304447296,0.4490924088928135,0.4540021452370879,0.4589267217483421,0.46386948406587364,0.4688297189338491,0.4738069993796045,0.47880391802811756,0.4838189959751383,0.4888528043494526,0.4939069187369674,0.4989795927548926,0.5040742179411092,0.5091882394741989,0.5143242934939541,0.519480841366499,0.5246605189695782,0.529860774144959,0.5350842629867526,0.5403293668305715,0.5455974579589155,0.5508899121603867,0.5562040647384632,0.5615426592562283,0.5669040105140666,0.5722919147836466,0.577701610932102,0.5831365219782882,0.5885970455193761,0.5940804364977416,0.5995905102669348,0.6051244958149699,0.6106842151090096,0.6162678665363043,0.6218768875923327,0.627512735634536,0.6331745970184661,0.6388623236974159,0.6445750605508336,0.6503137041961821,0.656078387432029,0.6618685638537222,0.6676857251951027,0.6735279399230053,0.6793971833342757,0.6852914897975773,0.6912138827190051,0.6971603311791722,0.7031323629811758,0.7091335377504808,0.7151587894894975,0.7212111928588061,0.727288688382748,0.7333943830627802,0.7395241540309135,0.7456795754640698,0.7518622221225512,0.7580699614625952,0.7643039444114984,0.7705640387323482,0.7768514277320878,0.7831628955092512,0.7894990412104851,0.7958624980435167,0.8022480152561944,0.8086588498956894,0.8150957971478479,0.8215591141377789,0.8280444859186529,0.8345532077754317,0.8410839962117271,0.8476425831870477,0.8542225467143049,0.8608205474618611,0.8674491093347713,0.8740937063216282,0.8807618316376675,0.8874501087536273,0.8941562294005477,0.900886023729324,0.9076310541062192,0.9143938321288543,0.9211760913964651,0.9279732997552812,0.9347791439177715,0.9415965873644457,0.9484236309855288,0.9552511244435956,0.9620729727476299,0.9688661897650951,0.975611979389701,0.9822535793047805};
    constexpr std::array<Real, 256> b {0.5300734481832133,0.5352444603525081,0.5401496873258278,0.5448207274312517,0.5492861078596931,0.5535696573805636,0.5576886141014704,0.5616609964897479,0.5655049239276769,0.5692314014828168,0.5728505620174764,0.5763702949099744,0.5797988437955895,0.5831429817760583,0.5864060347449428,0.5895948785151374,0.5927123254827708,0.5957607906567519,0.5987448538514638,0.6016649977484743,0.6045237142774258,0.6073215647568714,0.610060282216021,0.6127415605257925,0.6153640224697419,0.6179290152489123,0.6204342239271381,0.6228829350718792,0.625271898043606,0.6276013187112769,0.6298710242222681,0.6320785940805594,0.6342255823320162,0.6363078797194855,0.6383265475601053,0.6402785532144175,0.6421628771547095,0.6439785694946553,0.6457230933084593,0.6473958766593488,0.6489950538473048,0.6505184306872968,0.6519642297582939,0.653331171675656,0.6546165641245861,0.655818742117949,0.65693799619788,0.6579697282780315,0.658913466231862,0.6597667197389312,0.6605289119772367,0.6611976600287218,0.6617721467945837,0.6622494877377035,0.6626294254823041,0.662909148803712,0.6630875095523548,0.6631645872632188,0.66313835041097,0.6630068027529852,0.6627708908204921,0.6624287113415619,0.6619800968378087,0.6614242681940923,0.6607599298174441,0.6599874350032627,0.6591055981902831,0.6581161801158861,0.6570196633950978,0.6558144934883569,0.6545022813523073,0.653083356821114,0.651560454782398,0.6499322013289422,0.6482001697143078,0.6463672240776391,0.6444334413675159,0.6424018048025009,0.6402732833430833,0.6380499641554425,0.6357347102857865,0.6333297196772426,0.6308373994447983,0.6282600961681855,0.6256021241199032,0.6228641296151941,0.6200505178730409,0.6171638506919247,0.6142086153036438,0.6111870741109839,0.6081004450341899,0.604956300073173,0.6017530389476595,0.5984983035745908,0.5951924266541045,0.5918411070177483,0.5884465014386215,0.585011735976789,0.5815395577903398,0.5780332012739695,0.5744974521016776,0.570932506897979,0.5673441858884501,0.5637335767466599,0.5601017463095771,0.5564545500366853,0.552792120928892,0.5491173170199668,0.545432268052493,0.5417408483359869,0.5380431032601405,0.5343400883294095,0.5306366819504315,0.5269329923911534,0.5232298863541315,0.5195274755763871,0.5158306312483361,0.512137399653064,0.5084498297525021,0.50476878993514,0.501093374729771,0.49742746464864046,0.49377016285413927,0.4901233238126223,0.4864869893410944,0.48285922068624404,0.47924185046887885,0.4756350094982365,0.4720385174170281,0.4684525186264465,0.4648778196983237,0.4613125719281808,0.4577605762350459,0.45421788462774554,0.45068557847021135,0.4471634196988565,0.4436495949958711,0.44014585494792785,0.43665038594107897,0.43316392259011754,0.4296875651350407,0.42621737099556123,0.42275507747129093,0.41929886852490544,0.41584947056189187,0.4124050878709716,0.4089674228462569,0.405535554348405,0.4021065453700518,0.3986820908198349,0.39526040918192384,0.3918421707163013,0.38842457298679345,0.385011298212545,0.38159840245673277,0.3781859518708741,0.37477461745019225,0.37136257807629336,0.36795049364876287,0.3645365324577698,0.3611223530824753,0.35770601202584773,0.3542865122627056,0.35086450022751575,0.34744014111882704,0.34401206951565483,0.34058041413685886,0.33714480107622696,0.3337042204871271,0.33025867912310986,0.3268087806988829,0.32335261768909423,0.31989179322193145,0.31642542736448936,0.31295201716155624,0.30947360171764265,0.30598704730811505,0.3024949522990142,0.2989963432701857,0.29549074588061,0.29197919964809826,0.28845817867177614,0.28492969964810494,0.28139526496601874,0.2778537072388743,0.27430383727935276,0.2707481638305188,0.26718454740444175,0.2636154464355518,0.2600356289949842,0.2564496036120887,0.25285703768193385,0.24925763348723312,0.2456517661988547,0.24204016134332634,0.23842311536200755,0.2348003777648286,0.23117419840272618,0.22754413813760505,0.22390965714476324,0.2202708516770076,0.21662943633800064,0.21298749726033284,0.20934668932625894,0.2057069050257183,0.20207065907463867,0.1984392920478823,0.19481481252537644,0.19119869381163737,0.18759494102079383,0.18400288978976168,0.1804297764591602,0.17687797122612808,0.17334762525376426,0.16984714768059186,0.16637921636618158,0.16295000971437262,0.15956346939328955,0.15622826644396326,0.15294946537452883,0.14973455824735152,0.1465937236780026,0.143534307799379,0.1405644983355272,0.1376912476803779,0.1349304220039219,0.13229401208376287,0.129792103544935,0.12743980433391283,0.12524775464689872,0.12322300321471347,0.12137867908430083,0.11973149121600749,0.118292153356203,0.11706646197810774,0.11605402308057439,0.11527104051689083,0.11471851258263263,0.11437731754446762,0.1142616518784757,0.11434747323584868,0.1146183402870086,0.11505417484570082,0.11561141583641374,0.11625593038842527,0.1169286187796818,0.11754606974768209,0.11802695488214923,0.11824175235960144,0.11802059587857097,0.11712470175631184,0.1152182936223895,0.11179479561713918,0.10607343156231328,0.09668148446675398,0.0810687655704728};

    constexpr boost::math::tools::custom_color_map<Real> user_plasma(r, g, b);
    constexpr boost::math::tools::plasma_color_map<Real> plasma;

    static_assert(user_plasma(0.5)[0] == plasma(0.5)[0]);
}

#if !defined(BOOST_MATH_NO_CONSTEXPR_DETECTION) && !defined(BOOST_MATH_USING_BUILTIN_CONSTANT_P)
int main()
{
    test<float>();
    test<double>();
    test<long double>();

    return 0;
}
#else
int main()
{
    return 0;
}
#endif
