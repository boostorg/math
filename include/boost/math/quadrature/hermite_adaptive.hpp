//  Copyright Jacob Hass 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_QUADRATURE_HERMITE_ADAPTIVE_HPP
#define BOOST_MATH_QUADRATURE_HERMITE_ADAPTIVE_HPP

#ifndef BOOST_MATH_BUILD_MODULE
#include <array>
#endif
#include <boost/math/tools/config.hpp>
#include <boost/math/tools/big_constant.hpp>
#include <boost/math/quadrature/detail/quadrature_constant.hpp>
#include <boost/math/quadrature/gauss_hermite.hpp>
#include <boost/math/constants/constants.hpp>

BOOST_MATH_NAMESPACE_BEGIN namespace quadrature { namespace detail {

template <class T, unsigned N, unsigned Category>
class nested_hermite_detail;

#ifndef BOOST_HAS_FLOAT128
template <class T>
class nested_hermite_detail<T, 9, 0>
{
	using storage_type = typename quadrature_constant_category<T>::storage_type;
public:
	static std::array<storage_type, 5> const& abscissa()
	{
		static std::array<storage_type, 5> data = {
			static_cast<storage_type>(0.00000000000000000000000000000000000e+0L),
			static_cast<storage_type>(5.24033547486957629324422214267542586e-1L),
			static_cast<storage_type>(1.22474487139158894066781613219063729e+0L),
			static_cast<storage_type>(2.02323019110051571800568126491270959e+0L),
			static_cast<storage_type>(2.95921077906383755973251936666201800e+0L),
		};
		return data;
	}
	static std::array<storage_type, 5> const& weights()
	{
		static std::array<storage_type, 5> data = {
			static_cast<storage_type>(4.50147009753781857099852459214162081e-1L),
			static_cast<storage_type>(4.78694285491141513499968596079270355e-1L),
			static_cast<storage_type>(1.68118928947677659380133263766765594e-1L),
			static_cast<storage_type>(1.41731178739791068010855923375856946e-2L),
			static_cast<storage_type>(1.67088263068823530532042154206351370e-4L),
		};
		return data;
	}
};

template <class T>
class nested_hermite_detail<T, 11, 0>
{
	using storage_type = typename quadrature_constant_category<T>::storage_type;
public:
	static std::array<storage_type, 6> const& abscissa()
	{
		static std::array<storage_type, 6> data = {
			static_cast<storage_type>(0.00000000000000000000000000000000000e+0L),
			static_cast<storage_type>(5.98416053181475288258521106181433424e-1L),
			static_cast<storage_type>(1.22474487139158894066781613219063729e+0L),
			static_cast<storage_type>(1.89790257396108263954204176116036251e+0L),
			static_cast<storage_type>(2.64144797181307611211309449572581798e+0L),
			static_cast<storage_type>(3.52078826295824320880001323530450463e+0L),
		};
		return data;
	}
	static std::array<storage_type, 6> const& weights()
	{
		static std::array<storage_type, 6> data = {
			static_cast<storage_type>(5.93134883440277205224333556543570012e-1L),
			static_cast<storage_type>(4.25396823823452008284107250801753253e-1L),
			static_cast<storage_type>(1.44353425902717991125001617547241040e-1L),
			static_cast<storage_type>(1.91650634766381174067806369976096903e-2L),
			static_cast<storage_type>(7.40019801327461163725773829469289922e-4L),
			static_cast<storage_type>(4.15072848381223175908121686417118212e-6L),
		};
		return data;
	}
};

template <class T>
class nested_hermite_detail<T, 35, 0>
{
	using storage_type = typename quadrature_constant_category<T>::storage_type;
public:
	static std::array<storage_type, 18> const& abscissa()
	{
		static std::array<storage_type, 18> data = {
			static_cast<storage_type>(0.00000000000000000000000000000000000e+0L),
			static_cast<storage_type>(1.76064142082008934231396324321394786e-1L),
			static_cast<storage_type>(5.24033547486957629324422214267542586e-1L),
			static_cast<storage_type>(8.70040895352902854398280396708287299e-1L),
			static_cast<storage_type>(1.22474487139158894066781613219063729e+0L),
			static_cast<storage_type>(1.57941213484676712397458686609752476e+0L),
			static_cast<storage_type>(1.83570797517518680130876873590750620e+0L),
			static_cast<storage_type>(2.02323019110051571800568126491270959e+0L),
			static_cast<storage_type>(2.26651326205678804726062480767723173e+0L),
			static_cast<storage_type>(2.57055837658429675585125551151577383e+0L),
			static_cast<storage_type>(2.95921077906383755973251936666201800e+0L),
			static_cast<storage_type>(3.34916395371319497442641477391589433e+0L),
			static_cast<storage_type>(3.66777421594633779022842645645141602e+0L),
			static_cast<storage_type>(4.02922014050437127252735081128776073e+0L),
			static_cast<storage_type>(4.49959939831038902013915503630414605e+0L),
			static_cast<storage_type>(5.03608994447309399333789770025759935e+0L),
			static_cast<storage_type>(5.64325785788574485479784925701096654e+0L),
			static_cast<storage_type>(6.37593927098223556271250345162115991e+0L),
		};
		return data;
	}
	static std::array<storage_type, 18> const& weights()
	{
		static std::array<storage_type, 18> data = {
			static_cast<storage_type>(9.12626753636617995928126401139479640e-4L),
			static_cast<storage_type>(3.39885955855852350815382578730350360e-1L),
			static_cast<storage_type>(2.62448714887843048693127911974443123e-1L),
			static_cast<storage_type>(1.63712215557357987361086770761176012e-1L),
			static_cast<storage_type>(8.02455181473911149936384390457533299e-2L),
			static_cast<storage_type>(2.77805089085351002153867483457361232e-2L),
			static_cast<storage_type>(5.59288289114694137094119241737644188e-3L),
			static_cast<storage_type>(4.09675277203441073892165036340884399e-3L),
			static_cast<storage_type>(1.45155804251558624880236170184844013e-3L),
			static_cast<storage_type>(4.87853993044438190209571848399150440e-4L),
			static_cast<storage_type>(6.33286208056174849045580343975814230e-5L),
			static_cast<storage_type>(4.84627997370203590578492669482812971e-6L),
			static_cast<storage_type>(4.37378180409268828710955129279724929e-7L),
			static_cast<storage_type>(3.79202223923196175913039528918707655e-8L),
			static_cast<storage_type>(8.15537218169173833861604442789549133e-10L),
			static_cast<storage_type>(5.48968369484997613689306791136315866e-12L),
			static_cast<storage_type>(9.65994662785610811262528024476703647e-15L),
			static_cast<storage_type>(1.86840148945094141654219452337805717e-18L),
		};
		return data;
	}
};

#else
template <class T>
class nested_hermite_detail<T, 9, 0>
{
	using storage_type = typename quadrature_constant_category<T>::storage_type;
public:
	static std::array<storage_type, 5> const& abscissa()
	{
		static std::array<storage_type, 5> data = {
			static_cast<storage_type>(0.00000000000000000000000000000000000e+0Q),
			static_cast<storage_type>(5.24033547486957629324422214267542586e-1Q),
			static_cast<storage_type>(1.22474487139158894066781613219063729e+0Q),
			static_cast<storage_type>(2.02323019110051571800568126491270959e+0Q),
			static_cast<storage_type>(2.95921077906383755973251936666201800e+0Q),
		};
		return data;
	}
	static std::array<storage_type, 5> const& weights()
	{
		static std::array<storage_type, 5> data = {
			static_cast<storage_type>(4.50147009753781857099852459214162081e-1Q),
			static_cast<storage_type>(4.78694285491141513499968596079270355e-1Q),
			static_cast<storage_type>(1.68118928947677659380133263766765594e-1Q),
			static_cast<storage_type>(1.41731178739791068010855923375856946e-2Q),
			static_cast<storage_type>(1.67088263068823530532042154206351370e-4Q),
		};
		return data;
	}
};

template <class T>
class nested_hermite_detail<T, 11, 0>
{
	using storage_type = typename quadrature_constant_category<T>::storage_type;
public:
	static std::array<storage_type, 6> const& abscissa()
	{
		static std::array<storage_type, 6> data = {
			static_cast<storage_type>(0.00000000000000000000000000000000000e+0Q),
			static_cast<storage_type>(5.98416053181475288258521106181433424e-1Q),
			static_cast<storage_type>(1.22474487139158894066781613219063729e+0Q),
			static_cast<storage_type>(1.89790257396108263954204176116036251e+0Q),
			static_cast<storage_type>(2.64144797181307611211309449572581798e+0Q),
			static_cast<storage_type>(3.52078826295824320880001323530450463e+0Q),
		};
		return data;
	}
	static std::array<storage_type, 6> const& weights()
	{
		static std::array<storage_type, 6> data = {
			static_cast<storage_type>(5.93134883440277205224333556543570012e-1Q),
			static_cast<storage_type>(4.25396823823452008284107250801753253e-1Q),
			static_cast<storage_type>(1.44353425902717991125001617547241040e-1Q),
			static_cast<storage_type>(1.91650634766381174067806369976096903e-2Q),
			static_cast<storage_type>(7.40019801327461163725773829469289922e-4Q),
			static_cast<storage_type>(4.15072848381223175908121686417118212e-6Q),
		};
		return data;
	}
};

template <class T>
class nested_hermite_detail<T, 35, 0>
{
	using storage_type = typename quadrature_constant_category<T>::storage_type;
public:
	static std::array<storage_type, 18> const& abscissa()
	{
		static std::array<storage_type, 18> data = {
			static_cast<storage_type>(0.00000000000000000000000000000000000e+0Q),
			static_cast<storage_type>(1.76064142082008934231396324321394786e-1Q),
			static_cast<storage_type>(5.24033547486957629324422214267542586e-1Q),
			static_cast<storage_type>(8.70040895352902854398280396708287299e-1Q),
			static_cast<storage_type>(1.22474487139158894066781613219063729e+0Q),
			static_cast<storage_type>(1.57941213484676712397458686609752476e+0Q),
			static_cast<storage_type>(1.83570797517518680130876873590750620e+0Q),
			static_cast<storage_type>(2.02323019110051571800568126491270959e+0Q),
			static_cast<storage_type>(2.26651326205678804726062480767723173e+0Q),
			static_cast<storage_type>(2.57055837658429675585125551151577383e+0Q),
			static_cast<storage_type>(2.95921077906383755973251936666201800e+0Q),
			static_cast<storage_type>(3.34916395371319497442641477391589433e+0Q),
			static_cast<storage_type>(3.66777421594633779022842645645141602e+0Q),
			static_cast<storage_type>(4.02922014050437127252735081128776073e+0Q),
			static_cast<storage_type>(4.49959939831038902013915503630414605e+0Q),
			static_cast<storage_type>(5.03608994447309399333789770025759935e+0Q),
			static_cast<storage_type>(5.64325785788574485479784925701096654e+0Q),
			static_cast<storage_type>(6.37593927098223556271250345162115991e+0Q),
		};
		return data;
	}
	static std::array<storage_type, 18> const& weights()
	{
		static std::array<storage_type, 18> data = {
			static_cast<storage_type>(9.12626753636617995928126401139479640e-4Q),
			static_cast<storage_type>(3.39885955855852350815382578730350360e-1Q),
			static_cast<storage_type>(2.62448714887843048693127911974443123e-1Q),
			static_cast<storage_type>(1.63712215557357987361086770761176012e-1Q),
			static_cast<storage_type>(8.02455181473911149936384390457533299e-2Q),
			static_cast<storage_type>(2.77805089085351002153867483457361232e-2Q),
			static_cast<storage_type>(5.59288289114694137094119241737644188e-3Q),
			static_cast<storage_type>(4.09675277203441073892165036340884399e-3Q),
			static_cast<storage_type>(1.45155804251558624880236170184844013e-3Q),
			static_cast<storage_type>(4.87853993044438190209571848399150440e-4Q),
			static_cast<storage_type>(6.33286208056174849045580343975814230e-5Q),
			static_cast<storage_type>(4.84627997370203590578492669482812971e-6Q),
			static_cast<storage_type>(4.37378180409268828710955129279724929e-7Q),
			static_cast<storage_type>(3.79202223923196175913039528918707655e-8Q),
			static_cast<storage_type>(8.15537218169173833861604442789549133e-10Q),
			static_cast<storage_type>(5.48968369484997613689306791136315866e-12Q),
			static_cast<storage_type>(9.65994662785610811262528024476703647e-15Q),
			static_cast<storage_type>(1.86840148945094141654219452337805717e-18Q),
		};
		return data;
	}
};
#endif

template <class T>
class nested_hermite_detail<T, 9, 4>
{
public:
	static std::array<T, 5> const& abscissa()
	{
		static std::array<T, 5> data = {
			BOOST_MATH_HUGE_CONSTANT(T, 0, 0.0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000e+0),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 5.2403354748695762932442221426754258573055267332196519397839877359075841110961705151091995144770948515315409643180016E-1),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 1.2247448713915889406678161321906372904777526861047393192332256042572912104565073871209880155355121727067285395315485E+0),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 2.0232301911005157180056812649127095937728881849563827206653667792388738387649919136578645178601078914494691427908352E+0),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 2.9592107790638375597325193666620180010795593264837784875480642651376878185417139636185980020527293377041394245537310E+0),
		};
		return data;
	}

	static std::array<T, 5> const& weights()
	{
		static std::array<T, 5> data = {
			BOOST_MATH_HUGE_CONSTANT(T, 0, 4.5014700975378185709985245921416208148002624511901306612463793484992454865188242809103833921041988089907415578637021E-1),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 4.7869428549114151349996859607927035540342330935721469287499831892315841311494293609098852752874582783038240582505223E-1),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 1.6811892894767765938013326376676559448242187500000000000000000000000000000000000000000000000000000000000000000000000E-1),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 1.4173117873979106801085592337585694622248411180709815267179017609478151762470956417709987309141861024137850297253183E-2),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 1.6708826306882353053204215420635136979399248963265637272740133423352017483907588268627856635314984188549978153536916E-4),
		};
		return data;
	}
};

template <class T>
class nested_hermite_detail<T, 11, 4>
{
public:
	static std::array<T, 6> const& abscissa()
	{
		static std::array<T, 6> data = {
			BOOST_MATH_HUGE_CONSTANT(T, 0, 0.0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000e+0),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 5.9841605318147528825852110618143342435995483388313526091455219979467814749304959636313321357678815825865314492382376E-1),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 1.2247448713915889406678161321906372904777526861047393192332256042572912104565073871209880155355121727067285395315485E+0),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 1.8979025739610826395420417611603625118732452399848918461425108078621080030634939667826423003706935879920156169652756E+0),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 2.6414479718130761121130944957258179783821105963028599853711558148167114661612840553917987408266481384967908705810001E+0),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 3.5207882629582432088000132353045046329498291020945831371110027289716511981280924028531626039368094062831810461092306E+0),
		};
		return data;
	}

	static std::array<T, 6> const& weights()
	{
		static std::array<T, 6> data = {
			BOOST_MATH_HUGE_CONSTANT(T, 0, 5.9313488344027720522433355654357001185417175291559303964848267277186573244981732144868482235973799746019143149184394E-1),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 4.2539682382345200828410725080175325274467468262602949091836603167735504093363861661767042438790967576925513705340474E-1),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 1.4435342590271799112500161754724103957414627080044429902698882298243412895307884846114195479126091758103716556205002E-1),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 1.9165063476638117406780636997609690297394990924740626302139891228744592043586824711107535087710977556357301286982220E-2),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 7.4001980132746116372577382946928992168977856643375766906012011856380569381376820897987640613293292911664403834949730E-4),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 4.1507284838122317590812168641711821237549884242947062811305573332896570770117590646558731113126077004010761988246957E-6),
		};
		return data;
	}
};

template <class T>
class nested_hermite_detail<T, 35, 4>
{
public:
	static std::array<T, 18> const& abscissa()
	{
		static std::array<T, 18> data = {
			BOOST_MATH_HUGE_CONSTANT(T, 0, 0.0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000e+0),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 1.7606414208200893423139632432139478623867034910034896730686489962704620085778090339183538314211937450695205274173553E-1),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 5.2403354748695762932442221426754258573055267332196519397839877359075841110961705151091995144770948515315409643180016E-1),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 8.7004089535290285439828039670828729867935180668508299196826528536883954663420581529651736460388008731812886137206890E-1),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 1.2247448713915889406678161321906372904777526861047393192332256042572912104565073871209880155355121727067285395315485E+0),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 1.5794121348467671239745868660975247621536254884086614848243688726131228530398242112711305889534242104641542992260497E+0),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 1.8357079751751868013087687359075061976909637453116815465232335929253416758944938967281597298036509971364437687017007E+0),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 2.0232301911005157180056812649127095937728881849563827206653667792388738387649919136578645178601078914494691427908352E+0),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 2.2665132620567880472606248076772317290306091312725524326732753092243515437588310417694831173009663111960560737074899E+0),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 2.5705583765842967558512555115157738327980041501996314230893607349811657111776354198331042046148565575403557302432013E+0),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 2.9592107790638375597325193666620180010795593264837784875480642651376878185417139636185980020527293377041394245537310E+0),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 3.3491639537131949744264147739158943295478820804503226186002355977895529723428458278052678084974479394654704651058175E+0),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 3.6677742159463377902284264564514160156250000000000000000000000000000000000000000000000000000000000000000000000000000E+0),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 4.0292201405043712725273508112877607345581054687500000000000000000000000000000000000000000000000000000000000000000000E+0),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 4.4995993983103890201391550363041460514068603515625000000000000000000000000000000000000000000000000000000000000000000E+0),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 5.0360899444730939933378977002575993537902832031250000000000000000000000000000000000000000000000000000000000000000000E+0),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 5.6432578578857448547978492570109665393829345703125000000000000000000000000000000000000000000000000000000000000000000E+0),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 6.3759392709822355627125034516211599111557006835937500000000000000000000000000000000000000000000000000000000000000000E+0),
		};
		return data;
	}

	static std::array<T, 18> const& weights()
	{
		static std::array<T, 18> data = {
			BOOST_MATH_HUGE_CONSTANT(T, 0, 9.1262675363661799592812640113947963982354849583078867607179254607992076718493219032645008384315989625425637447308587E-4),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 3.3988595585585235081538257873035036027431488042306009193926878794995700800105746629094870263384856060331462983025518E-1),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 2.6244871488784304869312791197444312274456024176664587637018064358264383035332355153807758391125459360734445781194714E-1),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 1.6371221555735798736108677076117601245641708371463528098737521211018447441306034401858610949457863659556333152112197E-1),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 8.0245518147391114993638439045753329992294311527701869100256195952300158736341378410959293406446839926361035204522463E-2),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 2.7780508908535100215386748345736123155802488333902425922902502403231021116243643383336903709643283695874563015945686E-2),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 5.5928828911469413709411924173764418810606002812402756597784931802369248665121849776215779386217334059139218270359601E-3),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 4.0967527720344107389216503634088439866900444032971466818144327105571031295530999834734827312700700962897100480096951E-3),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 1.4515580425155862488023617018484401341993361715103903021539926563382916024124574468831023646957712685417444984918072E-3),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 4.8785399304443819020957184839915043994551524528377071185797695296424221502638237639631787286179612033995413366497870E-4),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 6.3328620805617484904558034397581423036172054712453293491131880271678378912721060971432700233328117725411081758633383E-5),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 4.8462799737020359057849266948281297118228394541001410439736774642511947544892507590833492193160142399219403711387527E-6),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 4.3737818040926882871095512927972492889239219952956150804595523432936794879497010284992082351406857988242649106642875E-7),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 3.7920222392319617591303952891870765462556391879817432604254114992618848516878985971086122514202756511163032049645622E-8),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 8.1553721816917383386160444278954913255308412090132784323430891333338798015479989086969757249027311624011451716421113E-10),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 5.4896836948499761368930679113631586644658777451857664928119141365956562574546864042636316355538198477077878446029291E-12),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 9.6599466278561081126252802447670364703320548829356083214761659972413573956085002492713563170623702323459479013716496E-15),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 1.8684014894509414165421945233780571719104960824678226110638260008836965693712045529172540496186337464451683871921787E-18),
		};
		return data;
	}
};

} // namespace detail

template <class Real, unsigned N, class Policy = BOOST_MATH_NAMESPACE::policies::policy<> >
class nested_hermite : public detail::nested_hermite_detail<Real, N, detail::quadrature_constant_category<Real>::value>
{
   using base = detail::nested_hermite_detail<Real, N, detail::quadrature_constant_category<Real>::value>;

   // Abscissas computed on demand are NaN when the recurrence overflows Real.
   static bool overflowed()
   {
      return !(BOOST_MATH_NAMESPACE::isfinite)(static_cast<Real>(base::abscissa().back()));
   }

   static Real overflow_error()
   {
      return policies::raise_evaluation_error(
         "boost::math::quadrature::nested_hermite<%1%>::integrate",
         "The Hermite recurrence overflows with this many points; use fewer points or a type with a wider exponent range.",
         std::numeric_limits<Real>::quiet_NaN(), Policy());
   }

public:
   // Integrates f(x) exp(-x^2) over the real line.
   template <class F>
   static auto integrate(F f, Real* error = nullptr, Real* pL1 = nullptr)->decltype(std::declval<F>()(std::declval<Real>()))
   {
      // In many math texts, K represents the field of real or complex numbers.
      // Too bad we can't put blackboard bold into C++ source!
      using K = decltype(f(Real(0)));
      static_assert(!std::is_integral<K>::value,
                    "The return type cannot be integral, it must be either a real or complex floating point type.");
      using std::abs;
      if (overflowed())
      {
         Real overflow = overflow_error();
         if (pL1)
            *pL1 = overflow;
		 if (error)
		 	*error = overflow;
         return static_cast<K>(overflow);
      }
      unsigned non_zero_start = 1;
      K result = Real(0);
      if (N & 1)
      {
         result = f(Real(0)) * static_cast<Real>(base::weights()[0]);
      }
      else
      {
         non_zero_start = 0;
      }
      Real L1 = abs(result);
      for (unsigned i = non_zero_start; i < base::abscissa().size(); ++i)
      {
         K fp = f(static_cast<Real>(base::abscissa()[i]));
         K fm = f(-static_cast<Real>(base::abscissa()[i]));
         result += (fp + fm) * static_cast<Real>(base::weights()[i]);
         L1 += (abs(fp) + abs(fm)) * static_cast<Real>(base::weights()[i]);
      }

	  if (pL1)
         *pL1 = L1;

	  if (error)
	  {
		unsigned hermite_order = (N - 1) / 2;
		unsigned hermite_non_zero_start = 1;
		K hermite_result = Real(0);
		if (hermite_order & 1)
		{
			hermite_result = f(Real(0)) * static_cast<Real>(gauss_hermite<Real, (N - 1) / 2>::weights()[0]);
		}
		else
		{
			hermite_non_zero_start = 0;
		}
		for (unsigned i = hermite_non_zero_start; i < gauss_hermite<Real, (N - 1) / 2>::abscissa().size(); ++i)
		{
			K fp = f(static_cast<Real>(gauss_hermite<Real, (N - 1) / 2>::abscissa()[i]));
			K fm = f(-static_cast<Real>(gauss_hermite<Real, (N - 1) / 2>::abscissa()[i]));
			hermite_result += (fp + fm) * static_cast<Real>(gauss_hermite<Real, (N - 1) / 2>::weights()[i]);
		}
		*error = (std::max)(static_cast<Real>(abs(hermite_result - result)), static_cast<Real>(abs(result * tools::epsilon<Real>() * Real(2))));
	  }
      return result;
   }

   // Explicit zero and norm for vector- and matrix-valued integrands.
   template <class F, class Norm>
   static auto integrate(F f, const decltype(std::declval<F>()(std::declval<Real>()))& zero, Norm norm, Real* error = nullptr, Real* pL1 = nullptr)
      ->decltype(static_cast<Real>(norm(f(Real(0)))), f(Real(0)))
   {
      using K = decltype(f(Real(0)));
      static_assert(!std::is_integral<K>::value, "The return type cannot be integral.");
      if (overflowed())
      {
         // Preserve the scalar error policy, including non-throwing policies,
         // and use scalar multiplication to produce a result with the right shape.
         Real overflow = overflow_error();
         if (error)
            *error = overflow;
         if (pL1)
            *pL1 = overflow;
         K result = zero * overflow;
         return result;
      }
      unsigned non_zero_start = 1;
      K result = zero;
      if (N & 1)
      {
         result = f(Real(0)) * static_cast<Real>(base::weights()[0]);
      }
      else
      {
         non_zero_start = 0;
      }
      Real L1 = static_cast<Real>(norm(result));
      for (unsigned i = non_zero_start; i < base::abscissa().size(); ++i)
      {
         K fp = f(static_cast<Real>(base::abscissa()[i]));
         K fm = f(-static_cast<Real>(base::abscissa()[i]));
         result += (fp + fm) * static_cast<Real>(base::weights()[i]);
         L1 += (static_cast<Real>(norm(fp)) + static_cast<Real>(norm(fm))) * static_cast<Real>(base::weights()[i]);
      }

	  if (pL1)
         *pL1 = L1;

	  if (error)
	  {
		unsigned hermite_non_zero_start = 1;
		unsigned hermite_order = (N - 1) / 2;
		K hermite_result = zero;
		if (hermite_order & 1)
		{
			hermite_result = f(Real(0)) * static_cast<Real>(gauss_hermite<Real, (N - 1) / 2>::weights()[0]);
		}
		else
		{
			hermite_non_zero_start = 0;
		}
		for (unsigned i = hermite_non_zero_start; i < gauss_hermite<Real, (N - 1) / 2>::abscissa().size(); ++i)
		{
			K fp = f(static_cast<Real>(gauss_hermite<Real, (N - 1) / 2>::abscissa()[i]));
			K fm = f(-static_cast<Real>(gauss_hermite<Real, (N - 1) / 2>::abscissa()[i]));
			hermite_result += (fp + fm) * static_cast<Real>(gauss_hermite<Real, (N - 1) / 2>::weights()[i]);
		}
		*error = (std::max)(static_cast<Real>(norm(result - hermite_result)), static_cast<Real>(norm(result * tools::epsilon<Real>() * Real(2))));
	  }
      return result;
   }
};

} // namespace quadrature
BOOST_MATH_NAMESPACE_END

#endif // BOOST_MATH_QUADRATURE_HERMITE_ADAPTIVE_HPP
