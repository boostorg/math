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
			static_cast<storage_type>(5.24033547486957645148383913594834662e-1L),
			static_cast<storage_type>(1.22474487139158904909864203735294570e+0L),
			static_cast<storage_type>(2.02323019110051565920832089517988354e+0L),
			static_cast<storage_type>(2.95921077906383772231113850053493001e+0L),
		};
		return data;
	}
	static std::array<storage_type, 5> const& weights()
	{
		static std::array<storage_type, 5> data = {
			static_cast<storage_type>(4.50147009753781848202709202118386078e-1L),
			static_cast<storage_type>(4.78694285491141488088899111869391248e-1L),
			static_cast<storage_type>(1.68118928947677671965950845303360654e-1L),
			static_cast<storage_type>(1.41731178739791059714177897487213826e-2L),
			static_cast<storage_type>(1.67088263068823521461393689906267630e-4L),
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
			static_cast<storage_type>(5.98416053181475309272966877992894546e-1L),
			static_cast<storage_type>(1.22474487139158904909864203735294570e+0L),
			static_cast<storage_type>(1.89790257396108263540210714118268341e+0L),
			static_cast<storage_type>(2.64144797181307620627699962863523315e+0L),
			static_cast<storage_type>(3.52078826295824305883026642048068890e+0L),
		};
		return data;
	}
	static std::array<storage_type, 6> const& weights()
	{
		static std::array<storage_type, 6> data = {
			static_cast<storage_type>(5.93134883440277258808275654555991068e-1L),
			static_cast<storage_type>(4.25396823823451990761884375950046286e-1L),
			static_cast<storage_type>(1.44353425902718001059769312197649986e-1L),
			static_cast<storage_type>(1.91650634766381190126395802089337198e-2L),
			static_cast<storage_type>(7.40019801327461178819978561117164009e-4L),
			static_cast<storage_type>(4.15072848381223183266747482990200902e-6L),
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
			static_cast<storage_type>(1.76064142082008935029745686285161957e-1L),
			static_cast<storage_type>(5.24033547486957645148383913594834662e-1L),
			static_cast<storage_type>(8.70040895352902900134956696280610044e-1L),
			static_cast<storage_type>(1.22474487139158904909864203735294570e+0L),
			static_cast<storage_type>(1.57941213484676708572036758299947500e+0L),
			static_cast<storage_type>(1.83570797517518687377303661427767510e+0L),
			static_cast<storage_type>(2.02323019110051565920832089517988354e+0L),
			static_cast<storage_type>(2.26651326205678802746598617543855944e+0L),
			static_cast<storage_type>(2.57055837658429670911306400444247515e+0L),
			static_cast<storage_type>(2.95921077906383772231113850053493001e+0L),
			static_cast<storage_type>(3.34916395371319497736738502677139580e+0L),
			static_cast<storage_type>(3.66777421594633786003793251745765051e+0L),
			static_cast<storage_type>(4.02922014050437136479350497812196521e+0L),
			static_cast<storage_type>(4.49959939831038880288429511940039540e+0L),
			static_cast<storage_type>(5.03608994447309396868596432263094693e+0L),
			static_cast<storage_type>(5.64325785788574506280375428303786623e+0L),
			static_cast<storage_type>(6.37593927098223595171270375072721786e+0L),
		};
		return data;
	}
	static std::array<storage_type, 18> const& weights()
	{
		static std::array<storage_type, 18> data = {
			static_cast<storage_type>(9.12626753636618015784617239051855639e-4L),
			static_cast<storage_type>(3.39885955855852376532651613014502461e-1L),
			static_cast<storage_type>(2.62448714887843068475168131366598856e-1L),
			static_cast<storage_type>(1.63712215557357978737966786225489413e-1L),
			static_cast<storage_type>(8.02455181473911204774905924216860908e-2L),
			static_cast<storage_type>(2.77805089085350998625061959108534475e-2L),
			static_cast<storage_type>(5.59288289114694100849109490561159433e-3L),
			static_cast<storage_type>(4.09675277203441034424345187225229409e-3L),
			static_cast<storage_type>(1.45155804251558622875835195134407927e-3L),
			static_cast<storage_type>(4.87853993044438190544096198508849397e-4L),
			static_cast<storage_type>(6.33286208056174798926698340752153159e-5L),
			static_cast<storage_type>(4.84627997370203600119434942907382286e-6L),
			static_cast<storage_type>(4.37378180409268808882839681384931432e-7L),
			static_cast<storage_type>(3.79202223923196151758999511288434096e-8L),
			static_cast<storage_type>(8.15537218169173874335697571164954813e-10L),
			static_cast<storage_type>(5.48968369484997636702024987294648685e-12L),
			static_cast<storage_type>(9.65994662785610740367890572257052316e-15L),
			static_cast<storage_type>(1.86840148945094127438034771893579257e-18L),
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
			static_cast<storage_type>(5.24033547486957645148383913594834662e-1Q),
			static_cast<storage_type>(1.22474487139158904909864203735294570e+0Q),
			static_cast<storage_type>(2.02323019110051565920832089517988354e+0Q),
			static_cast<storage_type>(2.95921077906383772231113850053493001e+0Q),
		};
		return data;
	}
	static std::array<storage_type, 5> const& weights()
	{
		static std::array<storage_type, 5> data = {
			static_cast<storage_type>(4.50147009753781848202709202118386078e-1Q),
			static_cast<storage_type>(4.78694285491141488088899111869391248e-1Q),
			static_cast<storage_type>(1.68118928947677671965950845303360654e-1Q),
			static_cast<storage_type>(1.41731178739791059714177897487213826e-2Q),
			static_cast<storage_type>(1.67088263068823521461393689906267630e-4Q),
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
			static_cast<storage_type>(5.98416053181475309272966877992894546e-1Q),
			static_cast<storage_type>(1.22474487139158904909864203735294570e+0Q),
			static_cast<storage_type>(1.89790257396108263540210714118268341e+0Q),
			static_cast<storage_type>(2.64144797181307620627699962863523315e+0Q),
			static_cast<storage_type>(3.52078826295824305883026642048068890e+0Q),
		};
		return data;
	}
	static std::array<storage_type, 6> const& weights()
	{
		static std::array<storage_type, 6> data = {
			static_cast<storage_type>(5.93134883440277258808275654555991068e-1Q),
			static_cast<storage_type>(4.25396823823451990761884375950046286e-1Q),
			static_cast<storage_type>(1.44353425902718001059769312197649986e-1Q),
			static_cast<storage_type>(1.91650634766381190126395802089337198e-2Q),
			static_cast<storage_type>(7.40019801327461178819978561117164009e-4Q),
			static_cast<storage_type>(4.15072848381223183266747482990200902e-6Q),
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
			static_cast<storage_type>(1.76064142082008935029745686285161957e-1Q),
			static_cast<storage_type>(5.24033547486957645148383913594834662e-1Q),
			static_cast<storage_type>(8.70040895352902900134956696280610044e-1Q),
			static_cast<storage_type>(1.22474487139158904909864203735294570e+0Q),
			static_cast<storage_type>(1.57941213484676708572036758299947500e+0Q),
			static_cast<storage_type>(1.83570797517518687377303661427767510e+0Q),
			static_cast<storage_type>(2.02323019110051565920832089517988354e+0Q),
			static_cast<storage_type>(2.26651326205678802746598617543855944e+0Q),
			static_cast<storage_type>(2.57055837658429670911306400444247515e+0Q),
			static_cast<storage_type>(2.95921077906383772231113850053493001e+0Q),
			static_cast<storage_type>(3.34916395371319497736738502677139580e+0Q),
			static_cast<storage_type>(3.66777421594633786003793251745765051e+0Q),
			static_cast<storage_type>(4.02922014050437136479350497812196521e+0Q),
			static_cast<storage_type>(4.49959939831038880288429511940039540e+0Q),
			static_cast<storage_type>(5.03608994447309396868596432263094693e+0Q),
			static_cast<storage_type>(5.64325785788574506280375428303786623e+0Q),
			static_cast<storage_type>(6.37593927098223595171270375072721786e+0Q),
		};
		return data;
	}
	static std::array<storage_type, 18> const& weights()
	{
		static std::array<storage_type, 18> data = {
			static_cast<storage_type>(9.12626753636618015784617239051855639e-4Q),
			static_cast<storage_type>(3.39885955855852376532651613014502461e-1Q),
			static_cast<storage_type>(2.62448714887843068475168131366598856e-1Q),
			static_cast<storage_type>(1.63712215557357978737966786225489413e-1Q),
			static_cast<storage_type>(8.02455181473911204774905924216860908e-2Q),
			static_cast<storage_type>(2.77805089085350998625061959108534475e-2Q),
			static_cast<storage_type>(5.59288289114694100849109490561159433e-3Q),
			static_cast<storage_type>(4.09675277203441034424345187225229409e-3Q),
			static_cast<storage_type>(1.45155804251558622875835195134407927e-3Q),
			static_cast<storage_type>(4.87853993044438190544096198508849397e-4Q),
			static_cast<storage_type>(6.33286208056174798926698340752153159e-5Q),
			static_cast<storage_type>(4.84627997370203600119434942907382286e-6Q),
			static_cast<storage_type>(4.37378180409268808882839681384931432e-7Q),
			static_cast<storage_type>(3.79202223923196151758999511288434096e-8Q),
			static_cast<storage_type>(8.15537218169173874335697571164954813e-10Q),
			static_cast<storage_type>(5.48968369484997636702024987294648685e-12Q),
			static_cast<storage_type>(9.65994662785610740367890572257052316e-15Q),
			static_cast<storage_type>(1.86840148945094127438034771893579257e-18Q),
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
			BOOST_MATH_HUGE_CONSTANT(T, 0, 5.2403354748695764514838391359483466233420304283870559387493001948377007156468753268890967128171616186758390430922461E-1),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 1.2247448713915890490986420373529456959829737403283350642163462836254801887286575132699297165523201174092973006133071E+0),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 2.0232301911005156592083208951798835405134106140896567635236468204669731308178641289345019964860533407648486493612883E+0),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 2.9592107790638377223111385005349300110208537494323252529925306245788685818262550907459729341537440162617289883378384E+0),
		};
		return data;
	}

	static std::array<T, 5> const& weights()
	{
		static std::array<T, 5> data = {
			BOOST_MATH_HUGE_CONSTANT(T, 0, 4.5014700975378184820270920211838607817080621107870149287969721647058064370565896669840016679074150232793062187283367E-1),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 4.7869428549114148808889911186939124791184538408822216991207114092308933495846905698619426642046228614813403845981908E-1),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 1.6811892894767767196595084530336065419488951480835650267610642993184795382299763237973786717032266474137652798604306E-1),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 1.4173117873979105971417789748721382576345635622234006469045082247136049406856569411227856966286976697623746604527599E-2),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 1.6708826306882352146139368990626763029108800303013860983263358909198225436334856111525441682965395701378782533471293E-4),
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
			BOOST_MATH_HUGE_CONSTANT(T, 0, 5.9841605318147530927296687799289454626094804015976255953985507975174274925803277728925294454179525110610490126867802E-1),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 1.2247448713915890490986420373529456959829737403283350642163462836254801887286575132699297165523201174092973006133071E+0),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 1.8979025739610826354021071411826834127738342522530396192617037125371091012698240116851970926303633820545540169918927E+0),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 2.6414479718130762062769996286352331528623576200143594645233715091970301483443729850311828753925194017174096186142413E+0),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 3.5207882629582430588302664204806889048012197501881388362535584016249489542364317847957314772150851753276449384376375E+0),
		};
		return data;
	}

	static std::array<T, 6> const& weights()
	{
		static std::array<T, 6> data = {
			BOOST_MATH_HUGE_CONSTANT(T, 0, 5.9313488344027725880827565455599106770741524283311255532383633229064743641216240317906845506544762659680270176185142E-1),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 4.2539682382345199076188437595004628587686158487056349070112202428161586044351614949678153922032669088937885294083687E-1),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 1.4435342590271800105976931219764998580084754151084013842346311261163908076948487806427754039755453935908769942257101E-1),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 1.9165063476638119012639580208933719848935973842809912509686439732217075435583074014959405216662795314500661337244696E-2),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 7.4001980132746117881997856111716400939931017028778121727523781107807407369706366291665838939200322077170464298150478E-4),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 4.15072848381223183266747482990200902269625013596359343891434458183336715372385900595761261249062597314258758149849358E-6),
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
			BOOST_MATH_HUGE_CONSTANT(T, 0, 1.7606414208200893502974568628516195659203083721422991116474933229734072607662272187691874235194827274588738085452995E-1),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 5.2403354748695764514838391359483466233420304283870559387493001948377007156468753268890967128171616186758390430922461E-1),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 8.7004089535290290013495669628061004431655365634722323554340958728472997415189310305929038109585701566500242505994830E-1),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 1.2247448713915890490986420373529456959829737403283350642163462836254801887286575132699297165523201174092973006133071E+0),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 1.5794121348467670857203675829994749984448118700625756787296604004731493061239993938451136940463781480598323608108401E+0),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 1.8357079751751868737730366142776751003918259570839660905580513675260134800736001213972598631692284359287132511850212E+0),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 2.0232301911005156592083208951798835405134106140896567635236468204669731308178641289345019964860533407648486493612883E+0),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 2.2665132620567880274659861754385594433599181621460437501490609843716507610505839787938539918024628268231996046687956E+0),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 2.5705583765842967091130640044424751491766277443737424225502826093623108862325020470243771638339593980806475639753774E+0),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 2.9592107790638377223111385005349300110208537494323252529925306245788685818262550907459729341537440162617289883378384E+0),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 3.3491639537131949773673850267713958037641111667378569069641962722950738114287593396657637853425131249301284146570422E+0),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 3.6677742159463378600379325174576505145204038021378896733581590017259651263308722020865752907432590392345136026676227E+0),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 4.0292201405043713647935049781219652052808124117728748328853310891203520167360128069802223919620500375768781325797622E+0),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 4.4995993983103888028842951194003954000654644313136439364782096648105799974306887733838358325319930243399005003605842E+0),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 5.0360899444730939686859643226309469280847963344196632589865439531888353441510823110613688615646906148676400038430046E+0),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 5.6432578578857450628037542830378662310207930130774420026046822897984429696032941122779610378343375406655097758066947E+0),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 6.3759392709822359517127037507272178645929499888343595857714810863749752732749225626479585736280501578486632216476035E+0),
		};
		return data;
	}

	static std::array<T, 18> const& weights()
	{
		static std::array<T, 18> data = {
			BOOST_MATH_HUGE_CONSTANT(T, 0, 9.1262675363661801578461723905185563880880151345297350917022924108362037176546562106575739551835308880382756093967042E-4),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 3.3988595585585237653265161301450246144989577070183248283562067996550317571605351840788789091213092212304799674354603E-1),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 2.6244871488784306847516813136659885637306886986545338158216889756372556548464317088973152036952199495371147690804131E-1),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 1.6371221555735797873796678622548941320380961499304191088311319082257827181269242023999779442599696899624650543284521E-1),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 8.0245518147391120477490592421686090750751096122532003004330920653839151556433964198763992121228957472010560890278502E-2),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 2.7780508908535099862506195910853447479138574778251008228836025760828861191649881655495894801437404905149676420619177E-2),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 5.5928828911469410084910949056115943321252152809866088650264308216348689106756032206345612877650220006642019911120146E-3),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 4.0967527720344103442434518722522940850108270680968930386421899912467606797802657136322071739134577301912242524873433E-3),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 1.4515580425155862287583519513440792711406992365420180929996595806948567712159155284816842729655920752924204853437943E-3),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 4.8785399304443819054409619850884939685916337186481530335444823605211286910071447113615340411859811866251300787480405E-4),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 6.3328620805617479892669834075215315897168532489327656610957559718436447835796232961269824465663471569700029010604403E-5),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 4.8462799737020360011943494290738228597210020202131038001548330415639589742319388326154739208476350600166188582185136E-6),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 4.3737818040926880888283968138493143242076190374984571811171061484012340126605814886896170788677664308345348630131283E-7),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 3.7920222392319615175899951128843409596122655410526463735474181579001124925199122079120980099016503613437843438096693E-8),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 8.1553721816917387433569757116495481316245014435776720480956254093977330759735585606035997718150763108383313354157374E-10),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 5.4896836948499763670202498729464868453902121295324571971910475154692756211606054022054793821619812930891852536658986E-12),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 9.6599466278561074036789057225705231643273971844543914426618198701828842127606087910329390857148492091406620537753803536126925434E-15),
			BOOST_MATH_HUGE_CONSTANT(T, 0, 1.8684014894509412743803477189357925705011492609340575373109389180193866954500032513862031811341932560722403487306967E-18),
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
