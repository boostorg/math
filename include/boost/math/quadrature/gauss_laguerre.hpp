//  Copyright Jacob Hass 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_QUADRATURE_GAUSS_LAGUERRE_HPP
#define BOOST_MATH_QUADRATURE_GAUSS_LAGUERRE_HPP

#ifdef _MSC_VER
#pragma once
#endif

#include <boost/math/tools/config.hpp>

#ifdef BOOST_MATH_NO_CXX17_IF_CONSTEXPR
#error "The header <boost/math/quadrature/gauss_laguerre.hpp> requires C++17 or later."
#endif

#ifndef BOOST_MATH_BUILD_MODULE
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <type_traits>
#include <utility>
#include <vector>
#endif
#include <boost/math/constants/constants.hpp>
#include <boost/math/policies/error_handling.hpp>
#include <boost/math/policies/policy.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/math/tools/big_constant.hpp>
#include <boost/math/tools/precision.hpp>
#include <boost/math/quadrature/detail/quadrature_constant.hpp>
#include <boost/math/special_functions/detail/orthogonal_polynomial.hpp>
#include <boost/math/special_functions/laguerre.hpp>

BOOST_MATH_NAMESPACE_BEGIN namespace quadrature { namespace detail {

#ifndef BOOST_MATH_GAUSS_NO_COMPUTE_ON_DEMAND

template <class Real, unsigned N, unsigned Category>
class gauss_laguerre_detail
{
   static_assert(N > 1, "gauss_laguerre requires N > 1");
   static const boost::math::detail::orthogonal_polynomial<Real, boost::math::detail::laguerre_family<Real>>& polynomial()
   {
      static const boost::math::detail::orthogonal_polynomial<Real, boost::math::detail::laguerre_family<Real>> value(N);
      return value;
   }

public:
   static const std::vector<Real>& abscissa()
   {
      return polynomial().abscissa();
   }

   static const std::vector<Real>& weights()
   {
      return polynomial().weights();
   }
};

#else

template <class Real, unsigned N, unsigned Category>
class gauss_laguerre_detail;

#endif

#ifndef BOOST_HAS_FLOAT128
template <class T>
class gauss_laguerre_detail<T, 7, 0>
{
using storage_type = typename quadrature_constant_category<T>::storage_type;
public:
    static std::array<storage_type, 7> const & abscissa()
    {
        static std::array<storage_type, 7> data = {
        static_cast<storage_type>(1.93043676560362413838247885003823725e-01L),
        static_cast<storage_type>(1.02666489533919195034519944317483354e+00L),
        static_cast<storage_type>(2.56787674495074620690778622665959020e+00L),
        static_cast<storage_type>(4.90035308452648456810171437810242722e+00L),
        static_cast<storage_type>(8.18215344456286079108182755123396163e+00L),
        static_cast<storage_type>(1.27341802917978137580126424581951091e+01L),
        static_cast<storage_type>(1.93957278622625403117125820576302546e+01L),
        };
        return data;
    }
    static std::array<storage_type, 7> const & weights()
    {
        static std::array<storage_type, 7> data = {
        static_cast<storage_type>(4.09318951701273902130432880017778795e-01L),
        static_cast<storage_type>(4.21831277861719779929281005417498934e-01L),
        static_cast<storage_type>(1.47126348657505278395374184636545717e-01L),
        static_cast<storage_type>(2.06335144687169398657056149642012957e-02L),
        static_cast<storage_type>(1.07401014328074552213195962843031379e-03L),
        static_cast<storage_type>(1.58654643485642012687326223234055954e-05L),
        static_cast<storage_type>(3.17031547899558056227132215385324202e-08L),
        };
        return data;
    }
};

#else
template <class T>
class gauss_laguerre_detail<T, 7, 0>
{
using storage_type = typename quadrature_constant_category<T>::storage_type;
public:
    static std::array<storage_type, 7> const & abscissa()
    {
        static std::array<storage_type, 7> data = {
        static_cast<storage_type>(1.93043676560362413838247885003823725e-01Q),
        static_cast<storage_type>(1.02666489533919195034519944317483354e+00Q),
        static_cast<storage_type>(2.56787674495074620690778622665959020e+00Q),
        static_cast<storage_type>(4.90035308452648456810171437810242722e+00Q),
        static_cast<storage_type>(8.18215344456286079108182755123396163e+00Q),
        static_cast<storage_type>(1.27341802917978137580126424581951091e+01Q),
        static_cast<storage_type>(1.93957278622625403117125820576302546e+01Q),
        };
        return data;
    }
    static std::array<storage_type, 7> const & weights()
    {
        static std::array<storage_type, 7> data = {
        static_cast<storage_type>(4.09318951701273902130432880017778795e-01Q),
        static_cast<storage_type>(4.21831277861719779929281005417498934e-01Q),
        static_cast<storage_type>(1.47126348657505278395374184636545717e-01Q),
        static_cast<storage_type>(2.06335144687169398657056149642012957e-02Q),
        static_cast<storage_type>(1.07401014328074552213195962843031379e-03Q),
        static_cast<storage_type>(1.58654643485642012687326223234055954e-05Q),
        static_cast<storage_type>(3.17031547899558056227132215385324202e-08Q),
        };
        return data;
    }
};

#endif
template <class T>
class gauss_laguerre_detail<T, 7, 4>
{
using storage_type = typename quadrature_constant_category<T>::storage_type;
public:
    static std::array<T, 7> const & abscissa()
    {
        static std::array<T, 7> data = {
        BOOST_MATH_HUGE_CONSTANT(T, 0, 1.9304367656036241383824788500382372459120984891356750620225136201006415345782559871118003823365115612373429871679738e-01),
        BOOST_MATH_HUGE_CONSTANT(T, 0, 1.0266648953391919503451994431748335351021966391136589333607107349300045456686682315951892100561068930429806129371487e+00),
        BOOST_MATH_HUGE_CONSTANT(T, 0, 2.5678767449507462069077862266595901968369581616444700992866063303240503321974487821383346200011905724664254743082485e+00),
        BOOST_MATH_HUGE_CONSTANT(T, 0, 4.9003530845264845681017143781024272182326822536691323824442733012926988370938012768358608940248312092759732409230278e+00),
        BOOST_MATH_HUGE_CONSTANT(T, 0, 8.1821534445628607910818275512339616296407668517524296643048342450041449025202649459386747106022315685557229645415476e+00),
        BOOST_MATH_HUGE_CONSTANT(T, 0, 1.2734180291797813758012642458195109090853336141446261917460200935582863612586646759328152143354140682014199320963541e+01),
        BOOST_MATH_HUGE_CONSTANT(T, 0, 1.9395727862262540311712582057630254604742850103460479496941123090856173616475344405452608383727847918520964087609689e+01),
        };
        return data;
    }
    static std::array<T, 7> const & weights()
    {
        static std::array<T, 7> data = {
        BOOST_MATH_HUGE_CONSTANT(T, 0, 4.0931895170127390213043288001777879476694318320743384936496071108439327612531393650526756975329973583687149037045411e-01),
        BOOST_MATH_HUGE_CONSTANT(T, 0, 4.2183127786171977992928100541749893425900270593200299919010224347881378577765054361914173904057093683899914455920747e-01),
        BOOST_MATH_HUGE_CONSTANT(T, 0, 1.4712634865750527839537418463654571735447659494559067156700082178865140336466631595376166029179058453234935126953389e-01),
        BOOST_MATH_HUGE_CONSTANT(T, 0, 2.0633514468716939865705614964201295701345516595854480283021401106613334093662834206930986586836486163163369652374680e-02),
        BOOST_MATH_HUGE_CONSTANT(T, 0, 1.0740101432807455221319596284303137904257521979552018605203061608685037526387650167316790636375375082365847554380883e-03),
        BOOST_MATH_HUGE_CONSTANT(T, 0, 1.5865464348564201268732622323405595386094807050091717842696789736074648974951001982867672006063759990800667569423195e-05),
        BOOST_MATH_HUGE_CONSTANT(T, 0, 3.1703154789955805622713221538532420152314112706016551819590923622237092653696183497591858655360389258725422336556501e-08),
        };
        return data;
    }
};

#ifndef BOOST_HAS_FLOAT128
template <class T>
class gauss_laguerre_detail<T, 10, 0>
{
using storage_type = typename quadrature_constant_category<T>::storage_type;
public:
    static std::array<storage_type, 10> const & abscissa()
    {
        static std::array<storage_type, 10> data = {
        static_cast<storage_type>(1.37793470540492430830772505652711188e-01L),
        static_cast<storage_type>(7.29454549503170498160373121676078781e-01L),
        static_cast<storage_type>(1.80834290174031604823292007575060883e+00L),
        static_cast<storage_type>(3.40143369785489951448253222140839068e+00L),
        static_cast<storage_type>(5.55249614006380363241755848686876286e+00L),
        static_cast<storage_type>(8.33015274676449670023876719727452218e+00L),
        static_cast<storage_type>(1.18437858379000655649185389191416140e+01L),
        static_cast<storage_type>(1.62792578313781020995326539358336223e+01L),
        static_cast<storage_type>(2.19965858119807619512770901955944940e+01L),
        static_cast<storage_type>(2.99206970122738915599087933407991952e+01L),
        };
        return data;
    }
    static std::array<storage_type, 10> const & weights()
    {
        static std::array<storage_type, 10> data = {
        static_cast<storage_type>(3.08441115765020141547470834677860696e-01L),
        static_cast<storage_type>(4.01119929155273551515780309912819515e-01L),
        static_cast<storage_type>(2.18068287611809421588648523474646727e-01L),
        static_cast<storage_type>(6.20874560986777473929021293135179537e-02L),
        static_cast<storage_type>(9.50151697518110055383907219417199123e-03L),
        static_cast<storage_type>(7.53008388587538775455964353675663902e-04L),
        static_cast<storage_type>(2.82592334959956556742256382685002128e-05L),
        static_cast<storage_type>(4.24931398496268637258657665974712355e-07L),
        static_cast<storage_type>(1.83956482397963078092153522435593825e-09L),
        static_cast<storage_type>(9.91182721960900855837754728324473606e-13L),
        };
        return data;
    }
};

#else
template <class T>
class gauss_laguerre_detail<T, 10, 0>
{
using storage_type = typename quadrature_constant_category<T>::storage_type;
public:
    static std::array<storage_type, 10> const & abscissa()
    {
        static std::array<storage_type, 10> data = {
        static_cast<storage_type>(1.37793470540492430830772505652711188e-01Q),
        static_cast<storage_type>(7.29454549503170498160373121676078781e-01Q),
        static_cast<storage_type>(1.80834290174031604823292007575060883e+00Q),
        static_cast<storage_type>(3.40143369785489951448253222140839068e+00Q),
        static_cast<storage_type>(5.55249614006380363241755848686876286e+00Q),
        static_cast<storage_type>(8.33015274676449670023876719727452218e+00Q),
        static_cast<storage_type>(1.18437858379000655649185389191416140e+01Q),
        static_cast<storage_type>(1.62792578313781020995326539358336223e+01Q),
        static_cast<storage_type>(2.19965858119807619512770901955944940e+01Q),
        static_cast<storage_type>(2.99206970122738915599087933407991952e+01Q),
        };
        return data;
    }
    static std::array<storage_type, 10> const & weights()
    {
        static std::array<storage_type, 10> data = {
        static_cast<storage_type>(3.08441115765020141547470834677860696e-01Q),
        static_cast<storage_type>(4.01119929155273551515780309912819515e-01Q),
        static_cast<storage_type>(2.18068287611809421588648523474646727e-01Q),
        static_cast<storage_type>(6.20874560986777473929021293135179537e-02Q),
        static_cast<storage_type>(9.50151697518110055383907219417199123e-03Q),
        static_cast<storage_type>(7.53008388587538775455964353675663902e-04Q),
        static_cast<storage_type>(2.82592334959956556742256382685002128e-05Q),
        static_cast<storage_type>(4.24931398496268637258657665974712355e-07Q),
        static_cast<storage_type>(1.83956482397963078092153522435593825e-09Q),
        static_cast<storage_type>(9.91182721960900855837754728324473606e-13Q),
        };
        return data;
    }
};

#endif
template <class T>
class gauss_laguerre_detail<T, 10, 4>
{
using storage_type = typename quadrature_constant_category<T>::storage_type;
public:
    static std::array<T, 10> const & abscissa()
    {
        static std::array<T, 10> data = {
        BOOST_MATH_HUGE_CONSTANT(T, 0, 1.3779347054049243083077250565271118810799168074578329416336651137344596447646208656254375241734116529355248965683739e-01),
        BOOST_MATH_HUGE_CONSTANT(T, 0, 7.2945454950317049816037312167607878107607273331224924600781104880626782307305157204652499178629631939583813469724208e-01),
        BOOST_MATH_HUGE_CONSTANT(T, 0, 1.8083429017403160482329200757506088332830602823714481516932097683732967664043428239153383606124012239270279810042947e+00),
        BOOST_MATH_HUGE_CONSTANT(T, 0, 3.4014336978548995144825322214083906792731566142034732294267661021823269228666946410395095870847780590663426380062014e+00),
        BOOST_MATH_HUGE_CONSTANT(T, 0, 5.5524961400638036324175584868687628579740642873178138284904788336252636084290435460571874259935064862774598817665864e+00),
        BOOST_MATH_HUGE_CONSTANT(T, 0, 8.3301527467644967002387671972745221827094389720309893447446294283573954320532476951665824491858289576075525949818506e+00),
        BOOST_MATH_HUGE_CONSTANT(T, 0, 1.1843785837900065564918538919141613985802816909465335873339591495276220150697108232717212527603241705173109616546015e+01),
        BOOST_MATH_HUGE_CONSTANT(T, 0, 1.6279257831378102099532653935833622335255995603033059077788392165002827350932597940225841507596519165747859331720868e+01),
        BOOST_MATH_HUGE_CONSTANT(T, 0, 2.1996585811980761951277090195594493976806732340001887788241503511693106375798998823067686533346980361167369020416848e+01),
        BOOST_MATH_HUGE_CONSTANT(T, 0, 2.9920697012273891559908793340799195179710670577517960166104251135309849605268452639201572864373106556343888311203255e+01),
        };
        return data;
    }
    static std::array<T, 10> const & weights()
    {
        static std::array<T, 10> data = {
        BOOST_MATH_HUGE_CONSTANT(T, 0, 3.0844111576502014154747083467786069562872888653833744211555718020018556202932557003040685460579742269270912496733327e-01),
        BOOST_MATH_HUGE_CONSTANT(T, 0, 4.0111992915527355151578030991281951479548361696211301757496775516766091500023676907078271471115370136487330141071994e-01),
        BOOST_MATH_HUGE_CONSTANT(T, 0, 2.1806828761180942158864852347464672674277853841218894056650398132078664598285192349000126323603209982063165104569642e-01),
        BOOST_MATH_HUGE_CONSTANT(T, 0, 6.2087456098677747392902129313517953695909065683802092368676760303923866462990276943718400176674956216762717600777500e-02),
        BOOST_MATH_HUGE_CONSTANT(T, 0, 9.5015169751811005538390721941719912258624504015797533646318970848680322451340750784240877901570876971156524042858250e-03),
        BOOST_MATH_HUGE_CONSTANT(T, 0, 7.5300838858753877545596435367566390179203914014362886506164014367761532732553587715637849496678837150424361526704126e-04),
        BOOST_MATH_HUGE_CONSTANT(T, 0, 2.8259233495995655674225638268500212828033164744374677117086795259322508664275187881955981278223738401217175272290547e-05),
        BOOST_MATH_HUGE_CONSTANT(T, 0, 4.2493139849626863725865766597471235464810801986441570081254604829021791344429695503912683040809284918604347723905036e-07),
        BOOST_MATH_HUGE_CONSTANT(T, 0, 1.8395648239796307809215352243559382479826127765906584095132767065586955845890346516805867002792829341205430551248512e-09),
        BOOST_MATH_HUGE_CONSTANT(T, 0, 9.9118272196090085583775472832447360645810946112447692431304366686254543563865419652261172586997161662742217301111110e-13),
        };
        return data;
    }
};

#ifndef BOOST_HAS_FLOAT128
template <class T>
class gauss_laguerre_detail<T, 15, 0>
{
using storage_type = typename quadrature_constant_category<T>::storage_type;
public:
    static std::array<storage_type, 15> const & abscissa()
    {
        static std::array<storage_type, 15> data = {
        static_cast<storage_type>(9.33078120172818047629030383672077003e-02L),
        static_cast<storage_type>(4.92691740301883908960101791412401731e-01L),
        static_cast<storage_type>(1.21559541207094946372992716487900297e+00L),
        static_cast<storage_type>(2.26994952620374320247421741374817219e+00L),
        static_cast<storage_type>(3.66762272175143727724905959436028773e+00L),
        static_cast<storage_type>(5.42533662741355316534358132596053452e+00L),
        static_cast<storage_type>(7.56591622661306786049739555811855455e+00L),
        static_cast<storage_type>(1.01202285680191127347927394568171697e+01L),
        static_cast<storage_type>(1.31302824821757235640991204175965313e+01L),
        static_cast<storage_type>(1.66544077083299578225202408429835997e+01L),
        static_cast<storage_type>(2.07764788994487667729157175675934362e+01L),
        static_cast<storage_type>(2.56238942267287801445868285976698414e+01L),
        static_cast<storage_type>(3.14075191697539385152432196202214042e+01L),
        static_cast<storage_type>(3.85306833064860094162515167595161938e+01L),
        static_cast<storage_type>(4.80260855726857943465734308507556623e+01L),
        };
        return data;
    }
    static std::array<storage_type, 15> const & weights()
    {
        static std::array<storage_type, 15> data = {
        static_cast<storage_type>(2.18234885940086889856413236448110800e-01L),
        static_cast<storage_type>(3.42210177922883329638948956806686032e-01L),
        static_cast<storage_type>(2.63027577941680097414812275021951666e-01L),
        static_cast<storage_type>(1.26425818105930535843030549378389129e-01L),
        static_cast<storage_type>(4.02068649210009148415854789871163157e-02L),
        static_cast<storage_type>(8.56387780361183836391575987649400017e-03L),
        static_cast<storage_type>(1.21243614721425207621920522466788622e-03L),
        static_cast<storage_type>(1.11674392344251941992578595518123419e-04L),
        static_cast<storage_type>(6.45992676202290092465319025312414635e-06L),
        static_cast<storage_type>(2.22631690709627263033182809179377262e-07L),
        static_cast<storage_type>(4.22743038497936500735127949331007487e-09L),
        static_cast<storage_type>(3.92189726704108929038460981948602319e-11L),
        static_cast<storage_type>(1.45651526407312640633273963455309983e-13L),
        static_cast<storage_type>(1.48302705111330133546164737187127769e-16L),
        static_cast<storage_type>(1.60059490621113323104997812369887952e-20L),
        };
        return data;
    }
};

#else
template <class T>
class gauss_laguerre_detail<T, 15, 0>
{
using storage_type = typename quadrature_constant_category<T>::storage_type;
public:
    static std::array<storage_type, 15> const & abscissa()
    {
        static std::array<storage_type, 15> data = {
        static_cast<storage_type>(9.33078120172818047629030383672077003e-02Q),
        static_cast<storage_type>(4.92691740301883908960101791412401731e-01Q),
        static_cast<storage_type>(1.21559541207094946372992716487900297e+00Q),
        static_cast<storage_type>(2.26994952620374320247421741374817219e+00Q),
        static_cast<storage_type>(3.66762272175143727724905959436028773e+00Q),
        static_cast<storage_type>(5.42533662741355316534358132596053452e+00Q),
        static_cast<storage_type>(7.56591622661306786049739555811855455e+00Q),
        static_cast<storage_type>(1.01202285680191127347927394568171697e+01Q),
        static_cast<storage_type>(1.31302824821757235640991204175965313e+01Q),
        static_cast<storage_type>(1.66544077083299578225202408429835997e+01Q),
        static_cast<storage_type>(2.07764788994487667729157175675934362e+01Q),
        static_cast<storage_type>(2.56238942267287801445868285976698414e+01Q),
        static_cast<storage_type>(3.14075191697539385152432196202214042e+01Q),
        static_cast<storage_type>(3.85306833064860094162515167595161938e+01Q),
        static_cast<storage_type>(4.80260855726857943465734308507556623e+01Q),
        };
        return data;
    }
    static std::array<storage_type, 15> const & weights()
    {
        static std::array<storage_type, 15> data = {
        static_cast<storage_type>(2.18234885940086889856413236448110800e-01Q),
        static_cast<storage_type>(3.42210177922883329638948956806686032e-01Q),
        static_cast<storage_type>(2.63027577941680097414812275021951666e-01Q),
        static_cast<storage_type>(1.26425818105930535843030549378389129e-01Q),
        static_cast<storage_type>(4.02068649210009148415854789871163157e-02Q),
        static_cast<storage_type>(8.56387780361183836391575987649400017e-03Q),
        static_cast<storage_type>(1.21243614721425207621920522466788622e-03Q),
        static_cast<storage_type>(1.11674392344251941992578595518123419e-04Q),
        static_cast<storage_type>(6.45992676202290092465319025312414635e-06Q),
        static_cast<storage_type>(2.22631690709627263033182809179377262e-07Q),
        static_cast<storage_type>(4.22743038497936500735127949331007487e-09Q),
        static_cast<storage_type>(3.92189726704108929038460981948602319e-11Q),
        static_cast<storage_type>(1.45651526407312640633273963455309983e-13Q),
        static_cast<storage_type>(1.48302705111330133546164737187127769e-16Q),
        static_cast<storage_type>(1.60059490621113323104997812369887952e-20Q),
        };
        return data;
    }
};

#endif
template <class T>
class gauss_laguerre_detail<T, 15, 4>
{
using storage_type = typename quadrature_constant_category<T>::storage_type;
public:
    static std::array<T, 15> const & abscissa()
    {
        static std::array<T, 15> data = {
        BOOST_MATH_HUGE_CONSTANT(T, 0, 9.3307812017281804762903038367207700278197076330453493009360137058347064015776121761668303610899058588907097107699871e-02),
        BOOST_MATH_HUGE_CONSTANT(T, 0, 4.9269174030188390896010179141240173099215246375371941601988508193506134664433011461933813698167861160364523544308066e-01),
        BOOST_MATH_HUGE_CONSTANT(T, 0, 1.2155954120709494637299271648790029749872681589377057427014003173216015106603771651065605856143208373528444729551797e+00),
        BOOST_MATH_HUGE_CONSTANT(T, 0, 2.2699495262037432024742174137481721885680048416478717117778472877455733945904832989178549502018189830165872698672666e+00),
        BOOST_MATH_HUGE_CONSTANT(T, 0, 3.6676227217514372772490595943602877290357430909410717845566428732064149379102861748582946304744040432223255333590429e+00),
        BOOST_MATH_HUGE_CONSTANT(T, 0, 5.4253366274135531653435813259605345208492623922927783354709435107483423740280199191227091479154617337151884199685196e+00),
        BOOST_MATH_HUGE_CONSTANT(T, 0, 7.5659162266130678604973955581185545482393108008999314491628694391139408365724029944489180550403629640393078390436542e+00),
        BOOST_MATH_HUGE_CONSTANT(T, 0, 1.0120228568019112734792739456817169743410911967058744917457454601072487825047671876216088650721887243391350486173212e+01),
        BOOST_MATH_HUGE_CONSTANT(T, 0, 1.3130282482175723564099120417596531310187270022548453696047974546863746229290362972280038996276387001030969282397023e+01),
        BOOST_MATH_HUGE_CONSTANT(T, 0, 1.6654407708329957822520240842983599653951575341030608785764448525487749797539465693536916048016362945551853788626249e+01),
        BOOST_MATH_HUGE_CONSTANT(T, 0, 2.0776478899448766772915717567593436247705141256327195360002221322945959767932746815884396953509961946757428840574146e+01),
        BOOST_MATH_HUGE_CONSTANT(T, 0, 2.5623894226728780144586828597669841397123940322411476950317672315339824238931352351127917412081280502785586672176735e+01),
        BOOST_MATH_HUGE_CONSTANT(T, 0, 3.1407519169753938515243219620221404196345622809889669925953438870665440073910962255366662433543120684101302171932893e+01),
        BOOST_MATH_HUGE_CONSTANT(T, 0, 3.8530683306486009416251516759516193798932119744594784284050364282888975230612867789684690042672596606317674205933510e+01),
        BOOST_MATH_HUGE_CONSTANT(T, 0, 4.8026085572685794346573430850755662259393479711335534147707476887606535372312894457067945653339456838525028684441788e+01),
        };
        return data;
    }
    static std::array<T, 15> const & weights()
    {
        static std::array<T, 15> data = {
        BOOST_MATH_HUGE_CONSTANT(T, 0, 2.1823488594008688985641323644811080021343324068445151647893418259657190410748792541872504691025036395044680815196202e-01),
        BOOST_MATH_HUGE_CONSTANT(T, 0, 3.4221017792288332963894895680668603182827472147567140620457276907638708199822947378530621505254574834996032558277546e-01),
        BOOST_MATH_HUGE_CONSTANT(T, 0, 2.6302757794168009741481227502195166616147477737942289338398907943112508686382961577835440203774163021418786314618248e-01),
        BOOST_MATH_HUGE_CONSTANT(T, 0, 1.2642581810593053584303054937838912877184803937971685684567321902219954509817127653766637180913751240893191067737467e-01),
        BOOST_MATH_HUGE_CONSTANT(T, 0, 4.0206864921000914841585478987116315744354255151588087542482770026083846642284268861374498795123424591421789047455108e-02),
        BOOST_MATH_HUGE_CONSTANT(T, 0, 8.5638778036118383639157598764940001725843892261220703970182442186170667853270532099641265249737915428749887172527299e-03),
        BOOST_MATH_HUGE_CONSTANT(T, 0, 1.2124361472142520762192052246678862218825774912216089252896089351714176877531295775371012761187312001078014473388155e-03),
        BOOST_MATH_HUGE_CONSTANT(T, 0, 1.1167439234425194199257859551812341890362625693267480181637850779341432852535113805714378050258665226088107162969510e-04),
        BOOST_MATH_HUGE_CONSTANT(T, 0, 6.4599267620229009246531902531241463548444527335439555949894778198068019741285302492718083846291056795408895968042879e-06),
        BOOST_MATH_HUGE_CONSTANT(T, 0, 2.2263169070962726303318280917937726215071887584934068919788512991159987070128917425821426819036682039019642914030592e-07),
        BOOST_MATH_HUGE_CONSTANT(T, 0, 4.2274303849793650073512794933100748675498159876301045574722416288222464623638254108884085708408168629608639640555277e-09),
        BOOST_MATH_HUGE_CONSTANT(T, 0, 3.9218972670410892903846098194860231885836676205596680875788417907215759905742825209015082901147543134852484805575028e-11),
        BOOST_MATH_HUGE_CONSTANT(T, 0, 1.4565152640731264063327396345530998285311241361150773694115515272365278527396127670028396473753138295863865512588769e-13),
        BOOST_MATH_HUGE_CONSTANT(T, 0, 1.4830270511133013354616473718712776892067670532770945176556036443926719007518190902345311172254087343813449801330193e-16),
        BOOST_MATH_HUGE_CONSTANT(T, 0, 1.6005949062111332310499781236988795197921487028677950643996761551559291234920430457494041052454552370257733201426544e-20),
        };
        return data;
    }
};

} // namespace detail

template <class Real, unsigned N, class Policy = BOOST_MATH_NAMESPACE::policies::policy<> >
class gauss_laguerre : public detail::gauss_laguerre_detail<Real, N, detail::quadrature_constant_category<Real>::value>
{
   using base = detail::gauss_laguerre_detail<Real, N, detail::quadrature_constant_category<Real>::value>;

   // Abscissas computed on demand are NaN when the recurrence overflows Real.
   static bool overflowed()
   {
      return !(BOOST_MATH_NAMESPACE::isfinite)(static_cast<Real>(base::abscissa().back()));
   }

   static Real overflow_error()
   {
      return policies::raise_evaluation_error(
         "boost::math::quadrature::gauss_laguerre<%1%>::integrate",
         "The Laguerre recurrence overflows with this many points; use fewer points or a type with a wider exponent range.",
         std::numeric_limits<Real>::quiet_NaN(), Policy());
   }

public:
   // Integrates f(x) exp(-x) over the real line.
   template <class F>
   static auto integrate(F f, Real* pL1 = nullptr)->decltype(std::declval<F>()(std::declval<Real>()))
   {
      // In many math texts, K represents the field of real or complex numbers.
      // Too bad we can't put blackboard bold into C++ source!
      using K = decltype(f(Real(0)));
      static_assert(!std::is_integral<K>::value,
                    "The return type cannot be integral, it must be either a real or complex floating point type.");
      using std::abs;
      if (overflowed())
      {
         Real error = overflow_error();
         if (pL1)
            *pL1 = error;
         return static_cast<K>(error);
      }
      K result = Real(0);
      Real L1 = abs(result);
      for (unsigned i = 0; i < base::abscissa().size(); ++i)
      {
         K fp = f(static_cast<Real>(base::abscissa()[i]));
         result += fp * static_cast<Real>(base::weights()[i]);
         L1 += abs(fp) * static_cast<Real>(base::weights()[i]);
      }
      if (pL1)
         *pL1 = L1;
      return result;
   }

   // Explicit zero and norm for vector- and matrix-valued integrands.
   template <class F, class Norm>
   static auto integrate(F f, const decltype(std::declval<F>()(std::declval<Real>()))& zero, Norm norm, Real* pL1 = nullptr)
      ->decltype(static_cast<Real>(norm(f(Real(0)))), f(Real(0)))
   {
      using K = decltype(f(Real(0)));
      static_assert(!std::is_integral<K>::value, "The return type cannot be integral.");
      if (overflowed())
      {
         // Preserve the scalar error policy, including non-throwing policies,
         // and use scalar multiplication to produce a result with the right shape.
         Real error = overflow_error();
         if (pL1)
            *pL1 = error;
         K result = zero * error;
         return result;
      }
      K result = zero;
      Real L1 = static_cast<Real>(norm(result));
      for (unsigned i = 0; i < base::abscissa().size(); ++i)
      {
         K fp = f(static_cast<Real>(base::abscissa()[i]));
         result += fp * static_cast<Real>(base::weights()[i]);
         L1 += static_cast<Real>(norm(fp)) * static_cast<Real>(base::weights()[i]);
      }
      if (pL1)
         *pL1 = L1;
      return result;
   }
};
} // namespace quadrature

BOOST_MATH_NAMESPACE_END

#endif // BOOST_MATH_QUADRATURE_GAUSS_LAGUERRE_HPP
