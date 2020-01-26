/*
 * Copyright Nick Thompson, 2019
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef BOOST_MATH_FILTERS_DAUBECHIES_HPP
#define BOOST_MATH_FILTERS_DAUBECHIES_HPP
#include <array>
#ifdef BOOST_HAS_FLOAT128
#include <boost/multiprecision/float128.hpp>
#endif
namespace boost::math::filters {

template <typename Real, unsigned p>
constexpr std::array<Real, 2*p> daubechies_scaling_filter()
{
    static_assert(sizeof(Real) <= 16, "Filter coefficients only computed up to 128 bits of precision.");
    static_assert(p < 25, "Filter coefficients only implemented up to 24.");
    if constexpr (p == 1) {
        if constexpr (std::is_same_v<Real, float>) {
            return {0x1.6a09e6p-1f, 0x1.6a09e6p-1f};
        }
        if constexpr (std::is_same_v<Real, double>) {
            return {0x1.6a09e667f3bcdp-1, 0x1.6a09e667f3bcdp-1};
        }
        if constexpr (std::is_same_v<Real, long double>) {
            return {0xb.504f333f9de6484p-4L, 0xb.504f333f9de6484p-4L};
        }
        #ifdef BOOST_HAS_FLOAT128
        if constexpr (std::is_same_v<Real, boost::multiprecision::float128>) {
            return {0x1.6a09e667f3bcc908b2fb1366ea95p-1Q, 0x1.6a09e667f3bcc908b2fb1366ea95p-1Q};
        }
        #endif
    }
    if constexpr (p == 2) {
        if constexpr (std::is_same_v<Real, float>) {
            return {0x1.ee8dd4p-2f, 0x1.ac4bdep-1f, 0x1.cb0bfp-3f, -0x1.0907dcp-3f};
        }
        if constexpr (std::is_same_v<Real, double>) {
            return {0x1.ee8dd4748bf15p-2, 0x1.ac4bdd6e3fd71p-1, 0x1.cb0bf0b6b7109p-3, -0x1.0907dc193069p-3};
        }
        if constexpr (std::is_same_v<Real, long double>) {
            return {0xf.746ea3a45f8a62ap-5L, 0xd.625eeb71feb8557p-4L, 0xe.585f85b5b8845bdp-6L, -0x8.483ee0c9834834cp-6L};
        }
        #ifdef BOOST_HAS_FLOAT128
        if constexpr (std::is_same_v<Real, boost::multiprecision::float128>) {
            return {0x1.ee8dd4748bf14c548b969de58fap-2Q, 0x1.ac4bdd6e3fd70aae9f48d8a63d1bp-1Q, 0x1.cb0bf0b6b7108b79b4bf11d08b16p-3Q, -0x1.0907dc1930690697b13714fd4a15p-3Q};
        }
        #endif
    }
    if constexpr (p == 3) {
        if constexpr (std::is_same_v<Real, float>) {
            return {0x1.54a796p-2f, 0x1.9d20e2p-1f, 0x1.d6ea2p-2f, -0x1.1480a8p-3f, -0x1.5df7acp-4f, 0x1.2092e4p-5f};
        }
        if constexpr (std::is_same_v<Real, double>) {
            return {0x1.54a796e50d264p-2, 0x1.9d20e247d28bbp-1, 0x1.d6ea20bf0f744p-2, -0x1.1480a85c59629p-3, -0x1.5df7ab50d483cp-4, 0x1.2092e373789b9p-5};
        }
        if constexpr (std::is_same_v<Real, long double>) {
            return {0xa.a53cb7286931db4p-5L, 0xc.e907123e945da19p-4L, 0xe.b75105f87ba23d6p-5L, -0x8.a40542e2cb14943p-6L, -0xa.efbd5a86a41e206p-7L, 0x9.04971b9bc4dcbc1p-8L};
        }
        #ifdef BOOST_HAS_FLOAT128
        if constexpr (std::is_same_v<Real, boost::multiprecision::float128>) {
            return {0x1.54a796e50d263b67e38e82a5584cp-2Q, 0x1.9d20e247d28bb431f1e6c634b34ep-1Q, 0x1.d6ea20bf0f7447ac92de97f0e152p-2Q, -0x1.1480a85c596292857548d060a171p-3Q, -0x1.5df7ab50d483c40c41dbcf2191cfp-4Q, 0x1.2092e373789b9781e66814a5fa37p-5Q};
        }
        #endif
    }
    if constexpr (p == 4) {
        if constexpr (std::is_same_v<Real, float>) {
            return {0x1.d7d052p-3f, 0x1.6e005ep-1f, 0x1.4302cep-1f, -0x1.ca7c7p-6f, -0x1.7f0c1cp-3f, 0x1.f94e22p-6f, 0x1.0d60acp-5f, -0x1.5b4174p-7f};
        }
        if constexpr (std::is_same_v<Real, double>) {
            return {0x1.d7d052af15ecp-3, 0x1.6e005ea45d748p-1, 0x1.4302cdd3de43ap-1, -0x1.ca7c6f9db5bfbp-6, -0x1.7f0c1b7c604d4p-3, 0x1.f94e2196383a9p-6, 0x1.0d60ac768117bp-5, -0x1.5b41730b72e29p-7};
        }
        if constexpr (std::is_same_v<Real, long double>) {
            return {0xe.be829578af5fe6dp-6L, 0xb.7002f522eba3e31p-4L, 0xa.18166e9ef21cf0dp-4L, -0xe.53e37cedadfd7d8p-9L, -0xb.f860dbe30269f7ap-6L, 0xf.ca710cb1c1d472dp-9L, 0x8.6b0563b408bdbabp-8L, -0xa.da0b985b97149dcp-10L};
        }
        #ifdef BOOST_HAS_FLOAT128
        if constexpr (std::is_same_v<Real, boost::multiprecision::float128>) {
            return {0x1.d7d052af15ebfcd98dc85ef4cc01p-3Q, 0x1.6e005ea45d747c624e43405b1919p-1Q, 0x1.4302cdd3de439e1a299a6a04b89fp-1Q, -0x1.ca7c6f9db5bfafb07b9dace22f7p-6Q, -0x1.7f0c1b7c604d3ef4ab50462dfcedp-3Q, 0x1.f94e2196383a8e592495a4baff62p-6Q, 0x1.0d60ac768117b7550c2a3307e30fp-5Q, -0x1.5b41730b72e293b823fb2cbd40c6p-7Q};
        }
        #endif
    }
    if constexpr (p == 5) {
        if constexpr (std::is_same_v<Real, float>) {
            return {0x1.47e3c4p-3f, 0x1.35291cp-1f, 0x1.72d892p-1f, 0x1.1b8038p-3f, -0x1.f0384ep-3f, -0x1.082664p-5f, 0x1.3dbb9cp-4f, -0x1.990ad4p-8f, -0x1.9c3fp-7f, 0x1.b5385ep-9f};
        }
        if constexpr (std::is_same_v<Real, double>) {
            return {0x1.47e3c41a7b911p-3, 0x1.35291c2c4b00cp-1, 0x1.72d89143b54f5p-1, 0x1.1b80373befcc6p-3, -0x1.f0384d3f81474p-3, -0x1.0826648a8dc74p-5, 0x1.3dbb9b52515aap-4, -0x1.990ad4579f2e8p-8, -0x1.9c3eff3294128p-7, 0x1.b5385e04e3c09p-9};
        }
        if constexpr (std::is_same_v<Real, long double>) {
            return {0xa.3f1e20d3dc8840dp-6L, 0x9.a948e1625806185p-4L, 0xb.96c48a1daa7a767p-4L, 0x8.dc01b9df7e62cc6p-6L, -0xf.81c269fc0a3a1e2p-6L, -0x8.413324546e39fcbp-8L, 0x9.eddcda928ad4cf6p-7L, -0xc.c856a2bcf973f17p-11L, -0xc.e1f7f994a094303p-10L, 0xd.a9c2f0271e04896p-12L};
        }
        #ifdef BOOST_HAS_FLOAT128
        if constexpr (std::is_same_v<Real, boost::multiprecision::float128>) {
            return {0x1.47e3c41a7b91081ae95ed12670b7p-3Q, 0x1.35291c2c4b00c30a87eaf4d05e71p-1Q, 0x1.72d89143b54f4ecd724f0b79d111p-1Q, 0x1.1b80373befcc598b25a1d5436079p-3Q, -0x1.f0384d3f814743c3d93da5cd8909p-3Q, -0x1.0826648a8dc73f96ef19214f403ep-5Q, 0x1.3dbb9b52515a99ecaa66a225e968p-4Q, -0x1.990ad4579f2e7e2dc54832641bd8p-8Q, -0x1.9c3eff32941286062a45cb776526p-7Q, 0x1.b5385e04e3c0912c23cbbf7041e8p-9Q};
        }
        #endif
    }
    if constexpr (p == 6) {
        if constexpr (std::is_same_v<Real, float>) {
            return {0x1.c8def2p-4f, 0x1.fa7ebp-2f, 0x1.8094ap-1f, 0x1.42d0fcp-2f, -0x1.cf63dep-3f, -0x1.09c336p-3f, 0x1.8f5dd8p-4f, 0x1.c2ef44p-6f, -0x1.02b856p-5f, 0x1.225f72p-11f, 0x1.391514p-8f, -0x1.1a6874p-10f};
        }
        if constexpr (std::is_same_v<Real, double>) {
            return {0x1.c8def24dc3952p-4, 0x1.fa7eaf64539a9p-2, 0x1.80949fa3bc0bbp-1, 0x1.42d0fcfa92f21p-2, -0x1.cf63dd26916f1p-3, -0x1.09c33622722ebp-3, 0x1.8f5dd7f4e1752p-4, 0x1.c2ef43d612549p-6, -0x1.02b856404e8cep-5, 0x1.225f71210a7c1p-11, 0x1.391514c62a31bp-8, -0x1.1a6873b7a6466p-10};
        }
        if constexpr (std::is_same_v<Real, long double>) {
            return {0xe.46f7926e1ca8e15p-7L, 0xf.d3f57b229cd4bdfp-5L, 0xc.04a4fd1de05d701p-4L, 0xa.1687e7d497907f9p-5L, -0xe.7b1ee9348b785acp-6L, -0x8.4e19b11391756c5p-6L, 0xc.7aeebfa70ba93cep-7L, 0xe.177a1eb092a4be4p-9L, -0x8.15c2b2027466fc3p-8L, 0x9.12fb890853e0901p-14L, 0x9.c8a8a631518d6e7p-11L, -0x8.d3439dbd3232f08p-13L};
        }
        #ifdef BOOST_HAS_FLOAT128
        if constexpr (std::is_same_v<Real, boost::multiprecision::float128>) {
            return {0x1.c8def24dc3951c2a6dae413b8ddbp-4Q, 0x1.fa7eaf64539a97bd50da008738e1p-2Q, 0x1.80949fa3bc0bae02bf1b494d3af8p-1Q, 0x1.42d0fcfa92f20ff1b5b87f2f698cp-2Q, -0x1.cf63dd26916f0b589f7743608e82p-3Q, -0x1.09c33622722ead8a0ff038a97bp-3Q, 0x1.8f5dd7f4e175279c5a356a3fbcc9p-4Q, 0x1.c2ef43d6125497c75e24a80c6bfap-6Q, -0x1.02b856404e8cdf85933c3736d91cp-5Q, 0x1.225f71210a7c1202de723c931ed6p-11Q, 0x1.391514c62a31adcdfa747aebbdeep-8Q, -0x1.1a6873b7a6465e0ff5c033201061p-10Q};
        }
        #endif
    }
    if constexpr (p == 7) {
        if constexpr (std::is_same_v<Real, float>) {
            return {0x1.3ee1ccp-4f, 0x1.960e68p-2f, 0x1.7550cep-1f, 0x1.e10e9cp-2f, -0x1.26b83p-3f, -0x1.cad37cp-3f, 0x1.241522p-4f, 0x1.4a3072p-4f, -0x1.378a8ep-5f, -0x1.0f8eaap-6f, 0x1.9b4568p-7f, 0x1.c271f6p-12f, -0x1.d84a1p-10f, 0x1.72e554p-12f};
        }
        if constexpr (std::is_same_v<Real, double>) {
            return {0x1.3ee1cba38b6b1p-4, 0x1.960e674303003p-2, 0x1.7550cd294c1fep-1, 0x1.e10e9ba294ddcp-2, -0x1.26b830e491e33p-3, -0x1.cad37bbd5ab97p-3, 0x1.241522ca7821cp-4, 0x1.4a30727f2fa53p-4, -0x1.378a8eecf45ccp-5, -0x1.0f8eaa8ffe709p-6, 0x1.9b45682a50d7p-7, 0x1.c271f584373d4p-12, -0x1.d84a0f9cb2f31p-10, 0x1.72e5533fa10d3p-12};
        }
        if constexpr (std::is_same_v<Real, long double>) {
            return {0x9.f70e5d1c5b5851p-7L, 0xc.b0733a181801b18p-5L, 0xb.aa86694a60ff1fap-4L, 0xf.0874dd14a6ee015p-5L, -0x9.35c187248f196d4p-6L, -0xe.569bddead5cb84ap-6L, 0x9.20a91653c10e31cp-7L, 0xa.518393f97d2944bp-7L, -0x9.bc547767a2e6031p-8L, -0x8.7c75547ff384611p-9L, 0xc.da2b415286b8088p-10L, 0xe.138fac21b9ea0e9p-15L, -0xe.c2507ce59798a73p-13L, 0xb.972a99fd0869924p-15L};
        }
        #ifdef BOOST_HAS_FLOAT128
        if constexpr (std::is_same_v<Real, boost::multiprecision::float128>) {
            return {0x1.3ee1cba38b6b0a205a4c9aac651fp-4Q, 0x1.960e67430300363014da61e0a27dp-2Q, 0x1.7550cd294c1fe3f40f374c82e23dp-1Q, 0x1.e10e9ba294ddc02a47644f227fffp-2Q, -0x1.26b830e491e32da8b525e8ab6aa1p-3Q, -0x1.cad37bbd5ab970931a91ed1dfe9cp-3Q, 0x1.241522ca7821c638e85eb5e52839p-4Q, 0x1.4a30727f2fa52896ef76d0300cfp-4Q, -0x1.378a8eecf45cc0627b1c4b5cace1p-5Q, -0x1.0f8eaa8ffe708c22b1145a6b8cc8p-6Q, 0x1.9b45682a50d70110df27e4af9e36p-7Q, 0x1.c271f584373d41d177397b25f7ebp-12Q, -0x1.d84a0f9cb2f314e62607c084aab4p-10Q, 0x1.72e5533fa10d324758be56783dp-12Q};
        }
        #endif
    }
    if constexpr (p == 8) {
        if constexpr (std::is_same_v<Real, float>) {
            return {0x1.bdc64ap-5f, 0x1.40616ap-2f, 0x1.59ec46p-1f, 0x1.2bb39cp-1f, -0x1.035814p-6f, -0x1.22d4f8p-2f, 0x1.ef6f9cp-12f, 0x1.07acbcp-3f, -0x1.1c942p-6f, -0x1.692bc6p-5f, 0x1.ca215cp-7f, 0x1.1e978ep-7f, -0x1.3f2ef6p-8f, -0x1.9ac502p-12f, 0x1.622148p-11f, -0x1.ecbbbcp-14f};
        }
        if constexpr (std::is_same_v<Real, double>) {
            return {0x1.bdc64ada308ddp-5, 0x1.4061690b4c31ep-2, 0x1.59ec45992376p-1, 0x1.2bb39bedb5e28p-1, -0x1.03581459a95c6p-6, -0x1.22d4f8724d56fp-2, 0x1.ef6f9caf662bp-12, 0x1.07acbb163ba09p-3, -0x1.1c9420f07509dp-6, -0x1.692bc518a7fe2p-5, 0x1.ca215cd5b85b4p-7, 0x1.1e978df35f5fcp-7, -0x1.3f2ef6d3ac74ap-8, -0x1.9ac501798e65dp-12, 0x1.622148e2ef341p-11, -0x1.ecbbbc88e3fc3p-14};
        }
        if constexpr (std::is_same_v<Real, long double>) {
            return {0xd.ee3256d1846ea19p-8L, 0xa.030b485a618f383p-5L, 0xa.cf622cc91bb0386p-4L, 0x9.5d9cdf6daf13eaep-4L, -0x8.1ac0a2cd4ae2f65p-9L, -0x9.16a7c3926ab7ac2p-5L, 0xf.7b7ce57b3157c6dp-15L, 0x8.3d65d8b1dd0442p-6L, -0x8.e4a10783a84e546p-9L, -0xb.495e28c53ff0c66p-8L, 0xe.510ae6adc2d9ca2p-10L, 0x8.f4bc6f9afafe075p-10L, -0x9.f977b69d63a5188p-11L, -0xc.d6280bcc732e77dp-15L, 0xb.110a471779a0a94p-14L, -0xf.65dde4471fe1763p-17L};
        }
        #ifdef BOOST_HAS_FLOAT128
        if constexpr (std::is_same_v<Real, boost::multiprecision::float128>) {
            return {0x1.bdc64ada308dd431e90c51caaddfp-5Q, 0x1.4061690b4c31e70689f16560d4fdp-2Q, 0x1.59ec45992376070ce88329890c41p-1Q, 0x1.2bb39bedb5e27d5c19af764e6b08p-1Q, -0x1.03581459a95c5ec9b8a25898e79fp-6Q, -0x1.22d4f8724d56f584680effd871d8p-2Q, 0x1.ef6f9caf662af8d967c439494864p-12Q, 0x1.07acbb163ba0883fad048482f33cp-3Q, -0x1.1c9420f07509ca8be1e93402904ap-6Q, -0x1.692bc518a7fe18cb861085169737p-5Q, 0x1.ca215cd5b85b39449ec69ab26159p-7Q, 0x1.1e978df35f5fc0eabddff5ee9da6p-7Q, -0x1.3f2ef6d3ac74a30f6f64186b6be2p-8Q, -0x1.9ac501798e65cefa54cc7686cbfep-12Q, 0x1.622148e2ef341527c8d993cd7e86p-11Q, -0x1.ecbbbc88e3fc2ec652987ddfb595p-14Q};
        }
        #endif
    }
    if constexpr (p == 9) {
        if constexpr (std::is_same_v<Real, float>) {
            return {0x1.37ef3ep-5f, 0x1.f35f98p-3f, 0x1.35ab6p-1f, 0x1.50881p-1f, 0x1.10c9cap-3f, -0x1.2c4ff6p-2f, -0x1.8ca8ecp-4f, 0x1.303622p-3f, 0x1.f768dap-6f, -0x1.15062ap-4f, 0x1.07231ap-12f, 0x1.6e5f9cp-6f, -0x1.358a3ap-8f, -0x1.1897b6p-8f, 0x1.e4597cp-10f, 0x1.e3276ap-13f, -0x1.0833dap-12f, 0x1.4a11bap-15f};
        }
        if constexpr (std::is_same_v<Real, double>) {
            return {0x1.37ef3e540da7cp-5, 0x1.f35f9808bc2ap-3, 0x1.35ab60603a288p-1, 0x1.5088101e8fe35p-1, 0x1.10c9ca803fb22p-3, -0x1.2c4ff66fd53efp-2, -0x1.8ca8ebcdc98fcp-4, 0x1.303621e43e771p-3, 0x1.f768d94677997p-6, -0x1.1506294f451a2p-4, 0x1.07231a6b6ca0dp-12, 0x1.6e5f9be058887p-6, -0x1.358a39f783bbfp-8, -0x1.1897b64b3bfb6p-8, 0x1.e4597bbfc711fp-10, 0x1.e3276a3bc510bp-13, -0x1.0833da803978ap-12, 0x1.4a11ba1ad31b5p-15};
        }
        if constexpr (std::is_same_v<Real, long double>) {
            return {0x9.bf79f2a06d3e1d5p-8L, 0xf.9afcc045e1503b9p-6L, 0x9.ad5b0301d143d0bp-4L, 0xa.844080f47f1a981p-4L, 0x8.864e5401fd91397p-6L, -0x9.627fb37ea9f7708p-5L, -0xc.65475e6e4c7e255p-7L, 0x9.81b10f21f3b8936p-6L, 0xf.bb46ca33bccb6bep-9L, -0x8.a8314a7a28d0fa4p-7L, 0x8.3918d35b65065cap-15L, 0xb.72fcdf02c443b9ap-9L, -0x9.ac51cfbc1ddf90ep-11L, -0x8.c4bdb259dfdac19p-11L, 0xf.22cbddfe388f56ep-13L, 0xf.193b51de288565ap-16L, -0x8.419ed401cbc51f3p-15L, 0xa.508dd0d698da458p-18L};
        }
        #ifdef BOOST_HAS_FLOAT128
        if constexpr (std::is_same_v<Real, boost::multiprecision::float128>) {
            return {0x1.37ef3e540da7c3aa0d1dcf5bb03fp-5Q, 0x1.f35f9808bc2a0772556b997598d3p-3Q, 0x1.35ab60603a287a1609387fb893b8p-1Q, 0x1.5088101e8fe3530163beb11528dcp-1Q, 0x1.10c9ca803fb2272d1e1cf6fa2d0cp-3Q, -0x1.2c4ff66fd53eee10d811f3bf5459p-2Q, -0x1.8ca8ebcdc98fc4aa06bb294a039cp-4Q, 0x1.303621e43e77126c475b311fee41p-3Q, 0x1.f768d94677996d7b05f99f12d599p-6Q, -0x1.1506294f451a1f48db3652c40168p-4Q, 0x1.07231a6b6ca0cb934ff78371a261p-12Q, 0x1.6e5f9be058887734b49665246c12p-6Q, -0x1.358a39f783bbf21cd3bb737b9bf7p-8Q, -0x1.1897b64b3bfb583210060e933d6p-8Q, 0x1.e4597bbfc711eadc2635b08cc11dp-10Q, 0x1.e3276a3bc510acb32c14bfd1c3d9p-13Q, -0x1.0833da803978a3e6248db817fff9p-12Q, 0x1.4a11ba1ad31b48b0c2e32412ef31p-15Q};
        }
        #endif
    }
    if constexpr (p == 10) {
        if constexpr (std::is_same_v<Real, float>) {
            return {0x1.b4f654p-6f, 0x1.8162d6p-3f, 0x1.0ded5p-1f, 0x1.607db4p-1f, 0x1.1feba4p-2f, -0x1.ffaf7cp-3f, -0x1.914c48p-3f, 0x1.04da38p-3f, 0x1.7d29b8p-4f, -0x1.246e3p-4f, -0x1.e2a1dep-6f, 0x1.101406p-5f, 0x1.d8b7dcp-9f, -0x1.5fb466p-7f, 0x1.6dc878p-10f, 0x1.052608p-9f, -0x1.67962p-11f, -0x1.e87f56p-14f, 0x1.888a12p-14f, -0x1.bd12a2p-17f};
        }
        if constexpr (std::is_same_v<Real, double>) {
            return {0x1.b4f6549dc7ae3p-6, 0x1.8162d69198cfep-3, 0x1.0ded5071bf874p-1, 0x1.607db4062d775p-1, 0x1.1feba4923f567p-2, -0x1.ffaf7b6c111e3p-3, -0x1.914c47c1ca802p-3, 0x1.04da377a0ae83p-3, 0x1.7d29b819fd18dp-4, -0x1.246e307349ac4p-4, -0x1.e2a1dd5152b25p-6, 0x1.1014069cb8f3cp-5, 0x1.d8b7db3e21714p-9, -0x1.5fb466d770edcp-7, 0x1.6dc8787ae38ddp-10, 0x1.0526072a98cd8p-9, -0x1.67962098c50fp-11, -0x1.e87f555dc50ddp-14, 0x1.888a11cfae433p-14, -0x1.bd12a2a1a43dbp-17};
        }
        if constexpr (std::is_same_v<Real, long double>) {
            return {0xd.a7b2a4ee3d71687p-9L, 0xc.0b16b48cc67ed9dp-6L, 0x8.6f6a838dfc39cfcp-4L, 0xb.03eda0316bba9f7p-4L, 0x8.ff5d2491fab3641p-5L, -0xf.fd7bdb6088f17a6p-6L, -0xc.8a623e0e5400e29p-6L, 0x8.26d1bbd05741a5ep-6L, 0xb.e94dc0cfe8c6997p-7L, -0x9.2371839a4d62324p-7L, -0xf.150eea8a95926cep-9L, 0x8.80a034e5c79dfcbp-8L, 0xe.c5bed9f10b89f27p-12L, -0xa.fda336bb876e26p-10L, 0xb.6e43c3d71c6ebc5p-13L, 0x8.29303954c66bfc4p-12L, -0xb.3cb104c62878097p-14L, -0xf.43faaaee286e92ap-17L, 0xc.44508e7d72197a7p-17L, -0xd.e895150d21eda8p-20L};
        }
        #ifdef BOOST_HAS_FLOAT128
        if constexpr (std::is_same_v<Real, boost::multiprecision::float128>) {
            return {0x1.b4f6549dc7ae2d0da711a289f0cap-6Q, 0x1.8162d69198cfdb39a90922d70e5fp-3Q, 0x1.0ded5071bf8739f750ffec152978p-1Q, 0x1.607db4062d7753edab576af3e0dcp-1Q, 0x1.1feba4923f566c829244a5b4f4c2p-2Q, -0x1.ffaf7b6c111e2f4b082d1396dfd9p-3Q, -0x1.914c47c1ca801c5176763b1002edp-3Q, 0x1.04da377a0ae834bb94a9fa4f0243p-3Q, 0x1.7d29b819fd18d32d43de1a7d3641p-4Q, -0x1.246e307349ac4648ebbb94e378bbp-4Q, -0x1.e2a1dd5152b24d9c384763b949a8p-6Q, 0x1.1014069cb8f3bf95bd45be617d84p-5Q, 0x1.d8b7db3e21713e4ea2df5435db4ap-9Q, -0x1.5fb466d770edc4c0f4a29d7e3176p-7Q, 0x1.6dc8787ae38dd789b768877a70c2p-10Q, 0x1.0526072a98cd7f8751cc002da048p-9Q, -0x1.67962098c50f012db290a44515a5p-11Q, -0x1.e87f555dc50dd253428e6067701p-14Q, 0x1.888a11cfae432f4e4bd0c734fcc9p-14Q, -0x1.bd12a2a1a43db4fff0519cd821dfp-17Q};
        }
        #endif
    }
    if constexpr (p == 11) {
        if constexpr (std::is_same_v<Real, float>) {
            return {0x1.324992p-6f, 0x1.270c9cp-3f, 0x1.ccb286p-2f, 0x1.5f1256p-1f, 0x1.a5d9fcp-2f, -0x1.4c56f6p-3f, -0x1.18cff8p-2f, 0x1.0e83b8p-4f, 0x1.32d0a4p-3f, -0x1.7cc388p-5f, -0x1.10221ep-4f, 0x1.00b272p-5f, 0x1.557516p-6f, -0x1.f77976p-7f, -0x1.b5e4ap-9f, 0x1.42fd2p-8f, -0x1.439544p-12f, -0x1.d4338ep-11f, 0x1.05416p-12f, 0x1.c8ab06p-15f, -0x1.228a1p-15f, 0x1.2d9b0cp-18f};
        }
        if constexpr (std::is_same_v<Real, double>) {
            return {0x1.3249916076091p-6, 0x1.270c9c42314cfp-3, 0x1.ccb286198c9dfp-2, 0x1.5f125643212ccp-1, 0x1.a5d9fccefe816p-2, -0x1.4c56f6b2bf66cp-3, -0x1.18cff89a8cc46p-2, 0x1.0e83b8b6a4d8fp-4, 0x1.32d0a3f0ba732p-3, -0x1.7cc387e4a9a09p-5, -0x1.10221dbbeff7cp-4, 0x1.00b27276a8099p-5, 0x1.557516a958be5p-6, -0x1.f77975a6883f7p-7, -0x1.b5e49f3346a8bp-9, 0x1.42fd20a75f9acp-8, -0x1.439543c841133p-12, -0x1.d4338d3fdae3ap-11, 0x1.05415f0bc6ea2p-12, 0x1.c8ab05d193c38p-15, -0x1.228a0febb3e8cp-15, 0x1.2d9b0b4e10d78p-18};
        }
        if constexpr (std::is_same_v<Real, long double>) {
            return {0x9.924c8b03b048b06p-9L, 0x9.3864e2118a67afdp-6L, 0xe.659430cc64ef74bp-5L, 0xa.f892b2190965c91p-4L, 0xd.2ecfe677f40b2aap-5L, -0xa.62b7b595fb361ep-6L, -0x8.c67fc4d4662301ap-5L, 0x8.741dc5b526c7499p-7L, 0x9.96851f85d39931p-6L, -0xb.e61c3f254d048e4p-8L, -0x8.8110eddf7fbe097p-7L, 0x8.059393b5404c82ap-8L, 0xa.aba8b54ac5f295dp-9L, -0xf.bbcbad3441fb82ap-10L, -0xd.af24f99a3545bf1p-12L, 0xa.17e9053afcd5e0bp-11L, -0xa.1caa1e42089945fp-15L, -0xe.a19c69fed71d098p-14L, 0x8.2a0af85e37512d9p-15L, 0xe.45582e8c9e1bc7dp-18L, -0x9.14507f5d9f45ff4p-18L, 0x9.6cd85a7086bbdf3p-21L};
        }
        #ifdef BOOST_HAS_FLOAT128
        if constexpr (std::is_same_v<Real, boost::multiprecision::float128>) {
            return {0x1.324991607609160c660c469077c4p-6Q, 0x1.270c9c42314cf5fa956419fe6ca8p-3Q, 0x1.ccb286198c9dee96d12668532b6cp-2Q, 0x1.5f125643212cb922d29e40df5ae5p-1Q, 0x1.a5d9fccefe816554cfd2c1b4c2eap-2Q, -0x1.4c56f6b2bf66c3c093c0529f8f7fp-3Q, -0x1.18cff89a8cc4603466567c5566b2p-2Q, 0x1.0e83b8b6a4d8e932eca95459a69cp-4Q, 0x1.32d0a3f0ba73261f3f385a0ab52ep-3Q, -0x1.7cc387e4a9a091c75d9de8faea58p-5Q, -0x1.10221dbbeff7c12d2228c7667179p-4Q, 0x1.00b27276a80990547cbd63f84547p-5Q, 0x1.557516a958be52b928260c4c6ad5p-6Q, -0x1.f77975a6883f7054b8b01c477981p-7Q, -0x1.b5e49f3346a8b7e121e98be0fcfbp-9Q, 0x1.42fd20a75f9abc16cec7fe76464bp-8Q, -0x1.439543c8411328be7c3dbb239962p-12Q, -0x1.d4338d3fdae3a1309e0d7105695cp-11Q, 0x1.05415f0bc6ea25b1e91fb2adce4fp-12Q, 0x1.c8ab05d193c378f97fff6933dc9bp-15Q, -0x1.228a0febb3e8bfe721948bedcf9cp-15Q, 0x1.2d9b0b4e10d77be6e48f3790bbdep-18Q};
        }
        #endif
    }
    if constexpr (p == 12) {
        if constexpr (std::is_same_v<Real, float>) {
            return {0x1.ada998p-7f, 0x1.c0c89p-4f, 0x1.826962p-2f, 0x1.507c5ap-1f, 0x1.082246p-1f, -0x1.6eb4acp-5f, -0x1.43c448p-2f, -0x1.85997p-6f, 0x1.75b758p-3f, 0x1.5f3ea8p-8f, -0x1.8afc68p-4f, 0x1.63811ap-7f, 0x1.5458dcp-5f, -0x1.906176p-7f, -0x1.a4c4a6p-7f, 0x1.b7d844p-8f, 0x1.26babep-9f, -0x1.1dac02p-9f, 0x1.b73c72p-18f, 0x1.978844p-12f, -0x1.73369p-14f, -0x1.96b4a6p-16f, 0x1.acb93p-17f, -0x1.9a7502p-20f};
        }
        if constexpr (std::is_same_v<Real, double>) {
            return {0x1.ada9978d2fa57p-7, 0x1.c0c890682f879p-4, 0x1.82696272e40bdp-2, 0x1.507c59d8e6f6bp-1, 0x1.082245c297978p-1, -0x1.6eb4ac1d9855dp-5, -0x1.43c448da45b22p-2, -0x1.85996f0f3b3b2p-6, 0x1.75b757e56dd3bp-3, 0x1.5f3ea878d368ep-8, -0x1.8afc682193383p-4, 0x1.638119d1c4362p-7, 0x1.5458dbe394ebp-5, -0x1.90617513f389bp-7, -0x1.a4c4a623e2a51p-7, 0x1.b7d844bffa4e9p-8, 0x1.26babd1f8d181p-9, -0x1.1dac0186ec142p-9, 0x1.b73c724cbcd25p-18, 0x1.9788431be0bfbp-12, -0x1.7336904b8b4e8p-14, -0x1.96b4a56f63fa9p-16, 0x1.acb92f10c423ap-17, -0x1.9a7502d7dc2f3p-20};
        }
        if constexpr (std::is_same_v<Real, long double>) {
            return {0xd.6d4cbc697d2b459p-10L, 0xe.064483417c3c86cp-7L, 0xc.134b1397205e82fp-5L, 0xa.83e2cec737b5622p-4L, 0x8.41122e14bcbbcc8p-4L, -0xb.75a560ecc2ae67fp-8L, -0xa.1e2246d22d9117cp-5L, -0xc.2ccb7879d9d93dbp-9L, 0xb.adbabf2b6e9db71p-6L, 0xa.f9f543c69b47041p-11L, -0xc.57e3410c99c1b5p-7L, 0xb.1c08ce8e21b0dd8p-10L, 0xa.a2c6df1ca7582e3p-8L, -0xc.830ba89f9c4dbcbp-10L, -0xd.2625311f15287e1p-10L, 0xd.bec225ffd274461p-11L, 0x9.35d5e8fc68c04b8p-12L, -0x8.ed600c3760a133p-12L, 0xd.b9e39265e692542p-21L, 0xc.bc4218df05fd969p-15L, -0xb.99b4825c5a7426p-17L, -0xc.b5a52b7b1fd454bp-19L, 0xd.65c97886211cfe7p-20L, -0xc.d3a816bee17970fp-23L};
        }
        #ifdef BOOST_HAS_FLOAT128
        if constexpr (std::is_same_v<Real, boost::multiprecision::float128>) {
            return {0x1.ada9978d2fa568b21229688156b8p-7Q, 0x1.c0c890682f8790d7f6bb1c1aa44ep-4Q, 0x1.82696272e40bd05e57c851b9cb36p-2Q, 0x1.507c59d8e6f6ac435d1e88e090b6p-1Q, 0x1.082245c29797798f324c3b86f65cp-1Q, -0x1.6eb4ac1d9855ccfed03d4133c6fbp-5Q, -0x1.43c448da45b222f8d287f72956dcp-2Q, -0x1.85996f0f3b3b27b5f0ec32953a95p-6Q, 0x1.75b757e56dd3b6e26eea4e0a43bap-3Q, 0x1.5f3ea878d368e081a5d3cb36d936p-8Q, -0x1.8afc6821933836a0611906df913cp-4Q, 0x1.638119d1c4361bb00a2f22d1f90bp-7Q, 0x1.5458dbe394eb05c5aab3c4ceb31ap-5Q, -0x1.90617513f389b7954493854c387ep-7Q, -0x1.a4c4a623e2a50fc10158bea8a45cp-7Q, 0x1.b7d844bffa4e88c2ef52650eb0e3p-8Q, 0x1.26babd1f8d1809707edb54236067p-9Q, -0x1.1dac0186ec14265ff3ab908c13a3p-9Q, 0x1.b73c724cbcd24a84c3c597249f8dp-18Q, 0x1.9788431be0bfb2d1938d25135045p-12Q, -0x1.7336904b8b4e84c0e0396f9346efp-14Q, -0x1.96b4a56f63fa8a95abb9ecaeceaep-16Q, 0x1.acb92f10c4239fce7d915af19a6cp-17Q, -0x1.9a7502d7dc2f2e1e665f9ba74bc1p-20Q};
        }
        #endif
    }
    if constexpr (p == 13) {
        if constexpr (std::is_same_v<Real, float>) {
            return {0x1.2d8918p-7f, 0x1.53665p-4f, 0x1.3f7bf6p-2f, 0x1.38dc5p-1f, 0x1.2d82fp-1f, 0x1.644b26p-4f, -0x1.428842p-2f, -0x1.fe442cp-4f, 0x1.6f9128p-3f, 0x1.2acc8p-4f, -0x1.b16354p-4f, -0x1.b1fc6ep-6f, 0x1.cbe504p-5f, 0x1.37f29ep-9f, -0x1.86743ap-6f, 0x1.0128dp-8f, 0x1.db8098p-8f, -0x1.6a025cp-9f, -0x1.58e562p-10f, 0x1.e8ceb2p-11f, 0x1.9d26d8p-15f, -0x1.5a4cf4p-13f, 0x1.0159a8p-15f, 0x1.5e5f8p-17f, -0x1.3b708ap-18f, 0x1.183f9ep-21f};
        }
        if constexpr (std::is_same_v<Real, double>) {
            return {0x1.2d89174d37e3bp-7, 0x1.53664fc8a3d4cp-4, 0x1.3f7bf6c09afe8p-2, 0x1.38dc5001834bbp-1, 0x1.2d82ef0bea472p-1, 0x1.644b251290fb5p-4, -0x1.42884206fc5aep-2, -0x1.fe442b86ad763p-4, 0x1.6f91279c81aa4p-3, 0x1.2acc804557c8bp-4, -0x1.b163543c8eaccp-4, -0x1.b1fc6de239706p-6, 0x1.cbe5044520d5ep-5, 0x1.37f29dfe3e92cp-9, -0x1.867439245b0c7p-6, 0x1.0128d031aa3bp-8, 0x1.db80973172631p-8, -0x1.6a025cdac7e1p-9, -0x1.58e561b1f2cf6p-10, 0x1.e8ceb1ee24cabp-11, 0x1.9d26d847f17cap-15, -0x1.5a4cf360064a4p-13, 0x1.0159a865542d3p-15, 0x1.5e5f8028cd835p-17, -0x1.3b708a4c3be34p-18, 0x1.183f9db8da3ap-21};
        }
        if constexpr (std::is_same_v<Real, long double>) {
            return {0x9.6c48ba69bf1db2dp-10L, 0xa.9b327e451ea606ap-7L, 0x9.fbdfb604d7f42bfp-5L, 0x9.c6e2800c1a5d56fp-4L, 0x9.6c17785f5238ef1p-4L, 0xb.2259289487da5f5p-7L, -0xa.14421037e2d6f28p-5L, -0xf.f2215c356bb1af3p-7L, 0xb.7c893ce40d51d77p-6L, 0x9.5664022abe455a4p-7L, -0xd.8b1aa1e47565fb8p-7L, -0xd.8fe36f11cb8316bp-9L, 0xe.5f28222906af18bp-8L, 0x9.bf94eff1f49636bp-12L, -0xc.33a1c922d8636b9p-9L, 0x8.0946818d51d8108p-11L, 0xe.dc04b98b931878p-11L, -0xb.5012e6d63f07f27p-12L, -0xa.c72b0d8f967b059p-13L, 0xf.46758f712655aeep-14L, 0xc.e936c23f8be4d46p-18L, -0xa.d2679b00325201p-16L, 0x8.0acd432aa169479p-18L, 0xa.f2fc01466c1a50ap-20L, -0x9.db845261df1a01dp-21L, 0x8.c1fcedc6d1d027ap-24L};
        }
        #ifdef BOOST_HAS_FLOAT128
        if constexpr (std::is_same_v<Real, boost::multiprecision::float128>) {
            return {0x1.2d89174d37e3b6596d64d168e49p-7Q, 0x1.53664fc8a3d4c0d33efd4b4bc6fap-4Q, 0x1.3f7bf6c09afe857e70d916613d78p-2Q, 0x1.38dc5001834baade66ebb6ec9fc1p-1Q, 0x1.2d82ef0bea471de1f0a3864a5cccp-1Q, 0x1.644b251290fb4bea9c417da1f71ep-4Q, -0x1.42884206fc5ade5015e8dd0e5832p-2Q, -0x1.fe442b86ad7635e6a44fa231478p-4Q, 0x1.6f91279c81aa3aed801eca5172e8p-3Q, 0x1.2acc804557c8ab4899f3d0916648p-4Q, -0x1.b163543c8eacbf70b0092010e485p-4Q, -0x1.b1fc6de2397062d67f76ff78186p-6Q, 0x1.cbe5044520d5e315c90c0ed1b10cp-5Q, 0x1.37f29dfe3e92c6d51003ce786db8p-9Q, -0x1.867439245b0c6d71509d7a14369ep-6Q, 0x1.0128d031aa3b020f68dbd131fc68p-8Q, 0x1.db80973172630eff11f1aa6bf206p-8Q, -0x1.6a025cdac7e0fe4d3e8433e38ffbp-9Q, -0x1.58e561b1f2cf60b1ce144cfab0ffp-10Q, 0x1.e8ceb1ee24cab5dbcf0c32d1591bp-11Q, 0x1.9d26d847f17c9a8ce68c968f105dp-15Q, -0x1.5a4cf360064a4020bfd19dc8e5fp-13Q, 0x1.0159a865542d28f1113c2ad46cb5p-15Q, 0x1.5e5f8028cd834a146bf59b16d932p-17Q, -0x1.3b708a4c3be3403a7da544157c5bp-18Q, 0x1.183f9db8da3a04f4c66e4459dad2p-21Q};
        }
        #endif
    }
    if constexpr (p == 14) {
        if constexpr (std::is_same_v<Real, float>) {
            return {0x1.a7702ap-8f, 0x1.fee46p-5f, 0x1.04f778p-2f, 0x1.1bcdf2p-1f, 0x1.432b0ep-1f, 0x1.bfd66ap-3f, -0x1.163586p-2f, -0x1.be885ep-3f, 0x1.1b6ef4p-3f, 0x1.1eb29p-3f, -0x1.63524ep-4f, -0x1.251084p-4f, 0x1.c480a6p-5f, 0x1.ba103ap-6f, -0x1.ee8e8ap-6f, -0x1.6ffce6p-8f, 0x1.a3160cp-7f, -0x1.873bd2p-11f, -0x1.f89472p-9f, 0x1.1650e2p-10f, 0x1.7334fep-11f, -0x1.959f6ap-12f, -0x1.5e73f4p-15f, 0x1.20612ap-14f, -0x1.5adbf4p-17f, -0x1.26968ep-18f, 0x1.cf0cbcp-20f, -0x1.7fc91p-23f};
        }
        if constexpr (std::is_same_v<Real, double>) {
            return {0x1.a7702ace363acp-8, 0x1.fee460f5cb877p-5, 0x1.04f777f62f423p-2, 0x1.1bcdf22a1725cp-1, 0x1.432b0dc3136d4p-1, 0x1.bfd66ae4258acp-3, -0x1.1635861af071dp-2, -0x1.be885d06053bfp-3, 0x1.1b6ef32bc7359p-3, 0x1.1eb28fc03c55bp-3, -0x1.63524d6aa4cf7p-4, -0x1.2510847f3ed26p-4, 0x1.c480a659de0cdp-5, 0x1.ba103a92149f4p-6, -0x1.ee8e8a6bca9acp-6, -0x1.6ffce6192a67bp-8, 0x1.a3160ba7d924ep-7, -0x1.873bd13c8af18p-11, -0x1.f894721441f67p-9, 0x1.1650e1f6ec4c3p-10, 0x1.7334fd9e58c6ap-11, -0x1.959f69010da01p-12, -0x1.5e73f3c020201p-15, 0x1.20612a2e814f8p-14, -0x1.5adbf364f5e1cp-17, -0x1.26968e53fb76dp-18, 0x1.cf0cbb4133dbbp-20, -0x1.7fc90f0c46da2p-23};
        }
        if constexpr (std::is_same_v<Real, long double>) {
            return {0xd.3b815671b1d61adp-11L, 0xf.f72307ae5c3b8acp-8L, 0x8.27bbbfb17a114bp-5L, 0x8.de6f9150b92e008p-4L, 0xa.19586e189b69d98p-4L, 0xd.feb357212c55f71p-6L, -0x8.b1ac30d7838e67p-5L, -0xd.f442e83029df6fp-6L, 0x8.db77995e39ac42bp-6L, 0x8.f5947e01e2ad986p-6L, -0xb.1a926b55267b7d4p-7L, -0x9.288423f9f692d35p-7L, 0xe.240532cef066a7ap-8L, 0xd.d081d490a4fa2c2p-9L, -0xf.7474535e54d6257p-9L, -0xb.7fe730c9533d677p-11L, 0xd.18b05d3ec926ef4p-10L, -0xc.39de89e4578bd58p-14L, -0xf.c4a390a20fb38fap-12L, 0x8.b2870fb76261486p-13L, 0xb.99a7ecf2c634fb5p-14L, -0xc.acfb48086d004d1p-15L, -0xa.f39f9e01010090ep-18L, 0x9.030951740a7c181p-17L, -0xa.d6df9b27af0deefp-20L, -0x9.34b4729fdbb6709p-21L, 0xe.7865da099eddb46p-23L, -0xb.fe48786236d0c8ep-26L};
        }
        #ifdef BOOST_HAS_FLOAT128
        if constexpr (std::is_same_v<Real, boost::multiprecision::float128>) {
            return {0x1.a7702ace363ac35a31e24ad7f18bp-8Q, 0x1.fee460f5cb8771588fb11eee4edep-5Q, 0x1.04f777f62f422960e164233d2842p-2Q, 0x1.1bcdf22a1725c010c6b543c39465p-1Q, 0x1.432b0dc3136d3b304d5a3bf4569bp-1Q, 0x1.bfd66ae4258abee147e56af809ep-3Q, -0x1.1635861af071cce08b7655b9575p-2Q, -0x1.be885d06053beddff8849acb8697p-3Q, 0x1.1b6ef32bc7358856a3c8661f7dbbp-3Q, 0x1.1eb28fc03c55b30ce48dc30b7f77p-3Q, -0x1.63524d6aa4cf6fa795512cfa3aa5p-4Q, -0x1.2510847f3ed25a6abc4fb2ebbea9p-4Q, 0x1.c480a659de0cd4f4123d770896e6p-5Q, 0x1.ba103a92149f4583a4c79c292018p-6Q, -0x1.ee8e8a6bca9ac4ad53dedccf86b9p-6Q, -0x1.6ffce6192a67aceee5260da8193dp-8Q, 0x1.a3160ba7d924dde8cc52f58fbb97p-7Q, -0x1.873bd13c8af17ab00e30f03b73e4p-11Q, -0x1.f894721441f671f4ac26746366e6p-9Q, 0x1.1650e1f6ec4c290c75eee26e8b4ap-10Q, 0x1.7334fd9e58c69f69b95d47a6075fp-11Q, -0x1.959f69010da009a2eacf096b1dafp-12Q, -0x1.5e73f3c02020121b8b1dbc91cc8ap-15Q, 0x1.20612a2e814f83021e4b38a90018p-14Q, -0x1.5adbf364f5e1bddd8115a86dd662p-17Q, -0x1.26968e53fb76ce12ec5e07348937p-18Q, 0x1.cf0cbb4133dbb68c44095fc163dep-20Q, -0x1.7fc90f0c46da191bc033464815fbp-23Q};
        }
        #endif
    }
    if constexpr (p == 15) {
        if constexpr (std::is_same_v<Real, float>) {
            return {0x1.297006p-8f, 0x1.7eec02p-5f, 0x1.a5efd6p-3f, 0x1.f87476p-2f, 0x1.4aa806p-1f, 0x1.5b237cp-2f, -0x1.8bae9cp-3f, -0x1.27d0d6p-2f, 0x1.0b6624p-4f, 0x1.856ba4p-3f, -0x1.44f2p-5f, -0x1.c726cp-4f, 0x1.158586p-5f, 0x1.c0c324p-5f, -0x1.a62aaap-6f, -0x1.54f3aep-6f, 0x1.ee4514p-7f, 0x1.4e4c96p-8f, -0x1.a92e2p-8f, -0x1.fb0008p-13f, 0x1.fd6e44p-10f, -0x1.879fe8p-12f, -0x1.79081p-12f, 0x1.46f04ap-13f, 0x1.b0baccp-16f, -0x1.d7ff96p-16f, 0x1.c35f5cp-19f, 0x1.e6358ep-20f, -0x1.532292p-21f, 0x1.076d02p-24f};
        }
        if constexpr (std::is_same_v<Real, double>) {
            return {0x1.2970058a26dcep-8, 0x1.7eec010bca4c4p-5, 0x1.a5efd6f02857ep-3, 0x1.f8747691d3f73p-2, 0x1.4aa8051a530b4p-1, 0x1.5b237b0c28179p-2, -0x1.8bae9ca895b45p-3, -0x1.27d0d6e28e01fp-2, 0x1.0b6623378c72fp-4, 0x1.856ba3f0d1d5bp-3, -0x1.44f200621040dp-5, -0x1.c726bf22e4a16p-4, 0x1.15858527779bcp-5, 0x1.c0c32426f4663p-5, -0x1.a62aa972a6fc4p-6, -0x1.54f3ad3a0b7c6p-6, 0x1.ee4513500347p-7, 0x1.4e4c95b98ec2dp-8, -0x1.a92e1fc2fdcaep-8, -0x1.fb000715d8e65p-13, 0x1.fd6e43c3dd5e8p-10, -0x1.879fe7f24ba0ap-12, -0x1.79080f71ee322p-12, 0x1.46f04a6cb59p-13, 0x1.b0bacca20c60fp-16, -0x1.d7ff965f47a09p-16, 0x1.c35f5c1543372p-19, 0x1.e6358dfe44b53p-20, -0x1.532291d1cbfa4p-21, 0x1.076d016337019p-24};
        }
        if constexpr (std::is_same_v<Real, long double>) {
            return {0x9.4b802c5136e6cb1p-11L, 0xb.f760085e52620cp-8L, 0xd.2f7eb78142bee86p-6L, 0xf.c3a3b48e9fb9b9bp-5L, 0xa.554028d2985a2bep-4L, 0xa.d91bd86140bc88ap-5L, -0xc.5d74e544ada28ecp-6L, -0x9.3e86b714700f9f8p-5L, 0x8.5b3119bc6397595p-7L, 0xc.2b5d1f868ead57ap-6L, -0xa.2790031082066ebp-8L, -0xe.3935f917250aff5p-7L, 0x8.ac2c293bbcde0e2p-8L, 0xe.06192137a331763p-8L, -0xd.31554b9537e1c3bp-9L, -0xa.a79d69d05be2dedp-9L, 0xf.72289a801a381d6p-10L, 0xa.7264adcc7616af1p-11L, -0xd.4970fe17ee56c45p-11L, -0xf.d80038aec732ac5p-16L, 0xf.eb721e1eeaf4368p-13L, -0xc.3cff3f925d05241p-15L, -0xb.c8407b8f7191299p-15L, 0xa.37825365ac803f5p-16L, 0xd.85d665106307a4ap-19L, -0xe.bffcb2fa3d0481cp-19L, 0xe.1afae0aa19b8c5bp-22L, 0xf.31ac6ff225a9413p-23L, -0xa.99148e8e5fd1f6dp-24L, 0x8.3b680b19b80c80ep-27L};
        }
        #ifdef BOOST_HAS_FLOAT128
        if constexpr (std::is_same_v<Real, boost::multiprecision::float128>) {
            return {0x1.2970058a26dcd961f2ddf52bca08p-8Q, 0x1.7eec010bca4c417f3cf45c166751p-5Q, 0x1.a5efd6f02857dd0be400d7c02e2cp-3Q, 0x1.f8747691d3f737355530dd281d86p-2Q, 0x1.4aa8051a530b457c39051251b937p-1Q, 0x1.5b237b0c281791146b95d35023c4p-2Q, -0x1.8bae9ca895b451d87fbb9a690248p-3Q, -0x1.27d0d6e28e01f3ef0606971828aap-2Q, 0x1.0b6623378c72eb295fa0f634a67bp-4Q, 0x1.856ba3f0d1d5aaf4f59a12cbb286p-3Q, -0x1.44f200621040cdd51fafc7035006p-5Q, -0x1.c726bf22e4a15fea85834dbb1052p-4Q, 0x1.15858527779bc1c4ec7f8b6fdd21p-5Q, 0x1.c0c32426f4662ec60f158af78296p-5Q, -0x1.a62aa972a6fc3876e9915f565de4p-6Q, -0x1.54f3ad3a0b7c5bd971d938015013p-6Q, 0x1.ee451350034703acaa0dbe95c1aep-7Q, 0x1.4e4c95b98ec2d5e117559a2b0da9p-8Q, -0x1.a92e1fc2fdcad88a36c4f23f938p-8Q, -0x1.fb000715d8e6558a8e3b0bc34d7fp-13Q, 0x1.fd6e43c3dd5e86d041c185aa9428p-10Q, -0x1.879fe7f24ba0a4818b2d76fbd7e6p-12Q, -0x1.79080f71ee322531c0337932b9a5p-12Q, 0x1.46f04a6cb59007eab6af7b114fa4p-13Q, 0x1.b0bacca20c60f4939edefbff1354p-16Q, -0x1.d7ff965f47a09038cb01ea5da5e9p-16Q, 0x1.c35f5c15433718b508f832bbf1b8p-19Q, 0x1.e6358dfe44b5282546827ce8042cp-20Q, -0x1.532291d1cbfa3eda1a3275ffd36dp-21Q, 0x1.076d01633701901bb8cc971717bp-24Q};
        }
        #endif
    }
    if constexpr (p == 16) {
        if constexpr (std::is_same_v<Real, float>) {
            return {0x1.a2048p-9f, 0x1.1df6c8p-5f, 0x1.520d3ap-3f, 0x1.b8a3e6p-2f, 0x1.465392p-1f, 0x1.c2db72p-2f, -0x1.6f9ed6p-4f, -0x1.4ee9bp-2f, -0x1.c96974p-6f, 0x1.b084cp-3f, 0x1.bff16p-6f, -0x1.0f219ap-3f, -0x1.98ed2ep-8f, 0x1.36fc54p-4f, -0x1.f159dcp-8f, -0x1.2e3094p-5f, 0x1.516f08p-7f, 0x1.ca8c3ep-7f, -0x1.ca18fcp-8f, -0x1.dda9bcp-9f, 0x1.99ff0cp-9f, 0x1.abb604p-12f, -0x1.ed5dcep-11f, 0x1.df29e6p-14f, 0x1.6e8e3p-13f, -0x1.0000dep-14f, -0x1.d3f062p-17f, 0x1.7c64bap-17f, -0x1.1821aep-20f, -0x1.8b5554p-21f, 0x1.efcecp-23f, -0x1.6a61bcp-26f};
        }
        if constexpr (std::is_same_v<Real, double>) {
            return {0x1.a2047f25eec86p-9, 0x1.1df6c86cbec6p-5, 0x1.520d391a94107p-3, 0x1.b8a3e5feb4879p-2, 0x1.465391b465934p-1, 0x1.c2db72f84c204p-2, -0x1.6f9ed69c399c7p-4, -0x1.4ee9af39c16ecp-2, -0x1.c96973b93bd86p-6, 0x1.b084bf1dc763bp-3, 0x1.bff160aba4bc3p-6, -0x1.0f21999627c5fp-3, -0x1.98ed2d2871d0ap-8, 0x1.36fc54ec8b106p-4, -0x1.f159dcb973ff4p-8, -0x1.2e3093b602a51p-5, 0x1.516f07b10e2cp-7, 0x1.ca8c3dcce54c3p-7, -0x1.ca18fbf2cc141p-8, -0x1.dda9bb9568c15p-9, 0x1.99ff0c3f2ed81p-9, 0x1.abb6031610513p-12, -0x1.ed5dcd1b49a3fp-11, 0x1.df29e5ea12845p-14, 0x1.6e8e301063e5ap-13, -0x1.0000dea4283e9p-14, -0x1.d3f0626b575afp-17, 0x1.7c64ba15d2346p-17, -0x1.1821ad345e7cp-20, -0x1.8b555407022a1p-21, 0x1.efcebf5a99b7bp-23, -0x1.6a61bcee28103p-26};
        }
        if constexpr (std::is_same_v<Real, long double>) {
            return {0xd.1023f92f7642ed3p-12L, 0x8.efb64365f630239p-8L, 0xa.9069c8d4a08344dp-6L, 0xd.c51f2ff5a43c97dp-5L, 0xa.329c8da32c99c9dp-4L, 0xe.16db97c26101ca7p-5L, -0xb.7cf6b4e1cce36cbp-7L, -0xa.774d79ce0b75fdfp-5L, -0xe.4b4b9dc9dec3204p-9L, 0xd.8425f8ee3b1d5bfp-6L, 0xd.ff8b055d25e15c1p-9L, -0x8.790cccb13e2f58ep-6L, -0xc.c76969438e84e18p-11L, 0x9.b7e2a76458832d7p-7L, -0xf.8acee5cb9ff9e4cp-11L, -0x9.71849db015288a6p-8L, 0xa.8b783d88715ff76p-10L, 0xe.5461ee672a6184ep-10L, -0xe.50c7df9660a0b0ep-11L, -0xe.ed4ddcab460a85bp-12L, 0xc.cff861f976c0535p-12L, 0xd.5db018b082895d2p-15L, -0xf.6aee68da4d1f9e1p-14L, 0xe.f94f2f509422907p-17L, 0xb.747180831f2d336p-16L, -0x8.0006f52141f4728p-17L, -0xe.9f83135abad7bfp-20L, 0xb.e325d0ae91a325cp-20L, -0x8.c10d69a2f3dfe1fp-23L, -0xc.5aaaa0381150b42p-24L, 0xf.7e75fad4cdbdbc5p-26L, -0xb.530de771408178bp-29L};
        }
        #ifdef BOOST_HAS_FLOAT128
        if constexpr (std::is_same_v<Real, boost::multiprecision::float128>) {
            return {0x1.a2047f25eec85da6e57cd7f20965p-9Q, 0x1.1df6c86cbec60472edd73314093ap-5Q, 0x1.520d391a9410689a221d2a40ed0cp-3Q, 0x1.b8a3e5feb48792faed491d08e637p-2Q, 0x1.465391b4659339393dc6b44dd0b9p-1Q, 0x1.c2db72f84c20394e90d8314fee8dp-2Q, -0x1.6f9ed69c399c6d957cef0caf36cp-4Q, -0x1.4ee9af39c16ebfbef0c302e6791dp-2Q, -0x1.c96973b93bd864088cf77c5363e8p-6Q, 0x1.b084bf1dc763ab7d34b03cf7014cp-3Q, 0x1.bff160aba4bc2b81b5cb21d8ed3dp-6Q, -0x1.0f21999627c5eb1bdc396cd5510bp-3Q, -0x1.98ed2d2871d09c2ff301236b3aa5p-8Q, 0x1.36fc54ec8b1065ad541d9561e805p-4Q, -0x1.f159dcb973ff3c971dea81777e8cp-8Q, -0x1.2e3093b602a5114cced31fb3d21p-5Q, 0x1.516f07b10e2bfeeb35e5d488972dp-7Q, 0x1.ca8c3dcce54c309c6bfdafe5b818p-7Q, -0x1.ca18fbf2cc14161cd462d6382e8dp-8Q, -0x1.dda9bb9568c150b5b1b563e2a704p-9Q, 0x1.99ff0c3f2ed80a696d3f49d1ae83p-9Q, 0x1.abb6031610512ba452ab2c15607bp-12Q, -0x1.ed5dcd1b49a3f3c23377158c2ae9p-11Q, 0x1.df29e5ea1284520ed9350c7d7e2fp-14Q, 0x1.6e8e301063e5a66bc2dcbb672407p-13Q, -0x1.0000dea4283e8e4f97e320b5edbp-14Q, -0x1.d3f0626b575af7e0f5c7b70618f9p-17Q, 0x1.7c64ba15d23464b7681824c41ca1p-17Q, -0x1.1821ad345e7bfc3ef33cc03c406ap-20Q, -0x1.8b555407022a1683907685f46c79p-21Q, 0x1.efcebf5a99b7b789a7792201843ep-23Q, -0x1.6a61bcee28102f153c8290cebfbbp-26Q};
        }
        #endif
    }
    if constexpr (p == 17) {
        if constexpr (std::is_same_v<Real, float>) {
            return {0x1.25d69p-9f, 0x1.a9bea4p-6f, 0x1.0cba66p-3f, 0x1.7b3d38p-2f, 0x1.38d48cp-1f, 0x1.0960bp-1f, 0x1.bf874ap-6f, -0x1.50335p-2f, -0x1.0346bcp-3f, 0x1.941794p-3f, 0x1.9e404p-4f, -0x1.03b7f2p-3f, -0x1.d3b162p-5f, 0x1.4c35cap-4f, 0x1.6d90b8p-6f, -0x1.80637cp-5f, -0x1.acbb0ep-9f, 0x1.7477f4p-6f, -0x1.8ed9ccp-9f, -0x1.19e68ap-7f, 0x1.850572p-9f, 0x1.2d9fa2p-9f, -0x1.78a90ep-10f, -0x1.581268p-12f, 0x1.ccd01ap-12f, -0x1.adaa96p-16f, -0x1.582268p-14f, 0x1.85029cp-16f, 0x1.d5219cp-18f, -0x1.2e638p-18f, 0x1.43e648p-22f, 0x1.3d94aap-22f, -0x1.69ce4ap-24f, 0x1.f36b16p-28f};
        }
        if constexpr (std::is_same_v<Real, double>) {
            return {0x1.25d68f81bd86cp-9, 0x1.a9bea4085cef4p-6, 0x1.0cba6633a3a6dp-3, 0x1.7b3d3861245dp-2, 0x1.38d48c64aa0dbp-1, 0x1.0960af0f4068ep-1, 0x1.bf874a24b1d56p-6, -0x1.50335073492d2p-2, -0x1.0346bb1b97634p-3, 0x1.9417930ab08a3p-3, 0x1.9e403f27b713fp-4, -0x1.03b7f2067576fp-3, -0x1.d3b162866891bp-5, 0x1.4c35ca8307eb6p-4, 0x1.6d90b88bab06p-6, -0x1.80637c4a935c8p-5, -0x1.acbb0e1b24a77p-9, 0x1.7477f35c89773p-6, -0x1.8ed9cce7691cbp-9, -0x1.19e6894a383d1p-7, 0x1.8505717311392p-9, 0x1.2d9fa2834b37bp-9, -0x1.78a90e5fc8eecp-10, -0x1.58126708ae49dp-12, 0x1.ccd01a3d65145p-12, -0x1.adaa9541cd341p-16, -0x1.5822674fd0b1ap-14, 0x1.85029b61ff0fcp-16, 0x1.d5219c45df32ep-18, -0x1.2e63809757f53p-18, 0x1.43e648c7ce72fp-22, 0x1.3d94aa0876b2dp-22, -0x1.69ce4aed2f58ep-24, 0x1.f36b16a008a29p-28};
        }
        if constexpr (std::is_same_v<Real, long double>) {
            return {0x9.2eb47c0dec36183p-12L, 0xd.4df52042e779eb7p-9L, 0x8.65d3319d1d36a45p-6L, 0xb.d9e9c30922e8138p-5L, 0x9.c6a46325506d79dp-4L, 0x8.4b05787a03472d6p-4L, 0xd.fc3a51258eaaed3p-9L, -0xa.819a839a4968eafp-5L, -0x8.1a35d8dcbb19dd2p-6L, 0xc.a0bc9855845188cp-6L, 0xc.f201f93db89f96ep-7L, -0x8.1dbf9033abb79d9p-6L, -0xe.9d8b1433448d469p-8L, 0xa.61ae54183f5b11cp-7L, 0xb.6c85c45d582fc12p-9L, -0xc.031be2549ae3edfp-8L, -0xd.65d870d9253b88dp-12L, 0xb.a3bf9ae44bb963dp-9L, -0xc.76ce673b48e5651p-12L, -0x8.cf344a51c1e8512p-10L, 0xc.282b8b9889c92c5p-12L, 0x9.6cfd141a59bd766p-12L, -0xb.c54872fe4775c65p-13L, -0xa.c0933845724e821p-15L, 0xe.6680d1eb28a2645p-15L, -0xd.6d54aa0e69a09cdp-19L, -0xa.c1133a7e858d207p-17L, 0xc.2814db0ff87dcc6p-19L, 0xe.a90ce22ef996e32p-21L, -0x9.731c04babfa9664p-21L, 0xa.1f32463e73979d5p-25L, 0x9.eca55043b596609p-25L, -0xb.4e7257697ac7376p-27L, 0xf.9b58b5004514458p-31L};
        }
        #ifdef BOOST_HAS_FLOAT128
        if constexpr (std::is_same_v<Real, boost::multiprecision::float128>) {
            return {0x1.25d68f81bd86c306bc9043f43f4p-9Q, 0x1.a9bea4085cef3d6d79a074d9709ep-6Q, 0x1.0cba6633a3a6d489f4c961323b24p-3Q, 0x1.7b3d3861245d02709eb04b0b2b39p-2Q, 0x1.38d48c64aa0daf396243783e95bep-1Q, 0x1.0960af0f4068e5ac2ea6e7c9135p-1Q, 0x1.bf874a24b1d55da571e39f2bf7e7p-6Q, -0x1.50335073492d1d5e6945d8cdee28p-2Q, -0x1.0346bb1b97633ba3dcc115a9814cp-3Q, 0x1.9417930ab08a311830b19e500034p-3Q, 0x1.9e403f27b713f2db7263ac22d15fp-4Q, -0x1.03b7f2067576f3b2927e9c5889d3p-3Q, -0x1.d3b162866891a8d294e6e9eb979p-5Q, 0x1.4c35ca8307eb62387847e8508ae5p-4Q, 0x1.6d90b88bab05f824bb4e728b7a1p-6Q, -0x1.80637c4a935c7dbdd09fd2c6d8fdp-5Q, -0x1.acbb0e1b24a7711aad0c2cc58d1fp-9Q, 0x1.7477f35c89772c79eec161a3cdd2p-6Q, -0x1.8ed9cce7691caca1362d6d7abefcp-9Q, -0x1.19e6894a383d0a23e4d18574889p-7Q, 0x1.85057173113925897d4027a40b6ep-9Q, 0x1.2d9fa2834b37aecc772ff5c2d7b1p-9Q, -0x1.78a90e5fc8eeb8c9dbcf2179c3bbp-10Q, -0x1.58126708ae49d0421a9c690b913dp-12Q, 0x1.ccd01a3d65144c8ab80302b0f5dbp-12Q, -0x1.adaa9541cd34139a00a6638e9bf6p-16Q, -0x1.5822674fd0b1a40e647125cde434p-14Q, 0x1.85029b61ff0fb98b3294c269c9e6p-16Q, 0x1.d5219c45df32dc63e76e89984d9dp-18Q, -0x1.2e63809757f52cc8ed09e415ace1p-18Q, 0x1.43e648c7ce72f3a94ea0d55898bap-22Q, 0x1.3d94aa0876b2cc113a8ea22e3eecp-22Q, -0x1.69ce4aed2f58e6ecb8a211f55ad1p-24Q, 0x1.f36b16a008a288af85eb1e1cec02p-28Q};
        }
        #endif
    }
    if constexpr (p == 18) {
        if constexpr (std::is_same_v<Real, float>) {
            return {0x1.9d3864p-10f, 0x1.3c05f8p-6f, 0x1.a84c62p-4f, 0x1.423b32p-2f, 0x1.24c67cp-1f, 0x1.24c33p-1f, 0x1.2d835p-3f, -0x1.2cb3a6p-2f, -0x1.bb5a5cp-3f, 0x1.323edep-3f, 0x1.562ebap-3f, -0x1.7a31p-4f, -0x1.b541d8p-4f, 0x1.09c72ep-4f, 0x1.d35d24p-5f, -0x1.6cc216p-5f, -0x1.84d84cp-6f, 0x1.b4f90cp-6f, 0x1.9a65bep-8f, -0x1.ababc2p-7f, 0x1.f19208p-14f, 0x1.43f78cp-8f, -0x1.2544ep-10f, -0x1.5f6de6p-10f, 0x1.497f3cp-11f, 0x1.bfe9bap-13f, -0x1.a098a2p-13f, -0x1.49d5fp-23f, 0x1.39d678p-15f, -0x1.1de76cp-17f, -0x1.bf4c72p-19f, 0x1.dac908p-20f, -0x1.4a5a66p-24f, -0x1.f9216ep-24f, 0x1.079c6ap-25f, -0x1.58b01ap-29f};
        }
        if constexpr (std::is_same_v<Real, double>) {
            return {0x1.9d386358a1875p-10, 0x1.3c05f73c5b367p-6, 0x1.a84c610f33b44p-4, 0x1.423b3247213b9p-2, 0x1.24c67bbcaad51p-1, 0x1.24c32fbf17436p-1, 0x1.2d834f99026eap-3, -0x1.2cb3a5153503ap-2, -0x1.bb5a5b1e13b57p-3, 0x1.323ede758da42p-3, 0x1.562eba30ba3dbp-3, -0x1.7a30ff6ce5272p-4, -0x1.b541d802db16ep-4, 0x1.09c72d78f0163p-4, 0x1.d35d236774d12p-5, -0x1.6cc216694aaa1p-5, -0x1.84d84c9eacf2cp-6, 0x1.b4f90c6dc9842p-6, 0x1.9a65be90d86a5p-8, -0x1.ababc1ea3f978p-7, 0x1.f19207685977fp-14, 0x1.43f78b1e4d853p-8, -0x1.2544e0deeca72p-10, -0x1.5f6de508b2306p-10, 0x1.497f3c70a07ccp-11, 0x1.bfe9ba59913cep-13, -0x1.a098a12371c0dp-13, -0x1.49d5efd42da9p-23, 0x1.39d678c2c06afp-15, -0x1.1de76b2193bedp-17, -0x1.bf4c72e2a2d94p-19, 0x1.dac907ddfd9dbp-20, -0x1.4a5a655290e68p-24, -0x1.f9216d4e4d83bp-24, 0x1.079c6a0fc4518p-25, -0x1.58b0195a8266dp-29};
        }
        if constexpr (std::is_same_v<Real, long double>) {
            return {0xc.e9c31ac50c3a483p-13L, 0x9.e02fb9e2d9b3a9bp-9L, 0xd.426308799da225cp-7L, 0xa.11d9923909dc42bp-5L, 0x9.2633dde556a88a3p-4L, 0x9.26197df8ba1b1b2p-4L, 0x9.6c1a7cc81374c03p-6L, -0x9.659d28a9a81cd3ap-5L, -0xd.dad2d8f09dab48bp-6L, 0x9.91f6f3ac6d21107p-6L, 0xa.b175d185d1ed60ap-6L, -0xb.d187fb672938eb2p-7L, -0xd.aa0ec016d8b6c08p-7L, 0x8.4e396bc780b1464p-7L, 0xe.9ae91b3ba68905ep-8L, -0xb.6610b34a5550677p-8L, -0xc.26c264f56795d72p-9L, 0xd.a7c8636e4c20f3dp-9L, 0xc.d32df486c3529e3p-11L, -0xd.5d5e0f51fcbc222p-10L, 0xf.8c903b42cbbfa1ep-17L, 0xa.1fbc58f26c29bbap-11L, -0x9.2a2706f76538d5fp-13L, -0xa.fb6f284591832p-13L, 0xa.4bf9e38503e5e52p-14L, 0xd.ff4dd2cc89e6cb3p-16L, -0xd.04c5091b8e0651fp-16L, -0xa.4eaf7ea16d480ecp-26L, 0x9.ceb3c6160357432p-18L, -0x8.ef3b590c9df6482p-20L, -0xd.fa63971516c9e77p-22L, 0xe.d6483eefeced63bp-23L, -0xa.52d32a948733d34p-27L, -0xf.c90b6a726c1d6d6p-27L, 0x8.3ce3507e228be57p-28L, -0xa.c580cad41336728p-32L};
        }
        #ifdef BOOST_HAS_FLOAT128
        if constexpr (std::is_same_v<Real, boost::multiprecision::float128>) {
            return {0x1.9d386358a187490511616014e787p-10Q, 0x1.3c05f73c5b367536d2be9893ed8cp-6Q, 0x1.a84c610f33b444b779a9b1758c69p-4Q, 0x1.423b3247213b88555f9f2ada7d8ap-2Q, 0x1.24c67bbcaad511453636985416f3p-1Q, 0x1.24c32fbf17436364ee6d86ce0b36p-1Q, 0x1.2d834f99026e9806c6f6c1aa6cb3p-3Q, -0x1.2cb3a51535039a741562ba69c7a8p-2Q, -0x1.bb5a5b1e13b56915364cb7d6d1bfp-3Q, 0x1.323ede758da4220dcba304aad9d9p-3Q, 0x1.562eba30ba3dac1311ece1a1884ap-3Q, -0x1.7a30ff6ce5271d64766823319bd4p-4Q, -0x1.b541d802db16d80f141c678b2993p-4Q, 0x1.09c72d78f01628c7ccd127058828p-4Q, 0x1.d35d236774d120bc0fddfe320959p-5Q, -0x1.6cc216694aaa0ced2f850c2b3eb7p-5Q, -0x1.84d84c9eacf2bae4fd091b9bdb4fp-6Q, 0x1.b4f90c6dc9841e7a702107e1aca3p-6Q, 0x1.9a65be90d86a53c665da90605941p-8Q, -0x1.ababc1ea3f978443c6511c976c7fp-7Q, 0x1.f19207685977f43b68af0530733ap-14Q, 0x1.43f78b1e4d853774466d94263915p-8Q, -0x1.2544e0deeca71abd544eac93fa45p-10Q, -0x1.5f6de508b23063ffd5522faa2f1ap-10Q, 0x1.497f3c70a07cbca34a7835790cbfp-11Q, 0x1.bfe9ba59913cd96563f593913af3p-13Q, -0x1.a098a12371c0ca3dcc2d68b00891p-13Q, -0x1.49d5efd42da901d7159610c4d152p-23Q, 0x1.39d678c2c06ae863ef3b97fa7d6p-15Q, -0x1.1de76b2193bec903a9d077f1666ep-17Q, -0x1.bf4c72e2a2d93cedea468c1c5543p-19Q, 0x1.dac907ddfd9dac76e76088418dc5p-20Q, -0x1.4a5a655290e67a67e8d8672b7916p-24Q, -0x1.f9216d4e4d83adac1cb13466b74bp-24Q, 0x1.079c6a0fc4517cae3d72b9596e6p-25Q, -0x1.58b0195a8266ce508d99fc62d4eep-29Q};
        }
        #endif
    }
    if constexpr (p == 19) {
        if constexpr (std::is_same_v<Real, float>) {
            return {0x1.22a192p-10f, 0x1.d3f68ap-7f, 0x1.4cea48p-4f, 0x1.0ebbd8p-2f, 0x1.0c82ecp-1f, 0x1.34129ep-1f, 0x1.0b280cp-2f, -0x1.d32194p-3f, -0x1.24b2e2p-2f, 0x1.31c694p-4f, 0x1.b2e46cp-3f, -0x1.12957ap-5f, -0x1.246cd4p-3f, 0x1.c3f126p-6f, 0x1.63f856p-4f, -0x1.b2323ep-6f, -0x1.7629ccp-5f, 0x1.6248a8p-6f, 0x1.3d72f2p-6f, -0x1.ca5f1cp-7f, -0x1.807eap-8f, 0x1.cd6c24p-8f, 0x1.93274ep-11f, -0x1.604346p-9f, 0x1.66699p-12f, 0x1.81c5bep-11f, -0x1.1156b8p-12f, -0x1.054e8ap-13f, 0x1.6d608ep-14f, 0x1.56a79p-18f, -0x1.172d04p-16f, 0x1.941ff2p-19f, 0x1.9b3988p-20f, -0x1.7070fep-21f, 0x1.f1373cp-27f, 0x1.8e4f58p-25f, -0x1.7f97c4p-27f, 0x1.dc770ep-31f};
        }
        if constexpr (std::is_same_v<Real, double>) {
            return {0x1.22a1917fed85ap-10, 0x1.d3f6896603bf4p-7, 0x1.4cea4765b114ap-4, 0x1.0ebbd74f124f2p-2, 0x1.0c82ecc416391p-1, 0x1.34129e60b6e3p-1, 0x1.0b280be54f8b7p-2, -0x1.d32194b53c804p-3, -0x1.24b2e1dd4c73ep-2, 0x1.31c6940f7edfap-4, 0x1.b2e46c1684042p-3, -0x1.12957a28f5bb3p-5, -0x1.246cd3943959bp-3, 0x1.c3f126f461d98p-6, 0x1.63f8568e9e6ecp-4, -0x1.b2323dbfdcbe7p-6, -0x1.7629cb842f99bp-5, 0x1.6248a775ca66p-6, 0x1.3d72f2476d004p-6, -0x1.ca5f1bf2f3244p-7, -0x1.807e9f533b217p-8, 0x1.cd6c23ae5d8f9p-8, 0x1.93274e99dc9b5p-11, -0x1.6043462e3e51ap-9, 0x1.66698fcf3c249p-12, 0x1.81c5be1857466p-11, -0x1.1156b7a432a66p-12, -0x1.054e8a6af00d7p-13, 0x1.6d608dedb3cf6p-14, 0x1.56a78fb9bd527p-18, -0x1.172d0353b3f12p-16, 0x1.941ff236730f1p-19, 0x1.9b39878b95663p-20, -0x1.7070fedc9cdedp-21, 0x1.f1373b9dbeafap-27, 0x1.8e4f570c878dep-25, -0x1.7f97c3aa5d9ep-27, 0x1.dc770dc2ed00bp-31};
        }
        if constexpr (std::is_same_v<Real, long double>) {
            return {0x9.150c8bff6c2cf0bp-13L, 0xe.9fb44b301dfa184p-10L, 0xa.67523b2d88a4cd2p-7L, 0x8.75deba789279273p-5L, 0x8.64176620b1c8a9ap-4L, 0x9.a094f305b717f7bp-4L, 0x8.59405f2a7c5b569p-5L, -0xe.990ca5a9e402365p-6L, -0x9.25970eea639ef9ap-5L, 0x9.8e34a07bf6fd17cp-7L, 0xd.972360b420211e3p-6L, -0x8.94abd147add99bfp-8L, -0x9.23669ca1cacd462p-6L, 0xe.1f8937a30ecc2f7p-9L, 0xb.1fc2b474f37606ap-7L, -0xd.9191edfee5f3727p-9L, -0xb.b14e5c217ccd592p-8L, 0xb.12453bae532fcd5p-9L, 0x9.eb97923b6802088p-9L, -0xe.52f8df979921d82p-10L, -0xc.03f4fa99d90b962p-11L, 0xe.6b611d72ec7ca7cp-11L, 0xc.993a74cee4da67fp-14L, -0xb.021a3171f28cdccp-12L, 0xb.334c7e79e124bccp-15L, 0xc.0e2df0c2ba32ef3p-14L, -0x8.8ab5bd21953315ep-15L, -0x8.2a745357806b795p-16L, 0xb.6b046f6d9e7b04p-17L, 0xa.b53c7dcdea9365p-21L, -0x8.b9681a9d9f88fe1p-19L, 0xc.a0ff91b39878443p-22L, 0xc.d9cc3c5cab31a21p-23L, -0xb.8387f6e4e6f68p-24L, 0xf.89b9dcedf57d09dp-30L, 0xc.727ab8643c6ec27p-28L, -0xb.fcbe1d52ecefc57p-30L, 0xe.e3b86e176805b46p-34L};
        }
        #ifdef BOOST_HAS_FLOAT128
        if constexpr (std::is_same_v<Real, boost::multiprecision::float128>) {
            return {0x1.22a1917fed859e160446bb7b167ap-10Q, 0x1.d3f6896603bf4308e8ccc5a30641p-7Q, 0x1.4cea4765b11499a3de1dc5fdc4b7p-4Q, 0x1.0ebbd74f124f24e597d3a5bc26dbp-2Q, 0x1.0c82ecc4163915332b81143283e8p-1Q, 0x1.34129e60b6e2fef58b66c45f7099p-1Q, 0x1.0b280be54f8b6ad2d41e0f33532ep-2Q, -0x1.d32194b53c8046c9a7b9b8714cebp-3Q, -0x1.24b2e1dd4c73df34e82deb7bb6c6p-2Q, 0x1.31c6940f7edfa2f8c6ddc9cd2aefp-4Q, 0x1.b2e46c16840423c52b33ff7476eep-3Q, -0x1.12957a28f5bb337e3e66ea6f23e6p-5Q, -0x1.246cd3943959a8c450bf5ca92b86p-3Q, 0x1.c3f126f461d985ee9de87f3424e5p-6Q, 0x1.63f8568e9e6ec0d34220dce2e186p-4Q, -0x1.b2323dbfdcbe6e4ed8b2320fb129p-6Q, -0x1.7629cb842f99ab236f8e908d21c2p-5Q, 0x1.6248a775ca65f9a99ecab76ed61p-6Q, 0x1.3d72f2476d00410fe80cb42bd72bp-6Q, -0x1.ca5f1bf2f3243b0308970f1f944p-7Q, -0x1.807e9f533b2172c38f0c60523009p-8Q, 0x1.cd6c23ae5d8f94f8e6c997306168p-8Q, 0x1.93274e99dc9b4cfea213c9163decp-11Q, -0x1.6043462e3e519b97ce417a7ac3b3p-9Q, 0x1.66698fcf3c2497984ab865b5b6acp-12Q, 0x1.81c5be1857465de5530984af43cbp-11Q, -0x1.1156b7a432a662bbb70ff712e064p-12Q, -0x1.054e8a6af00d6f29e88cd98738bep-13Q, 0x1.6d608dedb3cf6080fb9e6ef72893p-14Q, 0x1.56a78fb9bd526ca09fc616b1b8ap-18Q, -0x1.172d0353b3f11fc1e8818be3f7a7p-16Q, 0x1.941ff236730f0886ee4074951d94p-19Q, 0x1.9b39878b9566344204396825f773p-20Q, -0x1.7070fedc9cdecfffd4c6dc5d1fe5p-21Q, 0x1.f1373b9dbeafa1390279c531a649p-27Q, 0x1.8e4f570c878dd84d17b935d5f69ap-25Q, -0x1.7f97c3aa5d9df8ae90e909fae15p-27Q, 0x1.dc770dc2ed00b68b9c525506b211p-31Q};
        }
        #endif
    }
    if constexpr (p == 20) {
        if constexpr (std::is_same_v<Real, float>) {
            return {0x1.98eb9ap-11f, 0x1.59aebcp-7f, 0x1.03c8a8p-4f, 0x1.c27102p-3f, 0x1.e40a78p-2f, 0x1.389292p-1f, 0x1.722da8p-2f, -0x1.1d1b3ap-3f, -0x1.4ea132p-2f, -0x1.120e7ep-6f, 0x1.d38a42p-3f, 0x1.467406p-5f, -0x1.3e6128p-3f, -0x1.94f5e4p-6f, 0x1.a2fca4p-4f, 0x1.711d6cp-8f, -0x1.f9a24ep-5f, 0x1.8100cep-8f, 0x1.088e0ep-5f, -0x1.200234p-7f, -0x1.c48b18p-7f, 0x1.b88232p-8f, 0x1.21b464p-8f, -0x1.d56f02p-9f, -0x1.b3fa62p-11f, 0x1.6d0d18p-10f, -0x1.c0c538p-15f, -0x1.93cfc4p-12f, 0x1.a9dc1cp-14f, 0x1.1c224ap-14f, -0x1.37443cp-15f, -0x1.25ad94p-18f, 0x1.e5f3b2p-18f, -0x1.0fa7b2p-20f, -0x1.6f998cp-21f, 0x1.1ad0c2p-22f, 0x1.baf44p-33f, -0x1.37c9a6p-26f, 0x1.16bc24p-28f, -0x1.49b9bep-32f};
        }
        if constexpr (std::is_same_v<Real, double>) {
            return {0x1.98eb9a1ad223fp-11, 0x1.59aebc7410b77p-7, 0x1.03c8a76dc45bbp-4, 0x1.c27102c5bb1f8p-3, 0x1.e40a78036ac7fp-2, 0x1.389291df573b5p-1, 0x1.722da89a0e974p-2, -0x1.1d1b3a29bcd1cp-3, -0x1.4ea132f283bf2p-2, -0x1.120e7e4faaa26p-6, 0x1.d38a4229ed991p-3, 0x1.46740628c1e1ap-5, -0x1.3e612851515f6p-3, -0x1.94f5e34522642p-6, 0x1.a2fca44817f2fp-4, 0x1.711d6c0e1d46ap-8, -0x1.f9a24d69ac62ep-5, 0x1.8100ce41923f9p-8, 0x1.088e0e0a9eca2p-5, -0x1.2002338d10ca4p-7, -0x1.c48b170cc05f1p-7, 0x1.b88231706930bp-8, 0x1.21b464fba22e8p-8, -0x1.d56f016e81eefp-9, -0x1.b3fa62b1e588ap-11, 0x1.6d0d181744b82p-10, -0x1.c0c537c57bf91p-15, -0x1.93cfc4d1b467bp-12, 0x1.a9dc1c3780427p-14, 0x1.1c224959ff66dp-14, -0x1.37443b45c500dp-15, -0x1.25ad943ec34c5p-18, 0x1.e5f3b2c4841b3p-18, -0x1.0fa7b2e66c174p-20, -0x1.6f998ba053586p-21, 0x1.1ad0c20194a77p-22, 0x1.baf43fba18386p-33, -0x1.37c9a671235e5p-26, 0x1.16bc244bc36b8p-28, -0x1.49b9be3c4a795p-32};
        }
        if constexpr (std::is_same_v<Real, long double>) {
            return {0xc.c75cd0d6911fa1bp-14L, 0xa.cd75e3a085bb58ep-10L, 0x8.1e453b6e22dd40bp-7L, 0xe.1388162dd8fc33cp-6L, 0xf.2053c01b563f9cp-5L, 0x9.c4948efab9da4ebp-4L, 0xb.916d44d074ba33ap-5L, -0x8.e8d9d14de68dd92p-6L, -0xa.750997941df90a2p-5L, -0x8.9073f27d5512ec4p-9L, 0xe.9c52114f6cc8785p-6L, 0xa.33a031460f0d07bp-8L, -0x9.f309428a8afb31p-6L, -0xc.a7af1a291320cf8p-9L, 0xd.17e52240bf97653p-7L, 0xb.88eb6070ea34f45p-11L, -0xf.cd126b4d6317224p-8L, 0xc.0806720c91fc6d9p-11L, 0x8.44707054f65127ep-8L, -0x9.00119c6886523a3p-10L, -0xe.2458b86602f87fcp-10L, 0xd.c4118b834985633p-11L, 0x9.0da327dd11741b6p-11L, -0xe.ab780b740f77b4ep-12L, -0xd.9fd3158f2c44f42p-14L, 0xb.6868c0ba25c0c5ep-13L, -0xe.0629be2bdfc8a9ap-18L, -0xc.9e7e268da33d44p-15L, 0xd.4ee0e1bc0213577p-17L, 0x8.e1124acffb364cp-17L, -0x9.ba21da2e280684dp-18L, -0x9.2d6ca1f61a624d6p-21L, 0xf.2f9d962420d947dp-21L, -0x8.7d3d973360b9f2p-23L, -0xb.7ccc5d029ac2c2ep-24L, 0x8.d686100ca53b51dp-25L, 0xd.d7a1fdd0c1c3167p-36L, -0x9.be4d33891af2665p-29L, 0x8.b5e1225e1b5bf96p-31L, -0xa.4dcdf1e253caa43p-35L};
        }
        #ifdef BOOST_HAS_FLOAT128
        if constexpr (std::is_same_v<Real, boost::multiprecision::float128>) {
            return {0x1.98eb9a1ad223f435df1887363eb7p-11Q, 0x1.59aebc7410b76b1badcc12950a6dp-7Q, 0x1.03c8a76dc45ba8152978f10f7c69p-4Q, 0x1.c27102c5bb1f8677bf80bd502bd1p-3Q, 0x1.e40a78036ac7f37fedcb5a80e979p-2Q, 0x1.389291df573b49d69709e6a39a77p-1Q, 0x1.722da89a0e9746749e19dd59185cp-2Q, -0x1.1d1b3a29bcd1bb23d3a59b3ec213p-3Q, -0x1.4ea132f283bf214303a247980416p-2Q, -0x1.120e7e4faaa25d883a426ff76e5ep-6Q, 0x1.d38a4229ed990f09767be4223746p-3Q, 0x1.46740628c1e1a0f69af0510e97b1p-5Q, -0x1.3e612851515f661f1b982f038f42p-3Q, -0x1.94f5e345226419f0d083e13269c2p-6Q, 0x1.a2fca44817f2eca68758db1b8b3cp-4Q, 0x1.711d6c0e1d469e893974f1ea6999p-8Q, -0x1.f9a24d69ac62e448cc934c83e493p-5Q, 0x1.8100ce41923f8db158db362c64bfp-8Q, 0x1.088e0e0a9eca24fb78af66e8d177p-5Q, -0x1.2002338d10ca4746416546df60b9p-7Q, -0x1.c48b170cc05f0ff881bd3e1b39f6p-7Q, 0x1.b88231706930ac6660faf6c0f2dep-8Q, 0x1.21b464fba22e836bf597b5aa439ep-8Q, -0x1.d56f016e81eef69c49d53f5d0ebbp-9Q, -0x1.b3fa62b1e5889e83d8f6afb75bd8p-11Q, 0x1.6d0d181744b818bb6e246e9454cp-10Q, -0x1.c0c537c57bf91533c45c4dad234ap-15Q, -0x1.93cfc4d1b467a87fbd054b646597p-12Q, 0x1.a9dc1c3780426aeee785e4f91b69p-14Q, 0x1.1c224959ff66c97f21d37d9c05fbp-14Q, -0x1.37443b45c500d09a28a777a35885p-15Q, -0x1.25ad943ec34c49ac55547e1ee8bep-18Q, 0x1.e5f3b2c4841b28f9ade386d031fap-18Q, -0x1.0fa7b2e66c173e3f4204ea792843p-20Q, -0x1.6f998ba05358585b78628b55c242p-21Q, 0x1.1ad0c20194a76a3ad5100084b452p-22Q, 0x1.baf43fba183862ce0fdd2ce6713ep-33Q, -0x1.37c9a671235e4ccaa17336a7b289p-26Q, 0x1.16bc244bc36b7f2c806c15a56631p-28Q, -0x1.49b9be3c4a79548641570a58d68ep-32Q};
        }
        #endif
    }
    if constexpr (p == 21) {
        if constexpr (std::is_same_v<Real, float>) {
            return {0x1.1fbdb6p-11f, 0x1.fda65ap-8f, 0x1.93701p-5f, 0x1.736cacp-3f, 0x1.adc2aep-2f, 0x1.33f89cp-1f, 0x1.c742b8p-2f, -0x1.24a464p-5f, -0x1.57b854p-2f, -0x1.cc60ep-4f, 0x1.b148bep-3f, 0x1.d7feep-4f, -0x1.1e9916p-3f, -0x1.4ef44ap-4f, 0x1.8bacdap-4f, 0x1.7690eep-5f, -0x1.0a25dap-4f, -0x1.319ff4p-6f, 0x1.457136p-5f, 0x1.b81b9cp-9f, -0x1.564bap-6f, 0x1.3b071ap-9f, 0x1.268bbap-7f, -0x1.7af914p-9f, -0x1.83c29p-9f, 0x1.c1ff8cp-10f, 0x1.4f3d4cp-11f, -0x1.6a1c4ep-11f, -0x1.0c224ep-15f, 0x1.9624eap-13f, -0x1.30f834p-15f, -0x1.2592c2p-15f, 0x1.019c78p-16f, 0x1.768306p-19f, -0x1.9ebc2ep-19f, 0x1.53f4f8p-22f, 0x1.41473ap-22f, -0x1.adab3ep-24f, -0x1.35ca1ap-29f, 0x1.e5063ep-28f, -0x1.949b92p-30f, 0x1.c8df1cp-34f};
        }
        if constexpr (std::is_same_v<Real, double>) {
            return {0x1.1fbdb5d9588b2p-11, 0x1.fda65a6714bdap-8, 0x1.937010050905p-5, 0x1.736cace0a6657p-3, 0x1.adc2ad3828f8ep-2, 0x1.33f89b5c2e57bp-1, 0x1.c742b82c92555p-2, -0x1.24a4646f28f23p-5, -0x1.57b853bbe6a56p-2, -0x1.cc60df28d7badp-4, 0x1.b148be3b1073fp-3, 0x1.d7fedf0356bp-4, -0x1.1e9915e3ca004p-3, -0x1.4ef44a8b3ce94p-4, 0x1.8bacd9d38a28ep-4, 0x1.7690ee8b59819p-5, -0x1.0a25da27e87d7p-4, -0x1.319ff4cde80a9p-6, 0x1.4571365df5489p-5, 0x1.b81b9bc4d4f3cp-9, -0x1.564b9fd2c3bbfp-6, 0x1.3b071a00e4e5cp-9, 0x1.268bb96001dc4p-7, -0x1.7af914f01f8fcp-9, -0x1.83c28f7e3c0e5p-9, 0x1.c1ff8c65fb698p-10, 0x1.4f3d4c629fe3fp-11, -0x1.6a1c4ee19b423p-11, -0x1.0c224d58d6e81p-15, 0x1.9624eae0b53fp-13, -0x1.30f833f09ef02p-15, -0x1.2592c13818da8p-15, 0x1.019c78bb61865p-16, 0x1.768306fdacad5p-19, -0x1.9ebc2e6b68816p-19, 0x1.53f4f805c9123p-22, 0x1.414739b599eb7p-22, -0x1.adab3d69a3265p-24, -0x1.35ca1a0003737p-29, 0x1.e5063d39ebf9ap-28, -0x1.949b918111d77p-30, 0x1.c8df1c0c683efp-34};
        }
        if constexpr (std::is_same_v<Real, long double>) {
            return {0x8.fdedaecac458cd2p-14L, 0xf.ed32d338a5ecf5cp-11L, 0xc.9b8080284827e72p-8L, 0xb.9b656705332bae6p-6L, 0xd.6e1569c147c6dadp-5L, 0x9.9fc4dae172bdbe9p-4L, 0xe.3a15c16492aa9dbp-5L, -0x9.252323794791b92p-8L, -0xa.bdc29ddf352b036p-5L, -0xe.6306f946bdd67aep-7L, 0xd.8a45f1d8839f921p-6L, 0xe.bff6f81ab580014p-7L, -0x8.f4c8af1e5002333p-6L, -0xa.77a25459e74a242p-7L, 0xc.5d66ce9c5146f15p-7L, 0xb.b487745acc0c8adp-8L, -0x8.512ed13f43ebad1p-7L, -0x9.8cffa66f405497bp-9L, 0xa.2b89b2efaa444bdp-8L, 0xd.c0dcde26a79e22p-12L, -0xa.b25cfe961ddf733p-9L, 0x9.d838d007272dcc4p-12L, 0x9.345dcb000ee2014p-10L, -0xb.d7c8a780fc7dcf5p-12L, -0xc.1e147bf1e072bc1p-12L, 0xe.0ffc632fdb4c29ep-13L, 0xa.79ea6314ff1f508p-14L, -0xb.50e2770cda119f3p-14L, -0x8.61126ac6b7407d4p-18L, 0xc.b1275705a9f8219p-16L, -0x9.87c19f84f78120fp-18L, -0x9.2c9609c0c6d43f4p-18L, 0x8.0ce3c5db0c32b2dp-19L, 0xb.b41837ed656a6e4p-22L, -0xc.f5e1735b440ad22p-22L, 0xa.9fa7c02e489184ap-25L, 0xa.0a39cdaccf5ba4ap-25L, -0xd.6d59eb4d1932a9ep-27L, -0x9.ae50d0001b9b402p-32L, 0xf.2831e9cf5fcccbdp-31L, -0xc.a4dc8c088ebb653p-33L, 0xe.46f8e06341f76fbp-37L};
        }
        #ifdef BOOST_HAS_FLOAT128
        if constexpr (std::is_same_v<Real, boost::multiprecision::float128>) {
            return {0x1.1fbdb5d9588b19a3bff14d1b38fep-11Q, 0x1.fda65a6714bd9eb85e8d7e88dc47p-8Q, 0x1.937010050904fce3eac118cf708bp-5Q, 0x1.736cace0a66575cc5580323567efp-3Q, 0x1.adc2ad3828f8db5a1d5b82619893p-2Q, 0x1.33f89b5c2e57b7d298074c990addp-1Q, 0x1.c742b82c925553b5685081f7bfa5p-2Q, -0x1.24a4646f28f23724c7651fdfae24p-5Q, -0x1.57b853bbe6a5606b16b5f4878728p-2Q, -0x1.cc60df28d7bacf5cd487dfbe8b15p-4Q, 0x1.b148be3b1073f24190a78355f27p-3Q, 0x1.d7fedf0356b00027d18f8e92c824p-4Q, -0x1.1e9915e3ca004666a067ba487474p-3Q, -0x1.4ef44a8b3ce944830f2730d25eb9p-4Q, 0x1.8bacd9d38a28de29668579705c28p-4Q, 0x1.7690ee8b5981915a6dd616884c3fp-5Q, -0x1.0a25da27e87d75a179ff300fb83fp-4Q, -0x1.319ff4cde80a92f66d084237880cp-6Q, 0x1.4571365df548897971310c0822e6p-5Q, 0x1.b81b9bc4d4f3c43f2ed25e559ae8p-9Q, -0x1.564b9fd2c3bbee66fedad254cc1dp-6Q, 0x1.3b071a00e4e5b987ae06dc453687p-9Q, 0x1.268bb96001dc402712698f3f0494p-7Q, -0x1.7af914f01f8fb9ea842adaa7f3cdp-9Q, -0x1.83c28f7e3c0e5781b7d022113abfp-9Q, 0x1.c1ff8c65fb69853c310221a8efcfp-10Q, 0x1.4f3d4c629fe3ea0f5f4b5c45e52p-11Q, -0x1.6a1c4ee19b4233e5be9da5626ae4p-11Q, -0x1.0c224d58d6e80fa7cb5e318a3bafp-15Q, 0x1.9624eae0b53f0431ec2fc4ac4a39p-13Q, -0x1.30f833f09ef0241df6514d4f75d1p-15Q, -0x1.2592c13818da87e8e3f7d01fb451p-15Q, 0x1.019c78bb618656595fd443f75a4ep-16Q, 0x1.768306fdacad4dc8fab15811b58ap-19Q, -0x1.9ebc2e6b68815a4433974307ca8ep-19Q, 0x1.53f4f805c91230935990691c9b37p-22Q, 0x1.414739b599eb74930d1d52f4adeap-22Q, -0x1.adab3d69a326553c517a3161ce4cp-24Q, -0x1.35ca1a000373680368e6e76ba1cap-29Q, 0x1.e5063d39ebf99979a7db1201c73dp-28Q, -0x1.949b918111d76ca513fa5f9592d7p-30Q, 0x1.c8df1c0c683eedf6229624f6ad7fp-34Q};
        }
        #endif
    }
    if constexpr (p == 22) {
        if constexpr (std::is_same_v<Real, float>) {
            return {0x1.9506cp-12f, 0x1.76fccap-8f, 0x1.37de72p-5f, 0x1.2fdb52p-3f, 0x1.788ddep-2f, 0x1.282856p-1f, 0x1.040b9cp-1f, 0x1.2df9bep-4f, -0x1.403b66p-2f, -0x1.9ac39cp-3f, 0x1.501016p-3f, 0x1.7095c8p-3f, -0x1.8dc40ep-4f, -0x1.0ddc74p-3f, 0x1.16d73p-4f, 0x1.5a58d6p-4f, -0x1.a4c6a6p-5f, -0x1.7d2e3p-5f, 0x1.2edd7cp-5f, 0x1.514ae8p-6f, -0x1.80b244p-6f, -0x1.9739fap-8f, 0x1.9bb88ep-7f, 0x1.3ab784p-12f, -0x1.658b52p-8f, 0x1.11bf26p-10f, 0x1.def09ap-10f, -0x1.941064p-11f, -0x1.bc5f9p-12f, 0x1.58926ap-12f, 0x1.6c8f7ep-15f, -0x1.8a7bd2p-14f, 0x1.7da8ecp-17f, 0x1.237bb8p-16f, -0x1.9dd79ap-18f, -0x1.a4264ap-20f, 0x1.5bac3cp-20f, -0x1.7917ccp-24f, -0x1.139828p-23f, 0x1.43164cp-25f, 0x1.cdd78ep-30f, -0x1.772814p-29f, 0x1.2558aap-31f, -0x1.3cd862p-35f};
        }
        if constexpr (std::is_same_v<Real, double>) {
            return {0x1.9506bf26823bdp-12, 0x1.76fcca8398a61p-8, 0x1.37de71ec845e4p-5, 0x1.2fdb523f49ab5p-3, 0x1.788dde34ceb66p-2, 0x1.282855bd3e661p-1, 0x1.040b9c9ba65ccp-1, 0x1.2df9be17f6307p-4, -0x1.403b658e1e271p-2, -0x1.9ac39bc6a693fp-3, 0x1.501016e351cb2p-3, 0x1.7095c864b0cccp-3, -0x1.8dc40d7416ec6p-4, -0x1.0ddc7410272ddp-3, 0x1.16d7307f2f199p-4, 0x1.5a58d5e218e8bp-4, -0x1.a4c6a60c82a99p-5, -0x1.7d2e2f61c66bdp-5, 0x1.2edd7c2485c6cp-5, 0x1.514ae8ff8882dp-6, -0x1.80b243789dd91p-6, -0x1.9739fa1f59fdbp-8, 0x1.9bb88df2b07f2p-7, 0x1.3ab784f96ef95p-12, -0x1.658b52a857f98p-8, 0x1.11bf26e4a8f4ap-10, 0x1.def099502d06ep-10, -0x1.941064b6e88a6p-11, -0x1.bc5f905d400afp-12, 0x1.58926afecb291p-12, 0x1.6c8f7d64132f3p-15, -0x1.8a7bd1deebb29p-14, 0x1.7da8eb6c7dd39p-17, 0x1.237bb7cb019efp-16, -0x1.9dd79a5aef3dp-18, -0x1.a4264a7b9f195p-20, 0x1.5bac3bb6b7df3p-20, -0x1.7917ccd42a174p-24, -0x1.139827c14728ep-23, 0x1.43164b7e12198p-25, 0x1.cdd78d40a514dp-30, -0x1.772813927397bp-29, 0x1.2558a9abd0467p-31, -0x1.3cd8627d9b085p-35};
        }
        if constexpr (std::is_same_v<Real, long double>) {
            return {0xc.a835f93411de553p-15L, 0xb.b7e6541cc530467p-11L, 0x9.bef38f6422f21e9p-8L, 0x9.7eda91fa4d5aa1ap-6L, 0xb.c46ef1a675b2f1p-5L, 0x9.4142ade9f3308e8p-4L, 0x8.205ce4dd32e5de6p-4L, 0x9.6fcdf0bfb18389dp-7L, -0xa.01db2c70f1385a8p-5L, -0xc.d61cde35349f534p-6L, 0xa.8080b71a8e58ef5p-6L, 0xb.84ae43258666006p-6L, -0xc.6e206ba0b762e61p-7L, -0x8.6ee3a081396e7a5p-6L, 0x8.b6b983f978cc515p-7L, 0xa.d2c6af10c745924p-7L, -0xd.26353064154caf2p-8L, -0xb.e9717b0e335e71fp-8L, 0x9.76ebe1242e35de2p-8L, 0xa.8a5747fc44168f9p-9L, -0xc.05921bc4eec8a0fp-9L, -0xc.b9cfd0facfed6d3p-11L, 0xc.ddc46f9583f92e4p-10L, 0x9.d5bc27cb77ca414p-15L, -0xb.2c5a9542bfcbcecp-11L, 0x8.8df9372547a50bfp-13L, 0xe.f784ca816836eb8p-13L, -0xc.a08325b74452d5dp-14L, -0xd.e2fc82ea0057a9fp-15L, 0xa.c49357f65948be4p-15L, 0xb.647beb209979943p-18L, -0xc.53de8ef75d94bb9p-17L, 0xb.ed475b63ee9c598p-20L, 0x9.1bddbe580cf7997p-19L, -0xc.eebcd2d779e82c2p-21L, -0xd.213253dcf8ca54ep-23L, 0xa.dd61ddb5bef98f1p-23L, -0xb.c8be66a150b9e14p-27L, -0x8.9cc13e0a3946c9cp-26L, 0xa.18b25bf090cbe47p-28L, 0xe.6ebc6a0528a6ab3p-33L, -0xb.b9409c939cbdb9ep-32L, 0x9.2ac54d5e82338aap-34L, -0x9.e6c313ecd8427c2p-38L};
        }
        #ifdef BOOST_HAS_FLOAT128
        if constexpr (std::is_same_v<Real, boost::multiprecision::float128>) {
            return {0x1.9506bf26823bcaa678f24755b1bep-12Q, 0x1.76fcca8398a608cea25db5f87af6p-8Q, 0x1.37de71ec845e43d133d0c299f8ecp-5Q, 0x1.2fdb523f49ab54335a2eba9c8437p-3Q, 0x1.788dde34ceb65e20730d6ddac0aep-2Q, 0x1.282855bd3e6611cfef751a206622p-1Q, 0x1.040b9c9ba65cbbcc9c304a8323e2p-1Q, 0x1.2df9be17f630713a1ee1adb275ebp-4Q, -0x1.403b658e1e270b5094f497dac5a2p-2Q, -0x1.9ac39bc6a693ea6805616ea07bd7p-3Q, 0x1.501016e351cb1dea06a54dd049bp-3Q, 0x1.7095c864b0ccc00c06f07ad7df02p-3Q, -0x1.8dc40d7416ec5cc24d97996382b6p-4Q, -0x1.0ddc7410272dcf4ab45934d4cf13p-3Q, 0x1.16d7307f2f198a2aec11f34d1e9fp-4Q, 0x1.5a58d5e218e8b247e862b3c1695fp-4Q, -0x1.a4c6a60c82a995e3152af9622fcap-5Q, -0x1.7d2e2f61c66bce3ea7e2b78a28c2p-5Q, 0x1.2edd7c2485c6bbc4d1a2c115a8a7p-5Q, 0x1.514ae8ff8882d1f2c15a47cd6559p-6Q, -0x1.80b243789dd9141e3fff7c22096dp-6Q, -0x1.9739fa1f59fdada5471944d8c62dp-8Q, 0x1.9bb88df2b07f25c8944100438d61p-7Q, 0x1.3ab784f96ef948279489735c5a1cp-12Q, -0x1.658b52a857f979d8d4ec9b9f73a9p-8Q, 0x1.11bf26e4a8f4a17dea14a455c454p-10Q, 0x1.def099502d06dd70a82c80e4bd5ep-10Q, -0x1.941064b6e88a5abad07f7589e66dp-11Q, -0x1.bc5f905d400af53d3eed609ad591p-12Q, 0x1.58926afecb2917c8a7cf6285c131p-12Q, 0x1.6c8f7d64132f3286043bbcbfb14ep-15Q, -0x1.8a7bd1deebb29772b3b880563654p-14Q, 0x1.7da8eb6c7dd38b3094b5e1389c85p-17Q, 0x1.237bb7cb019ef32d902afce7d618p-16Q, -0x1.9dd79a5aef3d0583abb88bd3c2cp-18Q, -0x1.a4264a7b9f194a9cadbc3e265f85p-20Q, 0x1.5bac3bb6b7df31e28d50f39c71ffp-20Q, -0x1.7917ccd42a173c285afd81440ef2p-24Q, -0x1.139827c14728d9384a9bd9717ddp-23Q, 0x1.43164b7e12197c8ececb2f261f6bp-25Q, 0x1.cdd78d40a514d5664ee974c433ep-30Q, -0x1.772813927397b73b1b90436f0da9p-29Q, 0x1.2558a9abd04671533bb0e6556b43p-31Q, -0x1.3cd8627d9b084f8342fbe01d8e11p-35Q};
        }
        #endif
    }
    if constexpr (p == 23) {
        if constexpr (std::is_same_v<Real, float>) {
            return {0x1.1d1cbap-12f, 0x1.136e6ep-8f, 0x1.e0371p-6f, 0x1.eda1b2p-4f, 0x1.4617f8p-2f, 0x1.170138p-1f, 0x1.1a1f1ap-1f, 0x1.737dfap-3f, -0x1.0baa62p-2f, -0x1.15ea6ep-2f, 0x1.79587ep-4f, 0x1.c9e0fcp-3f, -0x1.0ea48cp-5f, -0x1.4fe52ap-3f, 0x1.4c5162p-6f, 0x1.cbf7fcp-4f, -0x1.5a21c2p-6f, -0x1.1f91cap-4f, 0x1.649c9ep-6f, 0x1.3b5a9p-5f, -0x1.2f7d4p-6f, -0x1.1f53eep-6f, 0x1.a1db1p-7f, 0x1.8b4d7ep-8f, -0x1.cfb028p-8f, -0x1.297f88p-10f, 0x1.995258p-9f, -0x1.0279b8p-12f, -0x1.163206p-10f, 0x1.4eefcap-12f, 0x1.0d3c26p-12f, -0x1.3a9e5ep-13f, -0x1.1b7136p-15f, 0x1.734928p-15f, -0x1.61b10ep-19f, -0x1.181bb4p-17f, 0x1.41cbdcp-19f, 0x1.b56b6ap-21f, -0x1.1ea2bcp-21f, 0x1.3e5bd6p-26f, 0x1.d15d2ap-25f, -0x1.e103ap-27f, -0x1.046382p-30f, 0x1.20be9ap-30f, -0x1.a8f0bcp-33f, 0x1.b7e04ep-37f};
        }
        if constexpr (std::is_same_v<Real, double>) {
            return {0x1.1d1cb9f15ca3p-12, 0x1.136e6d0cada7ap-8, 0x1.e03710ca000b2p-6, 0x1.eda1b1bac15d3p-4, 0x1.4617f85b87d2cp-2, 0x1.1701372560847p-1, 0x1.1a1f195e16c2fp-1, 0x1.737dfa0d62dcfp-3, -0x1.0baa621ce4279p-2, -0x1.15ea6e8658872p-2, 0x1.79587da8cb753p-4, 0x1.c9e0fc732c538p-3, -0x1.0ea48c5a06944p-5, -0x1.4fe52af139a58p-3, 0x1.4c51617dbe58bp-6, 0x1.cbf7fc1c6fa45p-4, -0x1.5a21c1ca225bep-6, -0x1.1f91c92c06d32p-4, 0x1.649c9ed18296ap-6, 0x1.3b5a90477563fp-5, -0x1.2f7d3f55674e5p-6, -0x1.1f53eed1947bfp-6, 0x1.a1db0f12bb42p-7, 0x1.8b4d7e5397f08p-8, -0x1.cfb028e3aae51p-8, -0x1.297f885faab6dp-10, 0x1.995258715da31p-9, -0x1.0279b74217beap-12, -0x1.163205afadf6dp-10, 0x1.4eefc91e26cfap-12, 0x1.0d3c2558f1652p-12, -0x1.3a9e5e1c22af4p-13, -0x1.1b7136cbb5aep-15, 0x1.734927c4130bfp-15, -0x1.61b10dba0880bp-19, -0x1.181bb48217362p-17, 0x1.41cbdcc2df4d9p-19, 0x1.b56b6a75fa292p-21, -0x1.1ea2bb46afa08p-21, 0x1.3e5bd6b02b53fp-26, 0x1.d15d2a9531ecbp-25, -0x1.e103a0e0c3b96p-27, -0x1.0463822ac124p-30, 0x1.20be99079e759p-30, -0x1.a8f0bc9b067e6p-33, 0x1.b7e04d72be8f7p-37};
        }
        if constexpr (std::is_same_v<Real, long double>) {
            return {0x8.e8e5cf8ae517d2ep-15L, 0x8.9b7368656d3d159p-11L, 0xf.01b886500058d29p-9L, 0xf.6d0d8dd60ae956bp-7L, 0xa.30bfc2dc3e9620cp-5L, 0x8.b809b92b0423412p-4L, 0x8.d0f8caf0b617657p-4L, 0xb.9befd06b16e7779p-6L, -0x8.5d5310e7213cb45p-5L, -0x8.af537432c439318p-5L, 0xb.cac3ed465ba98efp-7L, 0xe.4f07e399629c099p-6L, -0x8.752462d034a1f39p-8L, -0xa.7f295789cd2c17ap-6L, 0xa.628b0bedf2c55d4p-9L, 0xe.5fbfe0e37d2292p-7L, -0xa.d10e0e5112df11ap-9L, -0x8.fc8e49603699185p-7L, 0xb.24e4f68c14b4ebbp-9L, 0x9.dad4823bab1fbp-8L, -0x9.7be9faab3a727aap-9L, -0x8.fa9f768ca3df844p-9L, 0xd.0ed87895da0fea4p-10L, 0xc.5a6bf29cbf83d82p-11L, -0xe.7d81471d5728bb2p-11L, -0x9.4bfc42fd55b69dep-13L, 0xc.ca92c38aed18631p-12L, -0x8.13cdba10bdf515p-15L, -0x8.b1902d7d6fb67f6p-13L, 0xa.777e48f1367d0a5p-15L, 0x8.69e12ac78b28dbcp-15L, -0x9.d4f2f0e11579f3cp-16L, -0x8.db89b65dad6fc16p-18L, 0xb.9a493e20985f9bfp-18L, -0xb.0d886dd044054a4p-22L, -0x8.c0dda410b9b125ap-20L, 0xa.0e5ee616fa6c626p-22L, 0xd.ab5b53afd14926ap-24L, -0x8.f515da357d03ea3p-24L, 0x9.f2deb5815a9f567p-29L, 0xe.8ae954a98f65bb4p-28L, -0xf.081d07061dcaca2p-30L, -0x8.231c1156091fc3ap-33L, 0x9.05f4c83cf3acb2cp-33L, -0xd.4785e4d833f33e2p-36L, 0xd.bf026b95f47b713p-40L};
        }
        #ifdef BOOST_HAS_FLOAT128
        if constexpr (std::is_same_v<Real, boost::multiprecision::float128>) {
            return {0x1.1d1cb9f15ca2fa5bd9123de7bf02p-12Q, 0x1.136e6d0cada7a2b148b6e8855daap-8Q, 0x1.e03710ca000b1a51a7d181d489bp-6Q, 0x1.eda1b1bac15d2ad6ddbdb3147e36p-4Q, 0x1.4617f85b87d2c4188e472323232cp-2Q, 0x1.1701372560846824bce2d06e180cp-1Q, 0x1.1a1f195e16c2ecad981d672e54f9p-1Q, 0x1.737dfa0d62dceef2a6ae0dd9fae6p-3Q, -0x1.0baa621ce427968960a0c4171c1cp-2Q, -0x1.15ea6e865887262fe21b5ed288ddp-2Q, 0x1.79587da8cb7531de57b391307b2ep-4Q, 0x1.c9e0fc732c5381328729c7899645p-3Q, -0x1.0ea48c5a06943e72865d8b0a9865p-5Q, -0x1.4fe52af139a582f4b4bfddd5990ap-3Q, 0x1.4c51617dbe58aba744aaa9e2e84dp-6Q, 0x1.cbf7fc1c6fa4523fc35d6c595afp-4Q, -0x1.5a21c1ca225be23388c3683c5f1ap-6Q, -0x1.1f91c92c06d3230ad5ace9dd7e2ap-4Q, 0x1.649c9ed182969d7683dcc682acbcp-6Q, 0x1.3b5a90477563f5ff5ff6e6498f97p-5Q, -0x1.2f7d3f55674e4f5496795e1de4a5p-6Q, -0x1.1f53eed1947bf0871dd1fee02046p-6Q, 0x1.a1db0f12bb41fd473a2e8421b507p-7Q, 0x1.8b4d7e5397f07b04d41a59615574p-8Q, -0x1.cfb028e3aae51764ab65b4e360fbp-8Q, -0x1.297f885faab6d3bcc655cccba49p-10Q, 0x1.995258715da30c618d8b26e9df7bp-9Q, -0x1.0279b74217bea2a0847d1e93003bp-12Q, -0x1.163205afadf6cfebfb69e0b0f9f9p-10Q, 0x1.4eefc91e26cfa14abc8a3c8c31bep-12Q, 0x1.0d3c2558f1651b77d56d78be87a6p-12Q, -0x1.3a9e5e1c22af3e7828bfcf5d2f48p-13Q, -0x1.1b7136cbb5adf82b0d9aad06bd7cp-15Q, 0x1.734927c4130bf37e18cb5a04041ap-15Q, -0x1.61b10dba0880a94866a6f32c90e9p-19Q, -0x1.181bb482173624b3f677525a8993p-17Q, 0x1.41cbdcc2df4d8c4b6d915e9702b2p-19Q, 0x1.b56b6a75fa2924d359ae84bf18d8p-21Q, -0x1.1ea2bb46afa07d45b08cabbea3b6p-21Q, 0x1.3e5bd6b02b53eacdcf10b59b874fp-26Q, 0x1.d15d2a9531ecb7682cdf70f5d33bp-25Q, -0x1.e103a0e0c3b959430c8cac42a206p-27Q, -0x1.0463822ac123f874c86ec02d6961p-30Q, 0x1.20be99079e7596586119e6bc965cp-30Q, -0x1.a8f0bc9b067e67c3b0dec15337bdp-33Q, 0x1.b7e04d72be8f6e258b11c231fbfep-37Q};
        }
        #endif
    }
    if constexpr (p == 24) {
        if constexpr (std::is_same_v<Real, float>) {
            return {0x1.91785p-13f, 0x1.93f98p-9f, 0x1.7059c4p-6f, 0x1.8e62d8p-4f, 0x1.17757p-2f, 0x1.023cecp-1f, 0x1.265e6ep-1f, 0x1.1fbaaep-2f, -0x1.7f8826p-3f, -0x1.4592dep-2f, 0x1.390a48p-8f, 0x1.e9f54ap-3f, 0x1.5c6536p-5f, -0x1.5e912ep-3f, -0x1.3da9ap-5f, 0x1.efaecap-4f, 0x1.57bcfap-6f, -0x1.5088bep-4f, -0x1.2c0d6ap-8f, 0x1.a4434cp-5f, -0x1.440e76p-8f, -0x1.ce3e5ap-6f, 0x1.f61e5cp-8f, 0x1.ab9f18p-7f, -0x1.9c50c6p-8f, -0x1.371236p-8f, 0x1.e9b0e8p-9f, 0x1.2e73dap-10f, -0x1.bcb74ap-10f, -0x1.7274d8p-15f, 0x1.334ca4p-11f, -0x1.ef71f4p-14f, -0x1.32337p-13f, 0x1.131ee6p-14f, 0x1.6e4982p-16f, -0x1.536264p-16f, 0x1.ccd17cp-27f, 0x1.05cc64p-18f, -0x1.e21fa8p-21f, -0x1.b0fcbap-22f, 0x1.d137c8p-23f, -0x1.160c08p-31f, -0x1.838882p-26f, 0x1.627092p-28f, 0x1.050b72p-31f, -0x1.ba8412p-32f, 0x1.3380aep-34f, -0x1.31989ep-38f};
        }
        if constexpr (std::is_same_v<Real, double>) {
            return {0x1.91785023a9879p-13, 0x1.93f9805696b9cp-9, 0x1.7059c4b494de9p-6, 0x1.8e62d8a0a8ffcp-4, 0x1.17756f5530d69p-2, 0x1.023cebcae9f67p-1, 0x1.265e6eff960ccp-1, 0x1.1fbaad4726276p-2, -0x1.7f8826c061991p-3, -0x1.4592ded92acf9p-2, 0x1.390a478f54bbbp-8, 0x1.e9f54ac655c14p-3, 0x1.5c6535de6b4e5p-5, -0x1.5e912ec022b6p-3, -0x1.3da9a08b26c34p-5, 0x1.efaeca9aee192p-4, 0x1.57bcf98fd11cep-6, -0x1.5088bd1cf023ep-4, -0x1.2c0d69ece949fp-8, 0x1.a4434b89b6915p-5, -0x1.440e7548329f5p-8, -0x1.ce3e5912a2303p-6, 0x1.f61e5ce9a7ff1p-8, 0x1.ab9f170d72cc6p-7, -0x1.9c50c528bc778p-8, -0x1.371235b499d34p-8, 0x1.e9b0e84ddc047p-9, 0x1.2e73daaf19dbbp-10, -0x1.bcb74a39ee4ap-10, -0x1.7274d8ffa9744p-15, 0x1.334ca46898daep-11, -0x1.ef71f401aeeebp-14, -0x1.32337036c1617p-13, 0x1.131ee63336347p-14, 0x1.6e4981b862b8bp-16, -0x1.536263fe83a12p-16, 0x1.ccd17bbbe3d04p-27, 0x1.05cc64bbcf96p-18, -0x1.e21fa8ce8bf7ap-21, -0x1.b0fcbaac8055ep-22, 0x1.d137c82c35ed7p-23, -0x1.160c0805b95d3p-31, -0x1.838882edd8948p-26, 0x1.627091a0c660dp-28, 0x1.050b72dfcacecp-31, -0x1.ba8411e011af6p-32, 0x1.3380aec6c4cc7p-34, -0x1.31989d5a96be3p-38};
        }
        if constexpr (std::is_same_v<Real, long double>) {
            return {0xc.8bc2811d4c3c975p-16L, 0xc.9fcc02b4b5cdd4fp-12L, 0xb.82ce25a4a6f4946p-9L, 0xc.7316c50547fde96p-7L, 0x8.bbab7aa986b4822p-5L, 0x8.11e75e574fb369fp-4L, 0x9.32f377fcb066266p-4L, 0x8.fdd56a39313b02ap-5L, -0xb.fc4136030cc8af3p-6L, -0xa.2c96f6c9567c64ep-5L, 0x9.c8523c7aa5dd6c6p-11L, 0xf.4faa5632ae09cacp-6L, 0xa.e329aef35a72807p-8L, -0xa.f489760115afe5ap-6L, -0x9.ed4d0459361a0cbp-8L, 0xf.7d7654d770c90fap-7L, 0xa.bde7cc7e88e7341p-9L, -0xa.8445e8e7811ed17p-7L, -0x9.606b4f674a4f729p-11L, 0xd.221a5c4db48a85p-8L, -0xa.2073aa4194fa9f4p-11L, -0xe.71f2c8951181b37p-9L, 0xf.b0f2e74d3ff8483p-11L, 0xd.5cf8b86b966301bp-10L, -0xc.e2862945e3bc0bdp-11L, -0x9.b891ada4ce99f42p-11L, 0xf.4d87426ee023471p-12L, 0x9.739ed578cedd928p-13L, -0xd.e5ba51cf724ff04p-13L, -0xb.93a6c7fd4ba228cp-18L, 0x9.9a652344c6d708fp-14L, -0xf.7b8fa00d7775af1p-17L, -0x9.919b81b60b0b6ddp-16L, 0x8.98f73199b1a34dep-17L, 0xb.724c0dc315c5bebp-19L, -0xa.9b131ff41d09255p-19L, 0xe.668bdddf1e8204p-30L, 0x8.2e6325de7cafddbp-21L, -0xf.10fd46745fbd3a3p-24L, -0xd.87e5d56402aecbp-25L, 0xe.89be4161af6badap-26L, -0x8.b060402dcae9786p-34L, -0xc.1c44176ec4a3d56p-29L, 0xb.13848d0633067f9p-31L, 0x8.285b96fe5675d57p-34L, -0xd.d4208f008d7b156p-35L, 0x9.9c0576362663aafp-37L, -0x9.8cc4ead4b5f196ap-41L};
        }
        #ifdef BOOST_HAS_FLOAT128
        if constexpr (std::is_same_v<Real, boost::multiprecision::float128>) {
            return {0x1.91785023a98792eac41b5272c04bp-13Q, 0x1.93f9805696b9ba9d87c3df41a8d9p-9Q, 0x1.7059c4b494de928c16cf6cacf0a7p-6Q, 0x1.8e62d8a0a8ffbd2cc57d79f83e08p-4Q, 0x1.17756f5530d69044d0184f790b3ep-2Q, 0x1.023cebcae9f66d3d390049f78dcp-1Q, 0x1.265e6eff960cc4ccf24089549b66p-1Q, 0x1.1fbaad4726276054f5d700d1bfb6p-2Q, -0x1.7f8826c0619915e6978c42a74dfcp-3Q, -0x1.4592ded92acf8c9c5ecd6ba9faddp-2Q, 0x1.390a478f54bbad8c711b3d202a1ap-8Q, 0x1.e9f54ac655c13958344e31a88778p-3Q, 0x1.5c6535de6b4e500e1e067af4d00ep-5Q, -0x1.5e912ec022b5fcb42ec3601ecf9bp-3Q, -0x1.3da9a08b26c34196db8bd596e813p-5Q, 0x1.efaeca9aee1921f4bc36a577b638p-4Q, 0x1.57bcf98fd11ce6827e350101b72ep-6Q, -0x1.5088bd1cf023da2d878ac9ebafc5p-4Q, -0x1.2c0d69ece949ee5284dcfb042baap-8Q, 0x1.a4434b89b691509ff6541a080cecp-5Q, -0x1.440e7548329f53e7873db8128a63p-8Q, -0x1.ce3e5912a230366d14a862ea3672p-6Q, 0x1.f61e5ce9a7ff0906541711dfe2e7p-8Q, 0x1.ab9f170d72cc6036a5092b63c08ap-7Q, -0x1.9c50c528bc77817ae0788ffc773ap-8Q, -0x1.371235b499d33e83c0b1b543f6e8p-8Q, 0x1.e9b0e84ddc0468e17e4fa3636bf9p-9Q, 0x1.2e73daaf19dbb24f6ae6d9d95255p-10Q, -0x1.bcb74a39ee49fe08e2b9fd7aeef1p-10Q, -0x1.7274d8ffa9744517f2f1470cebeep-15Q, 0x1.334ca46898dae11da29a57e9f92ap-11Q, -0x1.ef71f401aeeeb5e13dd40efc0a0cp-14Q, -0x1.32337036c1616db93444903d828ep-13Q, 0x1.131ee633363469bb9bebca29a839p-14Q, 0x1.6e4981b862b8b7d5bb86d25cc3eap-16Q, -0x1.536263fe83a124aaa2092f47c711p-16Q, 0x1.ccd17bbbe3d0407ffce581eb4165p-27Q, 0x1.05cc64bbcf95fbb57637394fb838p-18Q, -0x1.e21fa8ce8bf7a7465a57d37cce9ap-21Q, -0x1.b0fcbaac8055d95fcb97dc94d8aap-22Q, 0x1.d137c82c35ed75b3cca32ab07fp-23Q, -0x1.160c0805b95d2f0cf1ffb7605bf2p-31Q, -0x1.838882edd8947aab05cd00305346p-26Q, 0x1.627091a0c660cff16c28fc227cd6p-28Q, 0x1.050b72dfcacebaae0203f358d56p-31Q, -0x1.ba8411e011af62ac1faf9b23aa51p-32Q, 0x1.3380aec6c4cc755d395c84cee933p-34Q, -0x1.31989d5a96be32d359a4b1032025p-38Q};
        }
        #endif

    }
    std::array<Real, 2*p> m{};
    for (auto & x : m) {
        x = std::numeric_limits<Real>::quiet_NaN();
    }
    return m;

}

template<class Real, size_t p>
std::array<Real, 2*p> daubechies_wavelet_filter() {
    std::array<Real, 2*p> g;
    auto h = daubechies_scaling_filter<Real, p>();
    for (size_t i = 0; i < g.size(); i += 2)
    {
        g[i] = h[g.size() - i - 1];
        g[i+1] = -h[g.size() - i - 2];
    }
    return g;
}

} // namespaces
#endif
