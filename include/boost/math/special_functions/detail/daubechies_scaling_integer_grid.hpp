/*
 * Copyright Nick Thompson, 2019
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef BOOST_MATH_DAUBECHIES_SCALING_INTEGER_GRID_HPP
#define BOOST_MATH_DAUBECHIES_SCALING_INTEGER_GRID_HPP
#include <array>
#ifdef BOOST_HAS_FLOAT128
#include <boost/multiprecision/float128.hpp>
#endif
namespace boost::math::detail {

template <typename Real, unsigned p, unsigned order>
constexpr std::array<Real, 2*p> daubechies_scaling_integer_grid()
{
    static_assert(sizeof(Real) <= 16, "Integer grids only computed up to 128 bits of precision.");
    static_assert(p < 25, "Integer grids only implemented up to 24.");
    if constexpr (p == 2) {
        if constexpr (std::is_same_v<Real, float>) {
            if constexpr (order == 0) {
                return {0x0p+0f, -0x1.76cf5ep-2f, 0x1.5db3d8p+0f, -0x0p+0f};
            }
            if constexpr (order == 1) {
                return {0x0p+0f, 0x1p+0f, -0x1p+0f, -0x0p+0f};
            }
         }
        if constexpr (std::is_same_v<Real, double>) {
            if constexpr (order == 0) {
                return {0x0p+0, -0x1.76cf5d0b09955p-2, 0x1.5db3d742c2655p+0, -0x0p+0};
            }
            if constexpr (order == 1) {
                return {0x0p+0, 0x1p+0, -0x1p+0, -0x0p+0};
            }
         }
        if constexpr (std::is_same_v<Real, long double>) {
            if constexpr (order == 0) {
                return {0x0p+0L, -0xb.b67ae8584caa73bp-5L, 0xa.ed9eba16132a9cfp-3L, -0x0p+0L};
            }
            if constexpr (order == 1) {
                return {0x0p+0L, 0x8p-3L, -0x8p-3L, -0x0p+0L};
            }
         }
#ifdef BOOST_HAS_FLOAT128
        if constexpr (std::is_same_v<Real, boost::multiprecision::float128>) {
            if constexpr (order == 0) {
                return {0x0p+0Q, -0x1.76cf5d0b09954e764ae85ae0f168p-2Q, 0x1.5db3d742c265539d92ba16b83c5ap+0Q, -0x0p+0Q};
            }
            if constexpr (order == 1) {
                return {0x0p+0Q, 0x1p+0Q, -0x1p+0Q, -0x0p+0Q};
            }
         }
#endif
     }
    if constexpr (p == 3) {
        if constexpr (std::is_same_v<Real, float>) {
            if constexpr (order == 0) {
                return {0x0p+0f, 0x1.158088p-8f, 0x1.863744p-4f, -0x1.8b18d8p-2f, 0x1.494d42p+0f, 0x0p+0f};
            }
            if constexpr (order == 1) {
                return {0x0p+0f, -0x1.69a5e8p-5f, -0x1.19ae7cp-1f, 0x1.1dcb06p+1f, -0x1.a3719cp+0f, -0x0p+0f};
            }
            if constexpr (order == 2) {
                return {0x0p+0f, 0x1.8831a2p-4f, 0x1.6ced64p-1f, -0x1.b676b2p+0f, 0x1.cef9ccp-1f, 0x0p+0f};
            }
         }
        if constexpr (std::is_same_v<Real, double>) {
            if constexpr (order == 0) {
                return {0x0p+0, 0x1.158087f14084dp-8, 0x1.863743274d78dp-4, -0x1.8b18d8251ec88p-2, 0x1.494d414ee19a1p+0, 0x0p+0};
            }
            if constexpr (order == 1) {
                return {0x0p+0, -0x1.69a5e70c6cb47p-5, -0x1.19ae7cc6c0212p-1, 0x1.1dcb0537f529p+1, -0x1.a3719cd426dbdp+0, -0x0p+0};
            }
            if constexpr (order == 2) {
                return {0x0p+0, 0x1.8831a237a06dfp-4, 0x1.6ced632b23d6cp-1, -0x1.b676b19591eb6p+0, 0x1.cef9cbb90bf24p-1, 0x0p+0};
            }
         }
        if constexpr (std::is_same_v<Real, long double>) {
            if constexpr (order == 0) {
                return {0x0p+0L, 0x8.ac043f8a042675ap-11L, 0xc.31ba193a6bc67f2p-7L, -0xc.58c6c128f64431dp-5L, 0xa.4a6a0a770cd07e1p-3L, 0x0p+0L};
            }
            if constexpr (order == 1) {
                return {0x0p+0L, -0xb.4d2f386365a363ap-8L, -0x8.cd73e6360108efbp-4L, 0x8.ee5829bfa948208p-2L, -0xd.1b8ce6a136deae1p-3L, -0x0p+0L};
            }
            if constexpr (order == 2) {
                return {0x0p+0L, 0xc.418d11bd036f5a2p-7L, 0xb.676b19591eb63e3p-4L, -0xd.b3b58cac8f5b1f2p-3L, 0xe.77ce5dc85f9214cp-4L, 0x0p+0L};
            }
         }
#ifdef BOOST_HAS_FLOAT128
        if constexpr (std::is_same_v<Real, boost::multiprecision::float128>) {
            if constexpr (order == 0) {
                return {0x0p+0Q, 0x1.158087f14084ceb4f650d2c43fd9p-8Q, 0x1.863743274d78cfe42daf1d262e81p-4Q, -0x1.8b18d8251ec886399fc357ab2707p-2Q, 0x1.494d414ee19a0fc1701f9345a29ap+0Q, 0x0p+0Q};
            }
            if constexpr (order == 1) {
                return {0x0p+0Q, -0x1.69a5e70c6cb46c73632e8c1bb203p-5Q, -0x1.19ae7cc6c0211df5e41745563c88p-1Q, 0x1.1dcb0537f52904105ab1d13c6a8ep+1Q, -0x1.a3719cd426dbd5c2283e8b6cd95p+0Q, -0x0p+0Q};
            }
            if constexpr (order == 2) {
                return {0x0p+0Q, 0x1.8831a237a06deb43265d99170ff6p-4Q, 0x1.6ced632b23d6c7c6d19ce6975a17p-1Q, -0x1.b676b19591eb63e368ce734bad06p+0Q, 0x1.cef9cbb90bf242979b344cdd1e02p-1Q, 0x0p+0Q};
            }
         }
#endif
     }
    if constexpr (p == 4) {
        if constexpr (std::is_same_v<Real, float>) {
            if constexpr (order == 0) {
                return {0x0p+0f, 0x1.3be7b6p-16f, -0x1.3a0992p-10f, -0x1.817e96p-7f, 0x1.447d2ap-5f, -0x1.15313cp-5f, 0x1.01d5e4p+0f, -0x0p+0f};
            }
            if constexpr (order == 1) {
                return {0x0p+0f, -0x1.d0194cp-10f, 0x1.b6da96p-5f, 0x1.0d0e14p-3f, -0x1.314492p+0f, 0x1.648654p+1f, -0x1.c6aca8p+0f, -0x0p+0f};
            }
            if constexpr (order == 2) {
                return {0x0p+0f, 0x1.3a82a2p-6f, -0x1.0ae7c2p-2f, -0x1.26d9a4p-3f, 0x1.66e46cp+1f, -0x1.0eaaap+2f, 0x1.cf8cc6p+0f, 0x0p+0f};
            }
            if constexpr (order == 3) {
                return {0x0p+0f, -0x1.bcd16ep-4f, 0x1.233972p-1f, -0x1.83800cp-3f, -0x1.c2597p+0f, 0x1.2d48b8p+1f, -0x1.bb2f44p-1f, 0x0p+0f};
            }
         }
        if constexpr (std::is_same_v<Real, double>) {
            if constexpr (order == 0) {
                return {0x0p+0, 0x1.3be7b6cb630ap-16, -0x1.3a0992ca51117p-10, -0x1.817e96e0425edp-7, 0x1.447d293e37265p-5, -0x1.15313c61b3acbp-5, 0x1.01d5e443d831dp+0, -0x0p+0};
            }
            if constexpr (order == 1) {
                return {0x0p+0, -0x1.d0194bee1a174p-10, 0x1.b6da96c702378p-5, 0x1.0d0e14f1f7dc6p-3, -0x1.314491c6de2d3p+0, 0x1.6486543c8460bp+1, -0x1.c6aca7b3a61bp+0, -0x0p+0};
            }
            if constexpr (order == 2) {
                return {0x0p+0, 0x1.3a82a1b96295bp-6, -0x1.0ae7c18dd841ep-2, -0x1.26d9a4de19b74p-3, 0x1.66e46b718b501p+1, -0x1.0eaaa0c1e1c35p+2, 0x1.cf8cc69cc42a4p+0, 0x0p+0};
            }
            if constexpr (order == 3) {
                return {0x0p+0, -0x1.bcd16d3945752p-4, 0x1.233972e072028p-1, -0x1.83800be8c4de9p-3, -0x1.c2596fe640caep+0, 0x1.2d48b852668c6p+1, -0x1.bb2f43bc30b81p-1, 0x0p+0};
            }
         }
        if constexpr (std::is_same_v<Real, long double>) {
            if constexpr (order == 0) {
                return {0x0p+0L, 0x9.df3db65b184ff7ep-19L, -0x9.d04c9652888ba27p-13L, -0xc.0bf4b70212f684ap-10L, 0xa.23e949f1b93272cp-8L, -0x8.a989e30d9d65a3ap-8L, 0x8.0eaf221ec18e913p-3L, -0x0p+0L};
            }
            if constexpr (order == 1) {
                return {0x0p+0L, -0xe.80ca5f70d0ba167p-13L, 0xd.b6d4b63811bbfa2p-8L, 0x8.6870a78fbee2d2ap-6L, -0x9.8a248e36f1695bdp-3L, 0xb.2432a1e4230589p-2L, -0xe.35653d9d30d7c1cp-3L, -0x0p+0L};
            }
            if constexpr (order == 2) {
                return {0x0p+0L, 0x9.d4150dcb14ada24p-9L, -0x8.573e0c6ec20ec34p-5L, -0x9.36cd26f0cdb9e9ap-6L, 0xb.37235b8c5a805d9p-2L, -0x8.7555060f0e1aa31p-1L, 0xe.7c6634e62152087p-3L, 0x0p+0L};
            }
            if constexpr (order == 3) {
                return {0x0p+0L, -0xd.e68b69ca2ba92a2p-7L, 0x9.19cb97039014287p-4L, -0xc.1c005f4626f4b41p-6L, -0xe.12cb7f320656fb6p-3L, 0x9.6a45c2933463235p-2L, -0xd.d97a1de185c06cap-4L, 0x0p+0L};
            }
         }
#ifdef BOOST_HAS_FLOAT128
        if constexpr (std::is_same_v<Real, boost::multiprecision::float128>) {
            if constexpr (order == 0) {
                return {0x0p+0Q, 0x1.3be7b6cb6309fefbda0c9fa0e6d1p-16Q, -0x1.3a0992ca5111744ed1b4ed1a692ap-10Q, -0x1.817e96e0425ed094fa1a3a5cf0e9p-7Q, 0x1.447d293e37264e578b8650812244p-5Q, -0x1.15313c61b3acb474231dab078175p-5Q, 0x1.01d5e443d831d2252369827793d5p+0Q, -0x0p+0Q};
            }
            if constexpr (order == 1) {
                return {0x0p+0Q, -0x1.d0194bee1a1742cd9ebfbba71356p-10Q, 0x1.b6da96c702377f43a25f74c94587p-5Q, 0x1.0d0e14f1f7dc5a54f21ce3ff41cdp-3Q, -0x1.314491c6de2d2b798153746ca796p+0Q, 0x1.6486543c8460b11f1a4b4a357f12p+1Q, -0x1.c6aca7b3a61af838bb3208359f46p+0Q, -0x0p+0Q};
            }
            if constexpr (order == 2) {
                return {0x0p+0Q, 0x1.3a82a1b96295b447f46661fe5aa6p-6Q, -0x1.0ae7c18dd841d8689ede2f5256aep-2Q, -0x1.26d9a4de19b73d34cbab42cdd884p-3Q, 0x1.66e46b718b500bb2be3ec05c7a3dp+1Q, -0x1.0eaaa0c1e1c354610b35262e9003p+2Q, 0x1.cf8cc69cc42a410e51b272a7a2f6p+0Q, 0x0p+0Q};
            }
            if constexpr (order == 3) {
                return {0x0p+0Q, -0x1.bcd16d394575254330d81aed457dp-4Q, 0x1.233972e07202850d1be1b9fc81e9p-1Q, -0x1.83800be8c4de9681c9e31925caedp-3Q, -0x1.c2596fe640cadf6ca968f3b30e4fp+0Q, 0x1.2d48b852668c646a630392a2b548p+1Q, -0x1.bb2f43bc30b80d947c8a537a1f1p-1Q, 0x0p+0Q};
            }
         }
#endif
     }
    if constexpr (p == 5) {
        if constexpr (std::is_same_v<Real, float>) {
            if constexpr (order == 0) {
                return {0x0p+0f, 0x1.75f2b2p-23f, 0x1.3b27d2p-15f, -0x1.c4ab56p-10f, 0x1.8dd2bep-10f, 0x1.3100cap-5f, -0x1.754196p-3f, 0x1.cbd5c2p-2f, 0x1.646bf2p-1f, 0x0p+0f};
            }
            if constexpr (order == 1) {
                return {0x0p+0f, 0x1.984da2p-17f, 0x1.5e1ffp-10f, -0x1.08fcfcp-5f, 0x1.64d96ap-6f, 0x1.784198p-2f, -0x1.3c6448p+0f, 0x1.37cf44p+1f, -0x1.8eee7ap+0f, 0x0p+0f};
            }
            if constexpr (order == 2) {
                return {0x0p+0f, -0x1.65cdcp-11f, -0x1.3d5d0ep-5f, 0x1.882acap-2f, -0x1.549f4cp-1f, -0x1.2d10dcp+0f, 0x1.524518p+2f, -0x1.83ccap+2f, 0x1.22c5ccp+1f, -0x0p+0f};
            }
            if constexpr (order == 3) {
                return {0x0p+0f, 0x1.13ee46p-4f, 0x1.050134p+1f, -0x1.4a4504p+3f, 0x1.15aaf4p+4f, -0x1.0b1eap+3f, -0x1.f7fdf6p+2f, 0x1.48c27ap+3f, -0x1.9877acp+1f, -0x0p+0f};
            }
            if constexpr (order == 4) {
                return {0x0p+0f, 0x1.1a66dcp-7f, 0x1.2c68dap-3f, -0x1.9ada2ap-2f, 0x1.0227f8p-1f, -0x1.3d6cc6p+0f, 0x1.29ab2ep+1f, -0x1.e1e5e6p+0f, 0x1.13b9e8p-1f, 0x0p+0f};
            }
         }
        if constexpr (std::is_same_v<Real, double>) {
            if constexpr (order == 0) {
                return {0x0p+0, 0x1.75f2b16626e98p-23, 0x1.3b27d2d798eedp-15, -0x1.c4ab558ff2dcfp-10, 0x1.8dd2be8c89c52p-10, 0x1.3100cab7c3f5fp-5, -0x1.754196833f707p-3, 0x1.cbd5c1bab5148p-2, 0x1.646bf1ec64a31p-1, 0x0p+0};
            }
            if constexpr (order == 1) {
                return {0x0p+0, 0x1.984da10897439p-17, 0x1.5e1fef6185af9p-10, -0x1.08fcfb4783fa8p-5, 0x1.64d969ffaac15p-6, 0x1.7841978e876dbp-2, -0x1.3c64475174b4ap+0, 0x1.37cf445237f1bp+1, -0x1.8eee7927087b2p+0, 0x0p+0};
            }
            if constexpr (order == 2) {
                return {0x0p+0, -0x1.65cdc06ecf654p-11, -0x1.3d5d0eda32d8bp-5, 0x1.882ac9029ddc7p-2, -0x1.549f4c5e04f3ap-1, -0x1.2d10dba463fc2p+0, 0x1.5245172406b78p+2, -0x1.83cca07a3ead8p+2, 0x1.22c5cb8d3f23bp+1, -0x0p+0};
            }
            if constexpr (order == 3) {
                return {0x0p+0, 0x1.13ee46114fd85p-4, 0x1.05013480812dfp+1, -0x1.4a45041b91e1cp+3, 0x1.15aaf38ada16ap+4, -0x1.0b1ea0ee5f9edp+3, -0x1.f7fdf60845cfcp+2, 0x1.48c27a70b93cep+3, -0x1.9877ac926fb3cp+1, -0x0p+0};
            }
            if constexpr (order == 4) {
                return {0x0p+0, 0x1.1a66dc6d2cfbfp-7, 0x1.2c68da3893cd6p-3, -0x1.9ada2a6a9a218p-2, 0x1.0227f72aa277ep-1, -0x1.3d6cc61e0f6a4p+0, 0x1.29ab2ec864d45p+1, -0x1.e1e5e5cfd0a8p+0, 0x1.13b9e8c4fdc5p-1, 0x0p+0};
            }
         }
        if constexpr (std::is_same_v<Real, long double>) {
            if constexpr (order == 0) {
                return {0x0p+0L, 0xb.af958b31374c3acp-26L, 0x9.d93e96bcc776751p-18L, -0xe.255aac7f96e75fp-13L, 0xc.6e95f4644e28f3ep-13L, 0x9.880655be1fafa7dp-8L, -0xb.aa0cb419fb8363cp-6L, 0xe.5eae0dd5a8a4298p-5L, 0xb.235f8f632518463p-4L, 0x0p+0L};
            }
            if constexpr (order == 1) {
                return {0x0p+0L, 0xc.c26d0844ba1cb28p-20L, 0xa.f0ff7b0c2d7cb32p-13L, -0x8.47e7da3c1fd3e92p-8L, 0xb.26cb4ffd560a9d3p-9L, 0xb.c20cbc743b6db2dp-5L, -0x9.e3223a8ba5a5065p-3L, 0x9.be7a2291bf8db95p-2L, -0xc.7773c93843d9207p-3L, 0x0p+0L};
            }
            if constexpr (order == 2) {
                return {0x0p+0L, -0xb.2e6e03767b2a39ap-14L, -0x9.eae876d196c5b0fp-8L, 0xc.41564814eee395ap-5L, -0xa.a4fa62f0279cc8ap-4L, -0x9.6886dd231fe0c2p-3L, 0xa.9228b92035bbe8fp-1L, -0xc.1e6503d1f56bc54p-1L, 0x9.162e5c69f91da28p-2L, -0x0p+0L};
            }
            if constexpr (order == 3) {
                return {0x0p+0L, 0x8.9f72308a7ec28cep-7L, 0x8.2809a404096fa64p-2L, -0xa.522820dc8f0dd92p+0L, 0x8.ad579c56d0b531bp+1L, -0x8.58f50772fcf6838p+0L, -0xf.bfefb0422e7dd11p-1L, 0xa.4613d385c9e6eeep+0L, -0xc.c3bd64937d9deedp-2L, -0x0p+0L};
            }
            if constexpr (order == 4) {
                return {0x0p+0L, 0x8.d336e36967df815p-10L, 0x9.6346d1c49e6ae26p-6L, -0xc.d6d15354d10bf47p-5L, 0x8.113fb95513bed76p-4L, -0x9.eb6630f07b521f1p-3L, 0x9.4d59764326a2533p-2L, -0xf.0f2f2e7e85400fbp-3L, 0x8.9dcf4627ee27fcfp-4L, 0x0p+0L};
            }
         }
#ifdef BOOST_HAS_FLOAT128
        if constexpr (std::is_same_v<Real, boost::multiprecision::float128>) {
            if constexpr (order == 0) {
                return {0x0p+0Q, 0x1.75f2b16626e9875840e69284f9adp-23Q, 0x1.3b27d2d798eecea1adaf099bbe69p-15Q, -0x1.c4ab558ff2dcebdfdc6142054c48p-10Q, 0x1.8dd2be8c89c51e7b28b0431533bdp-10Q, 0x1.3100cab7c3f5f4fa75b80a67ed16p-5Q, -0x1.754196833f706c7834855fa5579dp-3Q, 0x1.cbd5c1bab5148530b0e77ecbff84p-2Q, 0x1.646bf1ec64a308c5df07d89c8ea8p-1Q, 0x0p+0Q};
            }
            if constexpr (order == 1) {
                return {0x0p+0Q, 0x1.984da1089743964ff5ac0958e062p-17Q, 0x1.5e1fef6185af966378a125174c6fp-10Q, -0x1.08fcfb4783fa7d234eba360539a9p-5Q, 0x1.64d969ffaac153a6af06dc9cfca7p-6Q, 0x1.7841978e876db6599aa4ad10999dp-2Q, -0x1.3c64475174b4a0c9b8a3ee3dfc83p+0Q, 0x1.37cf445237f1b72a99ae1be261d4p+1Q, -0x1.8eee7927087b240ec5ade1ac8229p+0Q, 0x0p+0Q};
            }
            if constexpr (order == 2) {
                return {0x0p+0Q, -0x1.65cdc06ecf6547344a6a8d5255a7p-11Q, -0x1.3d5d0eda32d8b61e7ff27a395499p-5Q, 0x1.882ac9029ddc72b4c74d388da778p-2Q, -0x1.549f4c5e04f3991456e630b15374p-1Q, -0x1.2d10dba463fc1840fd45963c56bap+0Q, 0x1.5245172406b77d1d4fdfd2192d64p+2Q, -0x1.83cca07a3ead78a73330d7cf1a68p+2Q, 0x1.22c5cb8d3f23b44faf592c365fbep+1Q, -0x0p+0Q};
            }
            if constexpr (order == 3) {
                return {0x0p+0Q, 0x1.13ee46114fd8519c5d7f29b6f3d6p-4Q, 0x1.05013480812df4c7388d78d9d59cp+1Q, -0x1.4a45041b91e1bb24fe9a63ac1989p+3Q, 0x1.15aaf38ada16a6358cd773551d0ap+4Q, -0x1.0b1ea0ee5f9ed06f5107d7991a08p+3Q, -0x1.f7fdf60845cfba22e0c6bd03e97cp+2Q, 0x1.48c27a70b93cdddc043571f40d16p+3Q, -0x1.9877ac926fb3bdd992f46d84088p+1Q, -0x0p+0Q};
            }
            if constexpr (order == 4) {
                return {0x0p+0Q, 0x1.1a66dc6d2cfbf029033eb2e6fe0ep-7Q, 0x1.2c68da3893cd5c4c84ac50c60876p-3Q, -0x1.9ada2a6a9a217e8dbe2b55c6b7e8p-2Q, 0x1.0227f72aa277daecfc1fddb3c5b5p-1Q, -0x1.3d6cc61e0f6a43e21f6553d1005fp+0Q, 0x1.29ab2ec864d44a65c8f8abfda3p+1Q, -0x1.e1e5e5cfd0a801f531eaef0fd8b4p+0Q, 0x1.13b9e8c4fdc4ff9d9c7b93fd9e49p-1Q, 0x0p+0Q};
            }
         }
#endif
     }
    if constexpr (p == 6) {
        if constexpr (std::is_same_v<Real, float>) {
            if constexpr (order == 0) {
                return {0x0p+0f, 0x1.0fbc42p-28f, -0x1.5a00dp-19f, 0x1.05937cp-16f, 0x1.cd4fecp-10f, -0x1.cec09ap-9f, -0x1.a1e9c8p-6f, 0x1.24737ep-3f, -0x1.89fd1p-2f, 0x1.aa2466p-1f, 0x1.bf65e8p-2f, -0x0p+0f};
            }
            if constexpr (order == 1) {
                return {0x0p+0f, 0x1.018922p-24f, -0x1.45b126p-16f, 0x1.a69b86p-15f, 0x1.af662cp-8f, -0x1.721a1ep-7f, -0x1.0b11b2p-4f, 0x1.0ffe44p-2f, -0x1.050336p-1f, 0x1.801498p+0f, -0x1.2faeacp+0f, 0x0p+0f};
            }
            if constexpr (order == 2) {
                return {0x0p+0f, -0x1.0c5bep-16f, 0x1.4ebb22p-9f, 0x1.7eb866p-6f, -0x1.4137dp-2f, 0x1.da7724p-1f, -0x1.0310f6p-1f, -0x1.45eb5ep+1f, 0x1.a34deap+2f, -0x1.980b4p+2f, 0x1.1e62dep+1f, 0x0p+0f};
            }
            if constexpr (order == 3) {
                return {0x0p+0f, -0x1.816ef8p-8f, 0x1.d3687ap-2f, 0x1.77f884p+1f, -0x1.804948p+4f, 0x1.d54dd2p+5f, -0x1.10b7fap+6f, 0x1.3c98d8p+5f, -0x1.484b54p+3f, 0x1.255a76p+0f, -0x1.422fap-2f, -0x0p+0f};
            }
            if constexpr (order == 4) {
                return {0x0p+0f, 0x1.b96594p-10f, -0x1.f8b0a6p-5f, -0x1.c41888p-3f, 0x1.77a6f6p+0f, -0x1.46377ep+1f, 0x1.c7801ep+1f, -0x1.c3f2cep+2f, 0x1.33788p+3f, -0x1.93816ep+2f, 0x1.8faefp+0f, -0x0p+0f};
            }
            if constexpr (order == 5) {
                return {0x0p+0f, -0x1.0ac90ap-9f, 0x1.0c11f4p-5f, 0x1.c1a09cp-5f, -0x1.67c138p-2f, 0x1.8461e6p-2f, -0x1.29da94p-1f, 0x1.fa645ep+0f, -0x1.7b2e48p+1f, 0x1.e6b982p+0f, -0x1.cb9194p-2f, 0x0p+0f};
            }
         }
        if constexpr (std::is_same_v<Real, double>) {
            if constexpr (order == 0) {
                return {0x0p+0, 0x1.0fbc42c672231p-28, -0x1.5a00cfad970a8p-19, 0x1.05937b8388af2p-16, 0x1.cd4feb82d2494p-10, -0x1.cec09a7ccf05bp-9, -0x1.a1e9c758a51a1p-6, 0x1.24737e6f8ebf1p-3, -0x1.89fd104c6ff15p-2, 0x1.aa2466d50e4a5p-1, 0x1.bf65e79817d28p-2, -0x0p+0};
            }
            if constexpr (order == 1) {
                return {0x0p+0, 0x1.01892163486cbp-24, -0x1.45b12643cf82p-16, 0x1.a69b860b67a33p-15, 0x1.af662c52c3b7bp-8, -0x1.721a1d14b8323p-7, -0x1.0b11b28c7d67bp-4, 0x1.0ffe44cc8bc7fp-2, -0x1.05033635f5402p-1, 0x1.8014975295b21p+0, -0x1.2faeacbb8e754p+0, 0x0p+0};
            }
            if constexpr (order == 2) {
                return {0x0p+0, -0x1.0c5bdfb916a45p-16, 0x1.4ebb22d1ac076p-9, 0x1.7eb866ddca6bap-6, -0x1.4137cfe85b9e6p-2, 0x1.da7724c65802cp-1, -0x1.0310f5a5d86dep-1, -0x1.45eb5e5e7e3cap+1, 0x1.a34dea8cd59b4p+2, -0x1.980b3fe2011f8p+2, 0x1.1e62ddd540b06p+1, 0x0p+0};
            }
            if constexpr (order == 3) {
                return {0x0p+0, -0x1.816ef86832622p-8, 0x1.d3687a3bc3933p-2, 0x1.77f8847acd28bp+1, -0x1.8049478b6740dp+4, 0x1.d54dd14cd31afp+5, -0x1.10b7f9c607064p+6, 0x1.3c98d8a74cd77p+5, -0x1.484b537c56dedp+3, 0x1.255a76df553ccp+0, -0x1.422f9f5227e8p-2, -0x0p+0};
            }
            if constexpr (order == 4) {
                return {0x0p+0, 0x1.b965933f534b2p-10, -0x1.f8b0a606c7f79p-5, -0x1.c4188733407f2p-3, 0x1.77a6f5e8321cdp+0, -0x1.46377db27f625p+1, 0x1.c7801db767804p+1, -0x1.c3f2cd3d53133p+2, 0x1.337880b3a9a71p+3, -0x1.93816e52ba1c8p+2, 0x1.8faeef62b3a8p+0, -0x0p+0};
            }
            if constexpr (order == 5) {
                return {0x0p+0, -0x1.0ac90a45b7a47p-9, 0x1.0c11f4efd5c0bp-5, 0x1.c1a09b7773a5ap-5, -0x1.67c137342024dp-2, 0x1.8461e6acf9cb3p-2, -0x1.29da9349c32a2p-1, 0x1.fa645ee8ca231p+0, -0x1.7b2e4714b2703p+1, 0x1.e6b982120ca32p+0, -0x1.cb91942378a6bp-2, 0x0p+0};
            }
         }
        if constexpr (std::is_same_v<Real, long double>) {
            if constexpr (order == 0) {
                return {0x0p+0L, 0x8.7de216339118a43p-31L, -0xa.d0067d6cb853e5ep-22L, 0x8.2c9bdc1c4579256p-19L, 0xe.6a7f5c169249c7ap-13L, -0xe.7604d3e6782db77p-12L, -0xd.0f4e3ac528d0419p-9L, 0x9.239bf37c75f8764p-6L, -0xc.4fe882637f8a4ep-5L, 0xd.512336a872527e6p-4L, 0xd.fb2f3cc0be93e54p-5L, -0x0p+0L};
            }
            if constexpr (order == 1) {
                return {0x0p+0L, 0x8.0c490b1a4365a0ep-27L, -0xa.2d89321e7c101dfp-19L, 0xd.34dc305b3d198bdp-18L, 0xd.7b3162961dbd72bp-11L, -0xb.90d0e8a5c191592p-10L, -0x8.588d9463eb3d6b2p-7L, 0x8.7ff226645e3f525p-5L, -0x8.2819b1afaa00e4ep-4L, 0xc.00a4ba94ad909b7p-3L, -0x9.7d7565dc73a9dap-3L, 0x0p+0L};
            }
            if constexpr (order == 2) {
                return {0x0p+0L, -0x8.62defdc8b5224edp-19L, 0xa.75d9168d603b071p-12L, 0xb.f5c336ee535d213p-9L, -0xa.09be7f42dcf2d82p-5L, 0xe.d3b92632c015e4bp-4L, -0x8.1887ad2ec36f306p-4L, -0xa.2f5af2f3f1e4f4dp-2L, 0xd.1a6f5466acda0c4p-1L, -0xc.c059ff1008fbde4p-1L, 0x8.f316eeaa0582d86p-2L, 0x0p+0L};
            }
            if constexpr (order == 3) {
                return {0x0p+0L, -0xc.0b77c3419310f5ep-11L, 0xe.9b43d1de1c99addp-5L, 0xb.bfc423d6694569bp-2L, -0xc.024a3c5b3a06863p+1L, 0xe.aa6e8a6698d7b2fp+2L, -0x8.85bfce3038320c3p+3L, 0x9.e4c6c53a66bb9b6p+2L, -0xa.425a9be2b6f651fp+0L, 0x9.2ad3b6faa9e5fc2p-3L, -0xa.117cfa913f3fd65p-5L, -0x0p+0L};
            }
            if constexpr (order == 4) {
                return {0x0p+0L, 0xd.cb2c99fa9a58d79p-13L, -0xf.c58530363fbca91p-8L, -0xe.20c4399a03f920dp-6L, 0xb.bd37af4190e69cdp-3L, -0xa.31bbed93fb124adp-2L, 0xe.3c00edbb3c02374p-2L, -0xe.1f9669ea9899554p-1L, 0x9.9bc4059d4d38bf3p+0L, -0xc.9c0b7295d0e4242p-1L, 0xc.7d777b159d3fe93p-3L, -0x0p+0L};
            }
            if constexpr (order == 5) {
                return {0x0p+0L, -0x8.5648522dbd23badp-12L, 0x8.608fa77eae05beap-8L, 0xe.0d04dbbb9d2ce8dp-8L, -0xb.3e09b9a10126b64p-5L, 0xc.230f3567ce5957dp-5L, -0x9.4ed49a4e1950cdfp-4L, 0xf.d322f7465118802p-3L, -0xb.d97238a59381a17p-2L, 0xf.35cc10906518e34p-3L, -0xe.5c8ca11bc53594fp-5L, 0x0p+0L};
            }
         }
#ifdef BOOST_HAS_FLOAT128
        if constexpr (std::is_same_v<Real, boost::multiprecision::float128>) {
            if constexpr (order == 0) {
                return {0x0p+0Q, 0x1.0fbc42c672231485a0bcb5ce651p-28Q, -0x1.5a00cfad970a7cbc3091390140a9p-19Q, 0x1.05937b8388af24ac43ec83961206p-16Q, 0x1.cd4feb82d24938f30e41e1af2ed2p-10Q, -0x1.cec09a7ccf05b6ee1ea8f3c1832ep-9Q, -0x1.a1e9c758a51a08321aa63e9af1cp-6Q, 0x1.24737e6f8ebf0ec8bb9fcc85bb74p-3Q, -0x1.89fd104c6ff149bf7514a7ddf39fp-2Q, 0x1.aa2466d50e4a4fcbbebe54f94c47p-1Q, 0x1.bf65e79817d27ca7c78972c65f46p-2Q, -0x0p+0Q};
            }
            if constexpr (order == 1) {
                return {0x0p+0Q, 0x1.01892163486cb41cb64452f47fccp-24Q, -0x1.45b12643cf8203bdf796367eefcdp-16Q, 0x1.a69b860b67a33179a8c7a5f6dcf6p-15Q, 0x1.af662c52c3b7ae55063672df3f39p-8Q, -0x1.721a1d14b8322b2431432e8b8593p-7Q, -0x1.0b11b28c7d67ad63229fbc3aa8bfp-4Q, 0x1.0ffe44cc8bc7ea4a6ce4734d8e7p-2Q, -0x1.05033635f5401c9b60636a3253fbp-1Q, 0x1.8014975295b2136d5a1adb4544dap+0Q, -0x1.2faeacbb8e753b40eb856e00f5ddp+0Q, 0x0p+0Q};
            }
            if constexpr (order == 2) {
                return {0x0p+0Q, -0x1.0c5bdfb916a449da5acc3cc37bffp-16Q, 0x1.4ebb22d1ac0760e14618810576c3p-9Q, 0x1.7eb866ddca6ba425f9e4eb850119p-6Q, -0x1.4137cfe85b9e5b04e1023975514ap-2Q, 0x1.da7724c65802bc95482213f92dfbp-1Q, -0x1.0310f5a5d86de60ba877dea64942p-1Q, -0x1.45eb5e5e7e3c9e9ae26ca4a9f7acp+1Q, 0x1.a34dea8cd59b4187717c6efad071p+2Q, -0x1.980b3fe2011f7bc878be7d39672p+2Q, 0x1.1e62ddd540b05b0cc5ce586fe909p+1Q, 0x0p+0Q};
            }
            if constexpr (order == 3) {
                return {0x0p+0Q, -0x1.816ef86832621ebc63e47ca48658p-8Q, 0x1.d3687a3bc39335bac7f281aa3a0bp-2Q, 0x1.77f8847acd28ad3549099aaf68abp+1Q, -0x1.8049478b6740d0c6b01004455301p+4Q, 0x1.d54dd14cd31af65e5343fc4601c9p+5Q, -0x1.10b7f9c607064185ab2dfc32d362p+6Q, 0x1.3c98d8a74cd7736b60af4e3ce406p+5Q, -0x1.484b537c56deca3daefc65807dc2p+3Q, 0x1.255a76df553cbf84e8a0b1b9656cp+0Q, -0x1.422f9f5227e7facaa94dd9538c75p-2Q, -0x0p+0Q};
            }
            if constexpr (order == 4) {
                return {0x0p+0Q, 0x1.b965933f534b1af276b672235eaep-10Q, -0x1.f8b0a606c7f795229d59473e7f87p-5Q, -0x1.c4188733407f241acae534e90799p-3Q, 0x1.77a6f5e8321cd39a930db38d10c9p+0Q, -0x1.46377db27f624959a2b9b775569bp+1Q, 0x1.c7801db7678046e8ebdb22cc275fp+1Q, -0x1.c3f2cd3d53132aa78f87aedc0bacp+2Q, 0x1.337880b3a9a717e67b9cd7fb9268p+3Q, -0x1.93816e52ba1c8483ff81ea3a089cp+2Q, 0x1.8faeef62b3a7fd26eb98dd7b97c2p+0Q, -0x0p+0Q};
            }
            if constexpr (order == 5) {
                return {0x0p+0Q, -0x1.0ac90a45b7a47759a7a9273c0c64p-9Q, 0x1.0c11f4efd5c0b7d44ac1c3a203ffp-5Q, 0x1.c1a09b7773a59d1a65f356784661p-5Q, -0x1.67c137342024d6c8052e3d885daap-2Q, 0x1.8461e6acf9cb2af90934d76ca1e6p-2Q, -0x1.29da9349c32a19bed5b810cc82bap-1Q, 0x1.fa645ee8ca231003817fbe341bfep+0Q, -0x1.7b2e4714b270342ea6ba8998bb7bp+1Q, 0x1.e6b982120ca31c68b89329c21facp+0Q, -0x1.cb91942378a6b29e2dd51c5322fap-2Q, 0x0p+0Q};
            }
         }
#endif
     }
    if constexpr (p == 7) {
        if constexpr (std::is_same_v<Real, float>) {
            if constexpr (order == 0) {
                return {0x0p+0f, -0x1.57f8e6p-36f, -0x1.509d12p-25f, 0x1.310ffep-19f, -0x1.17d0b4p-14f, -0x1.1d4c6cp-11f, 0x1.d0243ap-10f, 0x1.50e2d6p-7f, -0x1.1142c4p-4f, 0x1.7a0caap-3f, -0x1.919d92p-2f, 0x1.027f6cp+0f, 0x1.033292p-2f, 0x0p+0f};
            }
            if constexpr (order == 1) {
                return {0x0p+0f, 0x1.54fd9p-33f, 0x1.4e8b2ap-23f, -0x1.6a87a6p-15f, -0x1.0efaecp-14f, 0x1.051eecp-8f, -0x1.5cb5dep-6f, 0x1.4e8938p-5f, 0x1.4ab3dp-6f, -0x1.145888p-2f, 0x1.2d988ap-1f, 0x1.ca0692p-2f, -0x1.9ec434p-1f, -0x0p+0f};
            }
            if constexpr (order == 2) {
                return {0x0p+0f, 0x1.cb0cccp-25f, 0x1.c4a78ap-16f, -0x1.4e09dp-9f, 0x1.eb4ccap-12f, 0x1.fe16d2p-4f, -0x1.18b244p-1f, 0x1.9d599ap-1f, 0x1.570932p-2f, -0x1.86bc52p+1f, 0x1.70aa18p+2f, -0x1.5362d8p+2f, 0x1.e0eaccp+0f, 0x0p+0f};
            }
            if constexpr (order == 3) {
                return {0x0p+0f, 0x1.64d198p-14f, 0x1.6365bp-6f, -0x1.13b5aap+0f, -0x1.833e7p-4f, 0x1.a43b42p+4f, -0x1.7ff792p+6f, 0x1.415156p+7f, -0x1.23bcc4p+7f, 0x1.1e704ep+6f, -0x1.171de6p+4f, 0x1.492e7p+1f, -0x1.4ca934p-1f, 0x0p+0f};
            }
            if constexpr (order == 4) {
                return {0x0p+0f, -0x1.1f66b2p-14f, -0x1.23f976p-7f, 0x1.7c490ep-3f, -0x1.3be428p-3f, -0x1.ffd4fcp+0f, 0x1.8c2cccp+2f, -0x1.22247ap+3f, 0x1.9815ccp+3f, -0x1.40c582p+4f, 0x1.578438p+4f, -0x1.7eda4cp+3f, 0x1.525f0cp+1f, -0x0p+0f};
            }
            if constexpr (order == 5) {
                return {0x0p+0f, 0x1.c9a574p-12f, 0x1.e323c6p-6f, -0x1.16c54ap-2f, 0x1.9df5b4p-2f, 0x1.97de12p-1f, -0x1.2d09cap+1f, 0x1.34da8p+1f, -0x1.4010d2p+2f, 0x1.74c81ep+3f, -0x1.ac8b3cp+3f, 0x1.ce7e68p+2f, -0x1.80a316p+0f, 0x0p+0f};
            }
            if constexpr (order == 6) {
                return {0x0p+0f, 0x1.24061ep-9f, 0x1.4b884p-4f, -0x1.656934p-2f, 0x1.256afep-1f, -0x1.6af582p-1f, 0x1.3544c2p+0f, -0x1.e6f03p+0f, 0x1.2a5764p+1f, -0x1.40ce38p+1f, 0x1.099a08p+1f, -0x1.01d8c8p+0f, 0x1.a096c8p-3f, -0x0p+0f};
            }
         }
        if constexpr (std::is_same_v<Real, double>) {
            if constexpr (order == 0) {
                return {0x0p+0, -0x1.57f8e69552537p-36, -0x1.509d126dab99dp-25, 0x1.310ffd8b56401p-19, -0x1.17d0b3032ba5cp-14, -0x1.1d4c6c249bf0fp-11, 0x1.d0243aab613f8p-10, 0x1.50e2d6421ee26p-7, -0x1.1142c3515ec43p-4, 0x1.7a0ca906acc48p-3, -0x1.919d92c6034c4p-2, 0x1.027f6bf7952cep+0, 0x1.033291a6c76d7p-2, 0x0p+0};
            }
            if constexpr (order == 1) {
                return {0x0p+0, 0x1.54fd8f0625f12p-33, 0x1.4e8b2a27fe065p-23, -0x1.6a87a6b931edcp-15, -0x1.0efaec65edbccp-14, 0x1.051eebbefaf7cp-8, -0x1.5cb5de436bfaep-6, 0x1.4e8938899ee37p-5, 0x1.4ab3cf0bae2a6p-6, -0x1.145887eacb11p-2, 0x1.2d988a986b0dfp-1, 0x1.ca0692ca7743ap-2, -0x1.9ec434341eaeep-1, -0x0p+0};
            }
            if constexpr (order == 2) {
                return {0x0p+0, 0x1.cb0ccbca08e79p-25, 0x1.c4a78957d2798p-16, -0x1.4e09cf12f0cddp-9, 0x1.eb4cc917dd1cdp-12, 0x1.fe16d1a3cfb3ep-4, -0x1.18b244e50df67p-1, 0x1.9d599937b31fbp-1, 0x1.5709310e8d1f6p-2, -0x1.86bc52429f3a4p+1, 0x1.70aa17b7adc0ep+2, -0x1.5362d7fe36c57p+2, 0x1.e0eacca617ebcp+0, 0x0p+0};
            }
            if constexpr (order == 3) {
                return {0x0p+0, 0x1.64d19861a8338p-14, 0x1.6365af2495a87p-6, -0x1.13b5aaaf6c71dp+0, -0x1.833e6f5afaca5p-4, 0x1.a43b4147c122cp+4, -0x1.7ff7929036dddp+6, 0x1.415156a34428ap+7, -0x1.23bcc34d6a52ap+7, 0x1.1e704d84a1075p+6, -0x1.171de5ca88cecp+4, 0x1.492e70f2635c1p+1, -0x1.4ca9343b32368p-1, 0x0p+0};
            }
            if constexpr (order == 4) {
                return {0x0p+0, -0x1.1f66b1e781529p-14, -0x1.23f9756392d0fp-7, 0x1.7c490d447fad3p-3, -0x1.3be427404fe7dp-3, -0x1.ffd4fb5b3918p+0, 0x1.8c2cccc6ee34dp+2, -0x1.222479cd96eeap+3, 0x1.9815cb7223c1p+3, -0x1.40c582473fa8dp+4, 0x1.578437542fe4fp+4, -0x1.7eda4b6a23329p+3, 0x1.525f0cb11c272p+1, -0x0p+0};
            }
            if constexpr (order == 5) {
                return {0x0p+0, 0x1.c9a573f22c19bp-12, 0x1.e323c5b92bd16p-6, -0x1.16c54a229aae1p-2, 0x1.9df5b4fe2e272p-2, 0x1.97de12ea03ca1p-1, -0x1.2d09c9846b582p+1, 0x1.34da807d2aaedp+1, -0x1.4010d1eedca9ap+2, 0x1.74c81d1a27df6p+3, -0x1.ac8b3ccd718c4p+3, 0x1.ce7e682714367p+2, -0x1.80a316d21a0efp+0, 0x0p+0};
            }
            if constexpr (order == 6) {
                return {0x0p+0, 0x1.24061e2653bdbp-9, 0x1.4b88409bdd117p-4, -0x1.6569331439298p-2, 0x1.256afda077268p-1, -0x1.6af581a622d41p-1, 0x1.3544c2755b1d8p+0, -0x1.e6f030cb23fcp+0, 0x1.2a57634989c52p+1, -0x1.40ce386ca776ep+1, 0x1.099a073a7c7b4p+1, -0x1.01d8c748e11fap+0, 0x1.a096c8f7fc8abp-3, -0x0p+0};
            }
         }
        if constexpr (std::is_same_v<Real, long double>) {
            if constexpr (order == 0) {
                return {0x0p+0L, -0xa.bfc734aa929b8c7p-39L, -0xa.84e8936d5cce758p-28L, 0x9.887fec5ab2008edp-22L, -0x8.be8598195d2e2f7p-17L, -0x8.ea636124df87704p-14L, 0xe.8121d55b09fc1bcp-13L, 0xa.8716b210f7132a1p-10L, -0x8.8a161a8af621982p-7L, 0xb.d0654835662420dp-6L, -0xc.8cec96301a62254p-5L, 0x8.13fb5fbca966de8p-3L, 0x8.19948d363b6b511p-5L, 0x0p+0L};
            }
            if constexpr (order == 1) {
                return {0x0p+0L, 0xa.a7ec78312f88c67p-36L, 0xa.7459513ff03276p-26L, -0xb.543d35c98f6e2d6p-18L, -0x8.77d7632f6de608ep-17L, 0x8.28f75df7d7bdeb9p-11L, -0xa.e5aef21b5fd6dcfp-9L, 0xa.7449c44cf71b8f2p-8L, 0xa.559e785d7152dc2p-9L, -0x8.a2c43f565887ceap-5L, 0x9.6cc454c3586f53bp-4L, 0xe.50349653ba1cf48p-5L, -0xc.f621a1a0f576c66p-4L, -0x0p+0L};
            }
            if constexpr (order == 2) {
                return {0x0p+0L, 0xe.58665e50473cb0fp-28L, 0xe.253c4abe93cc026p-19L, -0xa.704e7897866e7ecp-12L, 0xf.5a6648bee8e6ab7p-15L, 0xf.f0b68d1e7d9f3f6p-7L, -0x8.c59227286fb3ad2p-4L, 0xc.eaccc9bd98fd7eap-4L, 0xa.b849887468fb1bp-5L, -0xc.35e29214f9d1fbdp-2L, 0xb.8550bdbd6e073adp-1L, -0xa.9b16bff1b62b5b5p-1L, 0xf.07566530bf5dc3p-3L, 0x0p+0L};
            }
            if constexpr (order == 3) {
                return {0x0p+0L, 0xb.268cc30d419bc95p-17L, 0xb.1b2d7924ad4363p-9L, -0x8.9dad557b638eaddp-3L, -0xc.19f37ad7d65262bp-7L, 0xd.21da0a3e0916379p+1L, -0xb.ffbc9481b6eebc6p+3L, 0xa.0a8ab51a2145353p+4L, -0x9.1de61a6b5294d86p+4L, 0x8.f3826c25083a83ap+3L, -0x8.b8ef2e544676158p+1L, 0xa.497387931ae0b61p-2L, -0xa.6549a1d991b412ap-4L, 0x0p+0L};
            }
            if constexpr (order == 4) {
                return {0x0p+0L, -0x8.fb358f3c0a94afap-17L, -0x9.1fcbab1c968774dp-10L, 0xb.e2486a23fd69714p-6L, -0x9.df213a027f3e892p-6L, -0xf.fea7dad9c8c0222p-3L, 0xc.6166663771a6508p-1L, -0x9.1123ce6cb774e8cp+0L, 0xc.c0ae5b911e082dep+0L, -0xa.062c1239fd4649cp+1L, 0xa.bc21baa17f27afep+1L, -0xb.f6d25b5119947a7p+0L, 0xa.92f86588e139024p-2L, -0x0p+0L};
            }
            if constexpr (order == 5) {
                return {0x0p+0L, 0xe.4d2b9f9160cd482p-15L, 0xf.191e2dc95e8b372p-9L, -0x8.b62a5114d570977p-5L, 0xc.efada7f17138d13p-5L, 0xc.bef097501e50849p-4L, -0x9.684e4c235ac0d81p-2L, 0x9.a6d403e95576488p-2L, -0xa.00868f76e54ce56p-1L, 0xb.a640e8d13efae1ap+0L, -0xd.6459e66b8c620fp+0L, 0xe.73f34138a1b3bfep-1L, -0xc.0518b690d0776a4p-3L, 0x0p+0L};
            }
            if constexpr (order == 6) {
                return {0x0p+0L, 0x9.2030f1329ded9e3p-12L, 0xa.5c4204dee88b95dp-7L, -0xb.2b4998a1c94bea3p-5L, 0x9.2b57ed03b933f03p-4L, -0xb.57ac0d3116a0b56p-4L, 0x9.aa2613aad8ec175p-3L, -0xf.378186591fdfef3p-3L, 0x9.52bb1a4c4e28fdbp-2L, -0xa.0671c3653bb7036p-2L, 0x8.4cd039d3e3da0c5p-2L, -0x8.0ec63a4708fd39bp-3L, 0xd.04b647bfe4558acp-6L, -0x0p+0L};
            }
         }
#ifdef BOOST_HAS_FLOAT128
        if constexpr (std::is_same_v<Real, boost::multiprecision::float128>) {
            if constexpr (order == 0) {
                return {0x0p+0Q, -0x1.57f8e6955253718d03c58931bd7dp-36Q, -0x1.509d126dab99ceafee21c512e8c2p-25Q, 0x1.310ffd8b564011dafec30495b01bp-19Q, -0x1.17d0b3032ba5c5edbe3115296ab6p-14Q, -0x1.1d4c6c249bf0ee08cde9105ae26bp-11Q, 0x1.d0243aab613f83777184da54aabp-10Q, 0x1.50e2d6421ee2654230b30972367bp-7Q, -0x1.1142c3515ec433047c393c699626p-4Q, 0x1.7a0ca906acc4841a1191b023a479p-3Q, -0x1.919d92c6034c44a7ed45a578713p-2Q, 0x1.027f6bf7952cdbcf397ade5b839ep+0Q, 0x1.033291a6c76d6a2177301a6491e3p-2Q, 0x0p+0Q};
            }
            if constexpr (order == 1) {
                return {0x0p+0Q, 0x1.54fd8f0625f118cd9bccf4de7011p-33Q, 0x1.4e8b2a27fe064ec01f62710b1c1dp-23Q, -0x1.6a87a6b931edc5ab5423dbf179adp-15Q, -0x1.0efaec65edbcc11b181ba70774bcp-14Q, 0x1.051eebbefaf7bd71954e297a6312p-8Q, -0x1.5cb5de436bfadb9e555d33c557e2p-6Q, 0x1.4e8938899ee371e40cd7b7528b62p-5Q, 0x1.4ab3cf0bae2a5b8432a0a181629ap-6Q, -0x1.145887eacb10f9d392eed69a935ap-2Q, 0x1.2d988a986b0dea756341b6b8b7e5p-1Q, 0x1.ca0692ca77439e8f752ca148c243p-2Q, -0x1.9ec434341eaed8ccc6e8d7244ab3p-1Q, -0x0p+0Q};
            }
            if constexpr (order == 2) {
                return {0x0p+0Q, 0x1.cb0ccbca08e7961d7bc2c289a48ap-25Q, 0x1.c4a78957d279804ccd5115ffad85p-16Q, -0x1.4e09cf12f0cdcfd873bae55853dap-9Q, 0x1.eb4cc917dd1cd56dce6a927f91d9p-12Q, 0x1.fe16d1a3cfb3e7ecc6834f20ab6ap-4Q, -0x1.18b244e50df675a487d76601d8f3p-1Q, 0x1.9d599937b31fafd3d66db2bf5b07p-1Q, 0x1.5709310e8d1f635f77e733f3179cp-2Q, -0x1.86bc52429f3a3f79bb4f1cb7604p+1Q, 0x1.70aa17b7adc0e75a9ffe0b5a59a9p+2Q, -0x1.5362d7fe36c56b6a443e9d61610ap+2Q, 0x1.e0eacca617ebb860573873f7d92p+0Q, 0x0p+0Q};
            }
            if constexpr (order == 3) {
                return {0x0p+0Q, 0x1.64d19861a833792a3ce80e60f903p-14Q, 0x1.6365af2495a86c5f48d9bcf9d227p-6Q, -0x1.13b5aaaf6c71d5b914583e2decb4p+0Q, -0x1.833e6f5afaca4c566ce5982603f2p-4Q, 0x1.a43b4147c122c6f143d7b566fba5p+4Q, -0x1.7ff7929036ddd78bf54cfc10d281p+6Q, 0x1.415156a34428a6a6c08cde5d7d99p+7Q, -0x1.23bcc34d6a529b0b422faf4a87d6p+7Q, 0x1.1e704d84a10750732762c3ddc376p+6Q, -0x1.171de5ca88cec2afbea5499d7e79p+4Q, 0x1.492e70f2635c16c262e61993b979p+1Q, -0x1.4ca9343b3236825367445eb4699fp-1Q, 0x0p+0Q};
            }
            if constexpr (order == 4) {
                return {0x0p+0Q, -0x1.1f66b1e7815295f48755fd6ea2d6p-14Q, -0x1.23f9756392d0ee9915448c215bcp-7Q, 0x1.7c490d447fad2e27fa6b12a68ad2p-3Q, -0x1.3be427404fe7d12432d600fa5a8dp-3Q, -0x1.ffd4fb5b391804449296747a79cbp+0Q, 0x1.8c2cccc6ee34ca0f2445838dfc51p+2Q, -0x1.222479cd96ee9d17f59bff50c277p+3Q, 0x1.9815cb7223c105bc2b63ed061cb6p+3Q, -0x1.40c582473fa8c93703a121712545p+4Q, 0x1.578437542fe4f5fbe4cefecbb8f5p+4Q, -0x1.7eda4b6a23328f4e2cec0994502dp+3Q, 0x1.525f0cb11c2720485c66b9e5d785p+1Q, -0x0p+0Q};
            }
            if constexpr (order == 5) {
                return {0x0p+0Q, 0x1.c9a573f22c19a9032c5430d094d8p-12Q, 0x1.e323c5b92bd166e33ab5b4e1191p-6Q, -0x1.16c54a229aae12ed509d34005188p-2Q, 0x1.9df5b4fe2e271a26f87f7a0d3b0bp-2Q, 0x1.97de12ea03ca10913f6e6b159efdp-1Q, -0x1.2d09c9846b581b01a4323f449897p+1Q, 0x1.34da807d2aaec910bcf3aa745ce6p+1Q, -0x1.4010d1eedca99cac8148c112fa0fp+2Q, 0x1.74c81d1a27df5c33454de218b7aap+3Q, -0x1.ac8b3ccd718c41df56b9c15d3ff2p+3Q, 0x1.ce7e6827143677fbd12dabd4f225p+2Q, -0x1.80a316d21a0eed480c84ea67c1c5p+0Q, 0x0p+0Q};
            }
            if constexpr (order == 6) {
                return {0x0p+0Q, 0x1.24061e2653bdb3c5fd153e85b015p-9Q, 0x1.4b88409bdd1172b9682df89a8a49p-4Q, -0x1.6569331439297d46e21a0662eb81p-2Q, 0x1.256afda077267e050d8164e36595p-1Q, -0x1.6af581a622d416ac06a3c40e75d6p-1Q, 0x1.3544c2755b1d82eac792a5881e19p+0Q, -0x1.e6f030cb23fbfde6dfddcd38ea34p+0Q, 0x1.2a57634989c51fb583fc11fe33cfp+1Q, -0x1.40ce386ca776e06ba6f31ba659f6p+1Q, 0x1.099a073a7c7b41898468a5b9ce7dp+1Q, -0x1.01d8c748e11fa735da05e73bbc5ap+0Q, 0x1.a096c8f7fc8ab1575820ee7479dp-3Q, -0x0p+0Q};
            }
         }
#endif
     }
    if constexpr (p == 8) {
        if constexpr (std::is_same_v<Real, float>) {
            if constexpr (order == 0) {
                return {0x0p+0f, 0x1.569072p-45f, -0x1.f6ebd8p-33f, 0x1.66eb1ap-25f, 0x1.b237bcp-20f, 0x1.da6faap-16f, -0x1.92df02p-14f, 0x1.8430d2p-11f, -0x1.96b60cp-8f, 0x1.82b42ep-6f, -0x1.90a568p-5f, 0x1.234876p-4f, -0x1.596f5cp-3f, 0x1.fbaa04p-1f, 0x1.184a0cp-3f, 0x0p+0f};
            }
            if constexpr (order == 1) {
                return {0x0p+0f, 0x1.06b5eep-37f, -0x1.81519ep-26f, -0x1.c99c42p-23f, 0x1.1649ecp-14f, -0x1.e972ap-12f, -0x1.27c566p-9f, 0x1.594cfap-6f, -0x1.7e9e8ep-5f, -0x1.458036p-6f, 0x1.4be40ep-2f, -0x1.c3013ep-1f, 0x1.7d2f18p+0f, -0x1.8347dap-2f, -0x1.03034p-1f, 0x0p+0f};
            }
            if constexpr (order == 2) {
                return {0x0p+0f, 0x1.7bbd84p-32f, -0x1.15f344p-21f, -0x1.2c21cp-17f, 0x1.7e62aep-11f, -0x1.fc6302p-9f, -0x1.66849cp-7f, 0x1.b6728ep-4f, -0x1.d9c16ep-3f, 0x1.a23c1ap-6f, 0x1.891e36p-1f, -0x1.d62736p+0f, 0x1.a43ebap+1f, -0x1.c03746p+1f, 0x1.666d54p+0f, 0x0p+0f};
            }
            if constexpr (order == 3) {
                return {0x0p+0f, -0x1.364c0ap-20f, 0x1.c48022p-11f, 0x1.2fc65p-7f, -0x1.32a2c6p-1f, 0x1.d27f3cp+0f, 0x1.4b77b4p+2f, -0x1.34930ep+5f, 0x1.70d204p+6f, -0x1.d1869cp+6f, 0x1.38f35cp+6f, -0x1.e5bffcp+3f, -0x1.1601a8p+4f, 0x1.c1f7ep+3f, -0x1.b44c4ep+1f, -0x0p+0f};
            }
            if constexpr (order == 4) {
                return {0x0p+0f, -0x1.a6f024p-19f, 0x1.3201b6p-10f, 0x1.0b828ap-6f, -0x1.6566f6p-2f, 0x1.0ea6cap+0f, 0x1.050e3ep-1f, -0x1.ee4b3cp+2f, 0x1.0e28c8p+4f, -0x1.6a2c2p+4f, 0x1.da6956p+4f, -0x1.306ff2p+5f, 0x1.0f6c04p+5f, -0x1.0a8c44p+4f, 0x1.afe05cp+1f, 0x0p+0f};
            }
            if constexpr (order == 5) {
                return {0x0p+0f, 0x1.e6e9a6p-15f, -0x1.5ad37ap-7f, -0x1.049fbp-3f, 0x1.4ac5cap+0f, -0x1.d0d0a8p+1f, 0x1.9f39e8p+1f, 0x1.17ffcap+1f, -0x1.6dba5cp+2f, 0x1.af306ap+2f, -0x1.233374p+4f, 0x1.164dap+5f, -0x1.0a779ep+5f, 0x1.f5abb2p+3f, -0x1.77885ep+1f, 0x0p+0f};
            }
            if constexpr (order == 6) {
                return {0x0p+0f, 0x1.5c5b72p-13f, -0x1.e09e82p-7f, -0x1.ddc83ep-4f, 0x1.806826p-1f, -0x1.a6996ap+0f, 0x1.5e9b74p+1f, -0x1.3fe312p+2f, 0x1.dcba52p+2f, -0x1.20f516p+3f, 0x1.92a29ep+3f, -0x1.076d5p+4f, 0x1.b7d754p+3f, -0x1.85c00cp+2f, 0x1.1912b6p+0f, 0x0p+0f};
            }
            if constexpr (order == 7) {
                return {0x0p+0f, -0x1.5d1ed8p-14f, 0x1.c24ed2p-9f, 0x1.03531cp-6f, -0x1.41c46cp-4f, 0x1.c2f7cep-4f, -0x1.09dce6p-2f, 0x1.7a4f7p-1f, -0x1.ff103cp-1f, 0x1.025406p+0f, -0x1.ffefb8p+0f, 0x1.b139cp+1f, -0x1.851aep+1f, 0x1.5a965ap+0f, -0x1.eae802p-3f, -0x0p+0f};
            }
         }
        if constexpr (std::is_same_v<Real, double>) {
            if constexpr (order == 0) {
                return {0x0p+0, 0x1.5690727034a8ap-45, -0x1.f6ebd7398c4a4p-33, 0x1.66eb1a28068aep-25, 0x1.b237bc2c3b525p-20, 0x1.da6faa7679627p-16, -0x1.92df023007909p-14, 0x1.8430d19335685p-11, -0x1.96b60b10e5956p-8, 0x1.82b42ea0f8a16p-6, -0x1.90a568ef3b94bp-5, 0x1.234875f46000ep-4, -0x1.596f5b6b28cbp-3, 0x1.fbaa047a94a61p-1, 0x1.184a0c288d22bp-3, 0x0p+0};
            }
            if constexpr (order == 1) {
                return {0x0p+0, 0x1.06b5ee4f0712p-37, -0x1.81519dee82fp-26, -0x1.c99c42a59648bp-23, 0x1.1649ec5e54141p-14, -0x1.e972a09b0f4ddp-12, -0x1.27c56651be6eep-9, 0x1.594cf912461bcp-6, -0x1.7e9e8eefc30b4p-5, -0x1.4580356ca34ffp-6, 0x1.4be40ed86f2cdp-2, -0x1.c3013eaebde94p-1, 0x1.7d2f17f2a95afp+0, -0x1.8347da7a3030cp-2, -0x1.03033f41d9947p-1, 0x0p+0};
            }
            if constexpr (order == 2) {
                return {0x0p+0, 0x1.7bbd836b1bf85p-32, -0x1.15f344e8ed1b2p-21, -0x1.2c21c0130a3ebp-17, 0x1.7e62ae52149bdp-11, -0x1.fc6301aedca9p-9, -0x1.66849bebd2ebep-7, 0x1.b6728d569b4c7p-4, -0x1.d9c16d05d7787p-3, 0x1.a23c19ecf2debp-6, 0x1.891e36789053ap-1, -0x1.d627359f84d25p+0, 0x1.a43eb9f4500c3p+1, -0x1.c0374596ee4d9p+1, 0x1.666d532df9dfp+0, 0x0p+0};
            }
            if constexpr (order == 3) {
                return {0x0p+0, -0x1.364c0a7c4c64bp-20, 0x1.c48021ad9bc16p-11, 0x1.2fc64ff8d3421p-7, -0x1.32a2c5a3dfe81p-1, 0x1.d27f3b25c1fd7p+0, 0x1.4b77b4fb36749p+2, -0x1.34930dde67aeap+5, 0x1.70d2038145a54p+6, -0x1.d1869b3f8e2b4p+6, 0x1.38f35bf0e9ee1p+6, -0x1.e5bffb29b3c53p+3, -0x1.1601a8a979c79p+4, 0x1.c1f7df3ad7ab8p+3, -0x1.b44c4e267cbb5p+1, -0x0p+0};
            }
            if constexpr (order == 4) {
                return {0x0p+0, -0x1.a6f0242017b15p-19, 0x1.3201b5c3bc1c2p-10, 0x1.0b828ad29e43fp-6, -0x1.6566f5ac7d805p-2, 0x1.0ea6cad7a5a8p+0, 0x1.050e3e4d4732ap-1, -0x1.ee4b3cc46d89p+2, 0x1.0e28c771588c3p+4, -0x1.6a2c1f828fe2ep+4, 0x1.da69559cc9113p+4, -0x1.306ff1c62c833p+5, 0x1.0f6c0325c5e14p+5, -0x1.0a8c44cdd47bfp+4, 0x1.afe05c7a6ada2p+1, 0x0p+0};
            }
            if constexpr (order == 5) {
                return {0x0p+0, 0x1.e6e9a51e3663ap-15, -0x1.5ad3792bdac36p-7, -0x1.049faf6ae6899p-3, 0x1.4ac5cadf9fb15p+0, -0x1.d0d0a8f3b805ep+1, 0x1.9f39e7400cf5dp+1, 0x1.17ffc928dd1d9p+1, -0x1.6dba5c59397a4p+2, 0x1.af306993826a7p+2, -0x1.2333738e736f1p+4, 0x1.164d9f148bde5p+5, -0x1.0a779dee6c7bep+5, 0x1.f5abb22b1b6ecp+3, -0x1.77885e6e270a7p+1, 0x0p+0};
            }
            if constexpr (order == 6) {
                return {0x0p+0, 0x1.5c5b71613d71bp-13, -0x1.e09e81650c6b7p-7, -0x1.ddc83de4acfafp-4, 0x1.806826de1640bp-1, -0x1.a6996ae95f76fp+0, 0x1.5e9b74431b073p+1, -0x1.3fe31292805cp+2, 0x1.dcba519c323cfp+2, -0x1.20f516098425cp+3, 0x1.92a29d7ac0e1ep+3, -0x1.076d50e7f7851p+4, 0x1.b7d753e3c818ap+3, -0x1.85c00c80d7a19p+2, 0x1.1912b6ad90d43p+0, 0x0p+0};
            }
            if constexpr (order == 7) {
                return {0x0p+0, -0x1.5d1ed72a09587p-14, 0x1.c24ed1e998c51p-9, 0x1.03531c096423ap-6, -0x1.41c46b919c33fp-4, 0x1.c2f7ce5c8a37cp-4, -0x1.09dce627745a6p-2, 0x1.7a4f70cc2f5d8p-1, -0x1.ff103b26e087ep-1, 0x1.0254056f57104p+0, -0x1.ffefb70380437p+0, 0x1.b139c010d5d8fp+1, -0x1.851adfd49d1a7p+1, 0x1.5a965a7535366p+0, -0x1.eae80165a34adp-3, -0x0p+0};
            }
         }
        if constexpr (std::is_same_v<Real, long double>) {
            if constexpr (order == 0) {
                return {0x0p+0L, 0xa.b4839381a545259p-48L, -0xf.b75eb9cc6251d52p-36L, 0xb.3758d14034571f5p-28L, 0xd.91bde161da92bd9p-23L, 0xe.d37d53b3cb1391fp-19L, -0xc.96f811803c8446bp-17L, 0xc.21868c99ab4273bp-14L, -0xc.b5b058872cab3ebp-11L, 0xc.15a17507c50b09fp-9L, -0xc.852b4779dca57d6p-8L, 0x9.1a43afa30007147p-7L, -0xa.cb7adb594658341p-6L, 0xf.dd5023d4a530983p-4L, 0x8.c250614469156ecp-6L, 0x0p+0L};
            }
            if constexpr (order == 1) {
                return {0x0p+0L, 0x8.35af7278388fe4bp-40L, -0xc.0a8cef7417800e8p-29L, -0xe.4ce2152cb245585p-26L, 0x8.b24f62f2a0a04dcp-17L, -0xf.4b9504d87a6e47fp-15L, -0x9.3e2b328df376c3ap-12L, 0xa.ca67c89230ddcadp-9L, -0xb.f4f4777e185a2cbp-8L, -0xa.2c01ab651a7f8e8p-9L, 0xa.5f2076c37966861p-5L, -0xe.1809f575ef4a37dp-4L, 0xb.e978bf954ad7b4cp-3L, -0xc.1a3ed3d18185c5ap-5L, -0x8.1819fa0ecca3abdp-4L, 0x0p+0L};
            }
            if constexpr (order == 2) {
                return {0x0p+0L, 0xb.ddec1b58dfc26aap-35L, -0x8.af9a274768d9303p-24L, -0x9.610e009851f575p-20L, 0xb.f3157290a4debfp-14L, -0xf.e3180d76e5481d1p-12L, -0xb.3424df5e975efd8p-10L, 0xd.b3946ab4da635ap-7L, -0xe.ce0b682ebbc393bp-6L, 0xd.11e0cf6796f5a28p-9L, 0xc.48f1b3c4829d099p-4L, -0xe.b139acfc2692864p-3L, 0xd.21f5cfa280614c9p-2L, -0xe.01ba2cb7726c73bp-2L, 0xb.336a996fcef7c52p-3L, 0x0p+0L};
            }
            if constexpr (order == 3) {
                return {0x0p+0L, -0x9.b26053e2632577ap-23L, 0xe.24010d6cde0b23p-14L, 0x9.7e327fc69a10719p-10L, -0x9.95162d1eff408cfp-4L, 0xe.93f9d92e0feb674p-3L, 0xa.5bbda7d9b3a4ab2p-1L, -0x9.a4986ef33d75017p+2L, 0xb.86901c0a2d2a313p+3L, -0xe.8c34d9fc715a214p+3L, 0x9.c79adf874f70573p+3L, -0xf.2dffd94d9e295dcp+0L, -0x8.b00d454bce3c7c1p+1L, 0xe.0fbef9d6bd5bcc5p+0L, -0xd.a2627133e5da6a2p-2L, -0x0p+0L};
            }
            if constexpr (order == 4) {
                return {0x0p+0L, -0xd.37812100bd8a443p-22L, 0x9.900dae1de0e10ep-13L, 0x8.5c145694f21f58cp-9L, -0xb.2b37ad63ec0243dp-5L, 0x8.753656bd2d3fd71p-3L, 0x8.2871f26a39950e3p-4L, -0xf.7259e6236c483c5p-1L, 0x8.71463b8ac46182cp+1L, -0xb.5160fc147f16c08p+1L, 0xe.d34aace64889869p+1L, -0x9.837f8e316419b55p+2L, 0x8.7b60192e2f0a0a5p+2L, -0x8.5462266ea3dfa72p+1L, 0xd.7f02e3d356d0e5bp-2L, 0x0p+0L};
            }
            if constexpr (order == 5) {
                return {0x0p+0L, 0xf.374d28f1b31ceb4p-18L, -0xa.d69bc95ed61ae29p-10L, -0x8.24fd7b57344ca47p-6L, 0xa.562e56fcfd8a5e2p-3L, -0xe.8685479dc02ed7ep-2L, 0xc.f9cf3a0067ae91ap-2L, 0x8.bffe4946e8ec5f8p-2L, -0xb.6dd2e2c9cbd2213p-1L, 0xd.79834c9c135345p-1L, -0x9.199b9c739b78925p+1L, 0x8.b26cf8a45ef2551p+2L, -0x8.53bcef7363deeep+2L, 0xf.ad5d9158db75d9ap+0L, -0xb.bc42f3713853a26p-2L, 0x0p+0L};
            }
            if constexpr (order == 6) {
                return {0x0p+0L, 0xa.e2db8b09eb8dbb5p-16L, -0xf.04f40b28635ba7ep-10L, -0xe.ee41ef2567d7906p-7L, 0xc.034136f0b2057f9p-4L, -0xd.34cb574afbb77c5p-3L, 0xa.f4dba218d8396ddp-2L, -0x9.ff18949402e0231p-1L, 0xe.e5d28ce191e76fcp-1L, -0x9.07a8b04c212dd91p+0L, 0xc.9514ebd6070fp+0L, -0x8.3b6a873fbc285cap+1L, 0xd.beba9f1e40c4cbdp+0L, -0xc.2e006406bd0c7e6p-1L, 0x8.c895b56c86a1b02p-3L, 0x0p+0L};
            }
            if constexpr (order == 7) {
                return {0x0p+0L, -0xa.e8f6b9504ac38b3p-17L, 0xe.12768f4cc6289abp-12L, 0x8.1a98e04b211d3d4p-9L, -0xa.0e235c8ce19f70ep-7L, 0xe.17be72e451bdc7ep-7L, -0x8.4ee7313ba2d30a8p-5L, 0xb.d27b86617aebd37p-4L, -0xf.f881d937043f115p-4L, 0x8.12a02b7ab881d9bp-3L, -0xf.ff7db81c021b9bdp-3L, 0xd.89ce0086aec77acp-2L, -0xc.28d6fea4e8d3b6dp-2L, 0xa.d4b2d3a9a9b32a1p-3L, -0xf.57400b2d1a568f9p-6L, -0x0p+0L};
            }
         }
#ifdef BOOST_HAS_FLOAT128
        if constexpr (std::is_same_v<Real, boost::multiprecision::float128>) {
            if constexpr (order == 0) {
                return {0x0p+0Q, 0x1.5690727034a8a4b1397c221a594fp-45Q, -0x1.f6ebd7398c4a3aa441dfeeb3fd9ap-33Q, 0x1.66eb1a28068ae3ea1bd5d27e91aap-25Q, 0x1.b237bc2c3b5257b1d937078f68bp-20Q, 0x1.da6faa767962723e71f08f2e440bp-16Q, -0x1.92df0230079088d5ad074614b732p-14Q, 0x1.8430d19335684e7577fd3b699065p-11Q, -0x1.96b60b10e59567d6b2c6577bac45p-8Q, 0x1.82b42ea0f8a1613dc5ff9725bfffp-6Q, -0x1.90a568ef3b94afacf755834babcbp-5Q, 0x1.234875f46000e28dce7c3f342099p-4Q, -0x1.596f5b6b28cb0681bf78a3e54b9ap-3Q, 0x1.fbaa047a94a613069e6f550b4e14p-1Q, 0x1.184a0c288d22add7ea26104df301p-3Q, 0x0p+0Q};
            }
            if constexpr (order == 1) {
                return {0x0p+0Q, 0x1.06b5ee4f0711fc96e3467676dec9p-37Q, -0x1.81519dee82f001d0cfd2892e0a52p-26Q, -0x1.c99c42a59648ab096d412a77185dp-23Q, 0x1.1649ec5e541409b8483ac00a3d29p-14Q, -0x1.e972a09b0f4dc8fdb93114f9e2e9p-12Q, -0x1.27c56651be6ed873712071c18c64p-9Q, 0x1.594cf912461bb959724fb77ad207p-6Q, -0x1.7e9e8eefc30b45959dcb8ab52213p-5Q, -0x1.4580356ca34ff1cf3c18408b19e8p-6Q, 0x1.4be40ed86f2cd0c168ad2494909fp-2Q, -0x1.c3013eaebde946f9eaa42bce8a2dp-1Q, 0x1.7d2f17f2a95af6988be8eca123aep+0Q, -0x1.8347da7a3030b8b37cdc36ec3b2cp-2Q, -0x1.03033f41d994757a2b493442339p-1Q, 0x0p+0Q};
            }
            if constexpr (order == 2) {
                return {0x0p+0Q, 0x1.7bbd836b1bf84d53feb5483f20c3p-32Q, -0x1.15f344e8ed1b26069acbf2612d99p-21Q, -0x1.2c21c0130a3eaea0bd3721b16691p-17Q, 0x1.7e62ae52149bd7e06a18746f185p-11Q, -0x1.fc6301aedca903a123de4bcec794p-9Q, -0x1.66849bebd2ebdfafabb13137685ep-7Q, 0x1.b6728d569b4c6b3f03c9c20fdfacp-4Q, -0x1.d9c16d05d77872764189280b4337p-3Q, 0x1.a23c19ecf2deb44f4177f3e97885p-6Q, 0x1.891e36789053a132cb4397963d6ap-1Q, -0x1.d627359f84d250c8825852fe24b7p+0Q, 0x1.a43eb9f4500c29913bdae70f30b4p+1Q, -0x1.c0374596ee4d8e750b0b0f065cb5p+1Q, 0x1.666d532df9def8a39e69d4c9dc9p+0Q, 0x0p+0Q};
            }
            if constexpr (order == 3) {
                return {0x0p+0Q, -0x1.364c0a7c4c64aef451897574dab5p-20Q, 0x1.c48021ad9bc1646021c012b3c45bp-11Q, 0x1.2fc64ff8d3420e31f2a6dd144767p-7Q, -0x1.32a2c5a3dfe8119d379b71142e79p-1Q, 0x1.d27f3b25c1fd6ce83b2a8abfebc6p+0Q, 0x1.4b77b4fb3674956338abd14f8ae6p+2Q, -0x1.34930dde67aea02e0a703b473b02p+5Q, 0x1.70d2038145a54625eb25d83d60f8p+6Q, -0x1.d1869b3f8e2b44288688d2268c8p+6Q, 0x1.38f35bf0e9ee0ae5b177291e7059p+6Q, -0x1.e5bffb29b3c52bb8854979d60d4bp+3Q, -0x1.1601a8a979c78f82e398ce806f7ap+4Q, 0x1.c1f7df3ad7ab79897dce80bf1eccp+3Q, -0x1.b44c4e267cbb4d43745305e1a54p+1Q, -0x0p+0Q};
            }
            if constexpr (order == 4) {
                return {0x0p+0Q, -0x1.a6f0242017b1488503ceffcc93b6p-19Q, 0x1.3201b5c3bc1c21bf7086980a22e2p-10Q, 0x1.0b828ad29e43eb178f6b65a73c71p-6Q, -0x1.6566f5ac7d80487971d04a467cb1p-2Q, 0x1.0ea6cad7a5a7fae2eab7d73ce1d1p+0Q, 0x1.050e3e4d4732a1c55c0c97206d01p-1Q, -0x1.ee4b3cc46d890789fa8b15a770e5p+2Q, 0x1.0e28c771588c30581ce7070ea2d3p+4Q, -0x1.6a2c1f828fe2d81030f1ad46dd7p+4Q, 0x1.da69559cc91130d28dd4b6249ab5p+4Q, -0x1.306ff1c62c8336aad5c3f97b75abp+5Q, 0x1.0f6c0325c5e1414a9285f592915fp+5Q, -0x1.0a8c44cdd47bf4e495a53978c49dp+4Q, 0x1.afe05c7a6ada1cb6172616c12e8fp+1Q, 0x0p+0Q};
            }
            if constexpr (order == 5) {
                return {0x0p+0Q, 0x1.e6e9a51e36639d68c6db7d269d95p-15Q, -0x1.5ad3792bdac35c5163486dea174bp-7Q, -0x1.049faf6ae689948e66d5e5532967p-3Q, 0x1.4ac5cadf9fb14bc48382ff57eb1fp+0Q, -0x1.d0d0a8f3b805dafc2d2b54014ec9p+1Q, 0x1.9f39e7400cf5d234f227213340fbp+1Q, 0x1.17ffc928dd1d8bf0f209d0a7d892p+1Q, -0x1.6dba5c59397a4426e98586a614dap+2Q, 0x1.af306993826a689f75439948060dp+2Q, -0x1.2333738e736f124a0ec3de6cb61cp+4Q, 0x1.164d9f148bde4aa13b66fadbc689p+5Q, -0x1.0a779dee6c7bddbff2f0d8da401cp+5Q, 0x1.f5abb22b1b6ebb333fd9f15c1399p+3Q, -0x1.77885e6e270a744c86865505076ap+1Q, 0x0p+0Q};
            }
            if constexpr (order == 6) {
                return {0x0p+0Q, 0x1.5c5b71613d71b7692e7e5ecaf61ap-13Q, -0x1.e09e81650c6b74fbce41e496f9abp-7Q, -0x1.ddc83de4acfaf20be5fee4aab51ep-4Q, 0x1.806826de1640aff2fbfd8cd80665p-1Q, -0x1.a6996ae95f76ef890f12b58d56bp+0Q, 0x1.5e9b74431b072dba61e5268b88c5p+1Q, -0x1.3fe31292805c0462ca96adbfad48p+2Q, 0x1.dcba519c323cedf78dcb7708da67p+2Q, -0x1.20f516098425bb22664ce46332d4p+3Q, 0x1.92a29d7ac0e1e0002558aef201d1p+3Q, -0x1.076d50e7f7850b94ae780acea072p+4Q, 0x1.b7d753e3c818997a04ded8add4f3p+3Q, -0x1.85c00c80d7a18fcc693e469a28ffp+2Q, 0x1.1912b6ad90d43603284e2873fc2cp+0Q, 0x0p+0Q};
            }
            if constexpr (order == 7) {
                return {0x0p+0Q, -0x1.5d1ed72a09587166cbb2e5d97b8cp-14Q, 0x1.c24ed1e998c51356546e2489bdfp-9Q, 0x1.03531c096423a7a88f95e3f6a663p-6Q, -0x1.41c46b919c33ee1c6941ab004471p-4Q, 0x1.c2f7ce5c8a37b8fcabdaf25e3b67p-4Q, -0x1.09dce627745a61500e9471a5cdffp-2Q, 0x1.7a4f70cc2f5d7a6d2a4e62d3effp-1Q, -0x1.ff103b26e087e22aa6dbe27326ap-1Q, 0x1.0254056f57103b364cf4bcbd55a7p+0Q, -0x1.ffefb7038043737a18fc835fe278p+0Q, 0x1.b139c010d5d8ef58fb69c524b48p+1Q, -0x1.851adfd49d1a76d9d23c35621b66p+1Q, 0x1.5a965a7535366541f0a1b0109262p+0Q, -0x1.eae80165a34ad1f208010834839fp-3Q, -0x0p+0Q};
            }
         }
#endif
     }
    if constexpr (p == 9) {
        if constexpr (std::is_same_v<Real, float>) {
            if constexpr (order == 0) {
                return {0x0p+0f, -0x1.4fb4eep-52f, -0x1.705a7ap-38f, 0x1.278ba2p-29f, -0x1.58c1d8p-24f, -0x1.53f94ap-19f, 0x1.7e6498p-16f, 0x1.b33accp-15f, -0x1.291474p-10f, 0x1.3d3ef8p-8f, -0x1.097588p-7f, -0x1.7ee87p-8f, 0x1.fd1878p-5f, -0x1.4aef16p-3f, 0x1.8b08d6p-3f, 0x1.b1b652p-1f, 0x1.1d2fcp-4f, 0x0p+0f};
            }
            if constexpr (order == 1) {
                return {0x0p+0f, -0x1.08644ap-47f, -0x1.2234ep-34f, 0x1.13cc2ep-24f, -0x1.8562bp-23f, -0x1.9d3f0ap-16f, 0x1.b592fcp-12f, -0x1.0c494ap-11f, -0x1.96a364p-8f, 0x1.49841cp-6f, 0x1.540e9cp-7f, -0x1.82c926p-3f, 0x1.1f4f38p-1f, -0x1.169004p+0f, 0x1.d3c16cp+0f, -0x1.af96b8p-1f, -0x1.2b90f8p-2f, -0x0p+0f};
            }
            if constexpr (order == 2) {
                return {0x0p+0f, 0x1.1a9906p-40f, 0x1.36693ap-28f, -0x1.261d8ap-20f, 0x1.bca3f2p-17f, 0x1.b3d73ep-12f, -0x1.2bb318p-9f, -0x1.85c5f4p-10f, 0x1.723502p-5f, -0x1.66761ap-3f, 0x1.61aa48p-2f, -0x1.4f00eap-2f, -0x1.c92ba8p-4f, 0x1.5b4146p-1f, 0x1.135c58p-2f, -0x1.aad7p+0f, 0x1.e4c00ep-1f, 0x0p+0f};
            }
            if constexpr (order == 3) {
                return {0x0p+0f, 0x1.cc794ap-31f, 0x1.fa8378p-20f, -0x1.ffd29p-12f, 0x1.20ddep-10f, 0x1.73f12ep-4f, -0x1.fe4766p-2f, 0x1.9b8564p-2f, 0x1.e1a25p+1f, -0x1.d3cfacp+3f, 0x1.985f02p+4f, -0x1.7ce73p+4f, 0x1.effd26p+2f, 0x1.304a1ap+3f, -0x1.e79fd8p+3f, 0x1.2ded7ap+3f, -0x1.27d6ep+1f, 0x0p+0f};
            }
            if constexpr (order == 4) {
                return {0x0p+0f, 0x1.43e712p-26f, 0x1.654cep-16f, -0x1.22b4e2p-8f, -0x1.2e2904p-8f, 0x1.b42324p-2f, -0x1.023562p+1f, 0x1.69d6eap+1f, 0x1.d829e2p+1f, -0x1.35c8bp+4f, 0x1.12ea48p+5f, -0x1.541e9cp+5f, 0x1.8eeb18p+5f, -0x1.af3bfap+5f, 0x1.4df684p+5f, -0x1.2ab8cap+4f, 0x1.c7cedp+1f, -0x0p+0f};
            }
            if constexpr (order == 5) {
                return {0x0p+0f, -0x1.398368p-20f, -0x1.5bccdep-11f, 0x1.491f18p-4f, 0x1.2ade5ep-3f, -0x1.f4375ep+1f, 0x1.eccebap+3f, -0x1.b958d8p+4f, 0x1.88dca4p+4f, -0x1.669c08p+3f, 0x1.c94ed2p+2f, 0x1.016a68p+0f, -0x1.24c19cp+5f, 0x1.168c04p+6f, -0x1.dd0eb6p+5f, 0x1.94d9c6p+4f, -0x1.1617a4p+2f, 0x0p+0f};
            }
            if constexpr (order == 6) {
                return {0x0p+0f, -0x1.202646p-18f, -0x1.434438p-10f, 0x1.ebd4ap-5f, 0x1.75a7acp-5f, -0x1.6df696p+0f, 0x1.1eeca8p+2f, -0x1.050826p+3f, 0x1.ca4482p+3f, -0x1.815eep+4f, 0x1.f5b3e2p+4f, -0x1.25c844p+5f, 0x1.7d0e6cp+5f, -0x1.a86d9ep+5f, 0x1.30352ap+5f, -0x1.dd3c5ap+3f, 0x1.384148p+1f, -0x0p+0f};
            }
            if constexpr (order == 7) {
                return {0x0p+0f, 0x1.a81802p-17f, 0x1.e662fap-10f, -0x1.26fbfep-5f, 0x1.2e59eap-6f, 0x1.794e94p-2f, -0x1.d1d5cp-1f, 0x1.6bf238p+0f, -0x1.d31828p+1f, 0x1.d4ebfcp+2f, -0x1.12f0bcp+3f, 0x1.3ddf52p+3f, -0x1.0b593cp+4f, 0x1.5c5edep+4f, -0x1.04a7fep+4f, 0x1.96fe7ap+2f, -0x1.041f76p+0f, 0x0p+0f};
            }
            if constexpr (order == 8) {
                return {0x0p+0f, -0x1.e260bp-16f, -0x1.20ae94p-9f, 0x1.2926b8p-6f, -0x1.8ce86ep-6f, -0x1.481726p-5f, 0x1.054d18p-4f, -0x1.b7954cp-6f, 0x1.b5ef82p-2f, -0x1.1620dap+0f, 0x1.028beep+0f, -0x1.0bd63p+0f, 0x1.51b6cap+1f, -0x1.033a7ap+2f, 0x1.947e64p+1f, -0x1.3caa56p+0f, 0x1.9024ep-3f, 0x0p+0f};
            }
         }
        if constexpr (std::is_same_v<Real, double>) {
            if constexpr (order == 0) {
                return {0x0p+0, -0x1.4fb4ee678fa0ep-52, -0x1.705a7a4d8345p-38, 0x1.278ba12dc249ap-29, -0x1.58c1d88e09b31p-24, -0x1.53f9494edf128p-19, 0x1.7e6498cdbddd5p-16, 0x1.b33acccf47f1cp-15, -0x1.29147421d5711p-10, 0x1.3d3ef7388a3aep-8, -0x1.09758849af605p-7, -0x1.7ee87092fa243p-8, 0x1.fd1877f1138bep-5, -0x1.4aef16c7e666fp-3, 0x1.8b08d5cd4a2adp-3, 0x1.b1b6520101e0ep-1, 0x1.1d2fbfac2394ap-4, 0x0p+0};
            }
            if constexpr (order == 1) {
                return {0x0p+0, -0x1.08644a120d1efp-47, -0x1.2234e02ed8631p-34, 0x1.13cc2e0d34f25p-24, -0x1.8562af0388b7ep-23, -0x1.9d3f09be0cea7p-16, 0x1.b592fbb7be41bp-12, -0x1.0c49493de5a6fp-11, -0x1.96a364065f27dp-8, 0x1.49841c554fed2p-6, 0x1.540e9bea7c732p-7, -0x1.82c9256a70f8fp-3, 0x1.1f4f37268ff74p-1, -0x1.1690032de8cfbp+0, 0x1.d3c16b16ba87ap+0, -0x1.af96b810375aep-1, -0x1.2b90f77070093p-2, -0x0p+0};
            }
            if constexpr (order == 2) {
                return {0x0p+0, 0x1.1a9906065291fp-40, 0x1.36693a4632269p-28, -0x1.261d8947cb727p-20, 0x1.bca3f26f194dp-17, 0x1.b3d73e818e5c9p-12, -0x1.2bb318efa4d23p-9, -0x1.85c5f4a449273p-10, 0x1.7235017b58d32p-5, -0x1.6676194262896p-3, 0x1.61aa479f9dc68p-2, -0x1.4f00ea2c2d363p-2, -0x1.c92ba77b267a3p-4, 0x1.5b41454ff0957p-1, 0x1.135c57711e671p-2, -0x1.aad6ffd5d77f1p+0, 0x1.e4c00e361fc68p-1, 0x0p+0};
            }
            if constexpr (order == 3) {
                return {0x0p+0, 0x1.cc794a76e5f16p-31, 0x1.fa83788d3f78dp-20, -0x1.ffd28fb3400bfp-12, 0x1.20dddf659ee51p-10, 0x1.73f12d4155f36p-4, -0x1.fe4766e3faf54p-2, 0x1.9b8563bd7c7d8p-2, 0x1.e1a24f8472a68p+1, -0x1.d3cfac07a9d2bp+3, 0x1.985f0196c72e5p+4, -0x1.7ce72fdc2cb91p+4, 0x1.effd2573aa9abp+2, 0x1.304a1a210a5fap+3, -0x1.e79fd7c0b1c0fp+3, 0x1.2ded79e41fdb1p+3, -0x1.27d6df36685a8p+1, 0x0p+0};
            }
            if constexpr (order == 4) {
                return {0x0p+0, 0x1.43e711e678672p-26, 0x1.654cdf573874bp-16, -0x1.22b4e2598d44bp-8, -0x1.2e290398a3fe3p-8, 0x1.b423243b8b469p-2, -0x1.0235621854eefp+1, 0x1.69d6e9c20f73dp+1, 0x1.d829e1a762882p+1, -0x1.35c8b07214e14p+4, 0x1.12ea47c60decp+5, -0x1.541e9c27826e5p+5, 0x1.8eeb17bed133cp+5, -0x1.af3bf949d62bap+5, 0x1.4df684dede46dp+5, -0x1.2ab8c9c81252ep+4, 0x1.c7ced0f66aeadp+1, -0x0p+0};
            }
            if constexpr (order == 5) {
                return {0x0p+0, -0x1.398367062ce37p-20, -0x1.5bccddec65e54p-11, 0x1.491f178ce5646p-4, 0x1.2ade5d60c4581p-3, -0x1.f4375d664ce34p+1, 0x1.eccebab605f4bp+3, -0x1.b958d77e385a9p+4, 0x1.88dca456f3456p+4, -0x1.669c08c9522f3p+3, 0x1.c94ed28aba0f5p+2, 0x1.016a68491e4d8p+0, -0x1.24c19c047007fp+5, 0x1.168c0447c1b6fp+6, -0x1.dd0eb6aa76ffep+5, 0x1.94d9c534a28d4p+4, -0x1.1617a4f9d8857p+2, 0x0p+0};
            }
            if constexpr (order == 6) {
                return {0x0p+0, -0x1.20264586c6bb9p-18, -0x1.4344372e1b6fcp-10, 0x1.ebd49f3810c89p-5, 0x1.75a7ac677c4efp-5, -0x1.6df69657c972dp+0, 0x1.1eeca8acc3a79p+2, -0x1.050825dec75a9p+3, 0x1.ca4482f4236a5p+3, -0x1.815ee05c85e8bp+4, 0x1.f5b3e1ac99667p+4, -0x1.25c8439d2df8ap+5, 0x1.7d0e6c4e24171p+5, -0x1.a86d9ebb6edc4p+5, 0x1.303529a28e47p+5, -0x1.dd3c5aa19c52ep+3, 0x1.3841475e9a364p+1, -0x0p+0};
            }
            if constexpr (order == 7) {
                return {0x0p+0, 0x1.a81802ebf94f1p-17, 0x1.e662faa96110fp-10, -0x1.26fbfe966bb82p-5, 0x1.2e59e90e0a11dp-6, 0x1.794e93c6506ddp-2, -0x1.d1d5c0ac4dc0cp-1, 0x1.6bf237312f18dp+0, -0x1.d31827f76dc61p+1, 0x1.d4ebfbf56d996p+2, -0x1.12f0bcbe3fa44p+3, 0x1.3ddf51285055ap+3, -0x1.0b593b29c51bcp+4, 0x1.5c5ede521aae6p+4, -0x1.04a7fd5b5e592p+4, 0x1.96fe795439509p+2, -0x1.041f759e862aap+0, 0x0p+0};
            }
            if constexpr (order == 8) {
                return {0x0p+0, -0x1.e260b0ee2a093p-16, -0x1.20ae937e7a93p-9, 0x1.2926b8981e2a7p-6, -0x1.8ce86e181d425p-6, -0x1.4817258dd2847p-5, 0x1.054d18d1cc90cp-4, -0x1.b7954b06be8afp-6, 0x1.b5ef819dccd9ap-2, -0x1.1620d9bfc03ebp+0, 0x1.028bee5b6619fp+0, -0x1.0bd62fc1e0455p+0, 0x1.51b6c9e0d0443p+1, -0x1.033a79f260cefp+2, 0x1.947e636a9e374p+1, -0x1.3caa56e1c68fap+0, 0x1.9024e0f887aa3p-3, 0x0p+0};
            }
         }
        if constexpr (std::is_same_v<Real, long double>) {
            if constexpr (order == 0) {
                return {0x0p+0L, -0xa.7da7733c7d070cbp-55L, -0xb.82d3d26c1a28344p-41L, 0x9.3c5d096e124d267p-32L, -0xa.c60ec4704d98582p-27L, -0xa.9fca4a76f894114p-22L, 0xb.f324c66deeea6e7p-19L, 0xd.99d6667a3f8ded1p-18L, -0x9.48a3a10eab88924p-13L, 0x9.e9f7b9c451d6f5p-11L, -0x8.4bac424d7b02898p-10L, -0xb.f7438497d12158fp-11L, 0xf.e8c3bf889c5f1c6p-8L, -0xa.5778b63f3337b8dp-6L, 0xc.5846ae6a5156a59p-6L, 0xd.8db290080f07364p-4L, 0x8.e97dfd611ca4d0fp-7L, 0x0p+0L};
            }
            if constexpr (order == 1) {
                return {0x0p+0L, -0x8.4322509068f7b79p-50L, -0x9.11a70176c31860ap-37L, 0x8.9e617069a792bafp-27L, -0xc.2b15781c45bef81p-26L, -0xc.e9f84df067536c7p-19L, 0xd.ac97ddbdf20daf5p-15L, -0x8.624a49ef2d3764bp-14L, -0xc.b51b2032f93e669p-11L, 0xa.4c20e2aa7f68ebp-9L, 0xa.a074df53e3990b8p-10L, -0xc.16492b5387c74b6p-6L, 0x8.fa79b9347fb9e82p-4L, -0x8.b480196f467d6c9p-3L, 0xe.9e0b58b5d43cc82p-3L, -0xd.7cb5c081bad7186p-4L, -0x9.5c87bb838049b9dp-5L, -0x0p+0L};
            }
            if constexpr (order == 2) {
                return {0x0p+0L, 0x8.d4c83032948f76p-43L, 0x9.b349d2319134721p-31L, -0x9.30ec4a3e5b93bbap-23L, 0xd.e51f9378ca67c3cp-20L, 0xd.9eb9f40c72e4aacp-15L, -0x9.5d98c77d2691ac4p-12L, -0xc.2e2fa52249396c6p-13L, 0xb.91a80bdac698e63p-8L, -0xb.33b0ca13144b196p-6L, 0xb.0d523cfcee33daap-5L, -0xa.7807516169b18a5p-5L, -0xe.495d3bd933d14bbp-7L, 0xa.da0a2a7f84ab917p-4L, 0x8.9ae2bb88f3384cep-5L, -0xd.56b7feaebbf8a9ep-3L, 0xf.260071b0fe33d17p-4L, 0x0p+0L};
            }
            if constexpr (order == 3) {
                return {0x0p+0L, 0xe.63ca53b72f8b2fep-34L, 0xf.d41bc469fbc68cdp-23L, -0xf.fe947d9a005fa2fp-15L, 0x9.06eefb2cf7284e3p-13L, 0xb.9f896a0aaf9ad9fp-7L, -0xf.f23b371fd7a9d87p-5L, 0xc.dc2b1debe3ebee3p-5L, 0xf.0d127c23953405cp-2L, -0xe.9e7d603d4e95463p+0L, 0xc.c2f80cb639729f3p+1L, -0xb.e7397ee165c88c5p+1L, 0xf.7fe92b9d54d57e6p-1L, 0x9.8250d10852fce74p+0L, -0xf.3cfebe058e075ecp+0L, 0x9.6f6bcf20fed8622p+0L, -0x9.3eb6f9b342d3cadp-2L, 0x0p+0L};
            }
            if constexpr (order == 4) {
                return {0x0p+0L, 0xa.1f388f33c339298p-29L, 0xb.2a66fab9c3a5b39p-19L, -0x9.15a712cc6a2586fp-11L, -0x9.71481cc51ff17e4p-11L, 0xd.a11921dc5a34834p-5L, -0x8.11ab10c2a7777ccp-2L, 0xb.4eb74e107b9e597p-2L, 0xe.c14f0d3b14410dep-2L, -0x9.ae458390a70a077p+1L, 0x8.97523e306f6022p+2L, -0xa.a0f4e13c13724e9p+2L, 0xc.7758bdf6899e1b4p+2L, -0xd.79dfca4eb15ce04p+2L, 0xa.6fb426f6f2369dp+2L, -0x9.55c64e409297236p+1L, 0xe.3e7687b3575667ep-2L, -0x0p+0L};
            }
            if constexpr (order == 5) {
                return {0x0p+0L, -0x9.cc1b3831671b9dbp-23L, -0xa.de66ef632f29e0ep-14L, 0xa.48f8bc672b2306bp-7L, 0x9.56f2eb0622c089cp-6L, -0xf.a1baeb32671a172p-2L, 0xf.6675d5b02fa596cp+0L, -0xd.cac6bbf1c2d48b2p+1L, 0xc.46e522b79a2b2aep+1L, -0xb.34e0464a91798aep+0L, 0xe.4a769455d07a8e6p-1L, 0x8.0b534248f26c38ep-3L, -0x9.260ce023803fb4p+2L, 0x8.b460223e0db7486p+3L, -0xe.e875b553b7fed21p+2L, 0xc.a6ce29a51469f7p+1L, -0x8.b0bd27cec42baaep-1L, 0x0p+0L};
            }
            if constexpr (order == 6) {
                return {0x0p+0L, -0x9.01322c3635dc8d8p-21L, -0xa.1a21b970db7e3c3p-13L, 0xf.5ea4f9c08644463p-8L, 0xb.ad3d633be277714p-8L, -0xb.6fb4b2be4b96476p-3L, 0x8.f76545661d3c9p-1L, -0x8.28412ef63ad4b3p+0L, 0xe.522417a11b52a8bp+0L, -0xc.0af702e42f45538p+1L, 0xf.ad9f0d64cb33437p+1L, -0x9.2e421ce96fc4fd3p+2L, 0xb.e873627120b85ap+2L, -0xd.436cf5db76e1e52p+2L, 0x9.81a94d147237d62p+2L, -0xe.e9e2d50ce2972e4p+0L, 0x9.c20a3af4d1b1cfep-2L, -0x0p+0L};
            }
            if constexpr (order == 7) {
                return {0x0p+0L, 0xd.40c0175fca78546p-20L, 0xf.3317d54b0887a0cp-13L, -0x9.37dff4b35dc1001p-8L, 0x9.72cf4870508e45cp-9L, 0xb.ca749e32836e5f1p-5L, -0xe.8eae05626e0638ap-4L, 0xb.5f91b98978c679ep-3L, -0xe.98c13fbb6e30aacp-2L, 0xe.a75fdfab6ccb19fp-1L, -0x8.9785e5f1fd2220dp+0L, 0x9.eefa894282acd81p+0L, -0x8.5ac9d94e28ddf7fp+1L, 0xa.e2f6f290d572f2cp+1L, -0x8.253feadaf2c91fcp+1L, 0xc.b7f3caa1ca84594p-1L, -0x8.20fbacf4315510dp-3L, 0x0p+0L};
            }
            if constexpr (order == 8) {
                return {0x0p+0L, -0xf.1305877150497c2p-19L, -0x9.05749bf3d4982dep-12L, 0x9.4935c4c0f1539cfp-9L, -0xc.674370c0ea12a57p-9L, -0xa.40b92c6e94234bfp-8L, 0x8.2a68c68e64863cdp-7L, -0xd.bcaa5835f4575cbp-9L, 0xd.af7c0cee66cd13ep-5L, -0x8.b106cdfe01f5b31p-3L, 0x8.145f72db30cf502p-3L, -0x8.5eb17e0f022abcep-3L, 0xa.8db64f068221badp-2L, -0x8.19d3cf93067786bp-1L, 0xc.a3f31b54f1b9fa4p-2L, -0x9.e552b70e347cfbdp-3L, 0xc.812707c43d519efp-6L, 0x0p+0L};
            }
         }
#ifdef BOOST_HAS_FLOAT128
        if constexpr (std::is_same_v<Real, boost::multiprecision::float128>) {
            if constexpr (order == 0) {
                return {0x0p+0Q, -0x1.4fb4ee678fa0e19563ce4d37eb3fp-52Q, -0x1.705a7a4d83450688cba31db7817ap-38Q, 0x1.278ba12dc249a4cd49c8d7384aa5p-29Q, -0x1.58c1d88e09b30b0364784f04696p-24Q, -0x1.53f9494edf1282277c1de6de33a9p-19Q, 0x1.7e6498cdbddd4dce209effcda2aep-16Q, 0x1.b33acccf47f1bda17f46b81050d6p-15Q, -0x1.29147421d57112477f4b7b34942ap-10Q, 0x1.3d3ef7388a3ade9f5adf4b3f5817p-8Q, -0x1.09758849af60512fff4fa0607472p-7Q, -0x1.7ee87092fa242b1ef789aaf8df8fp-8Q, 0x1.fd1877f1138be38ba53889d724e3p-5Q, -0x1.4aef16c7e666f719f7698e178b01p-3Q, 0x1.8b08d5cd4a2ad4b181244736a702p-3Q, 0x1.b1b6520101e0e6c7ddaa9da59f33p-1Q, 0x1.1d2fbfac23949a1e129ec9f8392bp-4Q, 0x0p+0Q};
            }
            if constexpr (order == 1) {
                return {0x0p+0Q, -0x1.08644a120d1ef6f18093a15576bap-47Q, -0x1.2234e02ed8630c149656b2def8bp-34Q, 0x1.13cc2e0d34f2575dd283afe47ecbp-24Q, -0x1.8562af0388b7df01536ccb5c1371p-23Q, -0x1.9d3f09be0cea6d8d6f8ed457630cp-16Q, 0x1.b592fbb7be41b5e92f003f6147a5p-12Q, -0x1.0c49493de5a6ec95d0bb4c3ad7c7p-11Q, -0x1.96a364065f27ccd19657d67cde79p-8Q, 0x1.49841c554fed1d5f2443e15a58dap-6Q, 0x1.540e9bea7c732170245209a305e8p-7Q, -0x1.82c9256a70f8e96c68d16f2e32b6p-3Q, 0x1.1f4f37268ff73d0393fe0f6778a7p-1Q, -0x1.1690032de8cfad91fdf041ea0e4dp+0Q, 0x1.d3c16b16ba8799042d28d9c49795p+0Q, -0x1.af96b810375ae30c4c41bf8ddd41p-1Q, -0x1.2b90f770700937399f9f525fd0fcp-2Q, -0x0p+0Q};
            }
            if constexpr (order == 2) {
                return {0x0p+0Q, 0x1.1a9906065291eec0968b7cac4e42p-40Q, 0x1.36693a4632268e426c400b7938c5p-28Q, -0x1.261d8947cb7277733f476eb4b60fp-20Q, 0x1.bca3f26f194cf877e2358cae53a2p-17Q, 0x1.b3d73e818e5c9557fa6b9520a6e8p-12Q, -0x1.2bb318efa4d23587f50785ea3544p-9Q, -0x1.85c5f4a449272d8bf24d90db31f3p-10Q, 0x1.7235017b58d31cc6ca2ab84f79b3p-5Q, -0x1.667619426289632c2e8907dc1701p-3Q, 0x1.61aa479f9dc67b53226c77100d73p-2Q, -0x1.4f00ea2c2d363149504e045fb194p-2Q, -0x1.c92ba77b267a297500d5b5d4aa96p-4Q, 0x1.5b41454ff095722d0f803ea696abp-1Q, 0x1.135c57711e67099c1758b9c6ceb3p-2Q, -0x1.aad6ffd5d77f153c17f1c88228d5p+0Q, 0x1.e4c00e361fc67a2e12b083b2d4c3p-1Q, 0x0p+0Q};
            }
            if constexpr (order == 3) {
                return {0x0p+0Q, 0x1.cc794a76e5f165fc80f3eb984aa5p-31Q, 0x1.fa83788d3f78d19a901b74b2a48fp-20Q, -0x1.ffd28fb3400bf45efe5c3293f7b1p-12Q, 0x1.20dddf659ee509c6befdf1cc3517p-10Q, 0x1.73f12d4155f35b3e8452b5046574p-4Q, -0x1.fe4766e3faf53b0d1248280f8decp-2Q, 0x1.9b8563bd7c7d7dc5768ba3e50dffp-2Q, 0x1.e1a24f8472a680b7d24d0587f57cp+1Q, -0x1.d3cfac07a9d2a8c51198c07b128ap+3Q, 0x1.985f0196c72e53e6d1cfdeee0001p+4Q, -0x1.7ce72fdc2cb91189ae6d5e0e8417p+4Q, 0x1.effd2573aa9aafcb53e588a97fc7p+2Q, 0x1.304a1a210a5f9ce8bc8b8a0f5e8ap+3Q, -0x1.e79fd7c0b1c0ebd84f46088acc8p+3Q, 0x1.2ded79e41fdb0c43d1210e15dfebp+3Q, -0x1.27d6df36685a79598b6c86d12618p+1Q, 0x0p+0Q};
            }
            if constexpr (order == 4) {
                return {0x0p+0Q, 0x1.43e711e67867252f78e381e8e437p-26Q, 0x1.654cdf573874b672ebfe35905f4ap-16Q, -0x1.22b4e2598d44b0dd2f4a2967b02p-8Q, -0x1.2e290398a3fe2fc8cc3aebedad9dp-8Q, 0x1.b423243b8b46906701ef9314e55ap-2Q, -0x1.0235621854eeef97bfd3d50112ep+1Q, 0x1.69d6e9c20f73cb2d14ff738015b4p+1Q, 0x1.d829e1a7628821bc837255c378c5p+1Q, -0x1.35c8b07214e140ed7f1d44817eap+4Q, 0x1.12ea47c60dec043f14d3665cc419p+5Q, -0x1.541e9c27826e49d24f626bc79289p+5Q, 0x1.8eeb17bed133c367f65a681c8be9p+5Q, -0x1.af3bf949d62b9c07b05b4aec5471p+5Q, 0x1.4df684dede46d3a014591feedf25p+5Q, -0x1.2ab8c9c81252e46b55777dfee134p+4Q, 0x1.c7ced0f66aeaccfcda1658b5d776p+1Q, -0x0p+0Q};
            }
            if constexpr (order == 5) {
                return {0x0p+0Q, -0x1.398367062ce373b681067c794d7dp-20Q, -0x1.5bccddec65e53c1b0416c18acb16p-11Q, 0x1.491f178ce56460d544cbdffd257ap-4Q, 0x1.2ade5d60c45811382f7ff059fe3p-3Q, -0x1.f4375d664ce342e4e518659c49d4p+1Q, 0x1.eccebab605f4b2d725512ba0bf84p+3Q, -0x1.b958d77e385a9163a2eb3993efa8p+4Q, 0x1.88dca456f345655b7cb433cd85afp+4Q, -0x1.669c08c9522f315c7b041cb6d4dfp+3Q, 0x1.c94ed28aba0f51cc92c5ae118691p+2Q, 0x1.016a68491e4d871ca505a5327f18p+0Q, -0x1.24c19c047007f67f0d2fc475ceefp+5Q, 0x1.168c0447c1b6e90c3306023a0aap+6Q, -0x1.dd0eb6aa76ffda42c87705c1f26bp+5Q, 0x1.94d9c534a28d3ee0a8bf204242dp+4Q, -0x1.1617a4f9d885755cd347727e7722p+2Q, 0x0p+0Q};
            }
            if constexpr (order == 6) {
                return {0x0p+0Q, -0x1.20264586c6bb91b0aa0b6674418dp-18Q, -0x1.4344372e1b6fc786ddf9128f53fap-10Q, 0x1.ebd49f3810c888c6985578d06806p-5Q, 0x1.75a7ac677c4eee273d03fcdd2f5ap-5Q, -0x1.6df69657c972c8ebec75cae83705p+0Q, 0x1.1eeca8acc3a791ffedf80d0fadefp+2Q, -0x1.050825dec75a966028608e102521p+3Q, 0x1.ca4482f4236a5516d7f63b4c9594p+3Q, -0x1.815ee05c85e8aa704180637eb6b9p+4Q, 0x1.f5b3e1ac9966686ed844c9edc602p+4Q, -0x1.25c8439d2df89fa50079075b2c13p+5Q, 0x1.7d0e6c4e24170b40a4b959f30c55p+5Q, -0x1.a86d9ebb6edc3ca33aa7b2b81134p+5Q, 0x1.303529a28e46fac3a8f57111632ep+5Q, -0x1.dd3c5aa19c52e5c8a273faa4ddb2p+3Q, 0x1.3841475e9a3639fcb27c2cf7df2cp+1Q, -0x0p+0Q};
            }
            if constexpr (order == 7) {
                return {0x0p+0Q, 0x1.a81802ebf94f0a8bc8124ba5e036p-17Q, 0x1.e662faa96110f417469cde711495p-10Q, -0x1.26fbfe966bb82001d50ae4ea75bap-5Q, 0x1.2e59e90e0a11c8b7ce9f899484f3p-6Q, 0x1.794e93c6506dcbe26714d6565f73p-2Q, -0x1.d1d5c0ac4dc0c7147f916ba5d3b8p-1Q, 0x1.6bf237312f18cf3b573ff9d4ed94p+0Q, -0x1.d31827f76dc61557ae5c1a4ad21fp+1Q, 0x1.d4ebfbf56d99633e38e5707e4431p+2Q, -0x1.12f0bcbe3fa44419a6b5206f9071p+3Q, 0x1.3ddf512850559b028d83069f2d3dp+3Q, -0x1.0b593b29c51bbefd0a21a19d1bbdp+4Q, 0x1.5c5ede521aae5e58438f334329cdp+4Q, -0x1.04a7fd5b5e5923f7f627ba79cfd1p+4Q, 0x1.96fe795439508b28db4260b16a3ep+2Q, -0x1.041f759e862aa2196b99e0413bccp+0Q, 0x0p+0Q};
            }
            if constexpr (order == 8) {
                return {0x0p+0Q, -0x1.e260b0ee2a092f8355b8822b4d71p-16Q, -0x1.20ae937e7a9305bbde132925fe3cp-9Q, 0x1.2926b8981e2a739eefaa205a9026p-6Q, -0x1.8ce86e181d4254ae37719dc2b819p-6Q, -0x1.4817258dd284697ebc00b7126938p-5Q, 0x1.054d18d1cc90c79a1ef034bd3abbp-4Q, -0x1.b7954b06be8aeb9619601416a6bdp-6Q, 0x1.b5ef819dccd9a27cb2716f5076a6p-2Q, -0x1.1620d9bfc03eb662d8d6552f2e0dp+0Q, 0x1.028bee5b6619ea0406791eb40e8cp+0Q, -0x1.0bd62fc1e045579c86e940cb206bp+0Q, 0x1.51b6c9e0d0443759a67a3cd024e6p+1Q, -0x1.033a79f260cef0d616212bb8ecc2p+2Q, 0x1.947e636a9e373f489450faabc98bp+1Q, -0x1.3caa56e1c68f9f7ad07fd75f8dc7p+0Q, 0x1.9024e0f887aa33deb80195eab194p-3Q, 0x0p+0Q};
            }
         }
#endif
     }
    if constexpr (p == 10) {
        if constexpr (std::is_same_v<Real, float>) {
            if constexpr (order == 0) {
                return {0x0p+0f, -0x1.098874p-59f, 0x1.afee7cp-44f, 0x1.db3584p-37f, -0x1.58dbb8p-28f, 0x1.c9437p-24f, 0x1.403906p-24f, -0x1.b37dfap-16f, 0x1.44e94ep-14f, 0x1.94f202p-12f, -0x1.2885ecp-9f, 0x1.d4ed16p-10f, 0x1.1d23ep-6f, -0x1.494bcp-4f, 0x1.9e3daep-3f, -0x1.85d2eep-2f, 0x1.1c4636p-1f, 0x1.4e2c28p-1f, 0x1.12cb0ap-5f, 0x0p+0f};
            }
            if constexpr (order == 1) {
                return {0x0p+0f, -0x1.85b824p-55f, 0x1.3ced9ep-40f, 0x1.51db7ap-33f, -0x1.f125fap-26f, -0x1.2eea96p-21f, -0x1.06ba9ap-17f, -0x1.ce036p-15f, 0x1.f37546p-12f, -0x1.2c315ap-10f, 0x1.213456p-8f, -0x1.c84ec2p-6f, 0x1.ac223p-4f, -0x1.eef0eap-3f, 0x1.a2cb68p-2f, -0x1.61584ep-1f, 0x1.94bc28p+0f, -0x1.f5f9c2p-1f, -0x1.43cfccp-3f, -0x0p+0f};
            }
            if constexpr (order == 2) {
                return {0x0p+0f, 0x1.4b9c06p-46f, -0x1.0d99c2p-32f, 0x1.3f418cp-28f, 0x1.bb4b6cp-19f, -0x1.2b4a08p-15f, -0x1.67961ap-12f, 0x1.076728p-8f, -0x1.51b88cp-7f, -0x1.35a482p-6f, 0x1.55c762p-3f, -0x1.9dcbe6p-2f, 0x1.86e852p-2f, 0x1.6e6d86p-2f, -0x1.c8684ep+0f, 0x1.8bd77ap+1f, -0x1.0fd09ap+1f, -0x1.01b00cp-2f, 0x1.2d6cb8p-1f, -0x0p+0f};
            }
            if constexpr (order == 3) {
                return {0x0p+0f, 0x1.cd96fcp-39f, -0x1.7713ap-26f, 0x1.03c426p-22f, 0x1.351528p-13f, -0x1.dd5da6p-11f, -0x1.32422ap-7f, 0x1.60b8bp-4f, -0x1.c91578p-3f, -0x1.9cf2e8p-5f, 0x1.9984ccp+0f, -0x1.075a46p+2f, 0x1.41a436p+2f, -0x1.0dcc3cp+1f, -0x1.9e984cp+1f, 0x1.00fa12p+3f, -0x1.317cc4p+3f, 0x1.8d688ep+2f, -0x1.a8daep+0f, 0x0p+0f};
            }
            if constexpr (order == 4) {
                return {0x0p+0f, 0x1.962a6ep-31f, -0x1.49b0cap-19f, 0x1.a4010cp-19f, 0x1.0d9c8cp-7f, -0x1.fc20c8p-6f, -0x1.35466cp-2f, 0x1.244d46p+1f, -0x1.73b0a4p+2f, 0x1.11f34p+2f, 0x1.6502b8p+3f, -0x1.2299d6p+5f, 0x1.b92a7cp+5f, -0x1.eefc0cp+5f, 0x1.fdbc98p+5f, -0x1.dc75fap+5f, 0x1.4ce384p+5f, -0x1.18f0a4p+4f, 0x1.a04f8p+1f, -0x0p+0f};
            }
            if constexpr (order == 5) {
                return {0x0p+0f, -0x1.f5ed4cp-24f, 0x1.968e84p-13f, 0x1.38e018p-10f, -0x1.43d264p-2f, 0x1.9505e8p-1f, 0x1.959602p+2f, -0x1.3fa55cp+5f, 0x1.934736p+6f, -0x1.207a74p+7f, 0x1.12118cp+7f, -0x1.d17f22p+6f, 0x1.b35f4ap+6f, -0x1.fc1722p+5f, -0x1.17b23p+5f, 0x1.96ec36p+6f, -0x1.505f98p+6f, 0x1.09d96p+5f, -0x1.55886ep+2f, 0x0p+0f};
            }
            if constexpr (order == 6) {
                return {0x0p+0f, -0x1.16b49ap-22f, 0x1.c193f2p-13f, 0x1.ba72fap-9f, -0x1.40c0aep-3f, 0x1.a1f0a6p-2f, 0x1.4d0298p+0f, -0x1.fe186ep+2f, 0x1.259afcp+4f, -0x1.fb5cc4p+4f, 0x1.ae7316p+5f, -0x1.3f986p+6f, 0x1.7db87cp+6f, -0x1.addb2cp+6f, 0x1.f01b3p+6f, -0x1.d673d8p+6f, 0x1.260fdap+6f, -0x1.9eeda2p+4f, 0x1.f3cda2p+1f, -0x0p+0f};
            }
            if constexpr (order == 7) {
                return {0x0p+0f, 0x1.a79c62p-20f, -0x1.52be3cp-11f, -0x1.843d22p-7f, 0x1.87564p-3f, -0x1.181d4ep-1f, -0x1.6763c6p-4f, 0x1.4065c2p+1f, -0x1.38745p+2f, 0x1.2ddf68p+3f, -0x1.657834p+4f, 0x1.209b54p+5f, -0x1.4053c8p+5f, 0x1.83d54ap+5f, -0x1.14f1dp+6f, 0x1.276b64p+6f, -0x1.7c41fp+5f, 0x1.08fa64p+4f, -0x1.3586dap+1f, 0x0p+0f};
            }
            if constexpr (order == 8) {
                return {0x0p+0f, -0x1.d74fd6p-15f, 0x1.726478p-7f, 0x1.4b4084p-3f, -0x1.4f17d6p+0f, 0x1.b8a4e2p+1f, -0x1.247b58p+2f, 0x1.7e1a34p+2f, -0x1.7e0b78p+3f, 0x1.e686e4p+3f, -0x1.8ceb32p+2f, 0x1.9fa6aap-4f, -0x1.11ca3p+3f, 0x1.74ac2ep+2f, 0x1.4db4fcp+4f, -0x1.37885ep+5f, 0x1.cbd932p+4f, -0x1.48b25ep+3f, 0x1.7a3cf4p+0f, 0x0p+0f};
            }
            if constexpr (order == 9) {
                return {0x0p+0f, -0x1.30a436p-18f, 0x1.ce071ep-12f, 0x1.0952d6p-8f, -0x1.5730c4p-6f, 0x1.49980ep-5f, -0x1.8ec972p-4f, 0x1.faf556p-3f, -0x1.8145e4p-2f, 0x1.158768p-1f, -0x1.16aae2p+0f, 0x1.9bd984p+0f, -0x1.9413f4p+0f, 0x1.00d5dp+1f, -0x1.ae9244p+1f, 0x1.eac616p+1f, -0x1.3f509cp+1f, 0x1.b7ba92p-1f, -0x1.f6438ap-4f, -0x0p+0f};
            }
         }
        if constexpr (std::is_same_v<Real, double>) {
            if constexpr (order == 0) {
                return {0x0p+0, -0x1.0988745ca3ec4p-59, 0x1.afee7cd299cd6p-44, 0x1.db358382953c4p-37, -0x1.58dbb791a7fb7p-28, 0x1.c94370d7ecc92p-24, 0x1.403905df41d32p-24, -0x1.b37df99b75f09p-16, 0x1.44e94de53f6f6p-14, 0x1.94f202846521cp-12, -0x1.2885ec2406341p-9, 0x1.d4ed15a068488p-10, 0x1.1d23df75685c9p-6, -0x1.494bbf8e53bfcp-4, 0x1.9e3dada7662aap-3, -0x1.85d2edb33e46ep-2, 0x1.1c4636ef189edp-1, 0x1.4e2c27be21e0ep-1, 0x1.12cb0a33640b5p-5, 0x0p+0};
            }
            if constexpr (order == 1) {
                return {0x0p+0, -0x1.85b82452d9269p-55, 0x1.3ced9d37da5d3p-40, 0x1.51db791710beap-33, -0x1.f125f98e34b08p-26, -0x1.2eea964fb1e8p-21, -0x1.06ba9aa22e87cp-17, -0x1.ce035ff90a79bp-15, 0x1.f375455c5c9e6p-12, -0x1.2c315a17e4064p-10, 0x1.213456d8479d6p-8, -0x1.c84ec23c8cdbep-6, 0x1.ac22308c89e82p-4, -0x1.eef0eaaf4989ep-3, 0x1.a2cb67ff996a9p-2, -0x1.61584eb9fe3a6p-1, 0x1.94bc27297b9f9p+0, -0x1.f5f9c20090df3p-1, -0x1.43cfcb9d73923p-3, -0x0p+0};
            }
            if constexpr (order == 2) {
                return {0x0p+0, 0x1.4b9c0643323fbp-46, -0x1.0d99c100eff2ep-32, 0x1.3f418b80d3d45p-28, 0x1.bb4b6c6c9a6fbp-19, -0x1.2b4a0836ed60bp-15, -0x1.67961ae7bb9ap-12, 0x1.0767285c8f106p-8, -0x1.51b88ceed7adep-7, -0x1.35a481158ac34p-6, 0x1.55c762fd5a0dfp-3, -0x1.9dcbe6db37906p-2, 0x1.86e851ba85835p-2, 0x1.6e6d859b2e73p-2, -0x1.c8684eac579b4p+0, 0x1.8bd77a675e828p+1, -0x1.0fd0995df5aa3p+1, -0x1.01b00cc8d5d06p-2, 0x1.2d6cb779664d4p-1, -0x0p+0};
            }
            if constexpr (order == 3) {
                return {0x0p+0, 0x1.cd96fccfe9f0bp-39, -0x1.77139f17861ebp-26, 0x1.03c4250827b69p-22, 0x1.351528a4bd817p-13, -0x1.dd5da666efa59p-11, -0x1.324229742ef44p-7, 0x1.60b8b0c3b1139p-4, -0x1.c915787e3612ep-3, -0x1.9cf2e7196d1ccp-5, 0x1.9984cb75bad51p+0, -0x1.075a464fd9704p+2, 0x1.41a435a1df481p+2, -0x1.0dcc3bbac5e09p+1, -0x1.9e984b48af592p+1, 0x1.00fa12867f491p+3, -0x1.317cc4a10faf4p+3, 0x1.8d688d31d04a4p+2, -0x1.a8dadfc02dd04p+0, 0x0p+0};
            }
            if constexpr (order == 4) {
                return {0x0p+0, 0x1.962a6dd57b7dbp-31, -0x1.49b0c99eb9c65p-19, 0x1.a4010b80e5496p-19, 0x1.0d9c8b9d476f8p-7, -0x1.fc20c8dd92f1ap-6, -0x1.35466c5650c62p-2, 0x1.244d452e76301p+1, -0x1.73b0a319e688ap+2, 0x1.11f33f4ede648p+2, 0x1.6502b707871dfp+3, -0x1.2299d502734c4p+5, 0x1.b92a7c1b6a88ap+5, -0x1.eefc0cdda4efap+5, 0x1.fdbc97313e6edp+5, -0x1.dc75f9b0b9a1p+5, 0x1.4ce3831796a62p+5, -0x1.18f0a34d9f68ep+4, 0x1.a04f8069ea8b5p+1, -0x0p+0};
            }
            if constexpr (order == 5) {
                return {0x0p+0, -0x1.f5ed4b124df3bp-24, 0x1.968e83ac595b9p-13, 0x1.38e017336cd06p-10, -0x1.43d2632e32cd1p-2, 0x1.9505e84d98e61p-1, 0x1.959601005ec55p+2, -0x1.3fa55ba9206e4p+5, 0x1.934735f943482p+6, -0x1.207a74ca5ad8cp+7, 0x1.12118c5678b4dp+7, -0x1.d17f22f2480f9p+6, 0x1.b35f4a36b31c4p+6, -0x1.fc17211fa4d6fp+5, -0x1.17b23075f801p+5, 0x1.96ec354b30307p+6, -0x1.505f98e4bf277p+6, 0x1.09d95f1cede32p+5, -0x1.55886d2f64381p+2, 0x0p+0};
            }
            if constexpr (order == 6) {
                return {0x0p+0, -0x1.16b4995bb4f07p-22, 0x1.c193f16bb678ep-13, 0x1.ba72fac37f088p-9, -0x1.40c0addfd8554p-3, 0x1.a1f0a66406a12p-2, 0x1.4d0298c218927p+0, -0x1.fe186e46030ap+2, 0x1.259afcbb8804dp+4, -0x1.fb5cc3351ebc1p+4, 0x1.ae731643e9199p+5, -0x1.3f9860f6c4906p+6, 0x1.7db87bda6fdf9p+6, -0x1.addb2c4ad791ep+6, 0x1.f01b2fa154d03p+6, -0x1.d673d80eeade5p+6, 0x1.260fdaadd19ebp+6, -0x1.9eeda12641b62p+4, 0x1.f3cda2f9b1354p+1, -0x0p+0};
            }
            if constexpr (order == 7) {
                return {0x0p+0, 0x1.a79c626c8654fp-20, -0x1.52be3b0b8c7afp-11, -0x1.843d22529bb4ap-7, 0x1.87564028831fcp-3, -0x1.181d4d676ea1fp-1, -0x1.6763c5d11b775p-4, 0x1.4065c212ba5cp+1, -0x1.38744fdae2125p+2, 0x1.2ddf670e51dfdp+3, -0x1.65783340b3783p+4, 0x1.209b538b18df2p+5, -0x1.4053c8a2c8d3p+5, 0x1.83d54ac636b82p+5, -0x1.14f1cf74869cp+6, 0x1.276b64d22e135p+6, -0x1.7c41f01d6960dp+5, 0x1.08fa63342ecc2p+4, -0x1.3586d9a9e0e3ep+1, 0x0p+0};
            }
            if constexpr (order == 8) {
                return {0x0p+0, -0x1.d74fd6a021548p-15, 0x1.7264788c9473ep-7, 0x1.4b408438ac869p-3, -0x1.4f17d502240a9p+0, 0x1.b8a4e1c96768dp+1, -0x1.247b57131543cp+2, 0x1.7e1a3482e160cp+2, -0x1.7e0b780d60232p+3, 0x1.e686e3fbc3fbfp+3, -0x1.8ceb3288570cdp+2, 0x1.9fa6aaf0c8d98p-4, -0x1.11ca2f9f0b7dp+3, 0x1.74ac2dca21f56p+2, 0x1.4db4fb7a5ecdcp+4, -0x1.37885d5242031p+5, 0x1.cbd93276e7d59p+4, -0x1.48b25e52d274bp+3, 0x1.7a3cf471120ecp+0, 0x0p+0};
            }
            if constexpr (order == 9) {
                return {0x0p+0, -0x1.30a4361eb25b4p-18, 0x1.ce071dd8df1efp-12, 0x1.0952d6b855c64p-8, -0x1.5730c48bae56ap-6, 0x1.49980ed40883dp-5, -0x1.8ec972676b06bp-4, 0x1.faf556e68970ep-3, -0x1.8145e3cc3cfccp-2, 0x1.1587682f48216p-1, -0x1.16aae27d18a67p+0, 0x1.9bd9847a10dd5p+0, -0x1.9413f3bab26fdp+0, 0x1.00d5d03220f3cp+1, -0x1.ae9244fef585dp+1, 0x1.eac616f0d4f87p+1, -0x1.3f509b3842149p+1, 0x1.b7ba9255f9c78p-1, -0x1.f6438a3a8846p-4, -0x0p+0};
            }
         }
        if constexpr (std::is_same_v<Real, long double>) {
            if constexpr (order == 0) {
                return {0x0p+0L, -0x8.4c43a2e51f6200cp-62L, 0xd.7f73e694ce6b23p-47L, 0xe.d9ac1c14a9e22dfp-40L, -0xa.c6ddbc8d3fdb61fp-31L, 0xe.4a1b86bf6648ca3p-27L, 0xa.01c82efa0e991ep-27L, -0xd.9befccdbaf8457bp-19L, 0xa.274a6f29fb7b195p-17L, 0xc.a7901423290dc7ep-15L, -0x9.442f612031a0acp-12L, 0xe.a768ad03424408ap-13L, 0x8.e91efbab42e4488p-9L, -0xa.4a5dfc729dfe012p-7L, 0xc.f1ed6d3b3154c33p-6L, -0xc.2e976d99f236ef3p-5L, 0x8.e231b778c4f6662p-4L, 0xa.71613df10f0723dp-4L, 0x8.9658519b205ab72p-8L, 0x0p+0L};
            }
            if constexpr (order == 1) {
                return {0x0p+0L, -0xc.2dc12296c934a35p-58L, 0x9.e76ce9bed2e946cp-43L, 0xa.8edbc8b885f4c32p-36L, -0xf.892fcc71a583f65p-29L, -0x9.7754b27d8f401bap-24L, -0x8.35d4d511743e237p-20L, -0xe.701affc853cdb21p-18L, 0xf.9baa2ae2e4f30dfp-15L, -0x9.618ad0bf2031fa3p-13L, 0x9.09a2b6c23cead05p-11L, -0xe.427611e466df078p-9L, 0xd.611184644f4100cp-7L, -0xf.7787557a4c4f14ep-6L, 0xd.165b3ffccb5480ap-5L, -0xb.0ac275cff1d33d6p-4L, 0xc.a5e1394bdcfcb6dp-3L, -0xf.afce100486f966ap-4L, -0xa.1e7e5ceb9c91909p-6L, -0x0p+0L};
            }
            if constexpr (order == 2) {
                return {0x0p+0L, 0xa.5ce0321991fd8ep-49L, -0x8.6cce08077f97363p-35L, 0x9.fa0c5c069ea261fp-31L, 0xd.da5b6364d37d7c6p-22L, -0x9.5a5041b76b0571dp-18L, -0xb.3cb0d73ddccfdfep-15L, 0x8.3b3942e4788322ap-11L, -0xa.8dc46776bd6ee4dp-10L, -0x9.ad2408ac5619da4p-9L, 0xa.ae3b17ead06fb26p-6L, -0xc.ee5f36d9bc83386p-5L, 0xc.37428dd42c1a763p-5L, 0xb.736c2cd97397c01p-5L, -0xe.43427562bcd9f08p-3L, 0xc.5ebbd33af4141bap-2L, -0x8.7e84caefad51526p-2L, -0x8.0d806646ae83286p-5L, 0x9.6b65bbcb326a224p-4L, -0x0p+0L};
            }
            if constexpr (order == 3) {
                return {0x0p+0L, 0xe.6cb7e67f4f857ddp-42L, -0xb.b89cf8bc30f58d6p-29L, 0x8.1e2128413db47b5p-25L, 0x9.a8a94525ec0bb0bp-16L, -0xe.eaed33377d2c762p-14L, -0x9.92114ba177a20b1p-10L, 0xb.05c5861d889c6fap-7L, -0xe.48abc3f1b09716ap-6L, -0xc.e79738cb68e5e97p-8L, 0xc.cc265badd6a86efp-3L, -0x8.3ad2327ecb81caap-1L, 0xa.0d21ad0efa407b6p-1L, -0x8.6e61ddd62f0442fp-2L, -0xc.f4c25a457ac90c1p-2L, 0x8.07d09433fa48997p+0L, -0x9.8be625087d79ffdp+0L, 0xc.6b44698e825210ap-1L, -0xd.46d6fe016e82146p-3L, 0x0p+0L};
            }
            if constexpr (order == 4) {
                return {0x0p+0L, 0xc.b1536eabdbed537p-34L, -0xa.4d864cf5ce32429p-22L, 0xd.20085c072a4ae87p-22L, 0x8.6ce45cea3b7be89p-10L, -0xf.e10646ec978cdc7p-9L, -0x9.aa3362b28630f09p-5L, 0x9.226a2973b180497p-2L, -0xb.9d8518cf3444f4dp-1L, 0x8.8f99fa76f323eadp-1L, 0xb.2815b83c38efa28p+0L, -0x9.14cea8139a61e4ep+2L, 0xd.c953e0db5444e19p+2L, -0xf.77e066ed277d141p+2L, 0xf.ede4b989f376b43p+2L, -0xe.e3afcd85cd0820fp+2L, 0xa.671c18bcb530f7bp+2L, -0x8.c7851a6cfb47268p+1L, 0xd.027c034f545a722p-2L, -0x0p+0L};
            }
            if constexpr (order == 5) {
                return {0x0p+0L, -0xf.af6a58926f9d684p-27L, 0xc.b4741d62cadc637p-16L, 0x9.c700b99b66830a5p-13L, -0xa.1e93197196686f8p-5L, 0xc.a82f426cc7309e8p-4L, 0xc.acb00802f62a439p-1L, -0x9.fd2add490371c09p+2L, 0xc.9a39afca1a40e4p+3L, -0x9.03d3a652d6c5fe1p+4L, 0x8.908c62b3c5a6a21p+4L, -0xe.8bf91792407c46p+3L, 0xd.9afa51b598e1c25p+3L, -0xf.e0b908fd26b79fbp+2L, -0x8.bd9183afc008112p+2L, 0xc.b761aa5981834ffp+3L, -0xa.82fcc725f93bb38p+3L, 0x8.4ecaf8e76f18fe7p+2L, -0xa.ac43697b21c05d7p-1L, 0x0p+0L};
            }
            if constexpr (order == 6) {
                return {0x0p+0L, -0x8.b5a4cadda783993p-25L, 0xe.0c9f8b5db3c72d8p-16L, 0xd.d397d61bf843e3dp-12L, -0xa.06056efec2aa1d8p-6L, 0xd.0f8533203509157p-5L, 0xa.6814c610c493aa3p-3L, -0xf.f0c37230184fee7p-1L, 0x9.2cd7e5dc4026b75p+1L, -0xf.dae619a8f5e0bb5p+1L, 0xd.7398b21f48cc731p+2L, -0x9.fcc307b62482f39p+3L, 0xb.edc3ded37efc53ap+3L, -0xd.6ed96256bc8f10cp+3L, 0xf.80d97d0aa681a95p+3L, -0xe.b39ec07756f2504p+3L, 0x9.307ed56e8cf58f8p+3L, -0xc.f76d09320db0f5ep+1L, 0xf.9e6d17cd89aa012p-2L, -0x0p+0L};
            }
            if constexpr (order == 7) {
                return {0x0p+0L, 0xd.3ce3136432a75eep-23L, -0xa.95f1d85c63d7aadp-14L, -0xc.21e91294dda4d8cp-10L, 0xc.3ab2014418fdea7p-6L, -0x8.c0ea6b3b750f81fp-4L, -0xb.3b1e2e88dbba717p-7L, 0xa.032e1095d2e022bp-2L, -0x9.c3a27ed71092bccp-1L, 0x9.6efb38728efe53fp+0L, -0xb.2bc19a059bc163fp+1L, 0x9.04da9c58c6f9312p+2L, -0xa.029e45164697e35p+2L, 0xc.1eaa5631b5c1221p+2L, -0x8.a78e7ba434e03bep+3L, 0x9.3b5b2691709a7f2p+3L, -0xb.e20f80eb4b065edp+2L, 0x8.47d319a17660e21p+1L, -0x9.ac36cd4f071eccfp-2L, 0x0p+0L};
            }
            if constexpr (order == 8) {
                return {0x0p+0L, -0xe.ba7eb5010aa3e4p-18L, 0xb.9323c464a39ee53p-10L, 0xa.5a0421c564344a9p-6L, -0xa.78bea8112054738p-3L, 0xd.c5270e4b3b4644bp-2L, -0x9.23dab898aa1de5fp-1L, 0xb.f0d1a4170b05fe8p-1L, -0xb.f05bc06b0118d26p+0L, 0xf.34371fde1fdf9acp+0L, -0xc.67599442b866996p-1L, 0xc.fd35578646cc028p-7L, -0x8.8e517cf85be82f9p+0L, 0xb.a5616e510faadd9p-1L, 0xa.6da7dbd2f66dc91p+1L, -0x9.bc42ea921018403p+2L, 0xe.5ec993b73eacb21p+1L, -0xa.4592f29693a58f2p+0L, 0xb.d1e7a38890760ccp-3L, 0x0p+0L};
            }
            if constexpr (order == 9) {
                return {0x0p+0L, -0x9.8521b0f592d9c7fp-21L, 0xe.7038eec6f8f774fp-15L, 0x8.4a96b5c2ae31c5dp-11L, -0xa.b986245d72b5275p-9L, 0xa.4cc076a0441e56ep-8L, -0xc.764b933b583557ap-7L, 0xf.d7aab7344b86da8p-6L, -0xc.0a2f1e61e7e5e62p-5L, 0x8.ac3b417a410ae11p-4L, -0x8.b55713e8c533b1cp-3L, 0xc.decc23d086ea9cp-3L, -0xc.a09f9dd5937e42dp-3L, 0x8.06ae8191079dde4p-2L, -0xd.749227f7ac2e953p-2L, 0xf.5630b786a7c34dfp-2L, -0x9.fa84d9c210a45adp-2L, 0xd.bdd492afce3bd22p-4L, -0xf.b21c51d44230049p-7L, -0x0p+0L};
            }
         }
#ifdef BOOST_HAS_FLOAT128
        if constexpr (std::is_same_v<Real, boost::multiprecision::float128>) {
            if constexpr (order == 0) {
                return {0x0p+0Q, -0x1.0988745ca3ec40178fb79a68846ep-59Q, 0x1.afee7cd299cd64600a0752c73efp-44Q, 0x1.db358382953c45bea2586fcc9422p-37Q, -0x1.58dbb791a7fb6c3e33eb298f4d3ep-28Q, 0x1.c94370d7ecc9194547661e014d91p-24Q, 0x1.403905df41d323c0c33071ba892bp-24Q, -0x1.b37df99b75f08af621654dd9209cp-16Q, 0x1.44e94de53f6f632a67cbb58576c1p-14Q, 0x1.94f202846521b8fb6dee9352be82p-12Q, -0x1.2885ec240634157f42540c2ec4a5p-9Q, 0x1.d4ed15a06848811344377f87f482p-10Q, 0x1.1d23df75685c890f2df5ab229b92p-6Q, -0x1.494bbf8e53bfc023fe7e0915ab64p-4Q, 0x1.9e3dada7662a9865be6b99b2429cp-3Q, -0x1.85d2edb33e46dde58a5a504a9ac4p-2Q, 0x1.1c4636ef189eccc399682014d9cdp-1Q, 0x1.4e2c27be21e0e479b88ad4851617p-1Q, 0x1.12cb0a33640b56e423474a1328c7p-5Q, 0x0p+0Q};
            }
            if constexpr (order == 1) {
                return {0x0p+0Q, -0x1.85b82452d9269469aa073853ba4ep-55Q, 0x1.3ced9d37da5d28d70019b39cddf2p-40Q, 0x1.51db791710be9863d69b804f648bp-33Q, -0x1.f125f98e34b07ec90a3114694d5cp-26Q, -0x1.2eea964fb1e803735b7cbde59c5bp-21Q, -0x1.06ba9aa22e87c46e9b5b517ee7bfp-17Q, -0x1.ce035ff90a79b6418f9714e094d7p-15Q, 0x1.f375455c5c9e61be4fa670797afep-12Q, -0x1.2c315a17e4063f457d1b74aab03dp-10Q, 0x1.213456d8479d5a0a002e8a354438p-8Q, -0x1.c84ec23c8cdbe0f0d252ca4a6027p-6Q, 0x1.ac22308c89e820172f07afea78p-4Q, -0x1.eef0eaaf4989e29be4ba95e2f1f4p-3Q, 0x1.a2cb67ff996a90130726bb55194dp-2Q, -0x1.61584eb9fe3a67acec590408bd8fp-1Q, 0x1.94bc27297b9f96d98ca06b8b787cp+0Q, -0x1.f5f9c20090df2cd39669639ec81fp-1Q, -0x1.43cfcb9d739232114ccf7b42f2e2p-3Q, -0x0p+0Q};
            }
            if constexpr (order == 2) {
                return {0x0p+0Q, 0x1.4b9c0643323fb1c0ebf8a306241p-46Q, -0x1.0d99c100eff2e6c617789a9f70efp-32Q, 0x1.3f418b80d3d44c3e61feb87055afp-28Q, 0x1.bb4b6c6c9a6faf8c7dd84bcc5e06p-19Q, -0x1.2b4a0836ed60ae3a1490ad51d3e6p-15Q, -0x1.67961ae7bb99fbfcc0d26143de94p-12Q, 0x1.0767285c8f10645316c4c221ca25p-8Q, -0x1.51b88ceed7addc9a97168ddbb1f3p-7Q, -0x1.35a481158ac33b4716b6e91c1302p-6Q, 0x1.55c762fd5a0df64b9c06acd4f05ap-3Q, -0x1.9dcbe6db3790670c84be0eebddedp-2Q, 0x1.86e851ba85834ec58a127f16a3a3p-2Q, 0x1.6e6d859b2e72f8027855287734ebp-2Q, -0x1.c8684eac579b3e1019c502de48a5p+0Q, 0x1.8bd77a675e82837408b6453226d2p+1Q, -0x1.0fd0995df5aa2a4c5f39f9316485p+1Q, -0x1.01b00cc8d5d0650bb7d9e8751632p-2Q, 0x1.2d6cb779664d44484073db4511a4p-1Q, -0x0p+0Q};
            }
            if constexpr (order == 3) {
                return {0x0p+0Q, 0x1.cd96fccfe9f0afba331831b5a876p-39Q, -0x1.77139f17861eb1acd309511a0171p-26Q, 0x1.03c4250827b68f69b983fa4bde7p-22Q, 0x1.351528a4bd8176153593c8b22e24p-13Q, -0x1.dd5da666efa58ec47474db24205ep-11Q, -0x1.324229742ef4416185466da25dc1p-7Q, 0x1.60b8b0c3b1138df4716646361ef1p-4Q, -0x1.c915787e3612e2d3277f4d9a9293p-3Q, -0x1.9cf2e7196d1cbd2ef1ba36dde8bp-5Q, 0x1.9984cb75bad50ddd23b7c6bc0469p+0Q, -0x1.075a464fd9703954dc501fb568ddp+2Q, 0x1.41a435a1df480f6b501eca8095p+2Q, -0x1.0dcc3bbac5e0885ef946fb6acf6fp+1Q, -0x1.9e984b48af59218135976f5e8271p+1Q, 0x1.00fa12867f49132d26f610a5f757p+3Q, -0x1.317cc4a10faf3ff9c44aa0794aa9p+3Q, 0x1.8d688d31d04a4214094b65cddcbap+2Q, -0x1.a8dadfc02dd0428b41cd0c2467d5p+0Q, 0x0p+0Q};
            }
            if constexpr (order == 4) {
                return {0x0p+0Q, 0x1.962a6dd57b7daa6dc362e159413bp-31Q, -0x1.49b0c99eb9c648511ea73e827034p-19Q, 0x1.a4010b80e5495d0ed93317483b25p-19Q, 0x1.0d9c8b9d476f7d129de2a834c7aap-7Q, -0x1.fc20c8dd92f19b8d6ed30f2d1e03p-6Q, -0x1.35466c5650c61e120220d2c0eecdp-2Q, 0x1.244d452e7630092d6c490ecfe00ep+1Q, -0x1.73b0a319e6889e9aa2e8320b1af8p+2Q, 0x1.11f33f4ede647d59a412fde094acp+2Q, 0x1.6502b707871df45089990cd289cap+3Q, -0x1.2299d502734c3c9b7e4c9bf39012p+5Q, 0x1.b92a7c1b6a889c31c867c32e9611p+5Q, -0x1.eefc0cdda4efa281b21fcb804378p+5Q, 0x1.fdbc97313e6ed6865f909a084cb7p+5Q, -0x1.dc75f9b0b9a1041da91936d52bc5p+5Q, 0x1.4ce3831796a61ef6af32c0f62affp+5Q, -0x1.18f0a34d9f68e4d06a24b7b54699p+4Q, 0x1.a04f8069ea8b4e447d853d5c1328p+1Q, -0x0p+0Q};
            }
            if constexpr (order == 5) {
                return {0x0p+0Q, -0x1.f5ed4b124df3ad0700a16c248b74p-24Q, 0x1.968e83ac595b8c6d9c175d111059p-13Q, 0x1.38e017336cd0614a6ed32ba5fa09p-10Q, -0x1.43d2632e32cd0df0bb4abc9c579ep-2Q, 0x1.9505e84d98e613cfe4ec6f61f8d7p-1Q, 0x1.959601005ec5487295591034a7abp+2Q, -0x1.3fa55ba9206e3812ee32fe0f733p+5Q, 0x1.934735f943481c7f74f59af630cfp+6Q, -0x1.207a74ca5ad8bfc2ee105c93a5d8p+7Q, 0x1.12118c5678b4d442de8a317c0014p+7Q, -0x1.d17f22f2480f88bf1a42a2b503ep+6Q, 0x1.b35f4a36b31c3849a4f0aa4f41b2p+6Q, -0x1.fc17211fa4d6f3f6a326867d94b2p+5Q, -0x1.17b23075f8010224c3488cfe72dep+5Q, 0x1.96ec354b303069fd6646216403d4p+6Q, -0x1.505f98e4bf27766f2f4d4b86bac1p+6Q, 0x1.09d95f1cede31fcd15d2d55f183bp+5Q, -0x1.55886d2f64380baeb961187c1a0ap+2Q, 0x0p+0Q};
            }
            if constexpr (order == 6) {
                return {0x0p+0Q, -0x1.16b4995bb4f073263c24ab73bfbfp-22Q, 0x1.c193f16bb678e5af651147a8981bp-13Q, 0x1.ba72fac37f087c7aee3a6324de5cp-9Q, -0x1.40c0addfd85543af466d0e5b17a7p-3Q, 0x1.a1f0a66406a122ae435f38f9cb9p-2Q, 0x1.4d0298c21892754606c4092af2dp+0Q, -0x1.fe186e460309fdcedab63b0e1b1fp+2Q, 0x1.259afcbb8804d6eae552d647dd5ep+4Q, -0x1.fb5cc3351ebc176988bfb25cf9ffp+4Q, 0x1.ae731643e9198e621648a2b50d86p+5Q, -0x1.3f9860f6c4905e7132b1b1e5f348p+6Q, 0x1.7db87bda6fdf8a74bb08658d9d79p+6Q, -0x1.addb2c4ad791e217324cca49ed1dp+6Q, 0x1.f01b2fa154d0352990b1ce06f03ep+6Q, -0x1.d673d80eeade4a070ce149af8cb2p+6Q, 0x1.260fdaadd19eb1efc0791662953fp+6Q, -0x1.9eeda12641b61ebb7d846660b657p+4Q, 0x1.f3cda2f9b13540236201ea4e96d5p+1Q, -0x0p+0Q};
            }
            if constexpr (order == 7) {
                return {0x0p+0Q, 0x1.a79c626c8654ebdb591ab6934303p-20Q, -0x1.52be3b0b8c7af55af9d86785110ap-11Q, -0x1.843d22529bb49b17cf1513d3735dp-7Q, 0x1.87564028831fbd4eb29d0e5f8344p-3Q, -0x1.181d4d676ea1f03e7d8be530a435p-1Q, -0x1.6763c5d11b774e2d99bcbb38601bp-4Q, 0x1.4065c212ba5c04567835fae6b8f6p+1Q, -0x1.38744fdae212579791ea6bb69468p+2Q, 0x1.2ddf670e51dfca7d85e8b6fe072fp+3Q, -0x1.65783340b3782c7d703e3ec9694ap+4Q, 0x1.209b538b18df262448edf4b7d07bp+5Q, -0x1.4053c8a2c8d2fc6aef15d61af7eep+5Q, 0x1.83d54ac636b8244258f18503a19fp+5Q, -0x1.14f1cf74869c077b4ddd14fa8f81p+6Q, 0x1.276b64d22e134fe39031246250b4p+6Q, -0x1.7c41f01d6960cbda171cb5428a61p+5Q, 0x1.08fa63342ecc1c41a46aea3feb2p+4Q, -0x1.3586d9a9e0e3d99d2b9347eb730ap+1Q, 0x0p+0Q};
            }
            if constexpr (order == 8) {
                return {0x0p+0Q, -0x1.d74fd6a021547c7f9f7b587d2ec7p-15Q, 0x1.7264788c9473dca575618fd77c5ep-7Q, 0x1.4b408438ac8689521da09d55a0e1p-3Q, -0x1.4f17d502240a8e6f55c847ad377dp+0Q, 0x1.b8a4e1c96768c895a2f153e0232p+1Q, -0x1.247b57131543bcbee8aeed1657dep+2Q, 0x1.7e1a3482e160bfcf5baf7d4056c6p+2Q, -0x1.7e0b780d60231a4c786cba95bd49p+3Q, 0x1.e686e3fbc3fbf35811a3d925037ap+3Q, -0x1.8ceb3288570cd32c258dd6301b5cp+2Q, 0x1.9fa6aaf0c8d980502ec0fb5c781bp-4Q, -0x1.11ca2f9f0b7d05f1c15eae7277d6p+3Q, 0x1.74ac2dca21f55bb2488107bbadf4p+2Q, 0x1.4db4fb7a5ecdb9227e881ad0531p+4Q, -0x1.37885d5242030806d66dbf1a3c7p+5Q, 0x1.cbd93276e7d596419b901d50ac68p+4Q, -0x1.48b25e52d274b1e47c879b7b8ba3p+3Q, 0x1.7a3cf471120ec1987138adbab023p+0Q, 0x0p+0Q};
            }
            if constexpr (order == 9) {
                return {0x0p+0Q, -0x1.30a4361eb25b38fd9207567fafc2p-18Q, 0x1.ce071dd8df1eee9d16330b7911bep-12Q, 0x1.0952d6b855c638b97745763eaf9bp-8Q, -0x1.5730c48bae56a4ea4426af715c3dp-6Q, 0x1.49980ed40883cadc7b70359f9993p-5Q, -0x1.8ec972676b06aaf4a1909c599bbep-4Q, 0x1.faf556e68970db4f5325c83b3059p-3Q, -0x1.8145e3cc3cfcbcc409c892d12fcbp-2Q, 0x1.1587682f48215c21c42b35317ddap-1Q, -0x1.16aae27d18a67638fdda85894127p+0Q, 0x1.9bd9847a10dd537f2de86a0edf44p+0Q, -0x1.9413f3bab26fc85ab6068b638ed5p+0Q, 0x1.00d5d03220f3bbc8a51442e895a2p+1Q, -0x1.ae9244fef585d2a56f1ca6b2357p+1Q, 0x1.eac616f0d4f869be16cdf1a8cfc4p+1Q, -0x1.3f509b3842148b5acdb04f9867dp+1Q, 0x1.b7ba9255f9c77a442c3301f7dc97p-1Q, -0x1.f6438a3a88460091e0ea6e624f4dp-4Q, -0x0p+0Q};
            }
         }
#endif
     }
    if constexpr (p == 11) {
        if constexpr (std::is_same_v<Real, float>) {
            if constexpr (order == 0) {
                return {0x0p+0f, -0x1.65d844p-70f, -0x1.ad914cp-53f, -0x1.4f944ep-42f, -0x1.160128p-35f, 0x1.1409ccp-28f, -0x1.0af58p-26f, 0x1.63562ep-20f, 0x1.5c274ep-18f, -0x1.d7e9bap-15f, 0x1.f2bfaap-15f, 0x1.d91cd2p-14f, 0x1.0182ap-9f, -0x1.071062p-6f, 0x1.d853p-5f, -0x1.17e06ap-3f, 0x1.07dc8cp-2f, -0x1.d3dbc2p-2f, 0x1.a074eap-1f, 0x1.da87aep-2f, 0x1.f82d86p-7f, 0x0p+0f};
            }
            if constexpr (order == 1) {
                return {0x0p+0f, -0x1.29783p-62f, -0x1.651bd6p-46f, 0x1.d9433ap-36f, -0x1.739c16p-30f, -0x1.7db784p-25f, 0x1.e02848p-20f, 0x1.c4d82p-19f, -0x1.a2acf4p-14f, 0x1.693106p-12f, 0x1.185782p-11f, -0x1.092c0ap-7f, 0x1.fec808p-6f, -0x1.0f00fp-4f, 0x1.321c6cp-4f, 0x1.71d07ep-12f, -0x1.303436p-3f, 0x1.075e3cp-3f, 0x1.f03418p-1f, -0x1.cd9786p-1f, -0x1.4988f6p-4f, 0x0p+0f};
            }
            if constexpr (order == 2) {
                return {0x0p+0f, -0x1.32449ep-55f, -0x1.6fb512p-40f, 0x1.64323cp-30f, -0x1.6377bp-25f, -0x1.6928bp-19f, 0x1.0136c4p-15f, 0x1.af648p-17f, -0x1.cd1c36p-10f, 0x1.21e04ap-7f, -0x1.6c5862p-7f, -0x1.41156ep-5f, 0x1.4f999p-3f, -0x1.78ed1cp-3f, -0x1.150c78p-2f, 0x1.5f04aep+0f, -0x1.6a2a88p+1f, 0x1.08f33p+2f, -0x1.a41fd4p+1f, 0x1.307a2ap-1f, 0x1.5c07e6p-2f, 0x0p+0f};
            }
            if constexpr (order == 3) {
                return {0x0p+0f, 0x1.e2c778p-49f, 0x1.21deep-34f, -0x1.8b6cb6p-25f, 0x1.f90492p-21f, 0x1.950dccp-15f, -0x1.1c307ep-11f, -0x1.f33a86p-12f, 0x1.35b764p-6f, -0x1.2b8bd2p-4f, 0x1.3f8a26p-4f, 0x1.08b97cp-2f, -0x1.339d14p+0f, 0x1.38ad98p+1f, -0x1.80597ep+1f, 0x1.0fb22ap+1f, -0x1.02364ap-1f, 0x1.8f506ep-1f, -0x1.8e9818p+1f, 0x1.a72bc4p+1f, -0x1.1c297ep+0f, -0x0p+0f};
            }
            if constexpr (order == 4) {
                return {0x0p+0f, -0x1.8ec07p-39f, -0x1.df066ep-26f, 0x1.ba1b22p-17f, -0x1.695b64p-13f, -0x1.22b5dp-7f, 0x1.002e9ap-4f, 0x1.8f714p-5f, -0x1.9e20d2p+0f, 0x1.8e7526p+2f, -0x1.4b3d36p+3f, 0x1.64c18ap+1f, 0x1.6d1976p+4f, -0x1.aeeda8p+5f, 0x1.1abf8ep+6f, -0x1.17ce56p+6f, 0x1.f144c4p+5f, -0x1.9aba4p+5f, 0x1.0e7952p+5f, -0x1.c33b04p+3f, 0x1.51d8dcp+1f, -0x0p+0f};
            }
            if constexpr (order == 5) {
                return {0x0p+0f, 0x1.ca2f1p-31f, 0x1.136d24p-18f, -0x1.7e95e2p-10f, 0x1.230c94p-7f, 0x1.1851aap-1f, -0x1.7dbb7ep+1f, -0x1.66137p-1f, 0x1.6db5a8p+5f, -0x1.4eba06p+7f, 0x1.4159ep+8f, -0x1.90ad4cp+8f, 0x1.80bbfp+8f, -0x1.5f18c6p+8f, 0x1.2e04acp+8f, -0x1.4a282ep+7f, -0x1.36cbd8p+4f, 0x1.d2259cp+6f, -0x1.7429c4p+6f, 0x1.1648a2p+5f, -0x1.553842p+2f, 0x0p+0f};
            }
            if constexpr (order == 6) {
                return {0x0p+0f, 0x1.ebdb46p-30f, 0x1.28216ep-18f, -0x1.81c424p-10f, -0x1.1e668p-11f, 0x1.1ee7d4p-2f, -0x1.5e908ap+0f, 0x1.14b6dap+0f, 0x1.f4f5eap+2f, -0x1.cfad0ap+4f, 0x1.cafceap+5f, -0x1.7f74fcp+6f, 0x1.2da088p+7f, -0x1.9299f8p+7f, 0x1.c24b2cp+7f, -0x1.dd46b6p+7f, 0x1.e7564ap+7f, -0x1.917bfcp+7f, 0x1.be9ebp+6f, -0x1.201578p+5f, 0x1.43e57ap+2f, -0x0p+0f};
            }
            if constexpr (order == 7) {
                return {0x0p+0f, -0x1.5739aap-26f, -0x1.9e94aap-16f, 0x1.cc0f8cp-8f, 0x1.89a1eap-6f, -0x1.521b96p-1f, 0x1.6d63c6p+1f, -0x1.29abecp+2f, 0x1.0097bap+0f, 0x1.a5a3f6p+2f, -0x1.d23af8p+3f, 0x1.2bd7eep+5f, -0x1.4b192ep+6f, 0x1.d61ad6p+6f, -0x1.00a2fp+7f, 0x1.356154p+7f, -0x1.7fbdfcp+7f, 0x1.5d87f4p+7f, -0x1.8cc936p+6f, 0x1.f54b1ep+4f, -0x1.0e918ep+2f, -0x0p+0f};
            }
            if constexpr (order == 8) {
                return {0x0p+0f, -0x1.6c2d94p-22f, -0x1.baa15ap-13f, 0x1.b3efb4p-6f, 0x1.3f1c58p-4f, -0x1.3d9d88p+0f, 0x1.152962p+2f, -0x1.1522f8p+3f, 0x1.e00e48p+3f, -0x1.a72562p+4f, 0x1.483818p+5f, -0x1.b51d9ap+5f, 0x1.14413cp+6f, -0x1.493302p+6f, 0x1.678e1ep+6f, -0x1.84d65cp+6f, 0x1.91c80cp+6f, -0x1.452574p+6f, 0x1.5bbb8cp+5f, -0x1.a9c968p+3f, 0x1.c2a5p+0f, -0x0p+0f};
            }
            if constexpr (order == 9) {
                return {0x0p+0f, 0x1.b40e78p-22f, 0x1.0c4738p-13f, -0x1.7c528ap-8f, -0x1.c9b3fp-8f, 0x1.0b6248p-3f, -0x1.61b17ap-2f, 0x1.624786p-1f, -0x1.cb38bcp+0f, 0x1.c75ab6p+1f, -0x1.3ef834p+2f, 0x1.f62894p+2f, -0x1.af56b2p+3f, 0x1.0c0dcap+4f, -0x1.0953a8p+4f, 0x1.536d7p+4f, -0x1.d3d44ap+4f, 0x1.bc60c4p+4f, -0x1.f8f412p+3f, 0x1.38fe3p+2f, -0x1.480f32p-1f, 0x0p+0f};
            }
            if constexpr (order == 10) {
                return {0x0p+0f, -0x1.09aedap-21f, -0x1.4eea02p-14f, 0x1.71e7a4p-10f, -0x1.a3115p-12f, -0x1.b7227ep-7f, 0x1.8d3b22p-6f, -0x1.50a6f6p-5f, 0x1.68181ap-3f, -0x1.7340e8p-2f, 0x1.a30d2p-2f, -0x1.7019bcp-1f, 0x1.8867ecp+0f, -0x1.db8c8p+0f, 0x1.96ff88p+0f, -0x1.2e8a4ep+1f, 0x1.fbc728p+1f, -0x1.04a532p+2f, 0x1.2f9284p+1f, -0x1.79dddep-1f, 0x1.8a25fep-4f, 0x0p+0f};
            }
         }
        if constexpr (std::is_same_v<Real, double>) {
            if constexpr (order == 0) {
                return {0x0p+0, -0x1.65d8440cb56ep-70, -0x1.ad914b3b9c645p-53, -0x1.4f944d6c974a1p-42, -0x1.1601287f207ddp-35, 0x1.1409cb495cfabp-28, -0x1.0af5800900cccp-26, 0x1.63562e1b18da8p-20, 0x1.5c274e2e04a65p-18, -0x1.d7e9ba7b3ac47p-15, 0x1.f2bfaace30f59p-15, 0x1.d91cd19825f46p-14, 0x1.0182a0c9a44bbp-9, -0x1.0710616766f9ep-6, 0x1.d85300a8e9a5ap-5, -0x1.17e06a0deb07p-3, 0x1.07dc8c4b5880ap-2, -0x1.d3dbc22db3363p-2, 0x1.a074ea6831a18p-1, 0x1.da87adf545deap-2, 0x1.f82d860ed311dp-7, 0x0p+0};
            }
            if constexpr (order == 1) {
                return {0x0p+0, -0x1.29782f98b98dap-62, -0x1.651bd505e5aabp-46, 0x1.d9433a2520a2ep-36, -0x1.739c15dc80c55p-30, -0x1.7db784c0b0131p-25, 0x1.e0284772a21efp-20, 0x1.c4d81f25bbab6p-19, -0x1.a2acf3b1f9d07p-14, 0x1.693105a2bb707p-12, 0x1.1857810566103p-11, -0x1.092c0960e80ecp-7, 0x1.fec807d6bcde3p-6, -0x1.0f00f07475782p-4, 0x1.321c6b5970b99p-4, 0x1.71d07ef81053fp-12, -0x1.303436308c6p-3, 0x1.075e3bcf1d71fp-3, 0x1.f034187a5d753p-1, -0x1.cd9786f046dbdp-1, -0x1.4988f615d8fc2p-4, 0x0p+0};
            }
            if constexpr (order == 2) {
                return {0x0p+0, -0x1.32449da1b9609p-55, -0x1.6fb51174aa8e3p-40, 0x1.64323c063686ap-30, -0x1.6377af73884cp-25, -0x1.6928b0e495d08p-19, 0x1.0136c3fb4b19cp-15, 0x1.af6480e4adc47p-17, -0x1.cd1c35ab2ad93p-10, 0x1.21e04958552e7p-7, -0x1.6c5861062a836p-7, -0x1.41156e1c846c4p-5, 0x1.4f99902d11603p-3, -0x1.78ed1cb9b0ec2p-3, -0x1.150c77d31565dp-2, 0x1.5f04ae9c9034bp+0, -0x1.6a2a8758fd08cp+1, 0x1.08f3303fd9b5fp+2, -0x1.a41fd4c6e5304p+1, 0x1.307a2a8e783eap-1, 0x1.5c07e5083d709p-2, 0x0p+0};
            }
            if constexpr (order == 3) {
                return {0x0p+0, 0x1.e2c7778c519ap-49, 0x1.21dedfd7476c5p-34, -0x1.8b6cb6bf7586bp-25, 0x1.f904916865d4dp-21, 0x1.950dcb7c8e982p-15, -0x1.1c307da269deep-11, -0x1.f33a86db42466p-12, 0x1.35b764060a714p-6, -0x1.2b8bd1a4fd0ebp-4, 0x1.3f8a268debfb7p-4, 0x1.08b97b33caeccp-2, -0x1.339d144fd9551p+0, 0x1.38ad9781c7501p+1, -0x1.80597d49465bfp+1, 0x1.0fb2290536b32p+1, -0x1.02364935e2223p-1, 0x1.8f506e2a5e0ebp-1, -0x1.8e98178e99b32p+1, 0x1.a72bc3c00615fp+1, -0x1.1c297e1a72297p+0, -0x0p+0};
            }
            if constexpr (order == 4) {
                return {0x0p+0, -0x1.8ec06f9de2f7dp-39, -0x1.df066d7576ffdp-26, 0x1.ba1b21857af83p-17, -0x1.695b63166e6d6p-13, -0x1.22b5d07d79f25p-7, 0x1.002e9aac59197p-4, 0x1.8f713f264dac3p-5, -0x1.9e20d1654131ap+0, 0x1.8e75251a7f5bfp+2, -0x1.4b3d36d8f9771p+3, 0x1.64c18afe2844cp+1, 0x1.6d1975fc2f42ep+4, -0x1.aeeda88b7fe88p+5, 0x1.1abf8d1adcc3ep+6, -0x1.17ce56ad9d902p+6, 0x1.f144c4df7513dp+5, -0x1.9aba40cee29c1p+5, 0x1.0e79522ccf3ddp+5, -0x1.c33b049bfbce6p+3, 0x1.51d8dbccf9207p+1, -0x0p+0};
            }
            if constexpr (order == 5) {
                return {0x0p+0, 0x1.ca2f0f4b75006p-31, 0x1.136d2374e3712p-18, -0x1.7e95e283d8ce3p-10, 0x1.230c93cadc0b8p-7, 0x1.1851a9519e35ep-1, -0x1.7dbb7d4fa2b47p+1, -0x1.6613709f48646p-1, 0x1.6db5a82568a44p+5, -0x1.4eba06db33d53p+7, 0x1.4159e0e62036cp+8, -0x1.90ad4b61c673ep+8, 0x1.80bbefefb8959p+8, -0x1.5f18c6755ba05p+8, 0x1.2e04ab290e2ep+8, -0x1.4a282d1320ed7p+7, -0x1.36cbd86221f4cp+4, 0x1.d2259be1a7819p+6, -0x1.7429c31b8bcdcp+6, 0x1.1648a261cfac8p+5, -0x1.553841cb8504cp+2, 0x0p+0};
            }
            if constexpr (order == 6) {
                return {0x0p+0, 0x1.ebdb45ca159a6p-30, 0x1.28216d550397dp-18, -0x1.81c4231f27f24p-10, -0x1.1e667fbdaff18p-11, 0x1.1ee7d490234e2p-2, -0x1.5e908a00f74dap+0, 0x1.14b6d949957a5p+0, 0x1.f4f5eade96654p+2, -0x1.cfad0a751fefdp+4, 0x1.cafceab55c222p+5, -0x1.7f74fc0a8fbbcp+6, 0x1.2da088f392f17p+7, -0x1.9299f7098407ep+7, 0x1.c24b2cb18bb3p+7, -0x1.dd46b6e14bd99p+7, 0x1.e7564ac670465p+7, -0x1.917bfbcd50104p+7, 0x1.be9eb07640a6ep+6, -0x1.2015783789a13p+5, 0x1.43e579b56aec7p+2, -0x0p+0};
            }
            if constexpr (order == 7) {
                return {0x0p+0, -0x1.5739a9ee6eb1dp-26, -0x1.9e94a90cc2ee8p-16, 0x1.cc0f8b0b31976p-8, 0x1.89a1e9b6409a4p-6, -0x1.521b96f1f012p-1, 0x1.6d63c669ef72ap+1, -0x1.29abec5ca24e5p+2, 0x1.0097ba319005fp+0, 0x1.a5a3f6309aca8p+2, -0x1.d23af72a92131p+3, 0x1.2bd7eeccbbd8bp+5, -0x1.4b192d711636ap+6, 0x1.d61ad68084632p+6, -0x1.00a2f00b6f36p+7, 0x1.356154b8ce2cp+7, -0x1.7fbdfc1f5c92bp+7, 0x1.5d87f4fd9af88p+7, -0x1.8cc935a5962f3p+6, 0x1.f54b1dc71e94cp+4, -0x1.0e918e0751a5fp+2, -0x0p+0};
            }
            if constexpr (order == 8) {
                return {0x0p+0, -0x1.6c2d94cde3ae2p-22, -0x1.baa159e67dd4fp-13, 0x1.b3efb441a4e1bp-6, 0x1.3f1c58e4e2a68p-4, -0x1.3d9d88a64fae4p+0, 0x1.152961e6e0f57p+2, -0x1.1522f77734c2fp+3, 0x1.e00e48c9f7714p+3, -0x1.a72562d5d8065p+4, 0x1.483818775e0cbp+5, -0x1.b51d99f667842p+5, 0x1.14413c8e18076p+6, -0x1.4933025bdb689p+6, 0x1.678e1d36597aep+6, -0x1.84d65b83affb6p+6, 0x1.91c80b52c3671p+6, -0x1.452573aebb3bfp+6, 0x1.5bbb8c487ae5cp+5, -0x1.a9c9682497d43p+3, 0x1.c2a4ffe966119p+0, -0x0p+0};
            }
            if constexpr (order == 9) {
                return {0x0p+0, 0x1.b40e785cb8808p-22, 0x1.0c4738028aec9p-13, -0x1.7c528a3f2156cp-8, -0x1.c9b3ef0688e89p-8, 0x1.0b62472bd3d78p-3, -0x1.61b17ae5ff0c1p-2, 0x1.624785bf15dbep-1, -0x1.cb38bb2c5f4bfp+0, 0x1.c75ab5e1afc76p+1, -0x1.3ef834c3289aap+2, 0x1.f628948ceff74p+2, -0x1.af56b1534f29ep+3, 0x1.0c0dca303e8e2p+4, -0x1.0953a7a524acbp+4, 0x1.536d6f833b667p+4, -0x1.d3d44a6a16adbp+4, 0x1.bc60c3b54f90bp+4, -0x1.f8f411bacc9f9p+3, 0x1.38fe2f09aac58p+2, -0x1.480f31520e3f3p-1, 0x0p+0};
            }
            if constexpr (order == 10) {
                return {0x0p+0, -0x1.09aed93888304p-21, -0x1.4eea023db8bd5p-14, 0x1.71e7a3b4d88bfp-10, -0x1.a3114fdf2fbcp-12, -0x1.b7227d5381373p-7, 0x1.8d3b211eb839bp-6, -0x1.50a6f5e8f3b31p-5, 0x1.68181929e64c5p-3, -0x1.7340e7a399e64p-2, 0x1.a30d20b4e67c7p-2, -0x1.7019bc278ffcp-1, 0x1.8867eccd39d1bp+0, -0x1.db8c7f7bf2c3ap+0, 0x1.96ff87cfe4cc5p+0, -0x1.2e8a4ebac6b5ep+1, 0x1.fbc727fc23efp+1, -0x1.04a53246109b7p+2, 0x1.2f92842080f08p+1, -0x1.79dddd860d214p-1, 0x1.8a25fe79d1c4ep-4, 0x0p+0};
            }
         }
        if constexpr (std::is_same_v<Real, long double>) {
            if constexpr (order == 0) {
                return {0x0p+0L, -0xb.2ec22065ab6fd2ep-73L, -0xd.6c8a59dce3224ep-56L, -0xa.7ca26b64ba505d1p-45L, -0x8.b00943f903ee402p-38L, 0x8.a04e5a4ae7d54e6p-31L, -0x8.57ac00480665f53p-29L, 0xb.1ab170d8c6d3fd6p-23L, 0xa.e13a71702532bafp-21L, -0xe.bf4dd3d9d6239bp-18L, 0xf.95fd567187ac541p-18L, 0xe.c8e68cc12fa31cbp-17L, 0x8.0c15064d225d41p-12L, -0x8.38830b3b37cedbcp-9L, 0xe.c29805474d2cedap-8L, -0x8.bf03506f5837fe3p-6L, 0x8.3ee4625ac404c6dp-5L, -0xe.9ede116d99b1966p-5L, 0xd.03a753418d0bce6p-4L, 0xe.d43d6faa2ef51b8p-5L, 0xf.c16c3076988ebbcp-10L, 0x0p+0L};
            }
            if constexpr (order == 1) {
                return {0x0p+0L, -0x9.4bc17cc5cc6d1d6p-65L, -0xb.28dea82f2d5572ep-49L, 0xe.ca19d1290517188p-39L, -0xb.9ce0aee4062aafdp-33L, -0xb.edbc2605809849dp-28L, 0xf.01423b9510f7584p-23L, 0xe.26c0f92ddd5b332p-22L, -0xd.15679d8fce83a4ap-17L, 0xb.49882d15db83441p-15L, 0x8.c2bc082b308154p-14L, -0x8.49604b074075f3ap-10L, 0xf.f6403eb5e6f17b6p-9L, -0x8.780783a3abc1044p-7L, 0x9.90e35acb85cca72p-7L, 0xb.8e83f7c0829f529p-15L, -0x9.81a1b1846300118p-6L, 0x8.3af1de78eb8f98bp-6L, 0xf.81a0c3d2eba999fp-4L, -0xe.6cbc378236de837p-4L, -0xa.4c47b0aec7e1081p-7L, 0x0p+0L};
            }
            if constexpr (order == 2) {
                return {0x0p+0L, -0x9.9224ed0dcb04955p-58L, -0xb.7da88ba55471baap-43L, 0xb.2191e031b4353cfp-33L, -0xb.1bbd7b9c425fd7p-28L, -0xb.49458724ae83fdap-22L, 0x8.09b61fda58cdcd9p-18L, 0xd.7b2407256e2365dp-20L, -0xe.68e1ad5956c95f9p-13L, 0x9.0f024ac2a973782p-10L, -0xb.62c30831541b0b7p-10L, -0xa.08ab70e42362325p-8L, 0xa.7ccc81688b01b9cp-6L, -0xb.c768e5cd8760f3bp-6L, -0x8.a863be98ab2e9afp-5L, 0xa.f82574e481a558dp-3L, -0xb.51543ac7e845d9fp-2L, 0x8.479981fecdaf816p-1L, -0xd.20fea63729823dcp-2L, 0x9.83d15473c1f50fdp-4L, 0xa.e03f2841eb84819p-5L, 0x0p+0L};
            }
            if constexpr (order == 3) {
                return {0x0p+0L, 0xf.163bbc628cd01e9p-52L, 0x9.0ef6feba3b62aaap-37L, -0xc.5b65b5fbac3587ap-28L, 0xf.c8248b432ea6598p-24L, 0xc.a86e5be474c0e96p-18L, -0x8.e183ed134ef71abp-14L, -0xf.99d436da1232d03p-15L, 0x9.adbb20305389e61p-9L, -0x9.5c5e8d27e8757a5p-7L, 0x9.fc51346f5fdbbd1p-7L, 0x8.45cbd99e57663b5p-5L, -0x9.9ce8a27ecaa85dcp-3L, 0x9.c56cbc0e3a809b2p-2L, -0xc.02cbea4a32dfa62p-2L, 0x8.7d914829b598f95p-2L, -0x8.11b249af1111423p-4L, 0xc.7a837152f075af7p-4L, -0xc.74c0bc74cd98ea3p-2L, 0xd.395e1e0030afb9ep-2L, -0x8.e14bf0d3914bbd4p-3L, -0x0p+0L};
            }
            if constexpr (order == 4) {
                return {0x0p+0L, -0xc.76037cef17bea41p-42L, -0xe.f8336babb7fe566p-29L, 0xd.d0d90c2bd7c17dcp-20L, -0xb.4adb18b3736b1d7p-16L, -0x9.15ae83ebcf927bfp-10L, 0x8.0174d562c8cba6dp-7L, 0xc.7b89f9326d61459p-8L, -0xc.f1068b2a098cfb7p-3L, 0xc.73a928d3fadf73cp-1L, -0xa.59e9b6c7cbb8775p+0L, 0xb.260c57f14226301p-2L, 0xb.68cbafe17a170d9p+1L, -0xd.776d445bff43fe9p+2L, 0x8.d5fc68d6e61f284p+3L, -0x8.be72b56cec81364p+3L, 0xf.8a2626fba89eb9ep+2L, -0xc.d5d2067714e04fcp+2L, 0x8.73ca916679ee85cp+2L, -0xe.19d824dfde72f6dp+0L, 0xa.8ec6de67c9035e6p-2L, -0x0p+0L};
            }
            if constexpr (order == 5) {
                return {0x0p+0L, 0xe.51787a5ba8032eep-34L, 0x8.9b691ba71b88c22p-21L, -0xb.f4af141ec671a67p-13L, 0x9.18649e56e05bc3p-10L, 0x8.c28d4a8cf1aeca5p-4L, -0xb.eddbea7d15a390ap-2L, -0xb.309b84fa4322d85p-4L, 0xb.6dad412b4521e9bp+2L, -0xa.75d036d99ea9702p+4L, 0xa.0acf073101b6074p+5L, -0xc.856a5b0e339f38dp+5L, 0xc.05df7f7dc4ac954p+5L, -0xa.f8c633aadd02508p+5L, 0x9.702559487170317p+5L, -0xa.51416899076b9c3p+4L, -0x9.b65ec3110fa60e2p+1L, 0xe.912cdf0d3c0c767p+3L, -0xb.a14e18dc5e6e3bcp+3L, 0x8.b245130e7d6435p+2L, -0xa.a9c20e5c28262cfp-1L, 0x0p+0L};
            }
            if constexpr (order == 6) {
                return {0x0p+0L, 0xf.5eda2e50acd327dp-33L, 0x9.410b6aa81cbea18p-21L, -0xc.0e2118f93f91f36p-13L, -0x8.f333fded7f8c399p-14L, 0x8.f73ea4811a71163p-5L, -0xa.f4845007ba6cfc8p-3L, 0x8.a5b6ca4cabd2b23p-3L, 0xf.a7af56f4b329d5ep-1L, -0xe.7d6853a8ff7e50ap+1L, 0xe.57e755aae11112fp+2L, -0xb.fba7e0547ddde95p+3L, 0x9.6d04479c978badep+4L, -0xc.94cfb84c203ee2p+4L, 0xe.1259658c5d981b4p+4L, -0xe.ea35b70a5ecc9fep+4L, 0xf.3ab256338232753p+4L, -0xc.8bdfde6a8082091p+4L, 0xd.f4f583b2053726ep+3L, -0x9.00abc1bc4d0982ep+2L, 0xa.1f2bcdab57636dcp-1L, -0x0p+0L};
            }
            if constexpr (order == 7) {
                return {0x0p+0L, -0xa.b9cd4f73758e428p-29L, -0xc.f4a5486617740b1p-19L, 0xe.607c58598cbad8fp-11L, 0xc.4d0f4db204d1cc9p-9L, -0xa.90dcb78f808ffecp-4L, 0xb.6b1e334f7b950afp-2L, -0x9.4d5f62e5127240ep-1L, 0x8.04bdd18c802f85cp-3L, 0xd.2d1fb184d6542c6p-1L, -0xe.91d7b95490989a7p+0L, 0x9.5ebf7665dec59a9p+2L, -0xa.58c96b88b1b4fe8p+3L, 0xe.b0d6b4042318ee3p+3L, -0x8.0517805b79afde9p+4L, 0x9.ab0aa5c6715fdd9p+4L, -0xb.fdefe0fae49564ap+4L, 0xa.ec3fa7ecd7c41f4p+4L, -0xc.6649ad2cb179936p+3L, 0xf.aa58ee38f4a5c91p+1L, -0x8.748c703a8d2f937p-1L, -0x0p+0L};
            }
            if constexpr (order == 8) {
                return {0x0p+0L, -0xb.616ca66f1d70df8p-25L, -0xd.d50acf33eea750ap-16L, 0xd.9f7da20d270d536p-9L, 0x9.f8e2c727153415ep-7L, -0x9.ecec45327d721f6p-3L, 0x8.a94b0f3707ab7abp-1L, -0x8.a917bbb9a6176d9p+0L, 0xf.0072464fbb8a1e7p+0L, -0xd.392b16aec032585p+1L, 0xa.41c0c3baf0658a2p+2L, -0xd.a8eccfb33c210d9p+2L, 0x8.a209e470c03af8ap+3L, -0xa.499812dedb44ac4p+3L, 0xb.3c70e9b2cbd6e2ap+3L, -0xc.26b2dc1d7fdae55p+3L, 0xc.8e405a961b38a93p+3L, -0xa.292b9d75d9dfa3ap+3L, 0xa.dddc6243d72ddb1p+2L, -0xd.4e4b4124bea1ab6p+0L, 0xe.1527ff4b308c86p-3L, -0x0p+0L};
            }
            if constexpr (order == 9) {
                return {0x0p+0L, 0xd.a073c2e5c404378p-25L, 0x8.6239c0145764987p-16L, -0xb.e29451f90ab62c4p-11L, -0xe.4d9f7834474485cp-11L, 0x8.5b12395e9ebc248p-6L, -0xb.0d8bd72ff86061cp-5L, 0xb.123c2df8aedf04ap-4L, -0xe.59c5d962fa5f5d7p-3L, 0xe.3ad5af0d7e3ad17p-2L, -0x9.f7c1a61944d51aep-1L, 0xf.b144a4677fb9db7p-1L, -0xd.7ab58a9a794ee91p+0L, 0x8.606e5181f471318p+1L, -0x8.4a9d3d29256548ap+1L, 0xa.9b6b7c19db33b9dp+1L, -0xe.9ea25350b56d92ap+1L, 0xd.e3061daa7c85492p+1L, -0xf.c7a08dd664fcb0ep+0L, 0x9.c7f1784d562bdc9p-1L, -0xa.40798a9071f9886p-4L, 0x0p+0L};
            }
            if constexpr (order == 10) {
                return {0x0p+0L, -0x8.4d76c9c441821a3p-24L, -0xa.775011edc5ea48p-17L, 0xb.8f3d1da6c45f8b5p-13L, -0xd.188a7ef97de034bp-15L, -0xd.b913ea9c09b9ab5p-10L, 0xc.69d908f5c1cd609p-9L, -0xa.8537af479d9846bp-8L, 0xb.40c0c94f326298p-6L, -0xb.9a073d1ccf31d0dp-5L, 0xd.186905a733e3482p-5L, -0xb.80cde13c7fe01dap-4L, 0xc.433f6669ce8d77ap-3L, -0xe.dc63fbdf961d121p-3L, 0xc.b7fc3e7f2662817p-3L, -0x9.745275d635af24dp-2L, 0xf.de393fe11f77db2p-2L, -0x8.2529923084db5b6p-1L, 0x9.7c942104078413dp-2L, -0xb.ceeeec30690a3d4p-4L, 0xc.512ff3ce8e26ed6p-7L, 0x0p+0L};
            }
         }
#ifdef BOOST_HAS_FLOAT128
        if constexpr (std::is_same_v<Real, boost::multiprecision::float128>) {
            if constexpr (order == 0) {
                return {0x0p+0Q, -0x1.65d8440cb56dfa5bf9c15ab7a7b4p-70Q, -0x1.ad914b3b9c6449bf320bb309094p-53Q, -0x1.4f944d6c974a0ba102391a518b7cp-42Q, -0x1.1601287f207dc804b327d1ac7b7bp-35Q, 0x1.1409cb495cfaa9cc50f200c35599p-28Q, -0x1.0af5800900ccbea5a64700cabdf8p-26Q, 0x1.63562e1b18da7fac5905b8871f2bp-20Q, 0x1.5c274e2e04a6575e26e75abc9b19p-18Q, -0x1.d7e9ba7b3ac4736029b54851fd9fp-15Q, 0x1.f2bfaace30f58a825ed0b3a180adp-15Q, 0x1.d91cd19825f46395eba40bde3584p-14Q, 0x1.0182a0c9a44ba81f2b490da8801fp-9Q, -0x1.0710616766f9db787dd7848035fap-6Q, 0x1.d85300a8e9a59db3d8e3285446f6p-5Q, -0x1.17e06a0deb06ffc5a354422c28ebp-3Q, 0x1.07dc8c4b588098d99d74923feefdp-2Q, -0x1.d3dbc22db33632cb9d3eeea9305bp-2Q, 0x1.a074ea6831a179cbcefe13bdeb15p-1Q, 0x1.da87adf545dea36f29603c502041p-2Q, 0x1.f82d860ed311d7772f0a5fc02a37p-7Q, 0x0p+0Q};
            }
            if constexpr (order == 1) {
                return {0x0p+0Q, -0x1.29782f98b98da3ab04cf7e9484d3p-62Q, -0x1.651bd505e5aaae5cbe740745ad88p-46Q, 0x1.d9433a2520a2e30f13e84b94da1bp-36Q, -0x1.739c15dc80c555f970d15125701ap-30Q, -0x1.7db784c0b013093ad9874bfa593bp-25Q, 0x1.e0284772a21eeb07c87b76463408p-20Q, 0x1.c4d81f25bbab666371fdd703749p-19Q, -0x1.a2acf3b1f9d0749452e4dc44d399p-14Q, 0x1.693105a2bb706882bf9e39158cb4p-12Q, 0x1.1857810566102a80b03f6d573ebfp-11Q, -0x1.092c0960e80ebe73a0fb5f592b47p-7Q, 0x1.fec807d6bcde2f6b9da93df4ab99p-6Q, -0x1.0f00f0747578208709a35dbb57d5p-4Q, 0x1.321c6b5970b994e47d89204a797cp-4Q, 0x1.71d07ef81053ea5196f6a9c6bb88p-12Q, -0x1.303436308c600230256305532c1ep-3Q, 0x1.075e3bcf1d71f31524cf997af806p-3Q, 0x1.f034187a5d75333e9dbb5b104989p-1Q, -0x1.cd9786f046dbd06d9795ae53077p-1Q, -0x1.4988f615d8fc2101372678651801p-4Q, 0x0p+0Q};
            }
            if constexpr (order == 2) {
                return {0x0p+0Q, -0x1.32449da1b96092a92afa9c91b13ep-55Q, -0x1.6fb51174aa8e375496c43a8d2ba4p-40Q, 0x1.64323c063686a79e3f724a84a7b2p-30Q, -0x1.6377af73884bfae0ae854936cd2p-25Q, -0x1.6928b0e495d07fb474aad6c72b4bp-19Q, 0x1.0136c3fb4b19b9b10157d7ad132fp-15Q, 0x1.af6480e4adc46cb9e0c764772afcp-17Q, -0x1.cd1c35ab2ad92bf1561885236608p-10Q, 0x1.21e04958552e6f03502dd793e47dp-7Q, -0x1.6c5861062a83616dee9275c25f12p-7Q, -0x1.41156e1c846c464a303629e36dabp-5Q, 0x1.4f99902d11603738023acee12ff8p-3Q, -0x1.78ed1cb9b0ec1e760990f81cc249p-3Q, -0x1.150c77d31565d35edf80240e3b0ep-2Q, 0x1.5f04ae9c9034ab19fcfae00d6dbbp+0Q, -0x1.6a2a8758fd08bb3ed6b34742cbb9p+1Q, 0x1.08f3303fd9b5f02cbec87edc987fp+2Q, -0x1.a41fd4c6e53047b8904080871597p+1Q, 0x1.307a2a8e783ea1f99bb361a3c11cp-1Q, 0x1.5c07e5083d7090314f404490a308p-2Q, 0x0p+0Q};
            }
            if constexpr (order == 3) {
                return {0x0p+0Q, 0x1.e2c7778c519a03d1d1aa2da8d68bp-49Q, 0x1.21dedfd7476c555324241dc23af1p-34Q, -0x1.8b6cb6bf7586b0f4f1fce59263a7p-25Q, 0x1.f904916865d4cb30ef3840ef8eabp-21Q, 0x1.950dcb7c8e981d2c1ff88542b47bp-15Q, -0x1.1c307da269dee355c5c2d33e95b8p-11Q, -0x1.f33a86db42465a05a6ee3ed23fd3p-12Q, 0x1.35b764060a713cc242948b5b8db6p-6Q, -0x1.2b8bd1a4fd0eaf4aaa08a96c47f2p-4Q, 0x1.3f8a268debfb77a1cde76d546205p-4Q, 0x1.08b97b33caecc76a026ab832c03p-2Q, -0x1.339d144fd9550bb8463b5449b143p+0Q, 0x1.38ad9781c750136475272cd61fa1p+1Q, -0x1.80597d49465bf4c4d8869fb3fc1fp+1Q, 0x1.0fb2290536b31f2a2fbff33eddb1p+1Q, -0x1.02364935e2222846b9cdda15d684p-1Q, 0x1.8f506e2a5e0eb5edc0313f0b7341p-1Q, -0x1.8e98178e99b31d46ae26eebdf89bp+1Q, 0x1.a72bc3c00615f73ba93716ae2561p+1Q, -0x1.1c297e1a722977a767a1ff215d25p+0Q, -0x0p+0Q};
            }
            if constexpr (order == 4) {
                return {0x0p+0Q, -0x1.8ec06f9de2f7d48145d6ed1dbc3fp-39Q, -0x1.df066d7576ffcacbc6c975c4657bp-26Q, 0x1.ba1b21857af82fb76bf46f9a62f4p-17Q, -0x1.695b63166e6d63ad77759be760aap-13Q, -0x1.22b5d07d79f24f7de8397d0660dcp-7Q, 0x1.002e9aac591974da848a13bedefp-4Q, 0x1.8f713f264dac28b16bf918ff69c4p-5Q, -0x1.9e20d16541319f6d732cb19d28f8p+0Q, 0x1.8e75251a7f5bee7754c2a5097887p+2Q, -0x1.4b3d36d8f9770ee95e271e3f87c5p+3Q, 0x1.64c18afe2844c60255dd6c133594p+1Q, 0x1.6d1975fc2f42e1b2096261ca874dp+4Q, -0x1.aeeda88b7fe87fd16bcf8b24444bp+5Q, 0x1.1abf8d1adcc3e5071c7207af9e5bp+6Q, -0x1.17ce56ad9d9026c80fe902e6692p+6Q, 0x1.f144c4df7513d73b0fd718cace5fp+5Q, -0x1.9aba40cee29c09f71dc387a81452p+5Q, 0x1.0e79522ccf3dd0b86ea89759f6e5p+5Q, -0x1.c33b049bfbce5ed995e57b93ee75p+3Q, 0x1.51d8dbccf9206bcc4420e5965aa6p+1Q, -0x0p+0Q};
            }
            if constexpr (order == 5) {
                return {0x0p+0Q, 0x1.ca2f0f4b750065dc047ccc5e827ep-31Q, 0x1.136d2374e371184371f15ac8673ap-18Q, -0x1.7e95e283d8ce34ce8180a22d17f3p-10Q, 0x1.230c93cadc0b78607a6dc47e053p-7Q, 0x1.1851a9519e35d94904832cf7a4ccp-1Q, -0x1.7dbb7d4fa2b472146ed0ef453d1bp+1Q, -0x1.6613709f48645b09020ce0e5cf49p-1Q, 0x1.6db5a82568a43d359e12e8016c45p+5Q, -0x1.4eba06db33d52e049a495d989071p+7Q, 0x1.4159e0e62036c0e7d3a9017c60dep+8Q, -0x1.90ad4b61c673e71a7b1a1e950968p+8Q, 0x1.80bbefefb89592a773121cad9772p+8Q, -0x1.5f18c6755ba04a0fe58426b7577ap+8Q, 0x1.2e04ab290e2e062d601864f62d8dp+8Q, -0x1.4a282d1320ed73864b5c18dd9bb5p+7Q, -0x1.36cbd86221f4c1c381312344c114p+4Q, 0x1.d2259be1a7818ecd1bf97ff7d69cp+6Q, -0x1.7429c31b8bcdc777e0683b89fc4ep+6Q, 0x1.1648a261cfac86a07a575eca826p+5Q, -0x1.553841cb8504c59ec882711a9f8p+2Q, 0x0p+0Q};
            }
            if constexpr (order == 6) {
                return {0x0p+0Q, 0x1.ebdb45ca159a64fac22d479bba89p-30Q, 0x1.28216d550397d430f87787c4522bp-18Q, -0x1.81c4231f27f23e6bc64802e3f79cp-10Q, -0x1.1e667fbdaff187310af84e56a63fp-11Q, 0x1.1ee7d490234e22c547abf74a2774p-2Q, -0x1.5e908a00f74d9f8f6b638a72bdb4p+0Q, 0x1.14b6d949957a5646a54dc866fb8dp+0Q, 0x1.f4f5eade96653abb41cc870b3d96p+2Q, -0x1.cfad0a751fefca13c3b1245a80c1p+4Q, 0x1.cafceab55c22225d33cf4383b247p+5Q, -0x1.7f74fc0a8fbbbd2a66275a043fd8p+6Q, 0x1.2da088f392f175bbf14ea02b8328p+7Q, -0x1.9299f7098407dc40f69d2557cb5cp+7Q, 0x1.c24b2cb18bb303674bf40181cfffp+7Q, -0x1.dd46b6e14bd993fb7d25442ccc13p+7Q, 0x1.e7564ac670464ea6f5b828f486a4p+7Q, -0x1.917bfbcd50104121d7ccc0c5035cp+7Q, 0x1.be9eb07640a6e4dcf71bd794877fp+6Q, -0x1.2015783789a1305c4b0337d90efdp+5Q, 0x1.43e579b56aec6db8f44dc1cc046ep+2Q, -0x0p+0Q};
            }
            if constexpr (order == 7) {
                return {0x0p+0Q, -0x1.5739a9ee6eb1c850ac62e0d49d47p-26Q, -0x1.9e94a90cc2ee8162e5e960ae2735p-16Q, 0x1.cc0f8b0b31975b1da178761d2e23p-8Q, 0x1.89a1e9b6409a39917b3444a5488ep-6Q, -0x1.521b96f1f011ffd83508bc3cdf05p-1Q, 0x1.6d63c669ef72a15deadc3267c4d1p+1Q, -0x1.29abec5ca24e481b15102195037bp+2Q, 0x1.0097ba319005f0b8250a5e99956ap+0Q, 0x1.a5a3f6309aca858b554726858ecbp+2Q, -0x1.d23af72a9213134d138b60ca891ap+3Q, 0x1.2bd7eeccbbd8b3529ca232ebfc68p+5Q, -0x1.4b192d7116369fcf46c13ff5f7f6p+6Q, 0x1.d61ad68084631dc5b1d352c01432p+6Q, -0x1.00a2f00b6f35fbd2365fdf4a3308p+7Q, 0x1.356154b8ce2bfbb26ff7f681a29ep+7Q, -0x1.7fbdfc1f5c92ac9439fec54897d6p+7Q, 0x1.5d87f4fd9af883e7ca96bbbdc2dep+7Q, -0x1.8cc935a5962f326b8770312223fp+6Q, 0x1.f54b1dc71e94b922c8c248ee1d5ap+4Q, -0x1.0e918e0751a5f26d66d8bbb66469p+2Q, -0x0p+0Q};
            }
            if constexpr (order == 8) {
                return {0x0p+0Q, -0x1.6c2d94cde3ae1bef5e2de4ec940fp-22Q, -0x1.baa159e67dd4ea14a45d44a01585p-13Q, 0x1.b3efb441a4e1aa6b7bfa1244a249p-6Q, 0x1.3f1c58e4e2a682bc6ff11f32d654p-4Q, -0x1.3d9d88a64fae43ec46dd5a4c0059p+0Q, 0x1.152961e6e0f56f56720361babbafp+2Q, -0x1.1522f77734c2edb1934560bea09dp+3Q, 0x1.e00e48c9f77143cec0c3e3729a98p+3Q, -0x1.a72562d5d8064b092b7951100b12p+4Q, 0x1.483818775e0cb14480b189318289p+5Q, -0x1.b51d99f6678421b18134e01e3b73p+5Q, 0x1.14413c8e18075f1444c5de04e63cp+6Q, -0x1.4933025bdb6895888c0e5b9032e2p+6Q, 0x1.678e1d36597adc539936f5f8acb9p+6Q, -0x1.84d65b83affb5ca97cbeecbe2f35p+6Q, 0x1.91c80b52c3671525c325c7dd6caap+6Q, -0x1.452573aebb3bf474bf13b76bc15cp+6Q, 0x1.5bbb8c487ae5bb615e7ac55e1d0fp+5Q, -0x1.a9c9682497d4356cc963e6a8b598p+3Q, 0x1.c2a4ffe9661190c09a6113caab73p+0Q, -0x0p+0Q};
            }
            if constexpr (order == 9) {
                return {0x0p+0Q, 0x1.b40e785cb88086ef360bf1e352c1p-22Q, 0x1.0c4738028aec930d43185f6aae41p-13Q, -0x1.7c528a3f2156c58889a9031ef11ap-8Q, -0x1.c9b3ef0688e890b76da6575e8111p-8Q, 0x1.0b62472bd3d7849095fe37975dcap-3Q, -0x1.61b17ae5ff0c0c37a783c4f5d399p-2Q, 0x1.624785bf15dbe09349730f3e9f64p-1Q, -0x1.cb38bb2c5f4bebadf62f69736308p+0Q, 0x1.c75ab5e1afc75a2d9d356e0ab303p+1Q, -0x1.3ef834c3289aa35b262f7899ab76p+2Q, 0x1.f628948ceff73b6d1af4db511dfdp+2Q, -0x1.af56b1534f29dd21b9314098cba6p+3Q, 0x1.0c0dca303e8e262fc9705a8e32d8p+4Q, -0x1.0953a7a524aca9132832bbbd4497p+4Q, 0x1.536d6f833b66773a970c8daba54fp+4Q, -0x1.d3d44a6a16adb2540f7bf8d5e735p+4Q, 0x1.bc60c3b54f90a92419a8b4f6c112p+4Q, -0x1.f8f411bacc9f961bc610901d47bbp+3Q, 0x1.38fe2f09aac57b92aad8fc99fb12p+2Q, -0x1.480f31520e3f310c3077af46d25dp-1Q, 0x0p+0Q};
            }
            if constexpr (order == 10) {
                return {0x0p+0Q, -0x1.09aed93888304346e3035e51fe62p-21Q, -0x1.4eea023db8bd4900ec4f093701b5p-14Q, 0x1.71e7a3b4d88bf16a9f37815c3b2cp-10Q, -0x1.a3114fdf2fbc06953fb8d4ddd4a2p-12Q, -0x1.b7227d5381373569d61e5703ee64p-7Q, 0x1.8d3b211eb839ac111936e9c6a73ap-6Q, -0x1.50a6f5e8f3b308d50e58f5983b0bp-5Q, 0x1.68181929e64c530045dc3f092afap-3Q, -0x1.7340e7a399e63a19859a122eba5cp-2Q, 0x1.a30d20b4e67c690471d6c9051b85p-2Q, -0x1.7019bc278ffc03b406b1b1379edcp-1Q, 0x1.8867eccd39d1aef3728465048ae6p+0Q, -0x1.db8c7f7bf2c3a242351d1c9ade12p+0Q, 0x1.96ff87cfe4cc502ebda7f1f3316ep+0Q, -0x1.2e8a4ebac6b5e49a84f948f21e4ep+1Q, 0x1.fbc727fc23eefb638104ef681d42p+1Q, -0x1.04a53246109b6b6c0d3d0d0e6592p+2Q, 0x1.2f92842080f0827925add964c40bp+1Q, -0x1.79dddd860d2147a741626649523p-1Q, 0x1.8a25fe79d1c4ddac3358f05d288ep-4Q, 0x0p+0Q};
            }
         }
#endif
     }
    if constexpr (p == 12) {
        if constexpr (std::is_same_v<Real, float>) {
            if constexpr (order == 0) {
                return {0x0p+0f, -0x1.f67e04p-79f, 0x1.bb358cp-60f, 0x1.02e74cp-47f, -0x1.652334p-41f, 0x1.0e0d22p-33f, -0x1.3a5d96p-29f, -0x1.934a0ep-24f, -0x1.f0554cp-22f, 0x1.2d678cp-18f, -0x1.dafa2ep-18f, -0x1.43b7cep-15f, 0x1.f00fcp-12f, -0x1.8785b4p-9f, 0x1.86c906p-7f, -0x1.08554cp-5f, 0x1.042a86p-4f, -0x1.9d4fdcp-4f, 0x1.4840a4p-3f, -0x1.5c6f66p-2f, 0x1.da1522p-1f, 0x1.3af728p-2f, 0x1.ba5636p-8f, 0x0p+0f};
            }
            if constexpr (order == 1) {
                return {0x0p+0f, -0x1.13fd18p-69f, 0x1.e6d8e6p-52f, -0x1.71c42cp-44f, -0x1.b076eap-34f, 0x1.c87e12p-29f, 0x1.9ee91p-25f, -0x1.b2cff8p-20f, 0x1.13c0d6p-17f, 0x1.6a167p-15f, -0x1.c45b6ap-12f, 0x1.c5c0dp-11f, 0x1.3d712ap-9f, -0x1.025a88p-6f, 0x1.fa47fap-6f, -0x1.805958p-10f, -0x1.2402bcp-3f, 0x1.b8eefp-2f, -0x1.97243cp-1f, 0x1.f51bdap-1f, 0x1.16beeap-2f, -0x1.716f56p-1f, -0x1.3daa14p-5f, 0x0p+0f};
            }
            if constexpr (order == 2) {
                return {0x0p+0f, -0x1.a7d232p-65f, 0x1.75cc98p-48f, -0x1.84d7e4p-38f, -0x1.63226ep-31f, 0x1.0f785p-26f, 0x1.a3887p-21f, -0x1.da6da6p-19f, 0x1.09bd2p-15f, 0x1.90febp-15f, -0x1.8fe152p-10f, 0x1.935e46p-8f, -0x1.59c476p-7f, 0x1.091b1ep-6f, -0x1.4654eep-4f, 0x1.50c7e8p-2f, -0x1.abf5ep-1f, 0x1.80cee2p+0f, -0x1.240fe4p+1f, 0x1.bc153ap+1f, -0x1.9e9c7cp+1f, 0x1.e1d108p-1f, 0x1.77fed2p-3f, -0x0p+0f};
            }
            if constexpr (order == 3) {
                return {0x0p+0f, 0x1.051a2p-54f, -0x1.cc8984p-39f, 0x1.7e8cfp-32f, 0x1.9c90e4p-23f, -0x1.28c92p-18f, -0x1.f86e54p-15f, 0x1.dfc88ap-11f, -0x1.206bp-9f, -0x1.d58bf4p-7f, 0x1.a2cbbap-4f, -0x1.ffde0cp-3f, 0x1.ee149ap-4f, 0x1.c702ecp-1f, -0x1.5e6794p+1f, 0x1.f3ab98p+1f, -0x1.3621acp+1f, -0x1.c1feacp+0f, 0x1.80e0dep+2f, -0x1.8d30cp+2f, 0x1.f023b2p+0f, 0x1.1c49e2p+0f, -0x1.60087cp-1f, 0x0p+0f};
            }
            if constexpr (order == 4) {
                return {0x0p+0f, -0x1.1d1866p-45f, 0x1.f6c8c6p-31f, -0x1.d59d4cp-25f, -0x1.c361dcp-16f, 0x1.b6ad0ep-12f, 0x1.68d464p-8f, -0x1.003386p-4f, 0x1.17c48ap-3f, 0x1.335fd6p-1f, -0x1.0b934ap+2f, 0x1.52a34ep+3f, -0x1.822b1ap+3f, -0x1.90f39p+1f, 0x1.142ce2p+5f, -0x1.f8d9fp+5f, 0x1.1a6ca8p+6f, -0x1.d72d9cp+5f, 0x1.5a4b7p+5f, -0x1.014d7ap+5f, 0x1.596114p+4f, -0x1.34df3ep+3f, 0x1.efa02p+0f, -0x0p+0f};
            }
            if constexpr (order == 5) {
                return {0x0p+0f, 0x1.aa8cd8p-37f, -0x1.78041p-23f, 0x1.3e3df8p-18f, 0x1.5137c2p-9f, -0x1.951dp-6f, -0x1.43d702p-2f, 0x1.6e4fp+1f, -0x1.8ab9ccp+2f, -0x1.95561cp+3f, 0x1.8ddf94p+6f, -0x1.055f48p+8f, 0x1.9cacd4p+8f, -0x1.cb4de8p+8f, 0x1.a2714ep+8f, -0x1.628082p+8f, 0x1.fbf398p+7f, -0x1.677b64p+6f, -0x1.09c53p+6f, 0x1.e4421p+6f, -0x1.49fc06p+6f, 0x1.d05e84p+4f, -0x1.1672f8p+2f, -0x0p+0f};
            }
            if constexpr (order == 6) {
                return {0x0p+0f, 0x1.8b88d4p-34f, -0x1.5c7912p-21f, 0x1.260d96p-17f, 0x1.38d6e4p-8f, -0x1.6b5574p-6f, -0x1.4cc6cp-2f, 0x1.3f67e2p+1f, -0x1.7f694cp+2f, 0x1.8c237ap-1f, 0x1.bc14b2p+4f, -0x1.3b920ep+6f, 0x1.1f2d7ap+7f, -0x1.caa92ap+7f, 0x1.4bb5b2p+8f, -0x1.958a02p+8f, 0x1.aa319ep+8f, -0x1.a2f3c4p+8f, 0x1.7cdc34p+8f, -0x1.16422p+8f, 0x1.18d5fcp+7f, -0x1.50ba26p+5f, 0x1.6688e4p+2f, -0x0p+0f};
            }
            if constexpr (order == 7) {
                return {0x0p+0f, -0x1.869c12p-29f, 0x1.57bc52p-17f, -0x1.c53d4ap-15f, -0x1.343728p-5f, 0x1.39157p-4f, 0x1.6fa584p+0f, -0x1.21acdep+3f, 0x1.784654p+4f, -0x1.013c24p+5f, 0x1.93c8d6p+4f, -0x1.cc8ab2p+3f, -0x1.09d1ap+3f, 0x1.59e9b6p+6f, -0x1.969994p+7f, 0x1.12635ap+8f, -0x1.2e9626p+8f, 0x1.675a3ap+8f, -0x1.8f56fp+8f, 0x1.41512cp+8f, -0x1.48b0d8p+7f, 0x1.7dd222p+5f, -0x1.80f666p+2f, 0x0p+0f};
            }
            if constexpr (order == 8) {
                return {0x0p+0f, -0x1.6c7e6p-26f, 0x1.400292p-15f, 0x1.a1f2dcp-12f, -0x1.0f174cp-4f, 0x1.e57eccp-4f, 0x1.4e6abp+0f, -0x1.bab6f4p+2f, 0x1.14ef6ep+4f, -0x1.0a8b2ep+5f, 0x1.e8c646p+5f, -0x1.93d454p+6f, 0x1.20e874p+7f, -0x1.84734ap+7f, 0x1.ef2b6p+7f, -0x1.156986p+8f, 0x1.1c64b2p+8f, -0x1.29c6cep+8f, 0x1.23910ep+8f, -0x1.b2204ep+7f, 0x1.a873c2p+6f, -0x1.ddfb9cp+4f, 0x1.d573dcp+1f, 0x0p+0f};
            }
            if constexpr (order == 9) {
                return {0x0p+0f, 0x1.d792bap-25f, -0x1.9c18e6p-15f, -0x1.1fc084p-10f, 0x1.26b33ap-5f, -0x1.63e904p-4f, -0x1.04013ap-2f, 0x1.59e4e8p+0f, -0x1.780b78p+1f, 0x1.ab1578p+2f, -0x1.eced02p+3f, 0x1.9d3e5ep+4f, -0x1.217c84p+5f, 0x1.c4dfc6p+5f, -0x1.4d02e4p+6f, 0x1.73a4dcp+6f, -0x1.759cb8p+6f, 0x1.c308c4p+6f, -0x1.04d398p+7f, 0x1.a7b25ap+6f, -0x1.ab061ep+5f, 0x1.e196ep+3f, -0x1.d2f4ep+0f, -0x0p+0f};
            }
            if constexpr (order == 10) {
                return {0x0p+0f, -0x1.cf434ap-23f, 0x1.910e04p-14f, 0x1.17b65cp-9f, -0x1.b7f686p-6f, 0x1.2e7cdcp-4f, -0x1.4a01e6p-6f, -0x1.07ed62p-3f, 0x1.22ff2ap-3f, -0x1.75b64cp-1f, 0x1.64fab4p+1f, -0x1.21aae8p+2f, 0x1.68255cp+2f, -0x1.5bb8e2p+3f, 0x1.289832p+4f, -0x1.3946dp+4f, 0x1.239464p+4f, -0x1.a397ap+4f, 0x1.1efd1cp+5f, -0x1.f4e8dcp+4f, 0x1.0297dcp+4f, -0x1.249f96p+2f, 0x1.19eec2p-1f, -0x0p+0f};
            }
            if constexpr (order == 11) {
                return {0x0p+0f, -0x1.1b79fap-21f, 0x1.e19176p-14f, 0x1.e72d5p-10f, -0x1.8f4f48p-7f, 0x1.e6ef9ep-6f, -0x1.cec02p-5f, 0x1.07bd18p-3f, -0x1.00deep-2f, 0x1.8936aep-2f, -0x1.329c32p-1f, 0x1.daf36p-1f, -0x1.3cec22p+0f, 0x1.9cd0dcp+0f, -0x1.0e1264p+1f, 0x1.31bb74p+1f, -0x1.4e8702p+1f, 0x1.a9173ap+1f, -0x1.eb2b7ep+1f, 0x1.8aa498p+1f, -0x1.89e13cp+0f, 0x1.b83154p-2f, -0x1.a6c51cp-5f, -0x0p+0f};
            }
         }
        if constexpr (std::is_same_v<Real, double>) {
            if constexpr (order == 0) {
                return {0x0p+0, -0x1.f67e04f07f01ep-79, 0x1.bb358cf4adf4ep-60, 0x1.02e74c07e6081p-47, -0x1.652333b32760ap-41, 0x1.0e0d213886187p-33, -0x1.3a5d95f1928efp-29, -0x1.934a0de4da293p-24, -0x1.f0554b8e5712dp-22, 0x1.2d678c174b94dp-18, -0x1.dafa2eb6c7f5dp-18, -0x1.43b7ce1fa4e5bp-15, 0x1.f00fc00fce26fp-12, -0x1.8785b403c86f7p-9, 0x1.86c905bef28a7p-7, -0x1.08554bd18ab25p-5, 0x1.042a86a360592p-4, -0x1.9d4fdc34ae4eep-4, 0x1.4840a35a3e507p-3, -0x1.5c6f665240adep-2, 0x1.da1521352f488p-1, 0x1.3af72853b189p-2, 0x1.ba5636194c094p-8, 0x0p+0};
            }
            if constexpr (order == 1) {
                return {0x0p+0, -0x1.13fd176d0c1f4p-69, 0x1.e6d8e5785ac94p-52, -0x1.71c42c5b1e647p-44, -0x1.b076e96662045p-34, 0x1.c87e1185ca55bp-29, 0x1.9ee90facd3d13p-25, -0x1.b2cff79438463p-20, 0x1.13c0d643a27dcp-17, 0x1.6a1670435bfddp-15, -0x1.c45b69075cdb6p-12, 0x1.c5c0d04fc8743p-11, 0x1.3d712aaad7776p-9, -0x1.025a88022a414p-6, 0x1.fa47fa3f4b022p-6, -0x1.80595820660c9p-10, -0x1.2402bce395302p-3, 0x1.b8eef0110c9fap-2, -0x1.97243be7bf084p-1, 0x1.f51bd9fa1cf2ep-1, 0x1.16bee99870719p-2, -0x1.716f56bac1b11p-1, -0x1.3daa1461d0cf4p-5, 0x0p+0};
            }
            if constexpr (order == 2) {
                return {0x0p+0, -0x1.a7d23115de58ap-65, 0x1.75cc980993529p-48, -0x1.84d7e34db15d3p-38, -0x1.63226d50629cp-31, 0x1.0f7850b94c92bp-26, 0x1.a3886fa6d48c5p-21, -0x1.da6da6ad5e76cp-19, 0x1.09bd1f2d4e386p-15, 0x1.90feb0894b435p-15, -0x1.8fe1513629012p-10, 0x1.935e464b93aa4p-8, -0x1.59c475f85d7f4p-7, 0x1.091b1dd62c3e6p-6, -0x1.4654ed4660724p-4, 0x1.50c7e71a108bbp-2, -0x1.abf5df7fd7acp-1, 0x1.80cee17390c18p+0, -0x1.240fe3f248861p+1, 0x1.bc153af50c24cp+1, -0x1.9e9c7ce63835dp+1, 0x1.e1d10780fb061p-1, 0x1.77fed2e016936p-3, -0x0p+0};
            }
            if constexpr (order == 3) {
                return {0x0p+0, 0x1.051a1f819156cp-54, -0x1.cc8984cca607fp-39, 0x1.7e8cf0898f333p-32, 0x1.9c90e430666c9p-23, -0x1.28c920e4a3f6cp-18, -0x1.f86e53415f2ap-15, 0x1.dfc88961f6b15p-11, -0x1.206affc103ba1p-9, -0x1.d58bf4f6c76eap-7, 0x1.a2cbb989ecd0ap-4, -0x1.ffde0ca5fc1a1p-3, 0x1.ee149aa686bb3p-4, 0x1.c702eca93a33dp-1, -0x1.5e6794d526a25p+1, 0x1.f3ab98611d47ap+1, -0x1.3621ab0ea1b05p+1, -0x1.c1feac4e70c18p+0, 0x1.80e0de445b21ap+2, -0x1.8d30c07e51574p+2, 0x1.f023b1706ef0ep+0, 0x1.1c49e2797e56dp+0, -0x1.60087b087d2efp-1, 0x0p+0};
            }
            if constexpr (order == 4) {
                return {0x0p+0, -0x1.1d186572e3d31p-45, 0x1.f6c8c69c71b98p-31, -0x1.d59d4ba004a56p-25, -0x1.c361dbb11082cp-16, 0x1.b6ad0e03e206fp-12, 0x1.68d4635bbce68p-8, -0x1.0033850ae8807p-4, 0x1.17c48a3b6f0b8p-3, 0x1.335fd67def29fp-1, -0x1.0b934a511e089p+2, 0x1.52a34d125b1bep+3, -0x1.822b1a96e2d92p+3, -0x1.90f390767ce47p+1, 0x1.142ce1c1f06e2p+5, -0x1.f8d9ef0a2890cp+5, 0x1.1a6ca7d65632dp+6, -0x1.d72d9b409f2b2p+5, 0x1.5a4b6fa251ebp+5, -0x1.014d7a672544ap+5, 0x1.5961138d5be95p+4, -0x1.34df3dac7fcd4p+3, 0x1.efa01fe281c9p+0, -0x0p+0};
            }
            if constexpr (order == 5) {
                return {0x0p+0, 0x1.aa8cd8e6f8972p-37, -0x1.78040f51a9568p-23, 0x1.3e3df8b2da275p-18, 0x1.5137c147339bcp-9, -0x1.951cff21de4acp-6, -0x1.43d7023213a8p-2, 0x1.6e4f007cb3f6fp+1, -0x1.8ab9cc38b4949p+2, -0x1.95561c2d77ea8p+3, 0x1.8ddf94a470b87p+6, -0x1.055f4806a8346p+8, 0x1.9cacd3aee3a7ep+8, -0x1.cb4de777b4f36p+8, 0x1.a2714ddd5c22bp+8, -0x1.628081b6fdd81p+8, 0x1.fbf39786bea7bp+7, -0x1.677b64cbf326ep+6, -0x1.09c52f838d2fdp+6, 0x1.e442104e65336p+6, -0x1.49fc051fc99a8p+6, 0x1.d05e83fe81472p+4, -0x1.1672f76c48a3bp+2, -0x0p+0};
            }
            if constexpr (order == 6) {
                return {0x0p+0, 0x1.8b88d3fa551cp-34, -0x1.5c7911e0053bcp-21, 0x1.260d95329ea8ep-17, 0x1.38d6e371bb052p-8, -0x1.6b557354b8d85p-6, -0x1.4cc6c058f371dp-2, 0x1.3f67e104a689cp+1, -0x1.7f694c35d9d22p+2, 0x1.8c237a0fa26d7p-1, 0x1.bc14b21b62b93p+4, -0x1.3b920d0f7e0d1p+6, 0x1.1f2d7abf22a73p+7, -0x1.caa9291eb011fp+7, 0x1.4bb5b15c50a64p+8, -0x1.958a02a102449p+8, 0x1.aa319d6af6ec9p+8, -0x1.a2f3c3b8e6e4bp+8, 0x1.7cdc33f38c2d1p+8, -0x1.16422099d0e7p+8, 0x1.18d5fc503a177p+7, -0x1.50ba2688bc043p+5, 0x1.6688e31cf50ep+2, -0x0p+0};
            }
            if constexpr (order == 7) {
                return {0x0p+0, -0x1.869c11b6398bcp-29, 0x1.57bc52155abeap-17, -0x1.c53d49d197175p-15, -0x1.34372727b2155p-5, 0x1.39156fe5f0a3bp-4, 0x1.6fa583815aba2p+0, -0x1.21acdd98ca7e4p+3, 0x1.784654bbd7368p+4, -0x1.013c24badc627p+5, 0x1.93c8d50945095p+4, -0x1.cc8ab1005087cp+3, -0x1.09d1a0b7c5832p+3, 0x1.59e9b6e0b9c21p+6, -0x1.9699942af2221p+7, 0x1.1263590aebe89p+8, -0x1.2e9626a21624ep+8, 0x1.675a39773b248p+8, -0x1.8f56f07493efp+8, 0x1.41512ba753d05p+8, -0x1.48b0d88030936p+7, 0x1.7dd2228ef0449p+5, -0x1.80f665163dff6p+2, 0x0p+0};
            }
            if constexpr (order == 8) {
                return {0x0p+0, -0x1.6c7e600f14b1p-26, 0x1.4002914419fbbp-15, 0x1.a1f2dc12ec55dp-12, -0x1.0f174cea8c306p-4, 0x1.e57ecc0ee730ep-4, 0x1.4e6aaf12e910ep+0, -0x1.bab6f3f4dbb55p+2, 0x1.14ef6d11452bcp+4, -0x1.0a8b2e2c06173p+5, 0x1.e8c645d503448p+5, -0x1.93d454c0120c5p+6, 0x1.20e874beb6698p+7, -0x1.84734ae811783p+7, 0x1.ef2b60fc7ce31p+7, -0x1.156985f5c8bcep+8, 0x1.1c64b27b4c108p+8, -0x1.29c6cec686025p+8, 0x1.23910ebcbe60dp+8, -0x1.b2204de17a53ap+7, 0x1.a873c238f1ad9p+6, -0x1.ddfb9c5686a4fp+4, 0x1.d573db31c98a6p+1, 0x0p+0};
            }
            if constexpr (order == 9) {
                return {0x0p+0, 0x1.d792baac952f3p-25, -0x1.9c18e59aa15abp-15, -0x1.1fc083a0e2c36p-10, 0x1.26b339b06c46fp-5, -0x1.63e903ad632e8p-4, -0x1.04013ac7f90edp-2, 0x1.59e4e79f9d385p+0, -0x1.780b7892927c5p+1, 0x1.ab15772777a66p+2, -0x1.eced01c912f6cp+3, 0x1.9d3e5e5e09252p+4, -0x1.217c835a844bfp+5, 0x1.c4dfc682dd95ep+5, -0x1.4d02e405c664p+6, 0x1.73a4db524ebc3p+6, -0x1.759cb7054906p+6, 0x1.c308c46a51aa1p+6, -0x1.04d3977e22f61p+7, 0x1.a7b259e6037ecp+6, -0x1.ab061d20276fp+5, 0x1.e196df480fe72p+3, -0x1.d2f4e044749b5p+0, -0x0p+0};
            }
            if constexpr (order == 10) {
                return {0x0p+0, -0x1.cf434a44d45d4p-23, 0x1.910e0462410b6p-14, 0x1.17b65cad3d541p-9, -0x1.b7f6865bae123p-6, 0x1.2e7cdbbc90031p-4, -0x1.4a01e6ee0860ep-6, -0x1.07ed621dcab47p-3, 0x1.22ff2ac4d8b9cp-3, -0x1.75b64c97e0968p-1, 0x1.64fab317448fbp+1, -0x1.21aae7ae7df32p+2, 0x1.68255c99ba4c7p+2, -0x1.5bb8e2c464278p+3, 0x1.2898316e66a7fp+4, -0x1.3946d065739bfp+4, 0x1.239464bf9e82cp+4, -0x1.a397a05a5f775p+4, 0x1.1efd1b3c60e0bp+5, -0x1.f4e8dc0eb78aap+4, 0x1.0297dba66cdd3p+4, -0x1.249f95db3c6a8p+2, 0x1.19eec26520f72p-1, -0x0p+0};
            }
            if constexpr (order == 11) {
                return {0x0p+0, -0x1.1b79fa688c947p-21, 0x1.e19176d7ba691p-14, 0x1.e72d4f4950c8ap-10, -0x1.8f4f47250e3a5p-7, 0x1.e6ef9d5209f88p-6, -0x1.cec01fe9b05f6p-5, 0x1.07bd180d95969p-3, -0x1.00dee0f33922ap-2, 0x1.8936ae8cb907cp-2, -0x1.329c321d7ff6bp-1, 0x1.daf35f59b2442p-1, -0x1.3cec213a53552p+0, 0x1.9cd0db1b42b97p+0, -0x1.0e1263d4463c8p+1, 0x1.31bb74af8ae4ep+1, -0x1.4e87023b22be2p+1, 0x1.a91739cfcb437p+1, -0x1.eb2b7ee2bc88dp+1, 0x1.8aa497719f071p+1, -0x1.89e13c889a74dp+0, 0x1.b8315371d42a5p-2, -0x1.a6c51bc471a31p-5, -0x0p+0};
            }
         }
        if constexpr (std::is_same_v<Real, long double>) {
            if constexpr (order == 0) {
                return {0x0p+0L, -0xf.b3f02783f80f12bp-82L, 0xd.d9ac67a56fa6c6bp-63L, 0x8.173a603f30404d9p-50L, -0xb.29199d993b04f39p-44L, 0x8.706909c430c3afbp-36L, -0x9.d2ecaf8c9477601p-32L, -0xc.9a506f26d1499c5p-27L, -0xf.82aa5c72b896837p-25L, 0x9.6b3c60ba5ca6825p-21L, -0xe.d7d175b63fae5edp-21L, -0xa.1dbe70fd272d8b1p-18L, 0xf.807e007e7137774p-15L, -0xc.3c2da01e437b7bbp-12L, 0xc.36482df7945344bp-10L, -0x8.42aa5e8c5592b7bp-8L, 0x8.2154351b02c8eb6p-7L, -0xc.ea7ee1a57276e59p-7L, 0xa.42051ad1f283b58p-6L, -0xa.e37b3292056ef1ap-5L, 0xe.d0a909a97a43f8ep-4L, 0x9.d7b9429d8c4819ap-5L, 0xd.d2b1b0ca6049ce2p-11L, 0x0p+0L};
            }
            if constexpr (order == 1) {
                return {0x0p+0L, -0x8.9fe8bb6860fa129p-72L, 0xf.36c72bc2d649c0dp-55L, -0xb.8e2162d8f32359fp-47L, -0xd.83b74b331022454p-37L, 0xe.43f08c2e52ad9ddp-32L, 0xc.f7487d669e89999p-28L, -0xd.967fbca1c231ae8p-23L, 0x8.9e06b21d13ee1f3p-20L, 0xb.50b3821adfeea95p-18L, -0xe.22db483ae6db17p-15L, 0xe.2e06827e43a1813p-14L, 0x9.eb895556bbbac58p-12L, -0x8.12d440115209d8bp-9L, 0xf.d23fd1fa5810e5dp-9L, -0xc.02cac1033064418p-13L, -0x9.2015e71ca980e98p-6L, 0xd.c777808864fd26fp-5L, -0xc.b921df3df8422acp-4L, 0xf.a8decfd0e796cf1p-4L, 0x8.b5f74cc3838cb01p-5L, -0xb.8b7ab5d60d88817p-4L, -0x9.ed50a30e867a4p-8L, 0x0p+0L};
            }
            if constexpr (order == 2) {
                return {0x0p+0L, -0xd.3e9188aef2c4eadp-68L, 0xb.ae64c04c9a9495p-51L, -0xc.26bf1a6d8ae9884p-41L, -0xb.19136a8314dfe24p-34L, 0x8.7bc285ca6495b9cp-29L, 0xd.1c437d36a462b4ap-24L, -0xe.d36d356af3b5c2cp-22L, 0x8.4de8f96a71c3322p-18L, 0xc.87f5844a5a1a6fp-18L, -0xc.7f0a89b14808e6ap-13L, 0xc.9af2325c9d52281p-11L, -0xa.ce23afc2ebf9f61p-10L, 0x8.48d8eeb161f2f82p-9L, -0xa.32a76a330391e66p-7L, 0xa.863f38d0845d7acp-5L, -0xd.5faefbfebd5fdcdp-4L, 0xc.06770b9c860beebp-3L, -0x9.207f1f9244307c8p-2L, 0xd.e0a9d7a86125cacp-2L, -0xc.f4e3e731c1ae738p-2L, 0xf.0e883c07d830881p-4L, 0xb.bff69700b49b018p-6L, -0x0p+0L};
            }
            if constexpr (order == 3) {
                return {0x0p+0L, 0x8.28d0fc0c8ab5de3p-57L, -0xe.644c2665303f8d2p-42L, 0xb.f467844c7999b52p-35L, 0xc.e487218333646a4p-26L, -0x9.464907251fb5d88p-21L, -0xf.c3729a0af95006dp-18L, 0xe.fe444b0fb58abf3p-14L, -0x9.0357fe081dd0544p-12L, -0xe.ac5fa7b63b752dbp-10L, 0xd.165dcc4f668505cp-7L, -0xf.fef0652fe0d0a8dp-6L, 0xf.70a4d53435d97d5p-7L, 0xe.38176549d19e5a8p-4L, -0xa.f33ca6a935126f2p-2L, 0xf.9d5cc308ea3d3cbp-2L, -0x9.b10d58750d82957p-2L, -0xe.0ff56273860c0ebp-3L, 0xc.0706f222d90cfdfp-1L, -0xc.698603f28aba1a9p-1L, 0xf.811d8b837786d61p-3L, 0x8.e24f13cbf2b6563p-3L, -0xb.0043d843e977a5cp-4L, 0x0p+0L};
            }
            if constexpr (order == 4) {
                return {0x0p+0L, -0x8.e8c32b971e98576p-48L, 0xf.b64634e38dcc13ap-34L, -0xe.acea5d00252b3f7p-28L, -0xe.1b0edd888415d06p-19L, 0xd.b568701f1037966p-15L, 0xb.46a31adde7341a9p-11L, -0x8.019c285744039edp-7L, 0x8.be2451db785c08ap-6L, 0x9.9afeb3ef794f7dcp-4L, -0x8.5c9a5288f0446cfp-1L, 0xa.951a6892d8df2a8p+0L, -0xc.1158d4b716c8e65p+0L, -0xc.879c83b3e72381bp-2L, 0x8.a1670e0f8370cfp+2L, -0xf.c6cf785144862cap+2L, 0x8.d3653eb2b19664fp+3L, -0xe.b96cda04f958dddp+2L, 0xa.d25b7d128f57c87p+2L, -0x8.0a6bd3392a250e2p+2L, 0xa.cb089c6adf4abc4p+1L, -0x9.a6f9ed63fe69fc4p+0L, 0xf.7d00ff140e483d9p-3L, -0x0p+0L};
            }
            if constexpr (order == 5) {
                return {0x0p+0L, 0xd.5466c737c4b8e71p-40L, -0xb.c0207a8d4ab3f5dp-26L, 0x9.f1efc596d13a96ep-21L, 0xa.89be0a399cde158p-12L, -0xc.a8e7f90ef2563f2p-9L, -0xa.1eb811909d4016fp-5L, 0xb.727803e59fb7b2ap-2L, -0xc.55ce61c5a4a452bp-1L, -0xc.aab0e16bbf53e59p+0L, 0xc.6efca52385c3629p+3L, -0x8.2afa403541a338fp+5L, 0xc.e5669d771d3f14bp+5L, -0xe.5a6f3bbda79b16cp+5L, 0xd.138a6eeae115b02p+5L, -0xb.14040db7eec0855p+5L, 0xf.df9cbc35f53d89cp+4L, -0xb.3bdb265f99372aep+3L, -0x8.4e297c1c697e50ep+3L, 0xf.22108273299ad98p+3L, -0xa.4fe028fe4cd412cp+3L, 0xe.82f41ff40a38e2ap+1L, -0x8.b397bb62451d843p-1L, -0x0p+0L};
            }
            if constexpr (order == 6) {
                return {0x0p+0L, 0xc.5c469fd2a8dfc22p-37L, -0xa.e3c88f0029de0d1p-24L, 0x9.306ca994f547077p-20L, 0x9.c6b71b8dd828c64p-11L, -0xb.5aab9aa5c6c2553p-9L, -0xa.663602c79b8e8ecp-5L, 0x9.fb3f0825344de52p-2L, -0xb.fb4a61aece9122bp-1L, 0xc.611bd07d136ba02p-4L, 0xd.e0a590db15c99a5p+1L, -0x9.dc90687bf0684d6p+3L, 0x8.f96bd5f915394b1p+4L, -0xe.554948f5808f966p+4L, 0xa.5dad8ae285323a8p+5L, -0xc.ac50150812244abp+5L, 0xd.518ceb57b764b3fp+5L, -0xd.179e1dc73725588p+5L, 0xb.e6e19f9c6168b44p+5L, -0x8.b21104ce8738061p+5L, 0x8.c6afe281d0bb969p+4L, -0xa.85d13445e021855p+2L, 0xb.344718e7a8702b2p-1L, -0x0p+0L};
            }
            if constexpr (order == 7) {
                return {0x0p+0L, -0xc.34e08db1cc5dff7p-32L, 0xa.bde290aad5f5239p-20L, -0xe.29ea4e8cb8baa35p-18L, -0x9.a1b9393d90aa708p-8L, 0x9.c8ab7f2f851da2dp-7L, 0xb.7d2c1c0ad5d0d11p-3L, -0x9.0d66ecc653f1d47p+0L, 0xb.c232a5deb9b3ce5p+1L, -0x8.09e125d6e31381bp+2L, 0xc.9e46a84a284ab95p+1L, -0xe.64558802843df89p+0L, -0x8.4e8d05be2c19216p+0L, 0xa.cf4db705ce10894p+3L, -0xc.b4cca1579110a15p+4L, 0x8.931ac8575f4469p+5L, -0x9.74b13510b126f8ep+5L, 0xb.3ad1cbb9d9242edp+5L, -0xc.7ab783a49f782e2p+5L, 0xa.0a895d3a9e82ae4p+5L, -0xa.4586c401849b3a9p+4L, 0xb.ee9114778224972p+2L, -0xc.07b328b1effb03p-1L, 0x0p+0L};
            }
            if constexpr (order == 8) {
                return {0x0p+0L, -0xb.63f30078a588239p-29L, 0xa.00148a20cfdd93dp-18L, 0xd.0f96e09762ae4b7p-15L, -0x8.78ba67546182ecp-7L, 0xf.2bf66077398705p-7L, 0xa.735578974886db3p-3L, -0xd.d5b79fa6ddaa942p-1L, 0x8.a77b688a295dd83p+1L, -0x8.5459716030b981p+2L, 0xf.46322ea81a241fbp+2L, -0xc.9ea2a6009062863p+3L, 0x9.0743a5f5b34bf5fp+4L, -0xc.239a57408bc15f7p+4L, 0xf.795b07e3e718ad3p+4L, -0x8.ab4c2fae45e736fp+5L, 0x8.e32593da6083e54p+5L, -0x9.4e367634301286p+5L, 0x9.1c8875e5f3069d6p+5L, -0xd.91026f0bd29d199p+4L, 0xd.439e11c78d6c4d6p+3L, -0xe.efdce2b435277e9p+1L, 0xe.ab9ed98e4c52c95p-2L, 0x0p+0L};
            }
            if constexpr (order == 9) {
                return {0x0p+0L, 0xe.bc95d564a979727p-28L, -0xc.e0c72cd50ad59f6p-18L, -0x8.fe041d07161ad18p-13L, 0x9.3599cd8362377f9p-8L, -0xb.1f481d6b1973f98p-7L, -0x8.2009d63fc87650fp-5L, 0xa.cf273cfce9c2447p-3L, -0xb.c05bc49493e294p-2L, 0xd.58abb93bbd32d48p-1L, -0xf.67680e4897b60bep+0L, 0xc.e9f2f2f049291ep+1L, -0x9.0be41ad4225f98dp+2L, 0xe.26fe3416ecaefa6p+2L, -0xa.6817202e331ffc5p+3L, 0xb.9d26da9275e177fp+3L, -0xb.ace5b82a4830181p+3L, 0xe.184623528d5047ap+3L, -0x8.269cbbf117b07a3p+4L, 0xd.3d92cf301bf5fbep+3L, -0xd.5830e9013b783a8p+2L, 0xf.0cb6fa407f38df5p+0L, -0xe.97a70223a4da73dp-3L, -0x0p+0L};
            }
            if constexpr (order == 10) {
                return {0x0p+0L, -0xe.7a1a5226a2ea342p-26L, 0xc.88702312085adcbp-17L, 0x8.bdb2e569eaa0703p-12L, -0xd.bfb432dd7091673p-9L, 0x9.73e6dde48018975p-7L, -0xa.500f37704306c4bp-9L, -0x8.3f6b10ee55a3b6fp-6L, 0x9.17f95626c5ce28p-6L, -0xb.adb264bf04b4353p-4L, 0xb.27d598ba247d824p-2L, -0x9.0d573d73ef98e15p-1L, 0xb.412ae4cdd263584p-1L, -0xa.ddc71623213c14ap+0L, 0x9.44c18b73353f954p+1L, -0x9.ca36832b9cdf557p+1L, 0x9.1ca325fcf41627ep+1L, -0xd.1cbd02d2fbba7a7p+1L, 0x8.f7e8d9e307059b7p+2L, -0xf.a746e075bc553f8p+1L, 0x8.14bedd3366e9a26p+1L, -0x9.24fcaed9e3542f7p-1L, 0x8.cf76132907b923fp-4L, -0x0p+0L};
            }
            if constexpr (order == 11) {
                return {0x0p+0L, -0x8.dbcfd34464a37e2p-24L, 0xf.0c8bb6bdd348a37p-17L, 0xf.396a7a4a8644d53p-13L, -0xc.7a7a392871d2799p-10L, 0xf.377cea904fc3cf9p-9L, -0xe.7600ff4d82fb121p-8L, 0x8.3de8c06cacb4875p-6L, -0x8.06f70799c914d0ap-5L, 0xc.49b57465c83ddabp-5L, -0x9.94e190ebffb5bb9p-4L, 0xe.d79afacd9220eb9p-4L, -0x9.e76109d29aa8d6cp-3L, 0xc.e686d8da15cb65bp-3L, -0x8.70931ea231e404ep-2L, 0x9.8ddba57c5726fd3p-2L, -0xa.743811d915f0f83p-2L, 0xd.48b9ce7e5a1ba57p-2L, -0xf.595bf715e446977p-2L, 0xc.5524bb8cf8387b7p-2L, -0xc.4f09e444d3a6404p-3L, 0xd.c18a9b8ea152a17p-5L, -0xd.3628de238d184c6p-8L, -0x0p+0L};
            }
         }
#ifdef BOOST_HAS_FLOAT128
        if constexpr (std::is_same_v<Real, boost::multiprecision::float128>) {
            if constexpr (order == 0) {
                return {0x0p+0Q, -0x1.f67e04f07f01e25622e408c1fb5p-79Q, 0x1.bb358cf4adf4d8d516218fdb089fp-60Q, 0x1.02e74c07e60809b1c0066c0f0a8fp-47Q, -0x1.652333b327609e7169e3030b8ca1p-41Q, 0x1.0e0d2138861875f5d48c4e43d1dep-33Q, -0x1.3a5d95f1928eec0223d4ea86a954p-29Q, -0x1.934a0de4da29338aed7c38a1e1adp-24Q, -0x1.f0554b8e5712d06e0273c8db3f93p-22Q, 0x1.2d678c174b94d04aa232fc26ac09p-18Q, -0x1.dafa2eb6c7f5cbd9451ac8baf096p-18Q, -0x1.43b7ce1fa4e5b162d9e05ba18873p-15Q, 0x1.f00fc00fce26eee7a53c75fa8f29p-12Q, -0x1.8785b403c86f6f7548c23ebef8aap-9Q, 0x1.86c905bef28a68962be82a3bfb15p-7Q, -0x1.08554bd18ab256f6f2efbe750dfbp-5Q, 0x1.042a86a360591d6b2ba9327c6219p-4Q, -0x1.9d4fdc34ae4edcb2a06a461fd217p-4Q, 0x1.4840a35a3e5076b07b9bb98536c4p-3Q, -0x1.5c6f665240adde3380bdc0aea431p-2Q, 0x1.da1521352f487f1b71f5ba20803fp-1Q, 0x1.3af72853b1890334601f92a1284fp-2Q, 0x1.ba5636194c0939c427c9ec59581bp-8Q, 0x0p+0Q};
            }
            if constexpr (order == 1) {
                return {0x0p+0Q, -0x1.13fd176d0c1f425174c3b389d5ddp-69Q, 0x1.e6d8e5785ac93819db98735b4ceap-52Q, -0x1.71c42c5b1e646b3ea863aa0f1621p-44Q, -0x1.b076e966620448a8ccb970a07d5dp-34Q, 0x1.c87e1185ca55b3b93aa752e208f9p-29Q, 0x1.9ee90facd3d1333133d7e0787a21p-25Q, -0x1.b2cff794384635cf6c5bd857c344p-20Q, 0x1.13c0d643a27dc3e64f7aa5e99e04p-17Q, 0x1.6a1670435bfdd529f00d82c5a1c9p-15Q, -0x1.c45b69075cdb62e0f0501d72dc7dp-12Q, 0x1.c5c0d04fc87430265835a2d85971p-11Q, 0x1.3d712aaad77758afb71e5644c866p-9Q, -0x1.025a88022a413b155a9680077e65p-6Q, 0x1.fa47fa3f4b021cbace01025ddd1fp-6Q, -0x1.80595820660c8830537a9f5a3b9ap-10Q, -0x1.2402bce395301d2ffcfa2c7a4aa6p-3Q, 0x1.b8eef0110c9fa4de8e71d384473fp-2Q, -0x1.97243be7bf08455883f653e15b78p-1Q, 0x1.f51bd9fa1cf2d9e12262aabd168ep-1Q, 0x1.16bee9987071960277d4321df477p-2Q, -0x1.716f56bac1b1102dab376b1a427ap-1Q, -0x1.3daa1461d0cf47fff25f3df9666ep-5Q, 0x0p+0Q};
            }
            if constexpr (order == 2) {
                return {0x0p+0Q, -0x1.a7d23115de589d5af93186dc833dp-65Q, 0x1.75cc98099352929f72e59e34a35cp-48Q, -0x1.84d7e34db15d3108fe59165be447p-38Q, -0x1.63226d50629bfc487dd3f781635cp-31Q, 0x1.0f7850b94c92b737f20b8713819dp-26Q, 0x1.a3886fa6d48c5693a0b3eda3aa0ep-21Q, -0x1.da6da6ad5e76b8585347759d3bd2p-19Q, 0x1.09bd1f2d4e386644ae469ebc0daep-15Q, 0x1.90feb0894b434ddffdbc817d8b3bp-15Q, -0x1.8fe1513629011cd322fa618302fep-10Q, 0x1.935e464b93aa450281cf24f56eap-8Q, -0x1.59c475f85d7f3ec16b121b445e3dp-7Q, 0x1.091b1dd62c3e5f033d13bd912146p-6Q, -0x1.4654ed4660723ccc0f0acc157ee6p-4Q, 0x1.50c7e71a108baf578aa06cfbe437p-2Q, -0x1.abf5df7fd7abfb9aaab015c9ef57p-1Q, 0x1.80cee17390c17dd6afdc8b5de2ecp+0Q, -0x1.240fe3f248860f90966b383c4f86p+1Q, 0x1.bc153af50c24b9575bdc1b7e1e11p+1Q, -0x1.9e9c7ce63835ce70bd9f19490cc6p+1Q, 0x1.e1d10780fb06110288f8999d1c77p-1Q, 0x1.77fed2e01693602fb8497758cda5p-3Q, -0x0p+0Q};
            }
            if constexpr (order == 3) {
                return {0x0p+0Q, 0x1.051a1f819156bbc6f8ce464e2915p-54Q, -0x1.cc8984cca607f1a4cee533122be8p-39Q, 0x1.7e8cf0898f3336a3b6d3dcb4c532p-32Q, 0x1.9c90e430666c8d48daf02122d3bbp-23Q, -0x1.28c920e4a3f6bb0fc001a358c098p-18Q, -0x1.f86e53415f2a00da542f825467c2p-15Q, 0x1.dfc88961f6b157e54f77890a7f1p-11Q, -0x1.206affc103ba0a875625f4573ee9p-9Q, -0x1.d58bf4f6c76ea5b68a555417b03p-7Q, 0x1.a2cbb989ecd0a0b85c48b0ed9418p-4Q, -0x1.ffde0ca5fc1a1519c08ba08d56d2p-3Q, 0x1.ee149aa686bb2faa29ff6d4a357p-4Q, 0x1.c702eca93a33cb508d7b1cfeb33ap-1Q, -0x1.5e6794d526a24de36780aeb8d3a2p+1Q, 0x1.f3ab98611d47a796d9c54634845cp+1Q, -0x1.3621ab0ea1b052aecf48cff637a4p+1Q, -0x1.c1feac4e70c181d6b4fca6d0f0b8p+0Q, 0x1.80e0de445b219fbdea6ef7328fa3p+2Q, -0x1.8d30c07e51574352c1f1eb1eeca7p+2Q, 0x1.f023b1706ef0dac15d9852fd90fp+0Q, 0x1.1c49e2797e56cac67a206b13f47bp+0Q, -0x1.60087b087d2ef4b82953952c7fp-1Q, 0x0p+0Q};
            }
            if constexpr (order == 4) {
                return {0x0p+0Q, -0x1.1d186572e3d30aec4cf8700909f3p-45Q, 0x1.f6c8c69c71b982736f6617158cc1p-31Q, -0x1.d59d4ba004a567edb2272f7e4cb8p-25Q, -0x1.c361dbb11082ba0b86af830cbb87p-16Q, 0x1.b6ad0e03e206f2cc34c3b046448fp-12Q, 0x1.68d4635bbce683527c347aea7e82p-8Q, -0x1.0033850ae88073da0780e58be837p-4Q, 0x1.17c48a3b6f0b811400fc06ace46fp-3Q, 0x1.335fd67def29efb717b668f0220bp-1Q, -0x1.0b934a511e088d9de4577ff32c14p+2Q, 0x1.52a34d125b1be5508a31d64c3b1cp+3Q, -0x1.822b1a96e2d91cc911472ba87624p+3Q, -0x1.90f390767ce4703536c74c7086f9p+1Q, 0x1.142ce1c1f06e19e0c6f07bf53461p+5Q, -0x1.f8d9ef0a2890c5946bb990c0de7p+5Q, 0x1.1a6ca7d65632cc9e3b88e6689b3fp+6Q, -0x1.d72d9b409f2b1bb9491109cbb561p+5Q, 0x1.5a4b6fa251eaf90e3115b3506783p+5Q, -0x1.014d7a672544a1c4d79027708179p+5Q, 0x1.5961138d5be957877d8d0010a78bp+4Q, -0x1.34df3dac7fcd3f885ed5b3ae7a9bp+3Q, 0x1.efa01fe281c907b167ef0c507d2ap+0Q, -0x0p+0Q};
            }
            if constexpr (order == 5) {
                return {0x0p+0Q, 0x1.aa8cd8e6f8971ce2f7b981d30b62p-37Q, -0x1.78040f51a9567ebaa0ba4fc0fa47p-23Q, 0x1.3e3df8b2da2752dcb7f78a9b946ap-18Q, 0x1.5137c147339bc2af2a68f7c86934p-9Q, -0x1.951cff21de4ac7e4a1f08b551a66p-6Q, -0x1.43d7023213a802dd6433490a4d23p-2Q, 0x1.6e4f007cb3f6f6544bd7b03df217p+1Q, -0x1.8ab9cc38b4948a55522c145d0e8bp+2Q, -0x1.95561c2d77ea7cb1b418d8368cf7p+3Q, 0x1.8ddf94a470b86c5104b7334fb255p+6Q, -0x1.055f4806a834671e407050cd8491p+8Q, 0x1.9cacd3aee3a7e296cd7175d70d77p+8Q, -0x1.cb4de777b4f362d79e6ee730f44ap+8Q, 0x1.a2714ddd5c22b603c83ec9aa52d9p+8Q, -0x1.628081b6fdd810aa965a41523b09p+8Q, 0x1.fbf39786bea7b138ef97a86545fep+7Q, -0x1.677b64cbf326e55bcfc73720d51ep+6Q, -0x1.09c52f838d2fca1c275b6f360e68p+6Q, 0x1.e442104e65335b30f1932113632fp+6Q, -0x1.49fc051fc99a825862e5302e2f71p+6Q, 0x1.d05e83fe81471c548c375ec8ddcp+4Q, -0x1.1672f76c48a3b085a419c9982d9ep+2Q, -0x0p+0Q};
            }
            if constexpr (order == 6) {
                return {0x0p+0Q, 0x1.8b88d3fa551bf843cecd06395909p-34Q, -0x1.5c7911e0053bc1a15002ce79882dp-21Q, 0x1.260d95329ea8e0edce5a9f2ccdfbp-17Q, 0x1.38d6e371bb0518c7efe1c5194abdp-8Q, -0x1.6b557354b8d84aa675e942653d29p-6Q, -0x1.4cc6c058f371d1d7aaddd852cdf4p-2Q, 0x1.3f67e104a689bca4d8d8ba6c6d4ap+1Q, -0x1.7f694c35d9d2245677a9a51e547fp+2Q, 0x1.8c237a0fa26d7404976a5e332c0cp-1Q, 0x1.bc14b21b62b9334a0ec2194dee45p+4Q, -0x1.3b920d0f7e0d09abc1f3de638e46p+6Q, 0x1.1f2d7abf22a7296277ebf2580b71p+7Q, -0x1.caa9291eb011f2ccec54627c809p+7Q, 0x1.4bb5b15c50a647500af623803bf8p+8Q, -0x1.958a02a1024489567ee1c35d3d85p+8Q, 0x1.aa319d6af6ec967e38203aa3c33ap+8Q, -0x1.a2f3c3b8e6e4ab101398ccf1445p+8Q, 0x1.7cdc33f38c2d168782ad59775066p+8Q, -0x1.16422099d0e700c2cc8539851c78p+8Q, 0x1.18d5fc503a1772d2730cc681d1c6p+7Q, -0x1.50ba2688bc0430a9ee763156d70dp+5Q, 0x1.6688e31cf50e0564efec61e86e7fp+2Q, -0x0p+0Q};
            }
            if constexpr (order == 7) {
                return {0x0p+0Q, -0x1.869c11b6398bbfee9546a674e921p-29Q, 0x1.57bc52155abea4723bad8d6a408ap-17Q, -0x1.c53d49d19717546a332769322e23p-15Q, -0x1.34372727b2154e0f632d1590e6b7p-5Q, 0x1.39156fe5f0a3b45977ed8f707f76p-4Q, 0x1.6fa583815aba1a21d8ae91d95586p+0Q, -0x1.21acdd98ca7e3a8d979c5c932075p+3Q, 0x1.784654bbd73679c9f2d3052160d9p+4Q, -0x1.013c24badc62703530b4c2651bcfp+5Q, 0x1.93c8d5094509572acaec62e1e73ap+4Q, -0x1.cc8ab1005087bf119aa77a824edfp+3Q, -0x1.09d1a0b7c583242c5b04835d8f7ap+3Q, 0x1.59e9b6e0b9c21128db1702c4ca3fp+6Q, -0x1.9699942af2221429cba030693cf2p+7Q, 0x1.1263590aebe88d2006545352dbb2p+8Q, -0x1.2e9626a21624df1cdd8a3c6d9cb6p+8Q, 0x1.675a39773b2485da34f2a16d0b2fp+8Q, -0x1.8f56f07493ef05c3005ca7dde3acp+8Q, 0x1.41512ba753d055c75dc2be857e5dp+8Q, -0x1.48b0d88030936751bd00c84baed2p+7Q, 0x1.7dd2228ef04492e49ff5751ea432p+5Q, -0x1.80f665163dff606030fe03a98c35p+2Q, 0x0p+0Q};
            }
            if constexpr (order == 8) {
                return {0x0p+0Q, -0x1.6c7e600f14b10472b96cb740169ep-26Q, 0x1.4002914419fbb27ad678566766cbp-15Q, 0x1.a1f2dc12ec55c96d1325d8263eb5p-12Q, -0x1.0f174cea8c305d809bedec1060bbp-4Q, 0x1.e57ecc0ee730e09fd5abee0609cbp-4Q, 0x1.4e6aaf12e910db6611d8c77b6b67p+0Q, -0x1.bab6f3f4dbb55284f9a33da636a1p+2Q, 0x1.14ef6d11452bbb057169c0f225bp+4Q, -0x1.0a8b2e2c0617302078263e8b7b2bp+5Q, 0x1.e8c645d5034483f58d4a136da4dap+5Q, -0x1.93d454c0120c50c5e608611815d4p+6Q, 0x1.20e874beb6697ebef635bf1876d8p+7Q, -0x1.84734ae811782beee820481c96c9p+7Q, 0x1.ef2b60fc7ce315a6b8fcd6df6e53p+7Q, -0x1.156985f5c8bce6de73d0a173ad53p+8Q, 0x1.1c64b27b4c107ca71af0226ee84fp+8Q, -0x1.29c6cec6860250c0db50af163931p+8Q, 0x1.23910ebcbe60d3ac86a5c1c722fep+8Q, -0x1.b2204de17a53a3315386882bc37ap+7Q, 0x1.a873c238f1ad89ace06cfeb44c4dp+6Q, -0x1.ddfb9c5686a4efd2bd440ee03326p+4Q, 0x1.d573db31c98a592a2b8f7200fb6dp+1Q, 0x0p+0Q};
            }
            if constexpr (order == 9) {
                return {0x0p+0Q, 0x1.d792baac952f2e4df175e444ed1bp-25Q, -0x1.9c18e59aa15ab3ece42c2ee02612p-15Q, -0x1.1fc083a0e2c35a300d32795bf7f4p-10Q, 0x1.26b339b06c46eff131469b32afddp-5Q, -0x1.63e903ad632e7f30304780d7a81dp-4Q, -0x1.04013ac7f90eca1ecef43813a306p-2Q, 0x1.59e4e79f9d38488dc3af22ae7f61p+0Q, -0x1.780b7892927c5280dd461570147ap+1Q, 0x1.ab15772777a65a90ae1212c4eb48p+2Q, -0x1.eced01c912f6c17b1d2febdf7daep+3Q, 0x1.9d3e5e5e092523c01039d51bf71cp+4Q, -0x1.217c835a844bf319c8ff9d033a14p+5Q, 0x1.c4dfc682dd95df4b0a339a6f4e16p+5Q, -0x1.4d02e405c663ff8abf4a298a03edp+6Q, 0x1.73a4db524ebc2efdc4497c9e5037p+6Q, -0x1.759cb70549060301bbd5e18d7e45p+6Q, 0x1.c308c46a51aa08f3170dca57783bp+6Q, -0x1.04d3977e22f60f45b4ffdead7b28p+7Q, 0x1.a7b259e6037ebf7c8e0620dadfe4p+6Q, -0x1.ab061d20276f0750e281c91f0738p+5Q, 0x1.e196df480fe71bea5c91d4412cd9p+3Q, -0x1.d2f4e044749b4e7a1c48f9a5af7fp+0Q, -0x0p+0Q};
            }
            if constexpr (order == 10) {
                return {0x0p+0Q, -0x1.cf434a44d45d4684e001bf1f187p-23Q, 0x1.910e0462410b5b96ed7cfb7220f5p-14Q, 0x1.17b65cad3d540e05672585ab3067p-9Q, -0x1.b7f6865bae122ce58cfae562d0abp-6Q, 0x1.2e7cdbbc900312ead8652bff9ed1p-4Q, -0x1.4a01e6ee0860d896820cbbadf825p-6Q, -0x1.07ed621dcab476dea8f3423f5995p-3Q, 0x1.22ff2ac4d8b9c4ff137637f2c7d3p-3Q, -0x1.75b64c97e09686a5bb55cb243bd9p-1Q, 0x1.64fab317448fb047431d8e24a54p+1Q, -0x1.21aae7ae7df31c2992916da8c563p+2Q, 0x1.68255c99ba4c6b07d82aff3b0359p+2Q, -0x1.5bb8e2c4642782930893d0243ce9p+3Q, 0x1.2898316e66a7f2a70e272d3f11d2p+4Q, -0x1.3946d065739beaad7a0e2a40c507p+4Q, 0x1.239464bf9e82c4fc4fe84e572e24p+4Q, -0x1.a397a05a5f774f4ee219c8b14924p+4Q, 0x1.1efd1b3c60e0b36de50c2dcc8445p+5Q, -0x1.f4e8dc0eb78aa7eff213fd0b0496p+4Q, 0x1.0297dba66cdd344c60220a006936p+4Q, -0x1.249f95db3c6a85ee060f9f62e477p+2Q, 0x1.19eec26520f7247ef7d53934fea1p-1Q, -0x0p+0Q};
            }
            if constexpr (order == 11) {
                return {0x0p+0Q, -0x1.1b79fa688c946fc499418c780f3ap-21Q, 0x1.e19176d7ba69146eb054cee03b15p-14Q, 0x1.e72d4f4950c89aa59db59baf9ae5p-10Q, -0x1.8f4f47250e3a4f324685a323722dp-7Q, 0x1.e6ef9d5209f879f2bf8dfcb794e9p-6Q, -0x1.cec01fe9b05f624125cc39646ab6p-5Q, 0x1.07bd180d959690eaf4247b0a7651p-3Q, -0x1.00dee0f339229a14be03bd5af3b1p-2Q, 0x1.8936ae8cb907bb56ce7b22352259p-2Q, -0x1.329c321d7ff6b77198ead765845ep-1Q, 0x1.daf35f59b2441d717c1c5b34451bp-1Q, -0x1.3cec213a53551ad8400a57dd9116p+0Q, 0x1.9cd0db1b42b96cb69281b4281d17p+0Q, -0x1.0e1263d4463c809bd8425f14f04fp+1Q, 0x1.31bb74af8ae4dfa54f80bb87203ap+1Q, -0x1.4e87023b22be1f061886be12482dp+1Q, 0x1.a91739cfcb4374adc788307acde2p+1Q, -0x1.eb2b7ee2bc88d2eef9228d867b3bp+1Q, 0x1.8aa497719f070f6dc75bd776f2dcp+1Q, -0x1.89e13c889a74c808a9542df44197p+0Q, 0x1.b8315371d42a542d22a870736592p-2Q, -0x1.a6c51bc471a3098c7c2f6b4902fap-5Q, -0x0p+0Q};
            }
         }
#endif
     }
    if constexpr (p == 13) {
        if constexpr (std::is_same_v<Real, float>) {
            if constexpr (order == 0) {
                return {0x0p+0f, 0x1.aa25eep-86f, 0x1.1342e6p-65f, -0x1.d83da4p-53f, 0x1.52fae6p-45f, 0x1.2c5958p-38f, -0x1.07969cp-32f, 0x1.ddb5ep-29f, 0x1.74f33ap-24f, -0x1.2dc446p-21f, -0x1.5c34d2p-20f, 0x1.5f413ap-16f, -0x1.a08b66p-15f, -0x1.d4cc28p-13f, 0x1.0006aap-9f, -0x1.beac12p-8f, 0x1.c908bep-7f, -0x1.014518p-6f, -0x1.505c64p-12f, 0x1.2fca44p-5f, -0x1.f688bep-5f, -0x1.20e16ep-4f, 0x1.d0088ap-1f, 0x1.8ae6a8p-3f, 0x1.749078p-9f, 0x0p+0f};
            }
            if constexpr (order == 1) {
                return {0x0p+0f, 0x1.258ba4p-79f, 0x1.7b38aep-60f, -0x1.8616f6p-48f, 0x1.ce70acp-41f, 0x1.12414ap-33f, -0x1.5fce04p-29f, -0x1.6454c6p-26f, 0x1.bc40d8p-22f, -0x1.bebecp-18f, 0x1.17497ap-16f, 0x1.bbf1c4p-14f, -0x1.328522p-11f, 0x1.3482dcp-11f, 0x1.db22aap-10f, -0x1.288446p-11f, -0x1.1bc8f6p-5f, 0x1.3bc252p-3f, -0x1.8e439cp-2f, 0x1.72e38ep-1f, -0x1.22d30ap+0f, 0x1.856d7ep+0f, -0x1.3373bep-2f, -0x1.0c0768p-1f, -0x1.239236p-6f, 0x0p+0f};
            }
            if constexpr (order == 2) {
                return {0x0p+0f, -0x1.1de108p-72f, -0x1.715282p-54f, 0x1.401eb2p-42f, -0x1.b00a66p-36f, -0x1.287dbp-29f, 0x1.4c2776p-24f, -0x1.b40566p-23f, -0x1.8453f2p-17f, 0x1.f99c3cp-15f, 0x1.50ec24p-14f, -0x1.77d87p-10f, 0x1.1cbed4p-8f, -0x1.632f9ep-11f, -0x1.30a608p-5f, 0x1.30d718p-3f, -0x1.5831aep-2f, 0x1.017de8p-1f, -0x1.f8fd9cp-2f, 0x1.2f821p-2f, -0x1.58e6fap-2f, 0x1.a1c192p+0f, -0x1.368292p+1f, 0x1.e81e74p-1f, 0x1.7e9e9ap-4f, -0x0p+0f};
            }
            if constexpr (order == 3) {
                return {0x0p+0f, -0x1.15123cp-63f, -0x1.65f402p-46f, 0x1.417132p-35f, -0x1.a12416p-29f, -0x1.7e5cfep-23f, 0x1.379718p-18f, 0x1.dc01f8p-17f, -0x1.ec5d7p-12f, 0x1.2653ap-9f, 0x1.fddc3cp-14f, -0x1.35c4a6p-5f, 0x1.38a2dcp-3f, -0x1.05f57cp-2f, -0x1.3531bp-8f, 0x1.c857e6p-1f, -0x1.c6e272p+0f, 0x1.2378e8p+0f, 0x1.ef81d8p+0f, -0x1.9c4d14p+2f, 0x1.433598p+3f, -0x1.36cc46p+3f, 0x1.29379ap+2f, -0x1.e7908ep-3f, -0x1.95e034p-2f, -0x0p+0f};
            }
            if constexpr (order == 4) {
                return {0x0p+0f, 0x1.89a8b8p-54f, 0x1.fc9a2p-38f, -0x1.4bb43p-27f, 0x1.19c574p-21f, 0x1.fa353cp-16f, -0x1.23460cp-11f, -0x1.73814cp-10f, 0x1.5158e6p-5f, -0x1.660ee6p-3f, 0x1.ecd4c4p-5f, 0x1.d5a02cp+0f, -0x1.cda456p+2f, 0x1.a6237cp+3f, -0x1.35b5dcp+3f, -0x1.6d86aep+3f, 0x1.503102p+5f, -0x1.e69256p+5f, 0x1.b8b8bp+5f, -0x1.0ce3p+5f, 0x1.f006ccp+3f, -0x1.39e0f8p+3f, 0x1.24b0fp+3f, -0x1.5d2572p+2f, 0x1.4c5d52p+0f, 0x0p+0f};
            }
            if constexpr (order == 5) {
                return {0x0p+0f, -0x1.9990a4p-46f, -0x1.089ab2p-30f, 0x1.ea6ed6p-21f, -0x1.10aa3ep-15f, -0x1.dea5b4p-10f, 0x1.75f52cp-6f, 0x1.bffc3ap-5f, -0x1.3e3f5cp+0f, 0x1.36efc4p+2f, -0x1.1c0e0ep+2f, -0x1.85c7d2p+4f, 0x1.a875bap+6f, -0x1.b6219ap+7f, 0x1.2591e8p+8f, -0x1.200434p+8f, 0x1.cdbd7ep+7f, -0x1.3d735cp+7f, 0x1.14faaep+6f, 0x1.03810ep+5f, -0x1.98403ap+6f, 0x1.a68f66p+6f, -0x1.fa0a08p+5f, 0x1.5a1bc4p+4f, -0x1.a35158p+1f, -0x0p+0f};
            }
            if constexpr (order == 6) {
                return {0x0p+0f, -0x1.1b57acp-41f, -0x1.6e30fap-27f, 0x1.deb156p-18f, -0x1.51ea7ep-13f, -0x1.29374p-7f, 0x1.275514p-4f, 0x1.41600cp-3f, -0x1.7aa974p+1f, 0x1.66791ep+3f, -0x1.12ed1ep+4f, -0x1.08a1fep+2f, 0x1.1c029p+6f, -0x1.575e0cp+7f, 0x1.25dfd4p+8f, -0x1.bb4938p+8f, 0x1.2782eap+9f, -0x1.4e6b46p+9f, 0x1.48ee3ep+9f, -0x1.298ee6p+9f, 0x1.e67e2cp+8f, -0x1.416f5cp+8f, 0x1.2c245ep+7f, -0x1.54873ep+5f, 0x1.5d159cp+2f, 0x0p+0f};
            }
            if constexpr (order == 7) {
                return {0x0p+0f, 0x1.1fd836p-35f, 0x1.742b04p-22f, -0x1.599f88p-13f, 0x1.1d5d34p-9f, 0x1.f82518p-4f, -0x1.56d14ep-1f, -0x1.34c266p+0f, 0x1.30b382p+4f, -0x1.1609d2p+6f, 0x1.15bd06p+7f, -0x1.737018p+7f, 0x1.9e035ep+7f, -0x1.bce286p+7f, 0x1.48c0f8p+7f, 0x1.6a6f08p+5f, -0x1.2c2d3cp+8f, 0x1.bb4eb6p+8f, -0x1.0686p+9f, 0x1.3d934ap+9f, -0x1.495088p+9f, 0x1.e3622ep+8f, -0x1.c6da62p+7f, 0x1.ecb0f4p+5f, -0x1.d4fe28p+2f, -0x0p+0f};
            }
            if constexpr (order == 8) {
                return {0x0p+0f, 0x1.5b5402p-33f, 0x1.c17582p-21f, -0x1.7cd67ep-12f, 0x1.2b8db4p-10f, 0x1.1dedd2p-3f, -0x1.54f128p-1f, -0x1.1e3582p-2f, 0x1.214136p+3f, -0x1.f14982p+4f, 0x1.0ee55p+6f, -0x1.01b04cp+7f, 0x1.c4899p+7f, -0x1.59b2fp+8f, 0x1.d856p+8f, -0x1.32f292p+9f, 0x1.6f4bc2p+9f, -0x1.823eap+9f, 0x1.7af822p+9f, -0x1.7231bp+9f, 0x1.452582p+9f, -0x1.b10bfp+8f, 0x1.809c44p+7f, -0x1.8fef12p+5f, 0x1.6f91b2p+2f, -0x0p+0f};
            }
            if constexpr (order == 9) {
                return {0x0p+0f, -0x1.38bcd2p-31f, -0x1.95627cp-20f, 0x1.979caep-11f, 0x1.8ddc62p-9f, -0x1.267684p-3f, 0x1.4f8a82p-1f, -0x1.6165b8p-1f, -0x1.f43ec2p+0f, 0x1.d98372p+2f, -0x1.0b7b2ap+4f, 0x1.441338p+5f, -0x1.49a9a8p+6f, 0x1.f90c72p+6f, -0x1.658b66p+7f, 0x1.09589ap+8f, -0x1.592f3ap+8f, 0x1.685dfep+8f, -0x1.68e4a6p+8f, 0x1.907a86p+8f, -0x1.8e74c6p+8f, 0x1.1b0164p+8f, -0x1.0083b2p+7f, 0x1.09aad6p+5f, -0x1.e00f2p+1f, 0x0p+0f};
            }
            if constexpr (order == 10) {
                return {0x0p+0f, 0x1.4225e6p-27f, 0x1.a2ffd6p-17f, -0x1.25fbdep-8f, -0x1.8b4256p-6f, 0x1.99c706p-2f, -0x1.949fp+0f, 0x1.807b44p+1f, -0x1.ffc6e6p+1f, 0x1.b8e78ep+2f, -0x1.6174ccp+3f, 0x1.b5b93ep+2f, 0x1.325da8p+2f, -0x1.1183a6p+3f, 0x1.23fe36p+4f, -0x1.e50ce8p+5f, 0x1.92fe5p+6f, -0x1.7ea49ap+6f, 0x1.73ecdap+6f, -0x1.15e62ap+7f, 0x1.5b315p+7f, -0x1.1076bcp+7f, 0x1.ffd30ap+5f, -0x1.0aa6a6p+4f, 0x1.ddc992p+0f, -0x0p+0f};
            }
            if constexpr (order == 11) {
                return {0x0p+0f, 0x1.db60ap-27f, 0x1.373cf6p-17f, -0x1.2584f2p-10f, -0x1.0ee4f4p-8f, 0x1.904f34p-5f, -0x1.317678p-3f, 0x1.571ac2p-2f, -0x1.a5bdc8p-1f, 0x1.b27cf8p+0f, -0x1.65a332p+1f, 0x1.34062ap+2f, -0x1.05a0dap+3f, 0x1.64406ep+3f, -0x1.c6519cp+3f, 0x1.48587cp+4f, -0x1.a5a57ep+4f, 0x1.a382eep+4f, -0x1.af95c6p+4f, 0x1.154d2p+5f, -0x1.336156p+5f, 0x1.c92e2ap+4f, -0x1.a3afap+3f, 0x1.b108acp+1f, -0x1.82556p-2f, -0x0p+0f};
            }
            if constexpr (order == 12) {
                return {0x0p+0f, -0x1.7fa406p-27f, -0x1.fd199ap-19f, 0x1.4af176p-13f, 0x1.fdc736p-13f, -0x1.b2091cp-9f, 0x1.cd0926p-8f, -0x1.094f52p-6f, 0x1.d4c5bap-5f, -0x1.d5a152p-4f, 0x1.499f7p-3f, -0x1.577aa6p-2f, 0x1.543b96p-1f, -0x1.a643a4p-1f, 0x1.feadc8p-1f, -0x1.c40e7cp+0f, 0x1.3bab08p+1f, -0x1.1ce974p+1f, 0x1.1fd454p+1f, -0x1.cab7dp+1f, 0x1.240bf6p+2f, -0x1.cb9852p+1f, 0x1.aeb7a6p+0f, -0x1.bebdc6p-2f, 0x1.8dbb04p-5f, 0x0p+0f};
            }
         }
        if constexpr (std::is_same_v<Real, double>) {
            if constexpr (order == 0) {
                return {0x0p+0, 0x1.aa25ed9a2917ap-86, 0x1.1342e56840ec1p-65, -0x1.d83da475a639cp-53, 0x1.52fae6e83092p-45, 0x1.2c5958b74836cp-38, -0x1.07969b3a158bdp-32, 0x1.ddb5df0a412fdp-29, 0x1.74f33a2ccbc65p-24, -0x1.2dc4454eebde3p-21, -0x1.5c34d16af592p-20, 0x1.5f413ab3a7eb9p-16, -0x1.a08b65bbb5986p-15, -0x1.d4cc286b09e13p-13, 0x1.0006a916adee3p-9, -0x1.beac12be1e7dep-8, 0x1.c908becee5a6ap-7, -0x1.0145189baf014p-6, -0x1.505c64431d7a7p-12, 0x1.2fca4301354a5p-5, -0x1.f688bd1280c7p-5, -0x1.20e16d5539e3fp-4, 0x1.d008898983914p-1, 0x1.8ae6a79b4939p-3, 0x1.749078c0867aap-9, 0x0p+0};
            }
            if constexpr (order == 1) {
                return {0x0p+0, 0x1.258ba4722413cp-79, 0x1.7b38ad7c07a3ap-60, -0x1.8616f674ef406p-48, 0x1.ce70ab06ee7ebp-41, 0x1.1241493325c92p-33, -0x1.5fce040bf0fa2p-29, -0x1.6454c6cebdf98p-26, 0x1.bc40d833ef35ap-22, -0x1.bebec0cbd1f99p-18, 0x1.17497a4dcd45p-16, 0x1.bbf1c4306d1a3p-14, -0x1.3285223dae55ap-11, 0x1.3482dca87f01cp-11, 0x1.db22a9804aa11p-10, -0x1.2884461caa1e3p-11, -0x1.1bc8f6e410596p-5, 0x1.3bc2524066e84p-3, -0x1.8e439b194a1c1p-2, 0x1.72e38df054a24p-1, -0x1.22d30ac28dd9dp+0, 0x1.856d7e6b6468cp+0, -0x1.3373be839b47cp-2, -0x1.0c0767795d12bp-1, -0x1.23923550f457dp-6, 0x0p+0};
            }
            if constexpr (order == 2) {
                return {0x0p+0, -0x1.1de10823dfb11p-72, -0x1.71528144e6955p-54, 0x1.401eb1fad689ep-42, -0x1.b00a66e2af604p-36, -0x1.287daf168a21ap-29, 0x1.4c2776a5f116ep-24, -0x1.b4056601f4cp-23, -0x1.8453f2cbc2d99p-17, 0x1.f99c3b6dd2e63p-15, 0x1.50ec238da0096p-14, -0x1.77d86fa7164cp-10, 0x1.1cbed37ba4dbep-8, -0x1.632f9e9eb99fcp-11, -0x1.30a6075105b3bp-5, 0x1.30d718554c4cbp-3, -0x1.5831ad1464211p-2, 0x1.017de77081fc3p-1, -0x1.f8fd9be6cc911p-2, 0x1.2f82100ba4f28p-2, -0x1.58e6f9045a40ep-2, 0x1.a1c192a939b96p+0, -0x1.368292c7be10ep+1, 0x1.e81e73bf7d4e6p-1, 0x1.7e9e99d93265ep-4, -0x0p+0};
            }
            if constexpr (order == 3) {
                return {0x0p+0, -0x1.15123cc4efe5bp-63, -0x1.65f402b77f1abp-46, 0x1.4171311c3a7b2p-35, -0x1.a124166c57861p-29, -0x1.7e5cfdc9a1b18p-23, 0x1.379717b64ca38p-18, 0x1.dc01f8633f4c6p-17, -0x1.ec5d6f7ab50c1p-12, 0x1.2653a0321997dp-9, 0x1.fddc3c9f813d7p-14, -0x1.35c4a5041d9b9p-5, 0x1.38a2db882124ep-3, -0x1.05f57c1a6464dp-2, -0x1.3531afa42b77fp-8, 0x1.c857e67243154p-1, -0x1.c6e27151761e4p+0, 0x1.2378e73c6b10ap+0, 0x1.ef81d85b23a02p+0, -0x1.9c4d14b541486p+2, 0x1.4335983e13577p+3, -0x1.36cc467b08bf8p+3, 0x1.293799fcf6c8ep+2, -0x1.e7908dd5753c9p-3, -0x1.95e033565cd46p-2, -0x0p+0};
            }
            if constexpr (order == 4) {
                return {0x0p+0, 0x1.89a8b7b9d28dep-54, 0x1.fc9a1f4f09859p-38, -0x1.4bb4308ed52p-27, 0x1.19c5738584247p-21, 0x1.fa353cdebaaabp-16, -0x1.23460c48f637p-11, -0x1.73814bfd57b68p-10, 0x1.5158e51b27e5ep-5, -0x1.660ee59df078bp-3, 0x1.ecd4c4adc076p-5, 0x1.d5a02bb93708ap+0, -0x1.cda455c3ae393p+2, 0x1.a6237c764b4a3p+3, -0x1.35b5dc1025f2dp+3, -0x1.6d86add093f07p+3, 0x1.503102590c0b8p+5, -0x1.e6925570e3f51p+5, 0x1.b8b8afcf323dp+5, -0x1.0ce2ffa464611p+5, 0x1.f006cc927be53p+3, -0x1.39e0f8d13e8dbp+3, 0x1.24b0f0352262ap+3, -0x1.5d2571cda5544p+2, 0x1.4c5d516dd5341p+0, 0x0p+0};
            }
            if constexpr (order == 5) {
                return {0x0p+0, -0x1.9990a3735f2b8p-46, -0x1.089ab261bfe1dp-30, 0x1.ea6ed6c321e75p-21, -0x1.10aa3eff5f8abp-15, -0x1.dea5b49d2b29ap-10, 0x1.75f52c4bce394p-6, 0x1.bffc3922a3129p-5, -0x1.3e3f5c2c93c34p+0, 0x1.36efc345fd87ep+2, -0x1.1c0e0d6b929f9p+2, -0x1.85c7d166f26b7p+4, 0x1.a875ba068d845p+6, -0x1.b621996a63c18p+7, 0x1.2591e882d85cdp+8, -0x1.20043470ad457p+8, 0x1.cdbd7de801408p+7, -0x1.3d735b94c1054p+7, 0x1.14faaec768465p+6, 0x1.03810e60e0205p+5, -0x1.98403a971e631p+6, 0x1.a68f6523a827p+6, -0x1.fa0a0878d9a93p+5, 0x1.5a1bc35a1cc0ep+4, -0x1.a35157d8dba57p+1, -0x0p+0};
            }
            if constexpr (order == 6) {
                return {0x0p+0, -0x1.1b57abb7d7f1cp-41, -0x1.6e30fae8c34eep-27, 0x1.deb156b08b87fp-18, -0x1.51ea7d9f34504p-13, -0x1.29373f56ac109p-7, 0x1.275514f436a32p-4, 0x1.41600caa02956p-3, -0x1.7aa9733036516p+1, 0x1.66791e904ad3bp+3, -0x1.12ed1d6ca773p+4, -0x1.08a1fe0ac681ap+2, 0x1.1c028f2725fcep+6, -0x1.575e0c206b7eep+7, 0x1.25dfd4b9fa45bp+8, -0x1.bb493883e544dp+8, 0x1.2782e95d7bf9fp+9, -0x1.4e6b4521eb7a5p+9, 0x1.48ee3d90a4e35p+9, -0x1.298ee57db6ed6p+9, 0x1.e67e2b7f93304p+8, -0x1.416f5cc6c22c2p+8, 0x1.2c245e185ac8bp+7, -0x1.54873e7e1b3d2p+5, 0x1.5d159c98a090ap+2, 0x0p+0};
            }
            if constexpr (order == 7) {
                return {0x0p+0, 0x1.1fd835e2f3ac3p-35, 0x1.742b04eba227p-22, -0x1.599f88826e1dfp-13, 0x1.1d5d333a03734p-9, 0x1.f82518456ddfap-4, -0x1.56d14dbd1c0fp-1, -0x1.34c266d274d02p+0, 0x1.30b382548236p+4, -0x1.1609d2e86de1fp+6, 0x1.15bd065f2789cp+7, -0x1.737017c34f7dfp+7, 0x1.9e035ea535c3bp+7, -0x1.bce28601020e5p+7, 0x1.48c0f7b3bed9ep+7, 0x1.6a6f07f5d07ddp+5, -0x1.2c2d3c048caafp+8, 0x1.bb4eb52dbc23dp+8, -0x1.068600d74e1c9p+9, 0x1.3d934aacf56d1p+9, -0x1.49508717a90e8p+9, 0x1.e3622dd5ba05ep+8, -0x1.c6da62a20cd0bp+7, 0x1.ecb0f32ca76fbp+5, -0x1.d4fe27a0835aep+2, -0x0p+0};
            }
            if constexpr (order == 8) {
                return {0x0p+0, 0x1.5b5402fbe8db8p-33, 0x1.c1758243d6f69p-21, -0x1.7cd67d6099899p-12, 0x1.2b8db385d6eebp-10, 0x1.1dedd290f908p-3, -0x1.54f12799790c9p-1, -0x1.1e358290df39dp-2, 0x1.21413698e7cf2p+3, -0x1.f149829700499p+4, 0x1.0ee550f2b46fbp+6, -0x1.01b04bb5c7e4fp+7, 0x1.c4898ffe6b604p+7, -0x1.59b2f057a74e5p+8, 0x1.d856007cb8fap+8, -0x1.32f291524b4fdp+9, 0x1.6f4bc124eb7dcp+9, -0x1.823e9f8992b5ap+9, 0x1.7af82162ac35ap+9, -0x1.7231afc08bd27p+9, 0x1.45258188f7147p+9, -0x1.b10bf090c92b5p+8, 0x1.809c444472a8ep+7, -0x1.8fef11e8059cfp+5, 0x1.6f91b1fb16e46p+2, -0x0p+0};
            }
            if constexpr (order == 9) {
                return {0x0p+0, -0x1.38bcd286c35dep-31, -0x1.95627b187fdb2p-20, 0x1.979cae945ae1p-11, 0x1.8ddc624d415c9p-9, -0x1.267683e93d2b7p-3, 0x1.4f8a8171bb09p-1, -0x1.6165b77555fcbp-1, -0x1.f43ec1f2a328p+0, 0x1.d98372b1c68d4p+2, -0x1.0b7b2a297fde1p+4, 0x1.441338caa52dp+5, -0x1.49a9a7009f2b2p+6, 0x1.f90c71880a916p+6, -0x1.658b66c44d2dap+7, 0x1.09589acf597f5p+8, -0x1.592f39b72049ap+8, 0x1.685dfe59eb94fp+8, -0x1.68e4a68143eb6p+8, 0x1.907a858da198ep+8, -0x1.8e74c5b415b4dp+8, 0x1.1b01639f15a47p+8, -0x1.0083b1ed3add2p+7, 0x1.09aad588fc6e2p+5, -0x1.e00f201ca9b48p+1, 0x0p+0};
            }
            if constexpr (order == 10) {
                return {0x0p+0, 0x1.4225e5d67d6a9p-27, 0x1.a2ffd5a58a3e6p-17, -0x1.25fbde4ae0e29p-8, -0x1.8b4255eb2ab53p-6, 0x1.99c706ee33d33p-2, -0x1.949effe4fd5f1p+0, 0x1.807b43d233e13p+1, -0x1.ffc6e54f09b8ap+1, 0x1.b8e78e1a9ada5p+2, -0x1.6174cc292a978p+3, 0x1.b5b93d37eb78fp+2, 0x1.325da71616788p+2, -0x1.1183a5ac74637p+3, 0x1.23fe36ce5250dp+4, -0x1.e50ce86c3eca3p+5, 0x1.92fe5087093b6p+6, -0x1.7ea4991281a8p+6, 0x1.73ecdaa7dbcd1p+6, -0x1.15e629a2f07abp+7, 0x1.5b314fd5bd842p+7, -0x1.1076bbdb4cf75p+7, 0x1.ffd30a30f104ep+5, -0x1.0aa6a6fe474ap+4, 0x1.ddc9914b8fce2p+0, -0x0p+0};
            }
            if constexpr (order == 11) {
                return {0x0p+0, 0x1.db60a068c3378p-27, 0x1.373cf584c03efp-17, -0x1.2584f1f36b43ep-10, -0x1.0ee4f3163f55fp-8, 0x1.904f3491bf69cp-5, -0x1.31767817c6adfp-3, 0x1.571ac2983f85ep-2, -0x1.a5bdc7192b9c4p-1, 0x1.b27cf7ab49be3p+0, -0x1.65a3327bcf1acp+1, 0x1.34062a3bf7c11p+2, -0x1.05a0da2a87dbp+3, 0x1.64406dc5dece2p+3, -0x1.c6519b1cdb293p+3, 0x1.48587bdb4781p+4, -0x1.a5a57db59b67cp+4, 0x1.a382ed435f8edp+4, -0x1.af95c6a7f034fp+4, 0x1.154d20c532d38p+5, -0x1.336156f19a3ebp+5, 0x1.c92e2a69bebeap+4, -0x1.a3af9f5e28534p+3, 0x1.b108ac4bd9d5ep+1, -0x1.82555f0fca682p-2, -0x0p+0};
            }
            if constexpr (order == 12) {
                return {0x0p+0, -0x1.7fa406de421dcp-27, -0x1.fd199a3528256p-19, 0x1.4af1750eeda63p-13, 0x1.fdc735071589p-13, -0x1.b2091c25d2022p-9, 0x1.cd0926f4929adp-8, -0x1.094f52115b37cp-6, 0x1.d4c5bad41cb39p-5, -0x1.d5a1514a8be98p-4, 0x1.499f6f758e381p-3, -0x1.577aa571c0085p-2, 0x1.543b9648dc189p-1, -0x1.a643a47710e6p-1, 0x1.feadc8c000dcep-1, -0x1.c40e7b478efcp+0, 0x1.3bab08b8a360bp+1, -0x1.1ce97447c2138p+1, 0x1.1fd45435f9eb6p+1, -0x1.cab7cf17d2cf9p+1, 0x1.240bf593fb4cep+2, -0x1.cb98514b274a2p+1, 0x1.aeb7a5e38c253p+0, -0x1.bebdc5ccdeaa7p-2, 0x1.8dbb03f4beecap-5, 0x0p+0};
            }
         }
        if constexpr (std::is_same_v<Real, long double>) {
            if constexpr (order == 0) {
                return {0x0p+0L, 0xd.512f6cd148bd157p-89L, 0x8.9a172b420760a64p-68L, -0xe.c1ed23ad31ce183p-56L, 0xa.97d73741848fd34p-48L, 0x9.62cac5ba41b6085p-41L, -0x8.3cb4d9d0ac5e8bp-35L, 0xe.edaef852097ea04p-32L, 0xb.a799d1665e328c8p-27L, -0x9.6e222a775ef19bdp-24L, -0xa.e1a68b57ac901bdp-23L, 0xa.fa09d59d3f5ca6ep-19L, -0xd.045b2dddacc33b9p-18L, -0xe.a66143584f09aacp-16L, 0x8.003548b56f71a16p-12L, -0xd.f56095f0f3ef04ap-11L, 0xe.4845f6772d3537ep-10L, -0x8.0a28c4dd780a179p-9L, -0xa.82e32218ebd378fp-15L, 0x9.7e521809aa52714p-8L, -0xf.b445e89406380cap-8L, -0x9.070b6aa9cf1f929p-7L, 0xe.80444c4c1c89ecep-4L, 0xc.57353cda49c7dedp-6L, 0xb.a483c60433d4d39p-12L, 0x0p+0L};
            }
            if constexpr (order == 1) {
                return {0x0p+0L, 0x9.2c5d2391209dd6cp-82L, 0xb.d9c56be03d1d0d3p-63L, -0xc.30b7b3a77a02f76p-51L, 0xe.7385583773f5acap-44L, 0x8.920a49992e4937dp-36L, -0xa.fe70205f87d115p-32L, -0xb.22a63675efcbfa2p-29L, 0xd.e206c19f79acc2bp-25L, -0xd.f5f6065e8fcc9b9p-21L, 0x8.ba4bd26e6a27feep-19L, 0xd.df8e218368d1bf8p-17L, -0x9.942911ed72acc95p-14L, 0x9.a416e543f80e22fp-14L, 0xe.d9154c025508998p-13L, -0x9.442230e550f1904p-14L, -0x8.de47b72082cb31fp-8L, 0x9.de1292033741f59p-6L, -0xc.721cd8ca50e0b31p-5L, 0xb.971c6f82a511d45p-4L, -0x9.169856146ece4d3p-3L, 0xc.2b6bf35b234616fp-3L, -0x9.9b9df41cda3e206p-5L, -0x8.603b3bcae89572fp-4L, -0x9.1c91aa87a2be53cp-9L, 0x0p+0L};
            }
            if constexpr (order == 2) {
                return {0x0p+0L, -0x8.ef08411efd8849cp-75L, -0xb.8a940a2734aa54p-57L, 0xa.00f58fd6b44ef03p-45L, -0xd.805337157b01dc1p-39L, -0x9.43ed78b4510cd51p-32L, 0xa.613bb52f88b70ecp-27L, -0xd.a02b300fa5ffc5bp-26L, -0xc.229f965e16ccb5cp-20L, 0xf.cce1db6e973199cp-18L, 0xa.87611c6d004b37bp-17L, -0xb.bec37d38b25fe46p-13L, 0x8.e5f69bdd26decb6p-11L, -0xb.197cf4f5ccfe389p-14L, -0x9.85303a882d9d527p-8L, 0x9.86b8c2aa62655e8p-6L, -0xa.c18d68a32108919p-5L, 0x8.0bef3b840fe187dp-4L, -0xf.c7ecdf3664887bbp-5L, 0x9.7c10805d27943dp-5L, -0xa.c737c822d206e95p-5L, 0xd.0e0c9549cdcb3a8p-3L, -0x9.b414963df08728ap-2L, 0xf.40f39dfbea733adp-4L, 0xb.f4f4cec9932f0f6p-7L, -0x0p+0L};
            }
            if constexpr (order == 3) {
                return {0x0p+0L, -0x8.a891e6277f2d5e4p-66L, -0xb.2fa015bbf8d5914p-49L, 0xa.0b8988e1d3d9312p-38L, -0xd.0920b362bc30ba7p-32L, -0xb.f2e7ee4d0d8bee7p-26L, 0x9.bcb8bdb2651be88p-21L, 0xe.e00fc319fa6309dp-20L, -0xf.62eb7bd5a860bbfp-15L, 0x9.329d0190ccbe8fdp-12L, 0xf.eee1e4fc09eb578p-17L, -0x9.ae252820ecdca2bp-8L, 0x9.c516dc410926fe2p-6L, -0x8.2fabe0d32326ba1p-5L, -0x9.a98d7d215bbfa64p-11L, 0xe.42bf339218a9d4bp-4L, -0xe.37138a8bb0f1f2cp-3L, 0x9.1bc739e3588517cp-3L, 0xf.7c0ec2d91d00d64p-3L, -0xc.e268a5aa0a433dap-1L, 0xa.19acc1f09abba63p+0L, -0x9.b66233d845fbec5p+0L, 0x9.49bccfe7b646d1p-1L, -0xf.3c846eaba9e4babp-6L, -0xc.af019ab2e6a32f3p-5L, -0x0p+0L};
            }
            if constexpr (order == 4) {
                return {0x0p+0L, 0xc.4d45bdce946ee8cp-57L, 0xf.e4d0fa784c2c7aap-41L, -0xa.5da18476a8fff99p-30L, 0x8.ce2b9c2c2123732p-24L, 0xf.d1a9e6f5d555505p-19L, -0x9.1a306247b1b7d38p-14L, -0xb.9c0a5feabdb40fdp-13L, 0xa.8ac728d93f2edffp-8L, -0xb.30772cef83c5a19p-6L, 0xf.66a6256e03aff12p-8L, 0xe.ad015dc9b84519p-3L, -0xe.6d22ae1d71c94e9p-1L, 0xd.311be3b25a5141ap+0L, -0x9.adaee0812f969a5p+0L, -0xb.6c356e849f83a5ep+0L, 0xa.818812c8605bdbdp+2L, -0xf.3492ab871fa8811p+2L, 0xd.c5c57e7991e8074p+2L, -0x8.6717fd232308bfp+2L, 0xf.80366493df29948p+0L, -0x9.cf07c689f46d8f8p+0L, 0x9.258781a91314f19p+0L, -0xa.e92b8e6d2aa1c92p-1L, 0xa.62ea8b6ea9a0684p-3L, 0x0p+0L};
            }
            if constexpr (order == 5) {
                return {0x0p+0L, -0xc.cc851b9af95be1ep-49L, -0x8.44d5930dff0e45cp-33L, 0xf.5376b6190f3a9e1p-24L, -0x8.8551f7fafc55605p-18L, -0xe.f52da4e9594cf98p-13L, 0xb.afa9625e71ca263p-9L, 0xd.ffe1c91518944a2p-8L, -0x9.f1fae1649e1a3cap-3L, 0x9.b77e1a2fec3f3f9p-1L, -0x8.e0706b5c94fcaddp-1L, -0xc.2e3e8b37935bbc2p+1L, 0xd.43add0346c224eap+3L, -0xd.b10ccb531e0bf33p+4L, 0x9.2c8f4416c2e68c6p+5L, -0x9.0021a3856a2b57ap+5L, 0xe.6debef400a03eb9p+4L, -0x9.eb9adca6082a1b5p+4L, 0x8.a7d5763b4232a87p+3L, 0x8.1c08730701024ebp+2L, -0xc.c201d4b8f318b15p+3L, 0xd.347b291d41380b3p+3L, -0xf.d05043c6cd499c3p+2L, 0xa.d0de1ad0e606ddbp+1L, -0xd.1a8abec6dd2b59ap-2L, -0x0p+0L};
            }
            if constexpr (order == 6) {
                return {0x0p+0L, -0x8.dabd5dbebf8dd1p-44L, -0xb.7187d7461a76c1bp-30L, 0xe.f58ab5845c3f783p-21L, -0xa.8f53ecf9a28205dp-16L, -0x9.49b9fab56084ac4p-10L, 0x9.3aa8a7a1b518f2bp-7L, 0xa.0b00655014aaec4p-6L, -0xb.d54b9981b28ad7fp-2L, 0xb.33c8f482569da0fp+0L, -0x8.9768eb653b980edp+1L, -0x8.450ff056340cfcbp-1L, 0x8.e01479392fe733bp+3L, -0xa.baf061035bf7341p+4L, 0x9.2efea5cfd22d7f4p+5L, -0xd.da49c41f2a26a23p+5L, 0x9.3c174aebdfcfafep+6L, -0xa.735a290f5bd2b61p+6L, 0xa.4771ec85271a548p+6L, -0x9.4c772bedb76afeap+6L, 0xf.33f15bfc9982224p+5L, -0xa.0b7ae6361160fep+5L, 0x9.6122f0c2d6458a6p+4L, -0xa.a439f3f0d9e9083p+2L, 0xa.e8ace4c50485103p-1L, 0x0p+0L};
            }
            if constexpr (order == 7) {
                return {0x0p+0L, 0x8.fec1af179d61b9fp-38L, 0xb.a158275d1138117p-25L, -0xa.ccfc441370efba9p-16L, 0x8.eae999d01b9a034p-12L, 0xf.c128c22b6efd051p-7L, -0xa.b68a6de8e077c6ep-4L, -0x9.a6133693a681011p-3L, 0x9.859c12a411afe04p+1L, -0x8.b04e97436f0fa36p+3L, 0x8.ade832f93c4e129p+4L, -0xb.9b80be1a7bef779p+4L, 0xc.f01af529ae1d75fp+4L, -0xd.e71430081072aa3p+4L, 0xa.4607bd9df6ceea7p+4L, 0xb.53783fae83ee946p+2L, -0x9.6169e024655745ap+5L, 0xd.da75a96de11e455p+5L, -0x8.343006ba70e487p+6L, 0x9.ec9a5567ab68514p+6L, -0xa.4a8438bd48741b4p+6L, 0xf.1b116eadd02ec05p+5L, -0xe.36d31510668543p+4L, 0xf.658799653b7da4ap+2L, -0xe.a7f13d041ad70dp-1L, -0x0p+0L};
            }
            if constexpr (order == 8) {
                return {0x0p+0L, 0xa.daa017df46dbf37p-36L, 0xe.0bac121eb7b4789p-24L, -0xb.e6b3eb04cc4c723p-15L, 0x9.5c6d9c2eb77580cp-13L, 0x8.ef6e9487c83fc69p-6L, -0xa.a7893ccbc864ac9p-4L, -0x8.f1ac1486f9ce578p-5L, 0x9.0a09b4c73e78daep+0L, -0xf.8a4c14b8024c883p+1L, 0x8.772a8795a37d803p+3L, -0x8.0d825dae3f27671p+4L, 0xe.244c7ff35b01e7fp+4L, -0xa.cd9782bd3a72557p+5L, 0xe.c2b003e5c7d03a7p+5L, -0x9.97948a925a7e40dp+6L, 0xb.7a5e09275bee17p+6L, -0xc.11f4fc4c95acfe3p+6L, 0xb.d7c10b1561accb3p+6L, -0xb.918d7e045e93b18p+6L, 0xa.292c0c47b8a3474p+6L, -0xd.885f8486495a99bp+5L, 0xc.04e222239546c21p+4L, -0xc.7f788f402ce76dp+2L, 0xb.7c8d8fd8b7231c8p-1L, -0x0p+0L};
            }
            if constexpr (order == 9) {
                return {0x0p+0L, -0x9.c5e694361aef307p-34L, -0xc.ab13d8c3fed8cb7p-23L, 0xc.bce574a2d70819p-14L, 0xc.6ee3126a0ae46c7p-12L, -0x9.33b41f49e95bb8cp-6L, 0xa.7c540b8dd847e5dp-4L, -0xb.0b2dbbaaafe5407p-4L, -0xf.a1f60f951940304p-3L, 0xe.cc1b958e3469fa9p-1L, -0x8.5bd9514bfef0a5ep+1L, 0xa.2099c6552967d07p+2L, -0xa.4d4d3804f95933ep+3L, 0xf.c8638c40548b2c6p+3L, -0xb.2c5b3622696cf85p+4L, 0x8.4ac4d67acbfab3fp+5L, -0xa.c979cdb9024d1b7p+5L, 0xb.42eff2cf5ca79c7p+5L, -0xb.4725340a1f5b14p+5L, 0xc.83d42c6d0cc6daep+5L, -0xc.73a62da0ada64c1p+5L, 0x8.d80b1cf8ad238cbp+5L, -0x8.041d8f69d6e8efep+4L, 0x8.4d56ac47e370e97p+2L, -0xf.007900e54da3df6p-2L, 0x0p+0L};
            }
            if constexpr (order == 10) {
                return {0x0p+0L, 0xa.112f2eb3eb54818p-30L, 0xd.17fead2c51f312dp-20L, -0x9.2fdef25707148aap-11L, -0xc.5a12af5955a96fbp-9L, 0xc.ce3837719e99949p-5L, -0xc.a4f7ff27eaf885bp-3L, 0xc.03da1e919f0994ap-2L, -0xf.fe372a784dc531cp-2L, 0xd.c73c70d4d6d2496p-1L, -0xb.0ba6614954bbf44p+0L, 0xd.adc9e9bf5bc7971p-1L, 0x9.92ed38b0b3c3d8dp-1L, -0x8.8c1d2d63a31b528p+0L, 0x9.1ff1b672928682ep+1L, -0xf.28674361f651857p+2L, 0xc.97f2843849db0a7p+3L, -0xb.f524c8940d401f5p+3L, 0xb.9f66d53ede688c9p+3L, -0x8.af314d1783d59d9p+4L, 0xa.d98a7eadec21394p+4L, -0x8.83b5deda67ba93ep+4L, 0xf.fe98518788271d7p+2L, -0x8.553537f23a4ffcdp+1L, 0xe.ee4c8a5c7e71374p-3L, -0x0p+0L};
            }
            if constexpr (order == 11) {
                return {0x0p+0L, 0xe.db05034619bbe3ap-30L, 0x9.b9e7ac2601f7483p-20L, -0x9.2c278f9b5a1f3dep-13L, -0x8.772798b1faafa2dp-11L, 0xc.8279a48dfb4e393p-8L, -0x9.8bb3c0be356f62ap-6L, 0xa.b8d614c1fc2f1d3p-5L, -0xd.2dee38c95ce1d1p-4L, 0xd.93e7bd5a4df1be7p-3L, -0xb.2d1993de78d60e7p-2L, 0x9.a03151dfbe084ep-1L, -0x8.2d06d1543ed7f3cp+0L, 0xb.22036e2ef671076p+0L, -0xe.328cd8e6d949737p+0L, 0xa.42c3deda3c082c7p+1L, -0xd.2d2bedacdb3df28p+1L, 0xd.1c176a1afc76911p+1L, -0xd.7cae353f81a7b93p+1L, 0x8.aa690629969becdp+2L, -0x9.9b0ab78cd1f558dp+2L, 0xe.4971534df5f5093p+1L, -0xd.1d7cfaf14299ea2p+0L, 0xd.8845625eceaef48p-2L, -0xc.12aaf87e5340e4ep-5L, -0x0p+0L};
            }
            if constexpr (order == 12) {
                return {0x0p+0L, -0xb.fd2036f210ee33bp-30L, -0xf.e8ccd1a9412b0e4p-22L, 0xa.578ba8776d3173p-16L, 0xf.ee39a838ac47eep-16L, -0xd.9048e12e9010fcep-12L, 0xe.684937a494d69a2p-11L, -0x8.4a7a908ad9bddb3p-9L, 0xe.a62dd6a0e59c4b5p-8L, -0xe.ad0a8a545f4bd54p-7L, 0xa.4cfb7bac71c0bebp-6L, -0xa.bbd52b8e00425bfp-5L, 0xa.a1dcb246e0c46c2p-4L, -0xd.321d23b88730149p-4L, 0xf.f56e460006e6e42p-4L, -0xe.2073da3c77e022cp-3L, 0x9.dd5845c51b0587bp-2L, -0x8.e74ba23e109c11fp-2L, 0x8.fea2a1afcf5ac7cp-2L, -0xe.55be78be967ca6fp-2L, 0x9.205fac9fda66e7dp-1L, -0xe.5cc28a593a5127dp-2L, 0xd.75bd2f1c6129876p-3L, -0xd.f5ee2e66f553a64p-5L, 0xc.6dd81fa5f765105p-8L, 0x0p+0L};
            }
         }
#ifdef BOOST_HAS_FLOAT128
        if constexpr (std::is_same_v<Real, boost::multiprecision::float128>) {
            if constexpr (order == 0) {
                return {0x0p+0Q, 0x1.aa25ed9a2917a2ae55644835d21ap-86Q, 0x1.1342e56840ec14c80f3138e592fap-65Q, -0x1.d83da475a639c305d56b6abaa71dp-53Q, 0x1.52fae6e83091fa6866cddac840f4p-45Q, 0x1.2c5958b74836c10a90616f2a8cd3p-38Q, -0x1.07969b3a158bd160a5ce22c80f1p-32Q, 0x1.ddb5df0a412fd407d061c8680508p-29Q, 0x1.74f33a2ccbc6518fcdedc80ea4ddp-24Q, -0x1.2dc4454eebde3379c7645c57b31cp-21Q, -0x1.5c34d16af592037a11107a3c49b6p-20Q, 0x1.5f413ab3a7eb94dcb396b26439eep-16Q, -0x1.a08b65bbb598677114fddb7bc68cp-15Q, -0x1.d4cc286b09e13557f3e87a0f94abp-13Q, 0x1.0006a916adee342c6b2792c6bd9fp-9Q, -0x1.beac12be1e7de0939dd15f0315ffp-8Q, 0x1.c908becee5a6a6fcccc956053537p-7Q, -0x1.0145189baf0142f2a8862b8be423p-6Q, -0x1.505c64431d7a6f1e36954a7c7bbdp-12Q, 0x1.2fca4301354a4e285e82fc4b617dp-5Q, -0x1.f688bd1280c70194203f744266cbp-5Q, -0x1.20e16d5539e3f25296b1c6880b59p-4Q, 0x1.d008898983913d9bf1d194ac98fcp-1Q, 0x1.8ae6a79b4938fbd909b19934360cp-3Q, 0x1.749078c0867a9a724454443e75abp-9Q, 0x0p+0Q};
            }
            if constexpr (order == 1) {
                return {0x0p+0Q, 0x1.258ba4722413bad88d58509b825fp-79Q, 0x1.7b38ad7c07a3a1a5b76a01a7ffe3p-60Q, -0x1.8616f674ef405eeb4c6ec0a4eefp-48Q, 0x1.ce70ab06ee7eb593179514053762p-41Q, 0x1.1241493325c926f92e647132b6b8p-33Q, -0x1.5fce040bf0fa22a0bb7f561e2ddfp-29Q, -0x1.6454c6cebdf97f43adbbe0b2cd33p-26Q, 0x1.bc40d833ef3598561d72cd9d0703p-22Q, -0x1.bebec0cbd1f99371a7867cb4cb6ep-18Q, 0x1.17497a4dcd44ffdc92c6ebce1f08p-16Q, 0x1.bbf1c4306d1a37f02d18135ec0fep-14Q, -0x1.3285223dae55992a1d57ada13afbp-11Q, 0x1.3482dca87f01c45e946a258c3061p-11Q, 0x1.db22a9804aa1133069d05b96c6d8p-10Q, -0x1.2884461caa1e3207b6de1334ba2dp-11Q, -0x1.1bc8f6e41059663daa79bf4a8b19p-5Q, 0x1.3bc2524066e83eb2d5b3ff97e21ep-3Q, -0x1.8e439b194a1c16626ebcbea59c46p-2Q, 0x1.72e38df054a23a8a5ebabc498f9cp-1Q, -0x1.22d30ac28dd9c9a5de5413538c66p+0Q, 0x1.856d7e6b6468c2deff8342e75b74p+0Q, -0x1.3373be839b47c40b1d36ec3e48c4p-2Q, -0x1.0c0767795d12ae5d156ec5cf43dep-1Q, -0x1.23923550f457ca7814b0bb3fd5f3p-6Q, 0x0p+0Q};
            }
            if constexpr (order == 2) {
                return {0x0p+0Q, -0x1.1de10823dfb1093738d7e87643fp-72Q, -0x1.71528144e6954a7ffc62308a7d35p-54Q, 0x1.401eb1fad689de068b59fe468de8p-42Q, -0x1.b00a66e2af603b8210b62dc87672p-36Q, -0x1.287daf168a219aa1132b053c041bp-29Q, 0x1.4c2776a5f116e1d85dd049c72f2ep-24Q, -0x1.b4056601f4bff8b5e7d366f5d9afp-23Q, -0x1.8453f2cbc2d996b7775a12422d92p-17Q, 0x1.f99c3b6dd2e633388ada35d3b19p-15Q, 0x1.50ec238da00966f5e3d28fc7d1cep-14Q, -0x1.77d86fa7164bfc8c126a648166f9p-10Q, 0x1.1cbed37ba4dbd96c225d810586d2p-8Q, -0x1.632f9e9eb99fc711faf42ce3ac9ep-11Q, -0x1.30a6075105b3aa4d5336fc196ddp-5Q, 0x1.30d718554c4cabcfafc1434fa9e5p-3Q, -0x1.5831ad14642112326c714ae44775p-2Q, 0x1.017de77081fc30f9ad09eddb0ddbp-1Q, -0x1.f8fd9be6cc910f763e95aa8add16p-2Q, 0x1.2f82100ba4f2879faa0ab73173c2p-2Q, -0x1.58e6f9045a40dd2a472104dd1b9dp-2Q, 0x1.a1c192a939b9674ffcca9b2f45d8p+0Q, -0x1.368292c7be10e513188a104c7d63p+1Q, 0x1.e81e73bf7d4e6759d9f3db7648bfp-1Q, 0x1.7e9e99d93265e1ecb27f1a830c1fp-4Q, -0x0p+0Q};
            }
            if constexpr (order == 3) {
                return {0x0p+0Q, -0x1.15123cc4efe5abc719a7afe54ee3p-63Q, -0x1.65f402b77f1ab228cf374526611cp-46Q, 0x1.4171311c3a7b2623f568e5da2d31p-35Q, -0x1.a124166c5786174e1b9a6616834p-29Q, -0x1.7e5cfdc9a1b17dce328c7eac4902p-23Q, 0x1.379717b64ca37d0f86a74fce3867p-18Q, 0x1.dc01f8633f4c61396f200828148dp-17Q, -0x1.ec5d6f7ab50c177dc61e0a441837p-12Q, 0x1.2653a0321997d1fa4ad99294cae6p-9Q, 0x1.fddc3c9f813d6aef1f4f236ae26dp-14Q, -0x1.35c4a5041d9b9456d2520e84661cp-5Q, 0x1.38a2db882124dfc4e44179cb5386p-3Q, -0x1.05f57c1a6464d742d9d2e21ba65ep-2Q, -0x1.3531afa42b77f4c754018d53b69bp-8Q, 0x1.c857e67243153a965a648fba693bp-1Q, -0x1.c6e27151761e3e586f540385d17dp+0Q, 0x1.2378e73c6b10a2f7957b04749e9bp+0Q, 0x1.ef81d85b23a01ac83af86d6ad1bap+0Q, -0x1.9c4d14b5414867b3e2c99a375decp+2Q, 0x1.4335983e135774c5e5156a578a0cp+3Q, -0x1.36cc467b08bf7d8a89d5be9fffe2p+3Q, 0x1.293799fcf6c8da1f93fe2d122982p+2Q, -0x1.e7908dd5753c97550b1b5604f2eap-3Q, -0x1.95e033565cd465e59fc64a9eaba5p-2Q, -0x0p+0Q};
            }
            if constexpr (order == 4) {
                return {0x0p+0Q, 0x1.89a8b7b9d28ddd18e40bda1ebe6bp-54Q, 0x1.fc9a1f4f09858f53a1ad68a5af58p-38Q, -0x1.4bb4308ed51fff3108f194cb438bp-27Q, 0x1.19c5738584246e6431d3021e7dc6p-21Q, 0x1.fa353cdebaaaaa09ad0c35d40766p-16Q, -0x1.23460c48f636fa706d007644123dp-11Q, -0x1.73814bfd57b681fa89fa68091ab2p-10Q, 0x1.5158e51b27e5dbfe89f572c851f2p-5Q, -0x1.660ee59df078b432fb4f8324566ap-3Q, 0x1.ecd4c4adc075fe24577533d0672cp-5Q, 0x1.d5a02bb93708a320124ac6040ac8p+0Q, -0x1.cda455c3ae3929d25c3c535dc6a5p+2Q, 0x1.a6237c764b4a28340b84822e3294p+3Q, -0x1.35b5dc1025f2d34a60aa76eb49fep+3Q, -0x1.6d86add093f074bc432a4133b7f7p+3Q, 0x1.503102590c0b7b7ae2742afa26bp+5Q, -0x1.e6925570e3f510211d73253c9cacp+5Q, 0x1.b8b8afcf323d00e8ad39cbf76d21p+5Q, -0x1.0ce2ffa4646117df99ed8f422ee9p+5Q, 0x1.f006cc927be5328f1a33cf39b2e3p+3Q, -0x1.39e0f8d13e8db1f029ce5e78c5b8p+3Q, 0x1.24b0f03522629e31e5a0db947751p+3Q, -0x1.5d2571cda5543923016a33c5946p+2Q, 0x1.4c5d516dd5340d08798328e30c5bp+0Q, 0x0p+0Q};
            }
            if constexpr (order == 5) {
                return {0x0p+0Q, -0x1.9990a3735f2b7c3b4bcc9c5c6d94p-46Q, -0x1.089ab261bfe1c8b7eefd2af48d49p-30Q, 0x1.ea6ed6c321e753c1b3e83afca64ep-21Q, -0x1.10aa3eff5f8aac0902d77dc70b2ap-15Q, -0x1.dea5b49d2b299f2f340578d3117ap-10Q, 0x1.75f52c4bce3944c6d40a3262d92fp-6Q, 0x1.bffc3922a312894449dcabcfef0fp-5Q, -0x1.3e3f5c2c93c347939cabcb13bed9p+0Q, 0x1.36efc345fd87e7f12a7b84f216d3p+2Q, -0x1.1c0e0d6b929f95b9e5274f67c9c6p+2Q, -0x1.85c7d166f26b7783d79262788a0cp+4Q, 0x1.a875ba068d8449d487f680050a58p+6Q, -0x1.b621996a63c17e6687c94aa1d16ap+7Q, 0x1.2591e882d85cd18b0dc9613af4f8p+8Q, -0x1.20043470ad456af43bdd442cbbd3p+8Q, 0x1.cdbd7de801407d716cb9d8eda836p+7Q, -0x1.3d735b94c105436aa6c4b09ac11ap+7Q, 0x1.14faaec76846550eede2e14b17f6p+6Q, 0x1.03810e60e02049d57e6ea41d6f98p+5Q, -0x1.98403a971e63162a02dd30d22c56p+6Q, 0x1.a68f6523a82701663e042214f2c3p+6Q, -0x1.fa0a0878d9a933860f50245b2636p+5Q, 0x1.5a1bc35a1cc0dbb69c7105d46f75p+4Q, -0x1.a35157d8dba56b3359c49ce3ca7bp+1Q, -0x0p+0Q};
            }
            if constexpr (order == 6) {
                return {0x0p+0Q, -0x1.1b57abb7d7f1ba1f0926b76733bap-41Q, -0x1.6e30fae8c34ed8366c5fab4833d7p-27Q, 0x1.deb156b08b87ef0592d2a2f14595p-18Q, -0x1.51ea7d9f345040b9bbf370347fcep-13Q, -0x1.29373f56ac1095884f458a2ef40ep-7Q, 0x1.275514f436a31e564b7a89e74e72p-4Q, 0x1.41600caa02955d88aadc61b0b887p-3Q, -0x1.7aa9733036515afd48104df86f7ep+1Q, 0x1.66791e904ad3b41e47b16e9ebe0cp+3Q, -0x1.12ed1d6ca77301d9871a3004023ep+4Q, -0x1.08a1fe0ac6819f96e8282c250fefp+2Q, 0x1.1c028f2725fce675841d5a8c502ep+6Q, -0x1.575e0c206b7ee681892ff3aaebbbp+7Q, 0x1.25dfd4b9fa45afe8f2d3696cd273p+8Q, -0x1.bb493883e544d446bc28df66b84fp+8Q, 0x1.2782e95d7bf9f5fc360b57a94564p+9Q, -0x1.4e6b4521eb7a56c1ef734d3311afp+9Q, 0x1.48ee3d90a4e34a8f4d16c78e4d5bp+9Q, -0x1.298ee57db6ed5fd4bac8bd1360cep+9Q, 0x1.e67e2b7f93304447556e5c8f09e9p+8Q, -0x1.416f5cc6c22c1fc04eeb4968fc0dp+8Q, 0x1.2c245e185ac8b14c7ba88e92b596p+7Q, -0x1.54873e7e1b3d210625d001ab6dc3p+5Q, 0x1.5d159c98a090a206b1f32e2b3775p+2Q, 0x0p+0Q};
            }
            if constexpr (order == 7) {
                return {0x0p+0Q, 0x1.1fd835e2f3ac373d2ea72dba8e3p-35Q, 0x1.742b04eba227022dc26715ff5c1cp-22Q, -0x1.599f88826e1df7529e1dcd7c643ep-13Q, 0x1.1d5d333a03734067f69739a0ed86p-9Q, 0x1.f82518456ddfa0a1ea758d742e0dp-4Q, -0x1.56d14dbd1c0ef8dc07186c402e2fp-1Q, -0x1.34c266d274d02022671c56db6a3bp+0Q, 0x1.30b382548235fc08ce24d7d0f28fp+4Q, -0x1.1609d2e86de1f46bc028b6a7e95ep+6Q, 0x1.15bd065f2789c251ed4d5e04b6fep+7Q, -0x1.737017c34f7deef272252bf550c3p+7Q, 0x1.9e035ea535c3aebdf11dab668824p+7Q, -0x1.bce28601020e5546efd8a5d7daf3p+7Q, 0x1.48c0f7b3bed9dd4e2bc388185f3ep+7Q, 0x1.6a6f07f5d07dd28cc7321b0a4de5p+5Q, -0x1.2c2d3c048caae8b4e7ec9201a7dbp+8Q, 0x1.bb4eb52dbc23c8aa9dab0e6887ebp+8Q, -0x1.068600d74e1c90df9e5804aa0d8fp+9Q, 0x1.3d934aacf56d0a27027ca96526cep+9Q, -0x1.49508717a90e836846bcdacd328ap+9Q, 0x1.e3622dd5ba05d80a96357f5ac49cp+8Q, -0x1.c6da62a20cd0a85f326f76b51d4cp+7Q, 0x1.ecb0f32ca76fb494d56c96b31432p+5Q, -0x1.d4fe27a0835ae1a04e036356545fp+2Q, -0x0p+0Q};
            }
            if constexpr (order == 8) {
                return {0x0p+0Q, 0x1.5b5402fbe8db7e6eb19df50b45d1p-33Q, 0x1.c1758243d6f68f1116147eb57681p-21Q, -0x1.7cd67d6099898e46f31e243716b5p-12Q, 0x1.2b8db385d6eeb017f63c3d2cc4b4p-10Q, 0x1.1dedd290f907f8d2a830e53754bcp-3Q, -0x1.54f12799790c9591a34a2a400d4cp-1Q, -0x1.1e358290df39caefe5bc0f3a4095p-2Q, 0x1.21413698e7cf1b5b76b7b2b4a296p+3Q, -0x1.f149829700499105aba8989a6386p+4Q, 0x1.0ee550f2b46fb0067da9ad974911p+6Q, -0x1.01b04bb5c7e4ece2060d9f80313fp+7Q, 0x1.c4898ffe6b603cfd7a34a118a2f4p+7Q, -0x1.59b2f057a74e4aae3e99e4fbed0fp+8Q, 0x1.d856007cb8fa074ef05482f74bbcp+8Q, -0x1.32f291524b4fc819caac6200a944p+9Q, 0x1.6f4bc124eb7dc2e04be92f595dfbp+9Q, -0x1.823e9f8992b59fc66ded8fb98587p+9Q, 0x1.7af82162ac359965f673b98bf02bp+9Q, -0x1.7231afc08bd276307906f44234f2p+9Q, 0x1.45258188f71468e863d0e94993f4p+9Q, -0x1.b10bf090c92b53363545ffef6caap+8Q, 0x1.809c444472a8d84159fb3e880188p+7Q, -0x1.8fef11e8059ceda062ab926eba61p+5Q, 0x1.6f91b1fb16e463903ef56c841ef4p+2Q, -0x0p+0Q};
            }
            if constexpr (order == 9) {
                return {0x0p+0Q, -0x1.38bcd286c35de60d3e23dee22bb4p-31Q, -0x1.95627b187fdb196e7b629064b1eap-20Q, 0x1.979cae945ae1031f8b4e7263a405p-11Q, 0x1.8ddc624d415c8d8d987b2b2cbd9ep-9Q, -0x1.267683e93d2b771765d874c691d6p-3Q, 0x1.4f8a8171bb08fcb9cccfc5ebffd6p-1Q, -0x1.6165b77555fca80d9c5132d2189p-1Q, -0x1.f43ec1f2a328060887422c05e444p+0Q, 0x1.d98372b1c68d3f52be7c7c9ca575p+2Q, -0x1.0b7b2a297fde14bc6fabc0b44edp+4Q, 0x1.441338caa52cfa0e2a22db24d964p+5Q, -0x1.49a9a7009f2b267b39b548676c23p+6Q, 0x1.f90c71880a91658c0d91f6217ceap+6Q, -0x1.658b66c44d2d9f0911851306fe0ep+7Q, 0x1.09589acf597f567d376e5c50f028p+8Q, -0x1.592f39b72049a36e791edaea4a51p+8Q, 0x1.685dfe59eb94f38e35ee34e31246p+8Q, -0x1.68e4a68143eb6280f77bd8e601bdp+8Q, 0x1.907a858da198db5b4a7f9538990cp+8Q, -0x1.8e74c5b415b4c982c6228ffd195ap+8Q, 0x1.1b01639f15a47195db7ed19527ebp+8Q, -0x1.0083b1ed3add1dfc1895b602d516p+7Q, 0x1.09aad588fc6e1d2e5e2cdb128c9bp+5Q, -0x1.e00f201ca9b47bec13aa674a41b1p+1Q, 0x0p+0Q};
            }
            if constexpr (order == 10) {
                return {0x0p+0Q, 0x1.4225e5d67d6a903055bbacbfc331p-27Q, 0x1.a2ffd5a58a3e625ad94b3c26a0eap-17Q, -0x1.25fbde4ae0e29153aeb57176734ep-8Q, -0x1.8b4255eb2ab52df5075f139b4415p-6Q, 0x1.99c706ee33d33292c7b12a7534e2p-2Q, -0x1.949effe4fd5f10b6f36fb9a3c178p+0Q, 0x1.807b43d233e13293d03e0a7a9e08p+1Q, -0x1.ffc6e54f09b8a6376eed0359c272p+1Q, 0x1.b8e78e1a9ada492b51ca021308f2p+2Q, -0x1.6174cc292a977e876d4f6d6eb101p+3Q, 0x1.b5b93d37eb78f2e23cecae29ead8p+2Q, 0x1.325da71616787b1a91a66572e298p+2Q, -0x1.1183a5ac74636a4f4b4d339f3ea6p+3Q, 0x1.23fe36ce5250d05c2dd2bbdf6703p+4Q, -0x1.e50ce86c3eca30ad8c519882582ap+5Q, 0x1.92fe5087093b614dbab20c470227p+6Q, -0x1.7ea4991281a803e96e059fe4b67p+6Q, 0x1.73ecdaa7dbcd11927e9b6e5f041fp+6Q, -0x1.15e629a2f07ab3b28692b82bc7a2p+7Q, 0x1.5b314fd5bd8427278fa0d1efc894p+7Q, -0x1.1076bbdb4cf7527cd1169e255c3fp+7Q, 0x1.ffd30a30f104e3ad2c8abcee60e9p+5Q, -0x1.0aa6a6fe4749ff9a5422d9c4a3a9p+4Q, 0x1.ddc9914b8fce26e89bc7a2b0d628p+0Q, -0x0p+0Q};
            }
            if constexpr (order == 11) {
                return {0x0p+0Q, 0x1.db60a068c3377c744a93a4132504p-27Q, 0x1.373cf584c03ee9069f8f0a0190cp-17Q, -0x1.2584f1f36b43e7bbd7b609c5faf9p-10Q, -0x1.0ee4f3163f55f45905956286576cp-8Q, 0x1.904f3491bf69c7251cb44e046dbep-5Q, -0x1.31767817c6adec541d32e5a417f7p-3Q, 0x1.571ac2983f85e3a6dc1ec25443bcp-2Q, -0x1.a5bdc7192b9c3a20ec6d8bf75e52p-1Q, 0x1.b27cf7ab49be37cde074ac977d85p+0Q, -0x1.65a3327bcf1ac1cee0eedbfdbe74p+1Q, 0x1.34062a3bf7c109c055f0a9eaed07p+2Q, -0x1.05a0da2a87dafe77cb3821b24817p+3Q, 0x1.64406dc5dece20ec3bcccff17d1p+3Q, -0x1.c6519b1cdb292e6d288de07c138bp+3Q, 0x1.48587bdb4781058d97af1ab964p+4Q, -0x1.a5a57db59b67be5081a437e09db3p+4Q, 0x1.a382ed435f8ed22166b529a680fdp+4Q, -0x1.af95c6a7f034f726744d6082c503p+4Q, 0x1.154d20c532d37d99d4cd6220e55dp+5Q, -0x1.336156f19a3eab1917a5bc52b43p+5Q, 0x1.c92e2a69bebea12559e2af658dd3p+4Q, -0x1.a3af9f5e28533d44a132993cf76p+3Q, 0x1.b108ac4bd9d5de90a524ded5bfd6p+1Q, -0x1.82555f0fca681c9cec5be73239c9p-2Q, -0x0p+0Q};
            }
            if constexpr (order == 12) {
                return {0x0p+0Q, -0x1.7fa406de421dc676803df861dfeap-27Q, -0x1.fd199a35282561c7d7391492e277p-19Q, 0x1.4af1750eeda62e5f44851b7fd282p-13Q, 0x1.fdc735071588fdc0eb0810493aebp-13Q, -0x1.b2091c25d2021f9b77fa62403ef7p-9Q, 0x1.cd0926f4929ad34462a92ca3ac9ap-8Q, -0x1.094f52115b37bb655549976ff627p-6Q, 0x1.d4c5bad41cb389697ad26ff620ep-5Q, -0x1.d5a1514a8be97aa7ad9ab0a15385p-4Q, 0x1.499f6f758e3817d6d293038ad614p-3Q, -0x1.577aa571c0084b7d1102d7b713a1p-2Q, 0x1.543b9648dc188d831b43dde2b311p-1Q, -0x1.a643a47710e60292eadde5f1619cp-1Q, 0x1.feadc8c000dcdc8411b875b0972ep-1Q, -0x1.c40e7b478efc04571e6999136019p+0Q, 0x1.3bab08b8a360b0f528e7d7f84fd6p+1Q, -0x1.1ce97447c213823d4bcfb3e5d007p+1Q, 0x1.1fd45435f9eb58f769b49658b1b5p+1Q, -0x1.cab7cf17d2cf94dd9ebc0c7d4c22p+1Q, 0x1.240bf593fb4cdcf9b1d63d11a02p+2Q, -0x1.cb98514b274a24fad7f002adfa1cp+1Q, 0x1.aeb7a5e38c2530eb9809ef3f7f35p+0Q, -0x1.bebdc5ccdeaa74c83a6c1b000992p-2Q, 0x1.8dbb03f4beeca20a5c2c766dd4cfp-5Q, 0x0p+0Q};
            }
         }
#endif
     }
}

} // namespaces
#endif
