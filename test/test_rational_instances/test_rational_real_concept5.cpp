
#include <boost/detail/workaround.hpp>
#if !BOOST_WORKAROUND(__BORLANDC__, BOOST_TESTED_AT(0x582))

#include "test_rational.hpp"
#include <boost/math/concepts/real_concept.hpp>

template void do_test_spots<boost::math::concepts::real_concept, float>(boost::math::concepts::real_concept, float);

#endif
