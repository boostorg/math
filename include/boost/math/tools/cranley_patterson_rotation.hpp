/*
 *  (C) Copyright Nick Thompson 2018.
 *  Use, modification and distribution are subject to the
 *  Boost Software License, Version 1.0. (See accompanying file
 *  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 *
 *  Applies the Cranley-Patterson rotation to a sequence of quasi-random number generation.
 */
#ifndef BOOST_MATH_TOOLS_CRANLEY_PATTERSON_ROTATION_HPP
#define BOOST_MATH_TOOLS_CRANLEY_PATTERSON_ROTATION_HPP

namespace boost { namespace math {

template<class Real>
class cranley_patterson_rotation
{
public:
    cranley_patterson_rotation(size_t dimension)
    {
        std::random_device rd;
        std::mt19937_64 gen(rd());
        std::uniform_real_distribution<Real> dis(0, 1);
        m_rand_vec.resize(dimension);
        for (size_t i = 0; i < m_rand_vec.size(); ++i)
        {
            m_rand_vec[i] = dis(gen);
        }
    }

    template<class ForwardIterator>
    void operator()(ForwardIterator begin, ForwardIterator end)
    {
        auto it1 = begin;
        auto it2 = m_rand_vec.begin();
        size_t i = 0;
        while(it1 != end)
        {
            *it1 += *it2++;
            if (*it1 > 1)
            {
               *it1 -= 1;
            }
            ++it1;
            ++i;
        }
        assert(i == m_rand_vec.size());
    }

private:
    std::vector<Real> m_rand_vec;
};
}}
#endif
