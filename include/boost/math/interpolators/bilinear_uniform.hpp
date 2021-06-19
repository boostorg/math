// Copyright Nick Thompson, 2021
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

// This implements bilinear interpolation on a uniform grid.
// If dx and dy are both positive, then the (x,y) = (x0, y0) is associated with data index 0 (herein referred to as f[0])
// The point (x0 + dx, y0) is associated with f[1], and (x0 + i*dx, y0) is associated with f[i],
// i.e., we are assuming traditional C row major order.
// The y coordinate increases *downward*, as is traditional in 2D computer graphics.
// This is *not* how people generally think in numerical analysis (although it *is* how they lay out matrices).
// Providing the capability of a grid rotation is too expensive and not ergonomic; you'll need to perform any rotations at the call level.

// For clarity, the value f(x0 + i*dx, y0 + j*dy) must be stored in the f[j*cols + i] position.

#ifndef BOOST_MATH_INTERPOLATORS_BILINEAR_UNIFORM_HPP
#define BOOST_MATH_INTERPOLATORS_BILINEAR_UNIFORM_HPP
#include <stdexcept>

namespace boost::math::interpolators {

template <class RandomAccessContainer>
class bilinear_uniform
{
public:
    using Real = typename RandomAccessContainer::value_type;

    bilinear_uniform(RandomAccessContainer && fieldData, decltype(fieldData.size()) rows, decltype(fieldData.size()) cols, Real dx = 1, Real dy = 1, Real x0 = 0, Real y0 = 0)
    {
        if(fieldData.size() != rows*cols)
        {
            throw std::logic_error("The field data must have rows*cols elements.");
        }
        if (rows < 2) {
            throw std::logic_error("There must be at least two rows of data for bilinear interpolation to be well-defined.");
        }
        if (cols < 2) {
            throw std::logic_error("There must be at least two columns of data for bilinear interpolation to be well-defined.");
        }

        fieldData_ = std::move(fieldData);
        rows_ = rows;
        cols_ = cols;
        x0_ = x0;
        y0_ = y0;
        dx_ = dx;
        dy_ = dy;
        assert(dx_ != 0);
        assert(dy_ != 0);
    }

    Real operator()(Real x, Real y) const
    {
        using std::floor;
        if (x > x0_ + (cols_ - 1)*dx_ || x < x0_) {
            std::cerr << __FILE__ << ":" << __LINE__ << ":" << __func__ << "\n";
            std::cerr << "Querying the bilinear_uniform interpolator at (x,y) = (" << x << ", " << y << ") is not allowed.\n";
            std::cerr << "x must lie in the interval [" << x0_ << ", " << x0_ + (cols_ -1)*dx_ << "]\n";
            return std::numeric_limits<Real>::quiet_NaN();
        }
        if (y > y0_ + (rows_ - 1)*dy_ || y < y0_) {
            std::cerr << __FILE__ << ":" << __LINE__ << ":" << __func__ << "\n";
            std::cerr << "Querying the bilinear_uniform interpolator at (x,y) = (" << x << ", " << y << ") is not allowed.\n";
            std::cerr << "y must lie in the interval [" << y0_ << ", " << y0_ + (rows_ -1)*dy_ << "]\n";
            return std::numeric_limits<Real>::quiet_NaN(); 
        }

        Real s = (x - x0_)/dx_;
        Real s0 = floor(s);
        Real t = (y - y0_)/dy_;
        Real t0 = floor(t);
        decltype(fieldData_.size()) xidx = s0;
        decltype(fieldData_.size()) yidx = t0;
        decltype(fieldData_.size()) idx = yidx*cols_  + xidx;
        Real alpha = s - s0;
        Real beta = t - t0;

        Real fhi;
        // If alpha = 0, then we can segfault by reading fieldData_[idx+1]:
        if (alpha <= 2*s0*std::numeric_limits<Real>::epsilon())  {
            fhi = fieldData_[idx];
        } else {
            fhi = (1 - alpha)*fieldData_[idx] + alpha*fieldData_[idx + 1];
        }

        // Again, we can get OOB access without this check.
        // This corresponds to interpolation over a line segment aligned with the axes.
        if (beta <= 2*t0*std::numeric_limits<Real>::epsilon()) {
            return fhi;
        }

#ifndef NDEBUG
        if (idx + cols_ >= fieldData_.size()) {
            std::cerr << __FILE__ << ":" << __LINE__ << " About to segfault!!!\n";
            std::cerr << "idx + cols_ + 1 = " << idx + cols_ + 1 << ", but fieldData_.size() = " << fieldData_.size() << "\n";
            std::cerr << "alpha = " << alpha << ", beta = " << beta << ", s = " << s << ", t = " << t << "\n";
            std::cerr << *this;
        }
#endif
        auto bottom_left = fieldData_[idx + cols_];
        Real flo;
        if (alpha <= 2*s0*std::numeric_limits<Real>::epsilon())  {
            flo = bottom_left;
        } else {
#ifndef NDEBUG
            if (idx + cols_ + 1 >= fieldData_.size())  {
                std::cerr << __FILE__ << ":" << __LINE__ << ":" << __func__ << " Big trouble!\n";
                std::cerr << "idx + cols_ + 1 = " << idx + cols_ + 1 << ", but fieldData_.size() = " << fieldData_.size() << "\n";
                std::cerr << "alpha = " << alpha << ", beta = " << beta << "\n";
                std::cerr << *this;
            }
#endif
            flo = (1 - alpha)*bottom_left + alpha*fieldData_[idx + cols_ + 1];
        }
        // Convex combination over vertical to get the value:
        return (1 - beta)*fhi + beta*flo;
    }

    

    friend std::ostream& operator<<(std::ostream& out, bilinear_uniform<RandomAccessContainer> const & bu) {
        out << "(x₀, y₀) = (" << bu.x0_ << ", " << bu.y0_ << "), (Δx, Δy) = (" << bu.dx_ << ", " << bu.dy_ << "),";
        out << "(xf, yf) = (" << bu.x0_ + (bu.cols_ -1)*bu.dx_ << ", " << bu.y0_ + (bu.rows_ -1)*bu.dy_ << ")\n";
        for (decltype(bu.rows_) j = 0; j < bu.rows_; ++j) {
            out << "{";
            for (decltype(bu.cols_) i = 0; i < bu.cols_ - 1; ++i) {
                out << bu.fieldData_[j*bu.rows_ + i] << ", ";
            }
            out << bu.fieldData_[j*bu.rows_ + bu.cols_ - 1] << "}\n";
        }
        return out;
    }

private:
    RandomAccessContainer fieldData_;
    decltype(fieldData_.size()) rows_;
    decltype(fieldData_.size()) cols_;
    Real x0_;
    Real y0_;
    Real dx_;
    Real dy_;
};

}
#endif
