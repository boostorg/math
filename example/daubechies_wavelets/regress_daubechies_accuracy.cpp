#include <iostream>
#include <string>
#include <fstream>
#include <map>
#include <cmath>
#include <vector>
#include <iomanip>
#include <boost/algorithm/string.hpp>
#include <boost/math/statistics/linear_regression.hpp>


int main(int argc, char** argv)
{
    if (argc != 2)
    {
        std::cout << "Usage: ./regress_accuracy.x foo.csv\n";
        return 1;
    }
    std::string filename = std::string(argv[1]);
    std::ifstream ifs(filename.c_str());
    if (!ifs.good())
    {
        std::cerr << "Couldn't find file " << filename << "\n";
        return 1;
    }
    std::map<std::string, std::vector<double>> m;

    std::string header_line;
    std::getline(ifs, header_line);
    std::cout << "Header line = " << header_line << "\n";
    std::vector<std::string> header_strs;
    boost::split(header_strs, header_line, boost::is_any_of(","));
    for (auto & s : header_strs) {
        boost::algorithm::trim(s);
    }

    std::string line;
    std::vector<double> r;
    std::vector<double> matched_holder;
    std::vector<double> linear;
    std::vector<double> quadratic_b_spline;
    std::vector<double> cubic_b_spline;
    std::vector<double> quintic_b_spline;
    std::vector<double> cubic_hermite;
    std::vector<double> pchip;
    std::vector<double> makima;
    std::vector<double> fotaylor;
    std::vector<double> quintic_hermite;
    std::vector<double> sotaylor;
    std::vector<double> totaylor;
    std::vector<double> septic_hermite;
    while(std::getline(ifs, line))
    {
        std::vector<std::string> strs;
        boost::split(strs, line, boost::is_any_of(","));
        for (auto & s : strs)
        {
            boost::algorithm::trim(s);
        }
        std::vector<double> v(strs.size(), std::numeric_limits<double>::quiet_NaN());
        for (size_t i = 0; i < v.size(); ++i)
        {
            v[i] = std::stod(strs[i]);
        }
        r.push_back(v[0]);
        matched_holder.push_back(std::log2(v[1]));
        linear.push_back(std::log2(v[2]));
        quadratic_b_spline.push_back(std::log2(v[3]));
        cubic_b_spline.push_back(std::log2(v[4]));
        quintic_b_spline.push_back(std::log2(v[5]));
        cubic_hermite.push_back(std::log2(v[6]));
        pchip.push_back(std::log2(v[7]));
        makima.push_back(std::log2(v[8]));
        fotaylor.push_back(std::log2(v[9]));
        if (v.size() > 10) {
            quintic_hermite.push_back(std::log2(v[10]));
            sotaylor.push_back(std::log2(v[11]));
        }
        if (v.size() > 12) {
            totaylor.push_back(std::log2(v[12]));
            septic_hermite.push_back(std::log2(v[13]));
        }
    }

    std::cout << std::fixed << std::setprecision(16);
    auto p  = boost::math::statistics::simple_ordinary_least_squares(r, matched_holder);

    assert(p.second < 0);
    std::cout << "Matched Holder    : " << p.first << " - " << std::abs(p.second) << "r\n";

    p  = boost::math::statistics::simple_ordinary_least_squares(r, linear);
    assert(p.second < 0);
    std::cout << "Linear            : " << p.first << " - " << std::abs(p.second) << "r\n";

    p  = boost::math::statistics::simple_ordinary_least_squares(r, quadratic_b_spline);
    assert(p.second < 0);
    std::cout << "Quadratic B-spline: " << p.first << " - " << std::abs(p.second) << "r\n";

    p  = boost::math::statistics::simple_ordinary_least_squares(r, cubic_b_spline);
    assert(p.second < 0);
    std::cout << "Cubic B-spline    : " << p.first << " - " << std::abs(p.second) << "r\n";

    p  = boost::math::statistics::simple_ordinary_least_squares(r, quintic_b_spline);
    assert(p.second < 0);
    std::cout << "Quintic B-spline  : " << p.first << " - " << std::abs(p.second) << "r\n";

    p  = boost::math::statistics::simple_ordinary_least_squares(r, cubic_hermite);
    assert(p.second < 0);
    std::cout << "Cubic Hermite     : " << p.first << " - " << std::abs(p.second) << "r\n";

    p  = boost::math::statistics::simple_ordinary_least_squares(r, pchip);
    assert(p.second < 0);
    std::cout << "pchip             : " << p.first << " - " << std::abs(p.second) << "r\n";

    p  = boost::math::statistics::simple_ordinary_least_squares(r, makima);
    assert(p.second < 0);
    std::cout << "makima            : " << p.first << " - " << std::abs(p.second) << "r\n";

    p  = boost::math::statistics::simple_ordinary_least_squares(r, fotaylor);
    assert(p.second < 0);
    std::cout << "First-order Taylor: " << p.first << " - " << std::abs(p.second) << "r\n";

    if (sotaylor.size() > 0)
    {
    p  = boost::math::statistics::simple_ordinary_least_squares(r, quintic_hermite);
    assert(p.second < 0);
    std::cout << "Quintic Hermite   : " << p.first << " - " << std::abs(p.second) << "r\n";

    p  = boost::math::statistics::simple_ordinary_least_squares(r, sotaylor);
    assert(p.second < 0);
    std::cout << "Second-orde Taylor: " << p.first << " - " << std::abs(p.second) << "r\n";
    }

    if (totaylor.size() > 0)
    {
    p  = boost::math::statistics::simple_ordinary_least_squares(r, totaylor);
    assert(p.second < 0);
    std::cout << "Third-order Taylor: " << p.first << " - " << std::abs(p.second) << "r\n";

    p  = boost::math::statistics::simple_ordinary_least_squares(r, septic_hermite);
    assert(p.second < 0);
    std::cout << "septic hermite    : " << p.first << " - " << std::abs(p.second) << "r\n";
    }

}