/* Based on "Construction of higher order symplectic integrators" by Haruo Yoshida*/

#include <vector>
#include <cmath>
#include <valarray>

typedef std::valarray<double> derivative(std::valarray<double>);

template <class RealType, class ReturnType, class Func>
std::pair<ReturnType, ReturnType> second_order_yoshida(const ReturnType p0, 
                                                       const ReturnType q0, 
                                                       RealType dt, 
                                                       Func dHdp, 
                                                       Func dHdq)
{
    ReturnType p = p0;
    ReturnType q = q0;

    q = q + dt / 2 * dHdp(p);
    p = p - dt * dHdq(q);
    q = q + dt / 2 * dHdp(p);
    
    return std::make_pair(p, q);
}

template <class RealType, class ReturnType, class Func>
std::pair<ReturnType, ReturnType> fourth_order_yoshida(const ReturnType p0, 
                                                       const ReturnType q0, 
                                                       const RealType dt, 
                                                       Func dHdp, 
                                                       Func dHdq)
{
    ReturnType p = p0;
    ReturnType q = q0;

    RealType x0 = -(std::pow(2, 1/3) / (2 - std::pow(2, 1/3)));
    RealType x1 = 1 / (2 - std::pow(2, 1/3));

    std::vector<RealType> weights = { x1, x0, x1 };

    for (unsigned i=0; i < weights.size(); i++)
    {
        std::tie(p, q) = second_order_yoshida(p, q, weights[i] * dt, dHdp, dHdq);
    }
    
    return std::make_pair(p, q);
}

template <class RealType, class ReturnType, class Func>
std::pair<ReturnType, ReturnType> sixth_order_yoshida(const ReturnType p0, 
                                                       const ReturnType q0, 
                                                       RealType dt, 
                                                       Func dHdp, 
                                                       Func dHdq)
{
    ReturnType p = p0;
    ReturnType q = q0;
    
    // Choosing "System A" solution
    RealType w3 = 0.784513610477560;
    RealType w2 = 0.235573213359357;
    RealType w1 = -1.17767998417887;
    RealType w0 = 1 - 2 * (w1 + w2 + w3);
    std::vector<RealType> weights = { w3, w2, w1, w0, w1, w2, w3};

    for (unsigned i=0; i < weights.size(); i++)
    {
        std::tie(p, q) = second_order_yoshida(p, q, weights[i] * dt, dHdp, dHdq);
    }
    
    return std::make_pair(p, q);
}