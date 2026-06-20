/* Based on "Construction of higher order symplectic integrators" by Haruo Yoshida*/

#include <vector>
#include <cmath>
#include <valarray>

typedef std::valarray<double> derivative(std::valarray<double>);

template <class RealType>
std::pair<std::valarray<RealType>, std::valarray<RealType> > second_order_yoshida(const std::valarray<RealType> p0, 
                                                                                  const std::valarray<RealType> q0, 
                                                                                  RealType dt, 
                                                                                  derivative dHdp, 
                                                                                  derivative dHdq)
{
    std::valarray<RealType> p = p0;
    std::valarray<RealType> q = q0;

    q = q + dt / 2 * dHdp(p);
    p = p - dt * dHdq(q);
    q = q + dt / 2 * dHdp(p);
    
    return std::make_pair(p, q);
}

template <class RealType>
std::pair<std::valarray<RealType>, std::valarray<RealType> > fourth_order_yoshida(const std::valarray<RealType> p0, 
                                                                                  const std::valarray<RealType> q0, 
                                                                                  RealType dt, 
                                                                                  derivative dHdp, 
                                                                                  derivative dHdq)
{
    std::valarray<RealType> p = p0;
    std::valarray<RealType> q = q0;

    RealType x0 = -(std::pow(2, 1/3) / (2 - std::pow(2, 1/3)));
    RealType x1 = 1 / (2 - std::pow(2, 1/3));

    std::vector<RealType> weights = { x1, x0, x1 };

    for (unsigned i=0; i < weights.size(); i++)
    {
        std::tie(p, q) = second_order_yoshida(p, q, weights[i] * dt, dHdp, dHdq);
    }
    
    return std::make_pair(p, q);
}

template <class RealType>
std::pair<std::valarray<RealType>, std::valarray<RealType> > sixth_order_yoshida(const std::valarray<RealType> p0, 
                                                  const std::valarray<RealType> q0, 
                                                  RealType dt, 
                                                  derivative dHdp, 
                                                  derivative dHdq)
{
    std::valarray<RealType> p = p0;
    std::valarray<RealType> q = q0;
    
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