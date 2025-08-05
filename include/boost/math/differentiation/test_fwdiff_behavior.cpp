#include "autodiff.hpp"
#include<iostream>

using namespace boost::math::differentiation;

double myabs(double x)
{
	return (x>0) - (x<0);
}
int main()
{
	float y = -0.26647448539733887;
	float z = 0.8725948333740234;
	float w = abs(y);
	std::cout<<w<<std::endl;
	return 0;
}
