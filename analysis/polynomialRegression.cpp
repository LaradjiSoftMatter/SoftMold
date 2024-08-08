#include "polynomialRegressionInclude/polynomialRegression.h"

#include <vector>
#include <iostream>


int main()
{
	unsigned degree;
	std::cerr << "Degree? ";
	std::cin >> degree;
	
	std::cerr << "X and Y values" << std::endl;
	
	std::vector<double> X,Y;
	double x,y;
	while(std::cin >> x >> y)
	{
		X.push_back(x);
		Y.push_back(y);
	}
	
	auto coeff=mpd::polynomialRegression(X,Y,degree);
	std::cerr << "Coefficients from highest to lowest order:" << std::endl;
	for(auto c:coeff)
		std::cout << c << ' ';
	std::cout << std::endl;
	return 0;
}
