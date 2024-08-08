#include "polynomialRegressionInclude/polynomialRegression.h"

#include <vector>
#include <iostream>
#include <random>
#include <algorithm>

int main()
{
	unsigned degree;
	std::cerr << "Degree? ";
	std::cin >> degree;
	std::cerr << degree << std::endl;
	unsigned nBootstrap;
	std::cerr << "Number of bootstrap iterations? ";
	std::cin >> nBootstrap;
	std::cerr << nBootstrap << std::endl;
	double fBootstrap;
	std::cerr << "Fraction of bootstrap iterations? ";
	std::cin >> fBootstrap;
	std::cerr << fBootstrap << std::endl;
	int seed;
	std::cerr << "Seed? ";
	std::cin >> seed;
	std::cerr << seed << std::endl;
	
	std::cerr << "X and Y values" << std::endl;
	
	std::vector<double> X,Y;
	double x,y;
	while(std::cin >> x >> y)
	{
		X.push_back(x);
		Y.push_back(y);
	}
	
	unsigned lBootstrap=fBootstrap*X.size();
	
	std::vector<std::vector<double>> coeffs;
	
	std::mt19937_64 randGen(seed);
	std::vector<unsigned> ind(X.size());
	for(unsigned i=0;i<X.size();i++)
		ind[i]=i;
	
	std::vector<double> bX(lBootstrap),bY(lBootstrap);
	for(unsigned i=0;i<nBootstrap;i++)
	{
		std::shuffle(ind.begin(),ind.end(),randGen);
		for(unsigned j=0;j<lBootstrap;j++)
		{
			bX[j]=X[ind[j]];
			bY[j]=Y[ind[j]];
		}
		coeffs.push_back(mpd::polynomialRegression(bX,bY,degree));
		
	}
	std::vector<double> avg(degree+1,0.0);
	for(auto coeff:coeffs)
		for(unsigned i=0;i<coeff.size();i++)
			avg[i]+=coeff[i];
	for(auto &a:avg)
		a/=coeffs.size();
	std::vector<double> std(degree+1,0.0);
	for(auto coeff:coeffs)
		for(unsigned i=0;i<coeff.size();i++)
			std[i]+=std::pow(coeff[i]-avg[i],2);
	for(auto &s:std)
	{
		s/=coeffs.size();
		s=std::sqrt(s);
	}
	
	std::cerr << "Coefficients and standard deviations from highest to lowest order:" << std::endl;
	for(unsigned i=0;i<avg.size();i++)
		std::cout << avg[i] << ' ' << std[i] << ' ';
	std::cout << std::endl;
	
	return 0;
}
