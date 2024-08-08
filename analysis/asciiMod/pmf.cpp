#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <map>
#include <random>

template <typename T>
int getIndex(const T &input, const T &span, const T &min)
{
	return std::floor((input-min)/span);
}

template <typename T>
T getValue(const int &input, const T &span, const T &min)
{
	return input*span+span/2.0+min;
}

template<typename T, typename V>
std::map<int,T> getHistogram(V begin, V end, const T &span, const T &min)
{
	std::map<int,double> histogram;
	for(;begin!=end;++begin)
		histogram[getIndex(*begin,span,min)]++;
	double total=0;
	for(auto &h:histogram)
		total+=h.second*span;
	for(auto &h:histogram)
		h.second/=total;
	return histogram;
}

template <typename T>
struct data {
	T value,stdDev,count;
};

int main(int argc, char **argv)
{
	if(argc!=5)
	{
		std::cerr << "Usage: " << argv[0] << " temperature nPoints nBootstrap seed < inFile > outFile" << std::endl;
		std::cerr << " Data must be formated as time series (t) where y are histogram values:" << std::endl;
		std::cerr << "   t0 y0 t1 y1 t2 y2 t3 y3..." << std::endl;
		std::cerr << " Where any newline, space, or tab is whitespace!" << std::endl;
		std::cerr << " temperature is system temperature" << std::endl; 
		std::cerr << " nPoints is output/histogram points" << std::endl; 
		std::cerr << " nBootstrap is bootstrap iterations" << std::endl;
		std::cerr << " seed is for monte-carlo part of bootstrap" << std::endl;
		std::cerr << " output is -temperature*log(P(y))" << std::endl;
		return 0;
	}
	
	//std::cerr << argc << std::endl;
	double temperature=0;
	int nBootstrap=0,nPoints=0,seed=0;
	std::stringstream cmdArg;
	for(int i=1;i<argc;i++)
		cmdArg << argv[i] << ' ';
	cmdArg >> temperature >> nPoints >> nBootstrap >> seed;
	
	std::vector<double> v;
	std::string line;
	while(std::getline(std::cin,line,'\n'))
	{
		std::stringstream splitLine;
		splitLine << line;
		double t,y;
		splitLine >> t >> y;
		v.push_back(y);
	}
	
	std::sort(v.begin(), v.end());
	double xMin=v.front();
	double xMax=v.back();
	double xInc=(xMax-xMin)/static_cast<double>(nPoints-1);
	
	std::map<int,double> histogram=getHistogram(v.begin(),v.end(),xInc,xMin);
	std::map<int,data<double>> pmf;
	
	for(auto &h:histogram)
	{
		pmf[h.first].value=-temperature*std::log(h.second);
		pmf[h.first].stdDev=0;
		pmf[h.first].count=0;
	}
	
	std::mt19937_64 bsRand(seed);
	std::uniform_int_distribution<int> dist(0,v.size()-1);
	for(int i=0;i<nBootstrap;i++)
	{
		std::vector<double> bV;
		for(int i=0;i<v.size();i++)
			bV.push_back(v[dist(bsRand)]);
		std::map<int,double> bvHistogram=getHistogram(bV.begin(),bV.end(),xInc,xMin);
		for(auto &h:bvHistogram)
		{
			if(h.second>0)
			{
				double std=-temperature*std::log(h.second)-pmf[h.first].value;
				pmf[h.first].stdDev+=std*std;
				if(std!=std)
					std::cerr << "Error: " << h.second << ' ' << h.first << std::endl;
				pmf[h.first].count++;
			}
		}
	}
	
	for(auto &p:pmf)
	{
		std::cout << getValue(p.first,xInc,xMin) << ' ' << p.second.value;
		if(p.second.count!=0)
			std::cout << ' ' << std::sqrt(p.second.stdDev/p.second.count);
		else
			std::cout << ' ' << 0.0;
		std::cout << std::endl;
	}
	return 0;
}
