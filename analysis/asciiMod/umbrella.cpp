#include <iostream>
#include <sstream>
#include <fstream>
#include <map>
#include <vector>
#include <cmath>
#include <random>

struct twoVector {
	double x,y;
};
using histType=std::map<int,double>;

struct fileInfo {
	std::string name;
	double x,k;
};



template <typename STREAM>
histType getHistogram(STREAM &stream, double binWidth)
{
	histType histogram;
	std::string line;
	while(std::getline(stream,line))
	{
		std::stringstream input;
		input << line;
		twoVector v;
		input >> v.x >> v.y;
		int bin=std::floor(v.y/binWidth);
		histogram[bin]++;
	}
	double sum=0;
	for(auto h:histogram)
		sum+=h.second*binWidth;
	for(auto &h:histogram)
		h.second/=sum;
	return histogram;
}

//needs x y file data, binwidth, maximum samples
int main(int argc, char **argv)
{
	if(argc<6)
	{
		std::cerr << "Usage: " << argv[0] << " binWidth temperature file1 x1 k1 ..." << std::endl;
		return 0;
	}
	std::stringstream cmdArg;
	for(int i=1;i<argc;i++)
		cmdArg << argv[i] << ' ';
	double binWidth,temperature;
	cmdArg >> binWidth >> temperature;
	
	std::vector<histType> histograms;
	
	std::string fileName;
	std::vector<fileInfo> files;
	while(cmdArg >> fileName)
	{
		fileInfo file;
		file.name=fileName;
		cmdArg >> file.x >> file.k;
		files.push_back(file);
		std::cerr << "Using " << file.name << " with x=" << file.x << " and k=" << file.k << std::endl;
	}
	
	for(auto f:files)
	{
		std::fstream current(f.name,std::ios::in);
		histograms.push_back(getHistogram(current,binWidth));
		for(auto h:histograms.back())
		{
			twoVector point;
			point.x=h.first*binWidth+binWidth/2.0;
			point.y=-temperature*std::log(h.second)-(f.k/2.0)*(point.x-f.x)*(point.x-f.x);
			std::cout << point.x << ' ' << point.y << std::endl;
		}
		std::cout << std::endl;
	}
	
	return 0;
}
			
