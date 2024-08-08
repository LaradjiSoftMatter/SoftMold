#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <algorithm>
#include <cmath>

struct twoVector {
	double x,y;
};

bool sortByX(twoVector a, twoVector b)
{
	return a.x<b.x;	
}

int main(int argc, char **argv)
{
	if(argc!=3)
	{
		std::cerr << "Usage: " << argv[0] << " stdDev nPoints < inFile > outFile" << std::endl;
		std::cerr << " Data must be formated as:" << std::endl;
		std::cerr << "   x0 y0 x1 y1 x2 y2 x3 y3..." << std::endl;
		std::cerr << " Where any newline, space, or tab is whitespace!" << std::endl;
		std::cerr << "stdDev is standard deviation and nPoints is for interpolation" << std::endl;
		return 0;
	}
	
	//std::cerr << argc << std::endl;
	double stdDev=0;
	int nPoints=0;
	std::stringstream cmdArg;
	for(int i=1;i<argc;i++)
		cmdArg << argv[i] << ' ';
	cmdArg >> stdDev >> nPoints;
	
	std::vector<twoVector> v;
	std::string line;
	while(std::getline(std::cin,line,'\n'))
	{
		std::stringstream splitLine;
		splitLine << line;
		twoVector newV;
		splitLine >> newV.x >> newV.y;
		v.push_back(newV);
	}
	
	std::sort(v.begin(), v.end(), sortByX);
	double xMin=v.front().x;
	double xMax=v.back().x;
	double xInc=(xMax-xMin)/static_cast<double>(nPoints-1);
	
	std::vector<twoVector> gSum;
	for(double i=0;i<nPoints;i++)
	{
		twoVector g;
		g.x=xInc*i+xMin;
		g.y=0;
		gSum.push_back(g);
	}
	
	//double totalArea=0;
	
	//std::cerr << gSum.size() << std::endl;
	for(int i=1;i<v.size();i++)
	{
		twoVector mP;
		mP.x=(v[i].x+v[i-1].x)*0.5;
		mP.y=(v[i].y+v[i-1].y)*0.5;
		double inc=v[i].x-v[i-1].x;
		//double area=inc*(v[i].y-v[i-1].y)*0.5+inc*(v[i].y<v[i-1].y?v[i].y:v[i-1].y);
		double area=inc*mP.y;
		//totalArea+=area;
		for(auto& g:gSum)
			g.y+=area*exp(-pow((g.x-mP.x)/stdDev,2.0)/2.0)/(sqrt(2.0*M_PI)*stdDev);
	}
	//std::cerr << gSum.size() << std::endl;
	
	//double gArea=0;
	//for(int i=1;i<gSum.size();i++)
	//{
	//	twoVector mP;
	//	mP.x=(gSum[i].x+gSum[i-1].x)*0.5;
	//	mP.y=(gSum[i].y+gSum[i-1].y)*0.5;
	//	double inc=gSum[i].x-gSum[i-1].x;
	//	//double area=inc*(gSum[i].y-gSum[i-1].y)*0.5+inc*(gSum[i].y<gSum[i-1].y?gSum[i].y:gSum[i-1].y);
	//	double area=inc*mP.y;
	//	gArea+=area;
	//}
	//double correction=totalArea/gArea;
	//std::cerr << "totalArea: " << totalArea << " gArea: " << gArea << std::endl;
	for(auto& g:gSum)
		std::cout << g.x << ' ' << g.y << std::endl;//g.y*correction << std::endl;
	return 0;
}
