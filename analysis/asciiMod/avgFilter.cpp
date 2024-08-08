#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <algorithm>
#include <cmath>

struct twoVector {
	double x,y;
};

bool sortByY(twoVector a, twoVector b)
{
	return a.y<b.y;	
}

bool sortByX(twoVector a, twoVector b)
{
	return a.x<b.x;	
}

int main(int argc, char **argv)
{
	if(argc!=2)
	{
		std::cerr << "Usage: " << argv[0] << " window < inFile > outFile" << std::endl;
		std::cerr << " Data must be formated as:" << std::endl;
		std::cerr << "   x0 y0 x1 y1 x2 y2 x3 y3..." << std::endl;
		std::cerr << " Where any newline, space, or tab is whitespace!" << std::endl;
		std::cerr << "window is the length of the window to perform median filtering" << std::endl;
		return 0;
	}
	
	int window=0;
	std::stringstream cmdArg;
	for(int i=1;i<argc;i++)
		cmdArg << argv[i] << ' ';
	cmdArg >> window;
	int offset=window/2;
	if(window%2==0)
		offset+=1;
	
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
	std::vector<twoVector> out;
	for(int i=0;i<v.size();i++)
	{
		twoVector res=v[i];
		if(i-offset>=0 && i+offset<v.size())
		{
			double avg=0;
			for(int j=i-offset;j<=i+offset;j++)
				avg+=v[j].y;
			res.y=avg/(2*offset+1);
		}
		out.push_back(res);
	}
	
	for(auto& o:out)
		std::cout << o.x << ' ' << o.y << std::endl;
	return 0;
}
