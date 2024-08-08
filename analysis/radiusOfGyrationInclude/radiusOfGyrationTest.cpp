#include <iostream>
#include <vector>
#include <sstream>
#include <random>
#include <cmath>
#include "../../include/MD.h"
#include "../../include/system.h"
#include "radiusOfGyration.h"

int main(int argc, char **argv)
{
	if(argc!=4)
	{
		std::cerr << "Usage: " << argv[0] << " maxLength trials seed" << std::endl;
		std::cerr << "\tFinds radius of gyration for chains of monomers with " << std::endl;
		std::cerr << "a length of 1.0 between monomers. Pure random walk." << std::endl;
		return 0;
	}
	int trials,maxLength,seed;
	std::stringstream cmdArg;
	cmdArg << argv[1] << ' ' << argv[2] << ' ' << argv[3];
	cmdArg >> maxLength >> trials >> seed;
	std::mt19937 randNum(seed);
	std::vector<double> rog(maxLength,0.0);
	std::vector<int> indices;
	for(int i=0;i<maxLength;i++)
		indices.push_back(i);
	std::uniform_real_distribution<double> phiDist(0,M_PI);
	std::uniform_real_distribution<double> thetaDist(0,2.0*M_PI);
	threeVector<double> size=40000.0;
	for(int trial=0;trial<trials;++trial)
	{
		std::vector<position<double>> poly;
		for(int i=0;i<maxLength;i++)
		{
			position<double> p;
			if(i==0)
			{
				p.x=0;
				p.y=0;
				p.z=0;
				poly.push_back(p);
			}
			else
			{
				double theta=thetaDist(randNum);
				double phi=phiDist(randNum);
				p=poly.back();
				p.x+=std::sin(theta)*std::cos(phi);
				p.y+=std::sin(theta)*std::sin(phi);
				p.z+=std::cos(phi);
				poly.push_back(p);
				rog[i]+=radiusOfGyrationIndex(indices.begin(),indices.begin()+i,&poly[0],size);
			}
		}
	}
	int i=0;
	for(auto r:rog)
	{
		r/=trials;
		std::cout << i << ' ' << r << std::endl;
		i++;
	}
	return 0;
}
