#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include "../include/MD.h"

int main(int argc, char **argv)
{
	if(argc!=6)
	{
		std::cerr << "Usage: " << argv[0] << " size center[0,1.0] halfWidth k ds" << std::endl;
		return 0;
	}
	
	std::stringstream cmdArg;
	threeVector<double> size;
	double center, k, halfWidth,ds;
	for(int i=1;i<argc;i++)
		cmdArg << argv[i] << ' ';
	cmdArg >> size.x >> center >> halfWidth >> k >> ds; 

	std::vector<double> C;
	C.push_back(0);
	C.push_back(center);
	C.push_back(halfWidth);
	C.push_back(k/2.0);
	
	for(double pos=0;pos<size.x;pos+=ds)
	{
		position<double> p;
		p.x=pos;
		p.y=0;
		p.z=0;
		threeVector<double> a=0;
		boundaryF(p,a,&C[0],size);
		std::cout << pos << '\t' << boundaryP(p, &C[0],size) << '\t' << a.x << std::endl;
	}
	
	return 0;
}
