#include <iostream>
#include <sstream>
#include <vector>
#include <cmath>
#include "../include/MD.h"

int main(int argc, char ** argv)
{
	if(argc!=6)
	{
		std::cerr << "Usage: " << argv[0] << " offset ahalf kbond dr size" << std::endl;
		return 0;
	}
	double offset,ahalf,kbond,dr,size;
	std::stringstream cmdArg;
	for(int i=1;i<argc;i++)
		cmdArg << argv[i] << ' ';
	cmdArg >> offset >> ahalf >> kbond >> dr >> size;
	
	std::vector<double> U, r;
	threeVector<double> s=size;
	int nVal=2*floor((ahalf*2)/dr);
	
	//more consistent constants
	std::vector<double> C;
	C.push_back(2);//0, X direction
	C.push_back(offset);//1
	C.push_back(ahalf);//2
	C.push_back(kbond);//3
	
	for(int i=0;i<nVal;i++)
	{
		r.push_back(static_cast<double>(i-nVal/2.0)*dr);
		position<double> d;
		d.x=0;
		d.y=0;
		d.z=r[i];
		U.push_back(boundaryP(d,&C[0],s));
		std::cout << r[i] << '\t' << U[i] << std::endl;
	}
	
	
	return 0;
}
