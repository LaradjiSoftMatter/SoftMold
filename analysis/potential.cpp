#include <iostream>
#include <sstream>
#include <vector>
#include <cmath>
#include "../include/MD.h"

int main(int argc, char ** argv)
{
	if(argc!=6)
	{
		std::cerr << "Usage: " << argv[0] << " Umax Umin rc rmin dr" << std::endl;
		return 0;
	}
	double Umax, Umin, rc, rmin, dr;
	std::stringstream cmdArg;
	for(int i=1;i<argc;i++)
		cmdArg << argv[i] << ' ';
	cmdArg >> Umax >> Umin >> rc >> rmin >> dr;
	
	std::vector<double> U, r;
	int nVal=floor((rc)/dr);
	
	//more consistent constants
	std::vector<double> C;
	C.push_back(rmin);//C8
	C.push_back((Umax-Umin)/(rmin*rmin));//C4
	C.push_back(Umin);//C5,no index trick
	C.push_back(rc);//C7
	C.push_back((3.0*Umin)/((rc-rmin)*(rc-rmin)));//C1
	C.push_back((2.0*Umin)/((rc-rmin)*(rc-rmin)*(rc-rmin)));//C0

	double rcSquared=pow(rc,2.0);
	for(int i=0;i<nVal;i++)
	{
		r.push_back(static_cast<double>(i+1)*dr);
		threeVector<double> d;
		d.x=r[i];
		d.y=0;
		d.z=0;
		U.push_back(nonBondedP(d,&C[0],rcSquared));
		
		std::cout << r[i] << '\t' << U[i] << std::endl;
	}
	
	
	return 0;
}
