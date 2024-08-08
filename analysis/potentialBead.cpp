#include <iostream>
#include <sstream>
#include <vector>
#include <cmath>
#include "../include/MD.h"

int main(int argc, char ** argv)
{
	if(argc!=8)
	{
		std::cerr << "Usage: " << argv[0] << " Umax Umin rc rmin radius density dr" << std::endl;
		return 0;
	}
	double Umax, Umin, rc, rmin, radius, density, dr;
	std::stringstream cmdArg;
	for(int i=1;i<argc;i++)
		cmdArg << argv[i] << ' ';
	cmdArg >> Umax >> Umin >> rc >> rmin >> radius >> density >> dr;
	
	std::vector<double> U, r;
	int nVal=floor((radius+rc)/dr);
	
	//more consistent constants
	std::vector<double> C;
	C.push_back(rc+radius);//0
	C.push_back(-7.0/4.0*rmin);//1
	C.push_back(2.0*rmin*rmin);//2
	C.push_back(Umin*M_PI*radius*density/(rmin*rmin*rmin));//3
	C.push_back(radius);//4
	C.push_back(rmin);//5
	
	double D=density*M_PI*radius;
	double A=Umax-Umin;
	
	C.push_back(-D*A/(2.0*rmin*rmin));//6,0@B^4
	C.push_back(2.0*D*A/(3.0*rmin));//7,1,@B^3
	C.push_back(-D*Umin);//8,2,@B^2
	C.push_back(2.0*D*Umin*rmin);//9,3,@B^1
	C.push_back(D*1.3*Umin*rmin*rmin);//10,4
	
	double rcSquared=pow(radius+rc,2.0);
	for(int i=0;i<nVal;i++)
	{
		r.push_back(static_cast<double>(i+1)*dr);
		threeVector<double> d;
		d.x=r[i];
		d.y=0;
		d.z=0;
		U.push_back(beadPotential(d,&C[0],rcSquared));
		
		if(r[i]>radius)
			std::cout << r[i] << '\t' << U[i] << std::endl;
	}
	
	
	return 0;
}
