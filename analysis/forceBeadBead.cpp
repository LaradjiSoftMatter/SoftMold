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
	
	std::vector<double> F, r;
	int nVal=floor((2.0*radius+rc)/dr);
	
	std::vector<double> C;
	
	for(int i=0;i<nBEADCONST;i++)
		C.push_back(0);
	
	double D=pow(2.0*M_PI*density*radius,2.0)/(rmin*rmin*rmin);
	
	C.push_back(2.0*radius+rmin);//0
	C.push_back(2.0*radius+rc);//1
	C.push_back(Umin*D*2.0/30.0);//2, x^6
	C.push_back(-Umin*D*7.0*rmin/20.0);//3, x^5
	C.push_back(Umin*D*rmin*rmin/2.0);//4, x^4
	
	double A=Umax-Umin;
	
	C.push_back(-D*A*rmin/20.0);//5,0@B^5
	C.push_back(D*A*rmin*rmin/12.0);//6,1,@B^4
	C.push_back(-D*Umin*pow(rmin,3.0)/6.0);//7,2,@B^3
	C.push_back(D*Umin*pow(rmin,4.0)/2.0);//8,3,@B^2
	C.push_back(D*1.3*Umin*pow(rmin,5.0)/2.0);//9,@B
	C.push_back(D*13.0*Umin*pow(rmin,6.0)/60.0);//10
	
	double rcSquared=pow(radius*2.0+rc,2.0);
	
	for(int i=0;i<nVal;i++)
	{
		r.push_back(static_cast<double>(i+1)*dr);
		
		threeVector<double> d;
		d.x=r[i];
		d.y=0;
		d.z=0;
		threeVector<double> force=0, nothing=0;
		beadBeadForce(d,force,nothing,&C[0],rcSquared);
		F.push_back(force.x);
		
		if(r[i]>2.0*radius)
			std::cout << r[i] << '\t' << F[i] << std::endl;
	}
	
	
	return 0;
}
