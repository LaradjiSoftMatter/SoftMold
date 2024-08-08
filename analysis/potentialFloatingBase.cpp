#include <iostream>
#include <sstream>
#include <vector>
#include <cmath>
#include "../include/MD.h"

int main(int argc, char ** argv)
{
	if(argc!=8)
	{
		std::cerr << "Usage: " << argv[0] << " Umax Umin rc rm density dz size" << std::endl;
		return 0;
	}
	double Umax,Umin,rc,rm,density,dz,size;
	std::stringstream cmdArg;
	for(int i=1;i<argc;i++)
		cmdArg << argv[i] << ' ';
	cmdArg >> Umax >> Umin >> rc >> rm >>density >> dz >> size;
	
	std::vector<double> U, F, r;
	threeVector<double> s=size;
	int nVal=floor(size/dz);
	std::cerr << "Umax=" << Umax << std::endl;
	std::cerr << "Umin=" << Umin << std::endl;
	std::cerr << "rc=" << rc << std::endl;
	std::cerr << "rm=" << rm << std::endl;
	std::cerr << "density=" << density << std::endl;
	std::cerr << "dz=" << dz << std::endl;
	std::cerr << "size=" << size << std::endl;
	
	double rcm=rc-rm;
	//more consistent constants
	std::vector<double> C;
	C.push_back(rm);//0, X direction
	C.push_back(rc);//1
	C.push_back(M_PI*density*Umin);//2
	C.push_back(2*M_PI*density*(Umax-Umin)/(rm*rm));//3
	C.push_back(2*M_PI*density*Umin*((2.0/5.0)*rcm*rcm*rcm*rcm*rcm-
		rm*(7.0/4.0)*rcm*rcm*rcm*rcm+2.0*rm*rm*rcm*rcm*rcm)/(rm*rm*rm));//4
	C.push_back(2*M_PI*density*Umin/(rm*rm*rm));//5
	
	
	for(int i=0;i<nVal;i++)
	{
		r.push_back(static_cast<double>(i)*dz);
		threeVector<double> d,a=0;
		d.x=0;
		d.y=0;
		d.z=r[i];
		U.push_back(floatingBasePotential(d,&C[0],s));
		floatingBaseForce(d,a,&C[0],s);
		F.push_back(a.z);
		std::cout << r[i] << '\t' << U[i] << '\t' << F[i] << std::endl;
	}
	
	
	return 0;
}
