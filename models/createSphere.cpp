//makes a sphere and outputs in xyz format

#include "../include/algorithms/functions.h"

int main(int argc, char* argv[])
{
	if(argc!=3)
	{
		//so simple, this can't possibly mess it up
		std::cout << "usage: " << argv[0] << " diameter nPoints > filename.xyz\n";
		return 0;
	}
	
	double diameter;
	int nPoints;
	
	std::stringstream cmdArg;
	cmdArg << argv[1] << '\t' << argv[2];
	cmdArg >> diameter >> nPoints;

	double radius=diameter/2.0;
	
	double s=3.6/sqrt(static_cast<double>(nPoints));
	double length=0;
	double dz=2.0/static_cast<double>(nPoints);
	double z=1.0-dz/2.0;
	
	std::cout << nPoints << "\nSphere\n";
	for(int i=0;i<nPoints;i++)
	{
		position<double> p;
		
		double r=sqrt(1.0-z*z);
		//inner monomers to outer monomers
		p.x=r*cos(length)*radius;
		p.y=r*sin(length)*radius;
		p.z=z*radius;
		p.type=i/200;
		
		std::cout << p.type << '\t' << p.x << '\t' << p.y << '\t' << p.z << std::endl;
		
		z=z-dz;
		length=length+s/r;
	}
	
	return 0;
}

