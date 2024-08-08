//For molecular dynamics forces and potentials
#include "../include/MD.h"

//For the molecular dynamics variables
#include "../include/system.h"

int main(int argc, char* argv[])
{
	if(argc!=6)
	{
		//so simple, this can't possibly mess it up
		std::cout << "usage: " << argv[0] << " name type offsetX offsetY offsetZ\n";
		return 0;
	}
	
	///read in options
	char *name=argv[1];
	
	std::stringstream cmdArg;
	for(int i=2;i<argc;i++)
		cmdArg << argv[i] << ' ';
	threeVector<double> offset;
	int type;
	cmdArg >> type >> offset.x >> offset.y >> offset.z;
	
	///Configuration variables
	Blob<double> System;
	
	///load old variables, then initialize them
	Script<double, Blob <double> > fileIO(name,std::ios::in,&System);
	fileIO.read();
	fileIO.close();
	
	if(System.readNTypes()<=type)
	{
		std::cerr << "type out of range!" << std::endl;
		std::cerr << "\tnTypes: " << System.readNTypes() << '\t' << "type: " << type << std::endl;
		return -1;
	}
	
	position<double> *p=System.getPositions();
	threeVector<double> *v=System.getVelocities();
	threeVector<double> *a=System.getAccelerations();
	
	int nParticles=System.readNParticles();
	
	for(int i=0;i<nParticles;i++)
	{
		position<double> buf=p[i];
		if(buf.type==type)
		{
			buf.x+=offset.x;
			buf.y+=offset.y;
			buf.z+=offset.z;
			//this automatically puts particles within system bounds
			System.setParticle(i, buf, v[i], a[i]);
		}
	}
	
	//Save it
	fileIO.open(name,std::ios::out);
	fileIO.write();
	fileIO.close();
	
	return 0;
}

