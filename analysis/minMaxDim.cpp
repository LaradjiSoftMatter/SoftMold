#include <iostream>
#include <cstdlib>
#include <random>

//For molecular dynamics forces and potentials
#include "../include/MD.h"

//For the molecular dynamics variables
#include "../include/system.h"

int main(int argc, char **argv)
{
	if(argc!=4)
	{
		std::cout << "Usage: " << argv[0] << " name type dimension" << std::endl;
		return 0;
	}
	
	char *name=argv[1];
	std::stringstream cmdArg;
	for(int i=2;i<argc;i++)
		cmdArg << argv[i] << '\t';
	
	double minDistance;
	int type;
	std::string dimStr;
	cmdArg >> type >> dimStr;
	int dim=-1;
	if(dimStr=="x" || dimStr=="X" || dimStr=="0") dim=0;
	if(dimStr=="y" || dimStr=="Y" || dimStr=="1") dim=1;
	if(dimStr=="z" || dimStr=="Z" || dimStr=="2") dim=2;
	if(dim==-1)
	{
		std::cerr << "dimension must be one of {x,y,z,X,Y,Z,0,1,2}!" << std::endl;
		return -1;
	}
	
	//the variables for the simulation
	Blob<double> System;
	
	//file input-output script engine
	//name part (argv[1] in this case) is used with no extension (*.mpd)
	Script<double, Blob <double> > fileIO(argv[1],std::ios::in,&System);
	
	//read the variables in
	fileIO.read();
	
	//close the file that it is reading from
	fileIO.close();
	
	if(System.readNTypes()<=type)
	{
		std::cerr << "type out of range!" << std::endl;
		std::cerr << "\tnTypes: " << System.readNTypes() << '\t' << "type: " << type << std::endl;
		return -1;
	}
	
	auto p=System.getPositions();
	double min=-1,max=-1;
	for(int i=0;i<System.readNParticles();i++)
	{
		if(p[i].type==type)
		{
			if(min<0) min=p[i].s[dim];
			if(max<0) max=p[i].s[dim];
			if(min>p[i].s[dim]) min=p[i].s[dim];
			if(max<p[i].s[dim]) max=p[i].s[dim];
		}
	}
	
	std::cout << std::setprecision(std::numeric_limits<double>::digits10+1) 
		<< "min: " << min << "\nmax: " << max << std::endl;
	
	return 0;
}


