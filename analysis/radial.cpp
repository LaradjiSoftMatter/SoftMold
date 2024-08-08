//This is a very short and simple program that shows you how to modify a variable already present in an mpd script.
//The reason this is used is because some variables (forces, potentials, geometry, etc...) are not easily modified
// in the mpd script. Using a file that opens the mpd script, modifies, and then saves it (as below) is far easier
// than modifying the actual variable. Most of it looks like the run.cpp file, it just replaces the "run loop"
// with the modification.

//For molecular dynamics forces and potentials
#include "../include/MD.h"

//For the molecular dynamics variables
#include "../include/system.h"

class comDistance {
	public:
		comDistance(position<double> *positions, int nPositions, threeVector<double> centerOfMass)
		{
			p=positions;
			nP=nPositions;
			com=centerOfMass;
		};
		~comDistance(){};
		bool operator() (int i, int j)
		{
			if(i>nP || i<0)
			{
				std::cout << "Overflow: " << i << '\n';
				std::cout << "Upper: " << nP << '\n';
				std::cout << "Lower: 0\n";
				throw 0;
			}
			
			if(j>nP || j<0)
			{
				std::cout << "Overflow: " << j << '\n';
				std::cout << "Upper: " << nP << '\n';
				std::cout << "Lower: 0\n";
				throw 0;
			}
			
			return com.distanceSqr(p[i])<com.distanceSqr(p[j]);
		};
	private:
		position<double> *p;
		int nP;
		threeVector<double> com;
};

int main(int argc, char* argv[])
{
	if(argc<4)
	{
		//so simple, this can't possibly mess it up
		std::cout << "usage: " << argv[0] << " name nDivisions centerType1 centerType2 centerType3 ... centerTypeN\n";
		std::cout << "Output of radial distribution is only to the first encountered edge.\n";
		std::cout << "\"centerType\"s are the types to center about.\n";
		return 0;
	}
	
	///read in options
	char *name=argv[1];
	int nDivisions=atoi(argv[2]);
	int nCenterType=argc-3;
	
	int *centerType=new int[nCenterType];
	
	for(int i=3;i<argc;i++)
		centerType[i-3]=atoi(argv[i]);
	
	///Configuration variables
	Blob<double> System;
	
	///load old variables, then initialize them
	Script<double, Blob <double> > fileIO(name,std::ios::in,&System);

	fileIO.read();

	fileIO.close();

	///Locate center of mass for centerTypes
	threeVector<double> com;
	com.x=0;
	com.y=0;
	com.z=0;
	int nPCom=0;
	for(int i=0;i<System.readNParticles();i++)
	{
		position<double> p=System.getPositions()[i];
		for(int j=0;j<nCenterType;j++)
		{
			if(centerType[j]==p.type)
			{
				com.x+=p.x;
				com.y+=p.y;
				com.z+=p.z;
				nPCom++;
			}
		}
	}
	if(nPCom==0)
	{
		std::cout << "Types selected are not present!\n";
		delete centerType;
		return 0;
	}
	
	com.x/=(double)nPCom;
	com.y/=(double)nPCom;
	com.z/=(double)nPCom;
	
	//nearest edge
	double radius=System.readSize().x-com.x;
	radius=(System.readSize().y-com.y>radius)?System.readSize().y-com.y:radius;
	radius=(System.readSize().z-com.z>radius)?System.readSize().z-com.z:radius;
	
	//sortable list
	int *indexCom=new int[System.readNParticles()];
	for(int i=0;i<System.readNParticles();i++)
		indexCom[i]=i;
	comDistance comD(System.getPositions(),System.readNParticles(), com);
	
	//sort the list radially
	std::sort(indexCom, indexCom+System.readNParticles(), comD);
	
	//division counting list
	double *divisions=new double[System.readNTypes()];
	
	///generate distributions for nDivisions
	for(int i=0, index=0;i<nDivisions;i++)
	{
		position<double> *p=System.getPositions();
		for(int j=0;j<System.readNTypes();j++)
			divisions[j]=0;
		
		double r=sqrt(com.distanceSqr(p[indexCom[index]]));
		while(r<(double)(i+1)*radius/(double)nDivisions && index<System.readNParticles())
		{
			divisions[p[indexCom[index++]].type]++;
			r=sqrt(com.distanceSqr(p[indexCom[index]]));
		}
		
		double sliceVolume=(4.0/3.0)*M_PI*(pow((double)(i+1)*radius/(double)nDivisions,3.0)\
						-pow((double)i*radius/(double)nDivisions,3.0));
		
		std::cout << (radius/(double)nDivisions)*(double)i; 
		
		for(int j=0;j<System.readNTypes();j++)
		{
			divisions[j]/=sliceVolume;
			std::cout << '\t' << divisions[j];
		}
		std::cout << '\n';
	}
	
	delete divisions,centerType;
	
	return 0;
}

