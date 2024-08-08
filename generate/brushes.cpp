#include "include.h"

#define BB 0
#define BM 1
#define NANO 2

int main(int argc, char **argv)
{
	if(argc!=10)
	{
		//10 arguments
		std::cout << "Usage: " << argv[0] << " name seed nBrushes arealDensity UminBrushBrush UmaxBrushBrush nMonomers bondLength temp" << std::endl;
		return 0;
	}
	char *name=argv[1];
	std::stringstream cmdArg;
	
	for(int i=2;i<argc;i++)
		cmdArg << argv[i] << ' ';
	
	int seed,nBrushes,nMonomers;
	double arealDensity, UminBrushBrush,UmaxBrushBrush,bondLength,temp;
	
	cmdArg >> seed >> nBrushes >> arealDensity >> UminBrushBrush >> UmaxBrushBrush >> nMonomers >> bondLength >> temp;
	
	//the variables for the simulation
	Blob<double> System;
	
	System.setGamma(1.0);
	System.setNTypes(4);
	System.setSeed(seed);
	
	threeVector<bool> wrap;//periodic boundaries flag, be forwarned, this doesn't work yet, it is always on
	wrap.x=true;
	wrap.y=true;
	wrap.z=true;
	System.setPeriodic(wrap);
	System.setCutoff(2.0);
	
	System.setInitialTime(0);
	System.setFinalTime(100);
	System.setDeltaT(0.02);
	System.setStoreInterval(10);
	System.setMeasureInterval(1);
	//System.setDeltaLXY(0.01);//if you don't set it, it doesn't activate
	System.setInitialTemp(temp);
	System.setFinalTemp(temp);
	
	//initialize constants (force and potential constants)
	//generate force constants
	std::vector<double> Umax(System.readNTypes()*System.readNTypes());
	std::vector<double> Umin(System.readNTypes()*System.readNTypes());
	
	//could be one loop, but it's good for an example
	//of how to set constants
	for(int i=0;i<System.readNTypes();i++)
	{
		for(int j=0;j<System.readNTypes();j++)
		{
			//i=first type, j=second type, indexed grid
			int k=i+j*System.readNTypes();
			Umin[k]=0;
			Umax[k]=100;
		}
	}
	
	//Brush monomers
	Umin[BM+BM*System.readNTypes()]=UminBrushBrush;
	Umax[BM+BM*System.readNTypes()]=UmaxBrushBrush;
	
	//two body constants
	std::vector<double> twoBodyFconst=mpd::laradjiRevaleeFC(Umax,Umin,CUTOFF,RMIN);
	for(auto &v:twoBodyFconst)
		System.addTwoBodyFconst(v);
	
	std::vector<double> twoBodyUconst=mpd::laradjiRevaleePC(Umax,Umin,CUTOFF,RMIN);
	for(auto &v:twoBodyUconst)
		System.addTwoBodyUconst(v);
	
	//double radius=sqrt(((double)nBrushes/arealDensity)/(4.0*M_PI));//density adjusted
	threeVector<double>size;
	size=0;
	size.z=static_cast<double>(nMonomers)*bondLength*2.0;
	System.setSize(size);
	
	threeVector<double> pos;
	pos.x=0;
	pos.y=0;
	pos.z=0;
	
	double *constants=new double[4];
	constants[0]=bondLength;
	constants[1]=100;
	constants[2]=1;
	constants[3]=100;
	
	//simple brush
	brush<double>(System, nBrushes, nMonomers, pos, bondLength, arealDensity, constants, 4, 1.0, BB, BM);

	molecule< double, fourVector<int> > boundaryZ;
	boundaryZ.setType(BOUNDARY);
	for(int i=0;i<System.readNParticles();i++)
	{
		fourVector<int> buf;
		buf.s[0]=i;
		boundaryZ.addBond(buf);
	}
	boundaryZ.addConstant(2);//along Z
	boundaryZ.addConstant(0);//at z position
	boundaryZ.addConstant(2);//width of z
	boundaryZ.addConstant(20/2.0);//position of potential
	
	System.addMolecule(boundaryZ);

	std::cout << "Storing configuration...";
	
	//Write configuration
	Script<double,Blob <double> > output(name,std::ios::out,&System);
	output.write();
	
	//This stores the xyz file to check new configuration
	position<double> *p=System.getPositions();
	int nParticles=System.readNParticles();
	
	std::string newName("");
	newName+=name;
	newName+=".xyz";
	
	xyzFormat<double> xyzFile(p, nParticles);
	xyzFile.open(newName.c_str(), std::ios::out);
	xyzFile.store();
	xyzFile.close();
	
	std::cout << "Done.\nExiting...\n";
	
	delete constants;
	return 0;
}


