#include "include.h"

#define NANOTYPE 4
//#define SINGLEBEAD 5

int main(int argc, char **argv)
{
	if(argc!=9)
	{
		//9, bilayer, nanoparticle
		std::cout << "Usage: " << argv[0] << " name seed nLipids lipidArealDensity ";
		std::cout << "nanoRadius nanoHeadUmin nanoHeadUmax sigma" << std::endl;
		std::cout << "\tFor liposome and nanoparticles." << std::endl << std::endl;
		
		return 0;
	}
	
	char *name=argv[1];
	std::stringstream cmdArg;
	for(int i=2;i<argc;i++)
		cmdArg << argv[i] << '\t';
	
	int seed, nLipids=0;
	double lipidArealDensity=0, nanoRadius=0, nanoHeadUmin=0;
	double nanoNanoUmin=0, nanoHeadUmax=0, sigma=0;
	
	cmdArg >> seed >> nLipids >> lipidArealDensity >> nanoRadius >> nanoHeadUmin >> nanoHeadUmax >> sigma;
	
	//the variables for the simulation
	Blob<double> System;
	
	System.setGamma(1.0);
	System.setNTypes(8);
	System.setSeed(seed);
	
	threeVector<bool> wrap;//periodic boundaries flag, be forwarned, this doesn't work yet, it is always on
	wrap.x=true;
	wrap.y=true;
	wrap.z=true;
	System.setPeriodic(wrap);
	System.setCutoff(2.0);
	
	System.setInitialTime(0);
	System.setFinalTime(200000);
	System.setDeltaT(0.02);
	System.setStoreInterval(1000);
	System.setMeasureInterval(100);
	//System.setDeltaLXY(0.02);//if you don't set it, it doesn't activate
	System.setInitialTemp(3.0);
	System.setFinalTemp(3.0);
	
	//initialize constants (force and potential constants)
	//generate force constants
        std::vector<double> Umin(System.readNTypes()*System.readNTypes(),0.0);
        std::vector<double> Umax(System.readNTypes()*System.readNTypes(),0.0);
	
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
	
	//tail types hydrophobicity
	Umin[TAIL+TAIL*System.readNTypes()]=-6;
	Umax[TAIL+TAIL*System.readNTypes()]=200;
        
	Umin[TAIL2+TAIL2*System.readNTypes()]=-6;
	Umax[TAIL2+TAIL2*System.readNTypes()]=200;
	
	Umin[TAIL2+TAIL*System.readNTypes()]=-6;
	Umax[TAIL2+TAIL*System.readNTypes()]=200;
	
	Umin[TAIL+TAIL2*System.readNTypes()]=Umin[TAIL2+TAIL*System.readNTypes()];
	Umax[TAIL+TAIL2*System.readNTypes()]=Umax[TAIL2+TAIL*System.readNTypes()];
	
	Umin[HEAD2+HEAD*System.readNTypes()]=-6;
	Umax[HEAD2+HEAD*System.readNTypes()]=200;
	
	Umin[HEAD+HEAD2*System.readNTypes()]=Umin[HEAD2+HEAD*System.readNTypes()];
	Umax[HEAD+HEAD2*System.readNTypes()]=Umax[HEAD2+HEAD*System.readNTypes()];
	
	
	//two body constants
	std::vector<double> twoBodyFconst=mpd::laradjiRevaleeFC(Umax,Umin,CUTOFF,RMIN);
	for(auto &v:twoBodyFconst)
		System.addTwoBodyFconst(v);
	
	std::vector<double> twoBodyUconst=mpd::laradjiRevaleePC(Umax,Umin,CUTOFF,RMIN);
	for(auto &v:twoBodyUconst)
		System.addTwoBodyUconst(v);
	
	double lipoRadius=sqrt(((double)nLipids/lipidArealDensity)/(4.0*M_PI));
	threeVector<double>size;
	size=lipoRadius*4;
	System.setSize(size);
	
	threeVector<double> pos;
	pos.x=System.readSize().x/2.0;
	pos.y=System.readSize().y/2.0;
	pos.z=System.readSize().z/2.0;
	double bondLength=0.7;
	
	std::vector<double> lipidConstants(4);
	lipidConstants[0]=bondLength;
	lipidConstants[1]=100;
	lipidConstants[2]=1;
	lipidConstants[3]=100;
	
	lipoRadius=0;
	if(nLipids>0)
		lipoRadius=liposome(System, nLipids, 3, pos, 0.7, lipidArealDensity, &lipidConstants[0], 4);
	if(lipoRadius<=nanoRadius)
	{
		std::cerr << "Liposome, radius " << lipoRadius << ", too small compared to nanoparticle, radius " << nanoRadius << "!" << std::endl;
		return -1;
	}
	
	//Add nanoparticle constants
	std::vector<double> UminContinuum(System.readNTypes()*System.readNTypes(),0.0);
	std::vector<double> UmaxContinuum(System.readNTypes()*System.readNTypes(),10.0);
	
	//exceptions
	UminContinuum[HEAD+NANOTYPE*System.readNTypes()]=nanoHeadUmin;
	UminContinuum[NANOTYPE+HEAD*System.readNTypes()]=nanoHeadUmin;
	UmaxContinuum[HEAD+NANOTYPE*System.readNTypes()]=nanoHeadUmax;
	UmaxContinuum[NANOTYPE+HEAD*System.readNTypes()]=nanoHeadUmax;
	
	std::vector<double> C=mpd::laradjiSpanglerFC(UmaxContinuum, UminContinuum, CUTOFF, 
						     RMIN, nanoRadius-1, 5.88);
	//Add a nanoparticle
	std::vector< threeVector<double> > nanoPos;
	nanoPos.push_back(pos);
	
	continuumSphere<double>(System, &(nanoPos[0]), nanoPos.size(), nanoRadius, &C[0], 
				System.readNTypes()*System.readNTypes()*nBEADCONST, NANOTYPE);
	
	std::cerr << "Storing configuration...";
	
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
	
	std::cerr << "Done.\nExiting...\n";
	return 0;
}


