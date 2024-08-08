#include "include.h"

#define NANOTYPE 4
//#define SINGLEBEAD 5

template <typename T>
T sphericalToEuclidean(T input)//theta,phi,r
{
        T output;
        output.x=sin(input.x)*cos(input.y)*input.z;
        output.y=sin(input.x)*sin(input.y)*input.z;
        output.z=cos(input.x)*input.z;
        return output;
}

template <typename T>
T euclideanToSpherical(T input)//theta,phi,r
{
        T output;
	output.z=sqrt(input.x*input.x+input.y*input.y+input.z*input.z);
        output.x=acos(input.z/output.z);
	output.y=atan2(input.y,input.x);
        return output;
}

int main(int argc, char **argv)
{
	if(argc!=10)
	{
		std::cout << "nArgs: " << argc << std::endl;
		
		//9, bilayer, nanoparticle
		std::cerr << "Usage: " << argv[0] << " name seed nLipids lipidArealDensity ";
		std::cerr << "nanoRadius nanoHeadUmin nanoHeadUmax nTessellations nNano " << std::endl;
		std::cerr << "\tFor liposome and nanoparticles." << std::endl << std::endl;
		return 0;
	}
	
	char *name=argv[1];
	std::stringstream cmdArg;
	for(int i=2;i<argc;i++)
		cmdArg << argv[i] << '\t';
	
	int seed, nLipids=0,nTessellations=0;
	double lipidArealDensity=0, nanoRadius=0, nanoHeadUmin=0;
	double nanoNanoUmin=0, nanoHeadUmax=0, nNano=0;
	double cutoffNP=2.0;
	//boost::vector<T>
	cmdArg >> seed >> nLipids >> lipidArealDensity >> nanoRadius >> nanoHeadUmin >> nanoHeadUmax >> nTessellations >> nNano;
	double rminNP=cutoffNP/2.0;
	
	std::cout << "Read arguments:" << std::endl;
	std::cout << "seed=" << seed << std::endl;
	std::cout << "nLipids=" << nLipids << std::endl;
	std::cout << "lipidArealDensity=" << lipidArealDensity << std::endl;
	std::cout << "nanoRadius=" << nanoRadius << std::endl;
	std::cout << "nanoHeadUmin=" << nanoHeadUmin << std::endl;
	std::cout << "nanoHeadUmax=" << nanoHeadUmax << std::endl;
	std::cout << "nTessellations=" << nTessellations << std::endl;
	std::cout << "nNano=" << nNano << std::endl;
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
	System.setFinalTime(20000);
	System.setDeltaT(0.02);
	System.setStoreInterval(1000);
	System.setMeasureInterval(100);
	//System.setDeltaLXY(0.02);//if you don't set it, it doesn't activate
	//System.setTension(1.0);
	System.setInitialTemp(3.0);
	System.setFinalTemp(3.0);
	//System.setSolventGamma(5.0);
	
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
	
	Umin[HEAD+5*System.readNTypes()]=nanoHeadUmin;
	Umax[HEAD+5*System.readNTypes()]=nanoHeadUmax;
	
	Umin[5+HEAD*System.readNTypes()]=nanoHeadUmin;
	Umax[5+HEAD*System.readNTypes()]=nanoHeadUmax;
	
	//two body constants
	std::vector<double> twoBodyFconst=mpd::laradjiRevaleeFC(Umax,Umin,CUTOFF,RMIN);
	for(auto &v:twoBodyFconst)
		System.addTwoBodyFconst(v);
	
	std::vector<double> twoBodyUconst=mpd::laradjiRevaleePC(Umax,Umin,CUTOFF,RMIN);
	for(auto &v:twoBodyUconst)
		System.addTwoBodyUconst(v);
	
	threeVector<double>size;
	size=0;
	size.z=nanoRadius*8.0+40.0;
	System.setSize(size);
	
	threeVector<double> pos;
	pos.x=0;
	pos.y=0;
	pos.z=20.0;
	double bondLength=0.7;
	
	std::vector<double> lipidConstants(4);
	lipidConstants[0]=bondLength;
	lipidConstants[1]=100;
	lipidConstants[2]=1;
	lipidConstants[3]=100;
	
	bilayer<double>(System, nLipids, 3, pos, 0.7, lipidArealDensity, &lipidConstants[0], 4, 1.0);
	
	MTRand randNum(System.readSeed());
	
	std::vector< threeVector<double> > nanoPos;
	
	for(int i=0;i<nNano;i++)
	{
		std::cout << "Placing nanoparticle " << i << "!" << std::endl;
		threeVector<double> toPlace;
		bool overlap;
		do
		{
			overlap=false;
			toPlace.x=(System.readSize().x-2.0*nanoRadius)*randNum.rand53()+nanoRadius;
			toPlace.y=(System.readSize().y-2.0*nanoRadius)*randNum.rand53()+nanoRadius;
			toPlace.z=pos.z+nanoRadius+rminNP+3.0*2.0*0.7+2.0; //Janus particle is a little larger
		} while(overlap);
		std::cerr << i << '\t' << toPlace.x << '\t' << toPlace.y << '\t' << toPlace.z << std::endl;
		nanoPos.push_back(toPlace);
	}
	
	for(auto &posNano:nanoPos)
		janus<double>(System,1200,45,250,50,0,posNano,nanoRadius,nTessellations);
	
	//To get radius of NP:
	//System.getMolecules
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

