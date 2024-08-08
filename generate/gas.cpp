#include "include.h"

int main(int argc, char **argv)
{
	if(argc!=6)
	{
		std::cout << "Usage: " << argv[0] << " name density size.x size.y size.z" << std::endl;
		return 0;
	}

	char *name=argv[1];
	double density=atof(argv[2]);
	threeVector<double> size;
	size.x=atof(argv[3]);
	size.y=atof(argv[4]);
	size.z=atof(argv[5]);
	
	//the variables for the simulation
	Blob<double> System;
	
	System.setGamma(1.0);
	System.setNTypes(3);
	System.setSeed(1234);
	
	threeVector<bool> wrap;//periodic boundaries flag, be forwarned, this doesn't work yet, it is always on
	wrap.x=true;
	wrap.y=true;
	wrap.z=true;
	System.setPeriodic(wrap);
	System.setCutoff(2.0);
	
	System.setInitialTime(0);
	System.setFinalTime(100);
	System.setDeltaT(0.02);
	System.setStoreInterval(2);
	System.setMeasureInterval(2);
	System.setInitialTemp(3.0);
	System.setFinalTemp(3.0);
	
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
	
	//two body constants
	std::vector<double> twoBodyFconst=mpd::laradjiRevaleeFC(Umax,Umin,CUTOFF,RMIN);
	for(auto &v:twoBodyFconst)
		System.addTwoBodyFconst(v);
	
	std::vector<double> twoBodyUconst=mpd::laradjiRevaleePC(Umax,Umin,CUTOFF,RMIN);
	for(auto &v:twoBodyUconst)
		System.addTwoBodyUconst(v);
	
	System.setSize(size);
	
	int nParticles=density*size.x*size.y*size.z;
	threeVector<int> nIndices;
	nIndices=pow(nParticles,1.0/3.0);
	threeVector<double> latticeLength;
	latticeLength.x=size.x/static_cast<double>(nIndices.x);
	latticeLength.y=size.y/static_cast<double>(nIndices.y);
	latticeLength.z=size.z/static_cast<double>(nIndices.z);
	
	MTRand *randNum=new MTRand(System.readSeed());
	
	for(int i=0;i<nIndices.x;i++)
	{
		for(int j=0;j<nIndices.y;j++)
		{
			for(int k=0;k<nIndices.z;k++)
			{
				position<double> pIn;
				pIn.x=static_cast<double>(i)*latticeLength.x;
				pIn.y=static_cast<double>(j)*latticeLength.y;
				pIn.z=static_cast<double>(k)*latticeLength.z;
				pIn.type=1;
				
				threeVector<double> v, a;
				double theta,phi;
				double Vrms=sqrt(3.0*System.readInitialTemp());
				
				//velocity
				theta=M_PI*randNum->rand53();
				phi=M_PI*2.0*randNum->rand53();
				v.x=Vrms*cos(phi)*sin(theta);
				v.y=Vrms*sin(phi)*sin(theta);
				v.z=Vrms*cos(theta);
				
				System.addParticle(pIn,v,a);
			}
		}
	}
	
	//Write configuration
	Script<double,Blob <double> > output(name,std::ios::out,&System);
	output.write();
	
	//This stores the xyz file to check new configuration
	position<double> *p=System.getPositions();
	nParticles=System.readNParticles();
	
	std::string newName("");
	newName+=name;
	newName+=".xyz";
	
	xyzFormat<double> xyzFile(p, nParticles);
	xyzFile.open(newName.c_str(), std::ios::out);
	xyzFile.store();
	xyzFile.close();
	
	return 0;
}


