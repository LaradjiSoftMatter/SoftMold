//update from JanusP_two.cpp, change the number of types of particle from two to four, remove the interaction inside the particle
//Eric: update from JanusP.cpp, I've reduced this to 1 NP, and I've made it insert by position, rather than assuming
// it is resting on a vesicle/liposome. My thought is that we would add these 1 by 1, and if we
// need bonds between NPs, we can use the addBond program to do so.

#include "include.h"

int main(int argc, char* argv[])
{
	if(argc!=14)
	{
		//so simple, this can't possibly mess it up
		std::cout << "Usage: " << argv[0] << " oldName newName vBond cBond kbend particlebond ratio pX pY pZ particleradius tesselation Umin\n";
		return 0;
	}
	
	//read in options to parameters
	char *oldName=argv[1];
	char *newName=argv[2];
	
	threeVector<double> pos;
	double vBond, cBond, kbend, ratio, pbond, particleradius,Umin;
	int tesselation;
	std::stringstream cmdIn;
	for(int i=3;i<argc;i++)
		cmdIn << argv[i] << ' ';
	cmdIn >> vBond >> cBond >> kbend >> pbond >> ratio >> pos.x >> pos.y >> pos.z >> particleradius >> tesselation >> Umin;
	if(tesselation<=0)
	{
		std::cerr << "tesselation must be greater than 0!" << std::endl;
		return -1;
	}
	
	//Configuration variables
	Blob<double> System;
	
	//load old variables, then initialize them
	Script<double, Blob <double> > fileIO(oldName,std::ios::in,&System);
	fileIO.read();
	fileIO.close();
	
	int Nbefore = System.readNTypes();
	
	if(Nbefore<8)
		System.setNTypes(8);
	//set enough type of particles
	
	int Nafter = System.readNTypes();
	
	std::cout << Nbefore << " " << Nafter << std::endl;
	
	for(int i = 0; i < (Nafter*Nafter*6 - Nbefore*Nbefore*6); i++)
	{
		System.addTwoBodyFconst(0);
		System.addTwoBodyUconst(0);
	}
	
	//set constant for interaction between JanusParticle and head
	//head = 2; attract part on Janus particle is 4
	///initialize constants (force and potential constants)
	//generate force constants
	std::vector<double> UmaxS(System.readNTypes()*System.readNTypes());
	std::vector<double> UminS(System.readNTypes()*System.readNTypes());
	
	//could be one loop, but it's good for an example
	//of how to set constants
	for(int i=0;i<System.readNTypes();i++)
	{
		for(int j=0;j<System.readNTypes();j++)
		{
			//i=first type, j=second type, indexed grid
			int k=i+j*System.readNTypes();
			UminS[k]=0;
			UmaxS[k]=100;
		}
	}
	
	//umin and umax exceptions, tail types
	UminS[TAIL+TAIL*System.readNTypes()]=-6;
	UmaxS[TAIL+TAIL*System.readNTypes()]=200;
	
	UminS[HEAD+4*System.readNTypes()] = Umin;
	UmaxS[HEAD+4*System.readNTypes()] = 200;
	UminS[4+HEAD*System.readNTypes()] = UminS[HEAD+4*System.readNTypes()];
	UmaxS[4+HEAD*System.readNTypes()] = UmaxS[HEAD+4*System.readNTypes()];
	
	UminS[HEAD+5*System.readNTypes()] = 0;
	UmaxS[HEAD+5*System.readNTypes()] = 100;
	UminS[5+HEAD*System.readNTypes()] = UminS[HEAD+5*System.readNTypes()];
	UmaxS[5+HEAD*System.readNTypes()] = UmaxS[HEAD+5*System.readNTypes()];
	
	UminS[HEAD+6*System.readNTypes()] = Umin;
	UmaxS[HEAD+6*System.readNTypes()] = 200;
	UminS[6+HEAD*System.readNTypes()] = UminS[HEAD+6*System.readNTypes()];
	UmaxS[6+HEAD*System.readNTypes()] = UmaxS[HEAD+6*System.readNTypes()];
	
	UminS[HEAD+7*System.readNTypes()] = 0;
	UmaxS[HEAD+7*System.readNTypes()] = 100;
	UminS[7+HEAD*System.readNTypes()] = UminS[HEAD+7*System.readNTypes()];
	UmaxS[7+HEAD*System.readNTypes()] = UmaxS[HEAD+7*System.readNTypes()];
	
	//remove two body interaction inside the particle
	UminS[4+4*System.readNTypes()] = 0;
	UmaxS[4+4*System.readNTypes()] = 0;
	
	UminS[5+5*System.readNTypes()] = 0;
	UmaxS[5+5*System.readNTypes()] = 0;
	
	UminS[6+6*System.readNTypes()] = 0;
	UmaxS[6+6*System.readNTypes()] = 0;
	
	UminS[7+7*System.readNTypes()] = 0;
	UmaxS[7+7*System.readNTypes()] = 0;
	
	UminS[4+5*System.readNTypes()] = 0;
	UmaxS[4+5*System.readNTypes()] = 0;
	UminS[5+4*System.readNTypes()] = UminS[4+5*System.readNTypes()];
	UmaxS[5+4*System.readNTypes()] = UmaxS[4+5*System.readNTypes()];
	
	UminS[6+7*System.readNTypes()] = 0;
	UmaxS[6+7*System.readNTypes()] = 0;
	UminS[7+6*System.readNTypes()] = UminS[6+7*System.readNTypes()];
	UmaxS[7+6*System.readNTypes()] = UmaxS[6+7*System.readNTypes()];
	
	//increase the force between particles
	UminS[4+6*System.readNTypes()] = 0;
	UmaxS[4+6*System.readNTypes()] = 200;
	UminS[6+4*System.readNTypes()] = UminS[4+6*System.readNTypes()];
	UmaxS[6+4*System.readNTypes()] = UmaxS[4+6*System.readNTypes()];
	
	UminS[4+7*System.readNTypes()] = 0;
	UmaxS[4+7*System.readNTypes()] = 200;
	UminS[7+4*System.readNTypes()] = UminS[4+7*System.readNTypes()];
	UmaxS[7+4*System.readNTypes()] = UmaxS[4+7*System.readNTypes()];
	
	UminS[5+6*System.readNTypes()] = 0;
	UmaxS[5+6*System.readNTypes()] = 200;
	UminS[6+5*System.readNTypes()] = UminS[5+6*System.readNTypes()];
	UmaxS[6+5*System.readNTypes()] = UmaxS[5+6*System.readNTypes()];
	
	UminS[5+7*System.readNTypes()] = 0;
	UmaxS[5+7*System.readNTypes()] = 200;
	UminS[7+5*System.readNTypes()] = UminS[5+7*System.readNTypes()];
	UmaxS[7+5*System.readNTypes()] = UmaxS[5+7*System.readNTypes()];
	
	
	
	
	//two body constants
	std::vector<double> twoBodyFconst=mpd::laradjiRevaleeFC(UmaxS,UminS,CUTOFF,RMIN);
	for(auto &v:twoBodyFconst)
		System.addTwoBodyFconst(v);
	
	std::vector<double> twoBodyUconst=mpd::laradjiRevaleePC(UmaxS,UminS,CUTOFF,RMIN);
	for(auto &v:twoBodyUconst)
		System.addTwoBodyUconst(v);
	
	janus(System,vBond,cBond,kbend,pbond,ratio,pos,particleradius,tesselation);
	
	//Save system
	fileIO.open(newName,std::ios::out);
	fileIO.write();
	fileIO.close();
	
	return 0;
}

