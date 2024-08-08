#include "include.h"


//using type 1 for the polymer, 2 and 3 are for the HEAD and TAIL types respectively
//0 is a special type, it is stationary
#define POLY 1

int main(int argc, char **argv)
{
	if(argc!=18)
	{
		std::cerr << "Usage: " << argv[0] << " name seed nPoly nMonomers bondLength bondStrength bendStrength nLipids bilayerDensity UminPolyBilayer UmaxPolyBilayer UminPolyPoly UmaxPolyPoly UminPolySubstrate UmaxPolySubstrate polyGap temp" << std::endl;
		std::cerr << "name: Name of simulation (name.mpd)" << std::endl;
		std::cerr << "seed: Random number seed (any integer)." << std::endl;
		std::cerr << "nPoly: Number of polymers (0 to 1000)." << std::endl;
		std::cerr << "nMonomers: Number of monomers per polymer (3 to 100000)." << std::endl;
		std::cerr << "bondLength: Preferred bond length between polymer monomers (0.5 to 1.5, where 0.7 is same as lipids)." << std::endl;
		std::cerr << "bondStrength: Harmonic bond strength between polymer monomers (0 to 500)." << std::endl;
		std::cerr << "bendStrength: Harmonic strength of 3 body bending between polymer monomers (0 to 500)." << std::endl;
		std::cerr << "nLipids: Total number of lipids in flat membrane (10000 or higher)." << std::endl;
		std::cerr << "bilayerDensity: Lipids per membrane projected area (both leaflets, 3.11 for free bilayer)." << std::endl;
		std::cerr << "Umin/maxPolyBilayer: Attractive (Umin<0) and repulsive (Umax>0) components of Monomers and Heads (0 to -5 for Umin, and 100 to 200 for Umax)." << std::endl;
		std::cerr << "Umin/maxPolyPoly: Attractive (Umin<0) and repulsive (Umax>0) components of Monomers and Monomers (0 to -10 for Umin, and 100 to 200 for Umax)." << std::endl;
		std::cerr << "Umin/maxPolySubstrate: Attractive (Umin<0) and repulsive (Umax>0) components of Monomers and Substrate field at z=0 (0 to -2.0 for Umin, and 40 for Umax, substrate density is 5.88 particles per nm)." << std::endl;
		std::cerr << "polyGap: Gap between membrane and substrate that polymer inhabits." << std::endl;
		std::cerr << "temp: Temperature of simulation (3=~295K, 2.75=~256K)." << std::endl;
		return 0;
	}
	char *name=argv[1];
	std::stringstream cmdArg;
	
	for(int i=2;i<argc;i++)
		cmdArg << argv[i] << ' ';
	
	int seed, nMonomers, nLipids, nPoly;
	double bondLength, bondStrength, bendStrength, bilayerDensity, UminPolyBilayer, UmaxPolyBilayer, UminPolyPoly; 
	double UmaxPolyPoly, UminPolySubstrate, UmaxPolySubstrate, polyGap, temp;
	
	cmdArg >> seed >> nPoly >> nMonomers >> bondLength >> bondStrength >> bendStrength >> nLipids >> bilayerDensity >> UminPolyBilayer >> UmaxPolyBilayer >> UminPolyPoly >> UmaxPolyPoly >> UminPolySubstrate >> UmaxPolySubstrate >> polyGap >> temp;
	
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
	System.setFinalTime(10000);
	System.setDeltaT(0.02);
	System.setStoreInterval(1000);
	System.setMeasureInterval(100);
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
	
	//Polymer monomers
	Umin[POLY+POLY*System.readNTypes()]=UminPolyPoly;
	Umax[POLY+POLY*System.readNTypes()]=UmaxPolyPoly;
	
	Umin[POLY+HEAD*System.readNTypes()]=UminPolyBilayer;
	Umax[POLY+HEAD*System.readNTypes()]=UmaxPolyBilayer;
	
	Umin[HEAD+POLY*System.readNTypes()]=UminPolyBilayer;
	Umax[HEAD+POLY*System.readNTypes()]=UmaxPolyBilayer;
	
	//Tail monomers
	Umin[TAIL+TAIL*System.readNTypes()]=-6.0;
	Umax[TAIL+TAIL*System.readNTypes()]=200;
	
	//two body constants
	std::vector<double> twoBodyFconst=mpd::laradjiRevaleeFC(Umax,Umin,CUTOFF,RMIN);
	for(auto &v:twoBodyFconst)
		System.addTwoBodyFconst(v);
	
	std::vector<double> twoBodyUconst=mpd::laradjiRevaleePC(Umax,Umin,CUTOFF,RMIN);
	for(auto &v:twoBodyUconst)
		System.addTwoBodyUconst(v);
	
	double area=nLipids*bilayerDensity;
	threeVector<double>size;
	size.x=0.0;
	size.y=0.0;
	size.z=10.0*polyGap;
	System.setSize(size);
	
	threeVector<double> pos;
	pos.x=0;
	pos.y=0;
	pos.z=polyGap+0.7;
	
	//Bilayer constants
	std::vector<double> constants(4);
	constants[0]=0.7;//bondLength
	constants[1]=100;//kBond
	constants[2]=1;//-cos(theta_0)
	constants[3]=100;//kBend
	int particlesPerLipid=3;
	
	//simple bilayer, 1.0=square aspect ratio
	bilayer<double>(System, nLipids, particlesPerLipid, pos, constants[0], bilayerDensity, constants.data(), constants.size(), 1.0);
	
	size=System.readSize();
	
	
	//For polymer random velocities
	MTRand randNum(System.readSeed());
	
	//placing polymers in zigzag within gap
	std::vector< position<double> > poly(1);
	poly[0].type=POLY;
	poly[0].x=0;
	poly[0].y=0;
	poly[0].z=bondLength;
	
	int polyOffset=System.readNParticles();	
	
	//direction vector, becomes negative if at top of gap or edge of system
	threeVector<double> direction;
	direction.x=1.0;
	direction.y=1.0;
	direction.z=1.0;
	
	for(int iPoly=0;iPoly<nPoly;iPoly++)
	{
		for(int i=iPoly==0?1:0;i<nMonomers;i++)
		{
			threeVector<bool> go=false;
			go.z=true;
			position<double> next=poly.back();
			if(next.z+direction.z*bondLength>=polyGap || next.z+direction.z*bondLength<=0)
			{
				direction.z*=-1.0;
				go.z=false;
				go.x=true;
			}
			if(next.x+direction.x*bondLength>=size.x || next.x+direction.x*bondLength<=0)
			{
				direction.x*=-1.0;
				go.x=false;
				go.y=true;
			}
			if(next.y+direction.y*bondLength>=size.y || next.y+direction.y*bondLength<=0)
			{
				std::cerr << "No more path available!" << std::endl;
			}
			
			if(go.x)
				next.x+=direction.x*bondLength;
			if(go.y)
				next.y+=direction.y*bondLength;
			if(go.z)
				next.z+=direction.z*bondLength;
			poly.push_back(next);
		}
	}
	
	//Adding polymer particles to system
	for(auto p:poly)
	{
		threeVector<double> v,a;
		double Vrms=sqrt(3.0*System.readInitialTemp());
		
		//velocity has random direction and fixed magnitude
		double theta=M_PI*randNum.rand53();
		double phi=M_PI*2*randNum.rand53();
		v.x=Vrms*std::cos(phi)*std::sin(theta);
		v.y=Vrms*std::sin(phi)*std::sin(theta);
		v.z=Vrms*std::cos(theta);
		System.addParticle(p,v,a);
	}
	
	//Adding new polymer molecule to system
	fourVector<int> bond;
	molecule<double,fourVector<int> > chainMol;
	
	chainMol.setType(CHAIN);
	
	//constants for CHAIN
	chainMol.addConstant(bondLength);
	chainMol.addConstant(bondStrength);
	chainMol.addConstant(1); //-cos(theta_0)
	chainMol.addConstant(bendStrength);
	
	//bond for a chain type
	bond.s[START]=polyOffset;
	bond.s[NCHAINS]=nPoly;
	bond.s[CHAINLENGTH]=nMonomers;
	chainMol.addBond(bond);
	System.addMolecule(chainMol);
	
	
	
	molecule< double, fourVector<int> > boundaryZ;
	boundaryZ.setType(FLOATING_BASE);
	for(int i=0;i<System.readNParticles();i++)
	{
		fourVector<int> buf;
		buf.s[0]=i;
		boundaryZ.addBond(buf);
	}
	
	//Floating boundary substrate
	double fbDensity=5.88;
	std::vector<double> UminFB, UmaxFB;
	for(int i=0;i<System.readNTypes();i++)
	{
		UminFB.push_back(0);
		UmaxFB.push_back(40);
	}
	//Everything is repulsive to this substrate except the polymer monomers
	UminFB[POLY]=UminPolySubstrate;
	UmaxFB[POLY]=UmaxPolySubstrate;
	
	std::vector<double> fbConstants=mpd::laradjiPoursoroushFC(UmaxFB,UminFB,CUTOFF,RMIN,fbDensity);
	for(auto &fbC:fbConstants)
		boundaryZ.addConstant(fbC);//for each substrate constant
	
	if(boundaryZ.readNBond()>0)
		System.addMolecule(boundaryZ);
	
	size=System.readSize();
	size.z*=2.0;
	System.setSize(size);
	
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
	
	return 0;
}


