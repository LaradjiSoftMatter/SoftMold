#include "include.h"
#include <random>


//using type 3 for the brush base monomer to create a floating substrate
#define BB 3
//using type 0 for fixed brush base monomers
//#define BB 0
#define BM 1
#define NANOTYPE 2

//aspectRatio is x/y
template <typename T>
void variableBrush(Blob<T> &System, std::vector<int> polyLengths, threeVector<T> pos, T bondLength, T arealDensity, std::vector<T> constants, T aspectRatio, int bottomType, int chainType, T k, T temperature)
{
	//Error checking
	if(constants.size()!=2 && constants.size()!=4)
	{
		std::cerr << "Error: nConstants is not 2 or 4![bool brush()]\n";
		throw 0;
	}
	
	//Count the bonds
	int nParticles=0;
	for(auto len:polyLengths) nParticles+=len;
	
	//Count the bonds
	int nBonds=0;
	for(auto len:polyLengths) nBonds+=len-1;
	
	//Count the bends
	int nBends=0;
	if(constants.size()==4)
		for(auto len:polyLengths) nBends+=len-2;
	
	//Initialize molecule bonds
	molecule<T,fourVector<int> > mBond;
	
	mBond.setType(BOND);
	mBond.allocBonds(nBonds);
	
	for(int i=0;i<3;i++) mBond.addConstant(constants[i]);
	
	//Initialize molecule bends
	molecule<T,fourVector<int> > mBend;
	
	if(constants.size()==4)
	{
		mBend.setType(BEND);
		mBend.allocBonds(nBends);
		
		for(int i=0;i<3;i++) mBend.addConstant(constants[i+2]);
	}
	
	//Make sure we have enough room for the polymers
	int maxLength=0;
	for(auto len:polyLengths) maxLength=std::max(len,maxLength);
	T maxZHeight=maxLength*bondLength+1.0+pos.z;
	
	//adjust size of system if it is out of bounds
	threeVector<T> s;
	
	s.x=sqrt((T)polyLengths.size()/arealDensity)*sqrt(aspectRatio);
	s.y=sqrt((T)polyLengths.size()/arealDensity)/sqrt(aspectRatio);
	s.z=pos.z+maxZHeight;
	std::cerr << "Prefered system size is (brush): " << s.x << '\t' << s.y << '\t' << s.z << '\n';
	
	threeVector<T> size=System.readSize();
	
	std::cerr << "Current system size is (brush): " << size.x << '\t' << size.y << '\t' << size.z << '\n';
	
	//Adjust number of brushes to match grafting density
	size.x=(size.x<s.x)?s.x:size.x;
	size.y=(size.y<s.y)?s.y:size.y;
	size.z=(size.z<s.z)?s.z:size.z;
	std::cerr << "Adjusted system size is (brush): " << size.x << '\t' << size.y << '\t' << size.z << '\n';	
	System.setSize(size);
	
	//Using this random number generator for the velocities
	MTRand randNum(System.readSeed());
	
	//Velocities
	T Vrms=sqrt(3.0*System.readInitialTemp());
	
	//This is important to reduce excessive reallocations
	System.allocParticle(System.readNParticles()+nParticles);
	
	//lattice extent
	threeVector<int> latticeSize;
	//lattice components
	threeVector<T> latticeLength;
	
	latticeSize.x=static_cast<int>(sqrt((arealDensity*size.x*size.x)/aspectRatio))+1;
	latticeSize.y=static_cast<int>(sqrt((arealDensity*size.y*size.y)*aspectRatio))+1;
	
	latticeLength.x=size.x/static_cast<T>(latticeSize.x);
	latticeLength.y=size.y/static_cast<T>(latticeSize.y);
	
	std::cerr << "Adjusted density is (brush): " << static_cast<T>(polyLengths.size())/(size.x*size.y) << '=';
	std::cerr << polyLengths.size() << "/(" << size.x << '*' << size.y << ')' << std::endl;
	
	//initialize substrate positions in lattice
	std::vector<threeVector<T>> substrate;
	for(T x=0;x<latticeSize.x;x++)
	{
		for(T y=0;y<latticeSize.y;y++)
		{
			if(substrate.size()<polyLengths.size())
			{
				threeVector<T> temp;
				temp.x=x*latticeLength.x;
				temp.y=y*latticeLength.y;
				temp.z=pos.z;
				substrate.push_back(temp);
			}
		}
	}
	
	T rSub=(latticeLength.x<latticeLength.y)?latticeLength.y:latticeLength.x;
	//randomly space substrate using metropolis hasting monte carlo with k strength
	substrate=relaxSubstrate(substrate,rSub*1.0,size,1000,System.readSeed(),k, temperature);
	
	//set up polymers with substrate
	for(int i=0;i<substrate.size();i++)
	{
		T theta,phi;
		
		threeVector<T> a=0;
		for(int j=0;j<polyLengths[i];j++)
		{
			//position
			position<T> p;
			p.x=substrate[i].x;
			p.y=substrate[i].y;
			p.z=substrate[i].z+bondLength*j;
			p.type=(j==0)?bottomType:chainType;
			
			//random velocity vector direction
			threeVector<T> v;
			theta=M_PI*randNum.rand53();
			phi=M_PI*2*randNum.rand53();
			v.x=Vrms*cos(phi)*sin(theta);
			v.y=Vrms*sin(phi)*sin(theta);
			v.z=Vrms*cos(theta);
			
			//save its state
			System.addParticle(p,v,a);
			
			//if it is not the last monomer
			if(j!=polyLengths[i]-1)
			{
				fourVector<int> bond;
				bond.s[0]=System.readNParticles()-1;
				bond.s[1]=System.readNParticles();
				mBond.addBond(bond);
			}
			//if it is not the last bend segment
			if(j!=polyLengths[i]-2 && constants.size()==4)
			{
				fourVector<int> bend;
				bend.s[0]=System.readNParticles()-2;
				bend.s[1]=System.readNParticles()-1;
				bend.s[1]=System.readNParticles();
				mBend.addBond(bend);
			}
		}
	}
	//Add the molecules to the system blob.
	System.addMolecule(mBond);
	if(constants.size()==4)
		System.addMolecule(mBend);
}

template <typename T>
void distributeNanoparticles(Blob<T> &System, std::vector<T> radii, std::vector<int> nTessellations, 
			     threeVector<T> lowerBound, threeVector<T> upperBound)
{
	// need this for the random placement
	auto dBound=upperBound-lowerBound;
	
	// Position of nanoparticles
	std::vector< threeVector<T> > nanoPos;
	
	// Using the system seed for the random number generator
	MTRand randNum(System.readSeed());
	
	for(int i=0;i<radii.size();i++)
	{
		std::cerr << "Placing nanoparticle " << i << "!" << std::endl;
		threeVector<T> toPlace;
		bool overlap;
		int attempts=0;
		do
		{
			overlap=false;
			toPlace.x=lowerBound.x+(dBound.x-2.0*radii[i])*randNum.rand53()+radii[i];
			toPlace.y=lowerBound.y+(dBound.y-2.0*radii[i])*randNum.rand53()+radii[i];
			toPlace.z=lowerBound.z+(dBound.z-2.0*radii[i])*randNum.rand53()+radii[i];
			for(int j=0;j<nanoPos.size();j++)
			{
				//difference vector
				threeVector<T> d;
				d.x=nanoPos[j].x-toPlace.x;
				d.y=nanoPos[j].y-toPlace.y;
				d.z=nanoPos[j].z-toPlace.z;
				
				//minimum image
				d.x-=(d.x>System.readSize().x/2.0)?System.readSize().x:0;
				d.x+=(d.x<-System.readSize().x/2.0)?System.readSize().x:0;
				
				d.y-=(d.y>System.readSize().y/2.0)?System.readSize().y:0;
				d.y+=(d.y<-System.readSize().y/2.0)?System.readSize().y:0;
				
				d.z-=(d.z>System.readSize().z/2.0)?System.readSize().z:0;
				d.z+=(d.z<-System.readSize().z/2.0)?System.readSize().z:0;
				
				if(sqrt(d.x*d.x+d.y*d.y+d.z*d.z)<radii[i]+radii[j]+0.5)
					overlap=true;
			}
			attempts++;
			if(attempts>10000)
			{
				std::cerr << "10,000 nanoparticle placement attempts exceeded. Please check system dimensions!" << std::endl;
				throw 0;
			}
		} while(overlap);
		std::cerr << i << '\t' << toPlace.x << '\t' << toPlace.y << '\t' << toPlace.z << std::endl;
		nanoPos.push_back(toPlace);
		
	}
	
	//To get the nanoparticle core (in order to expulse any polymers) we can add the 
	//NANOCORE type
	molecule<T,fourVector<int> > mCore;
	mCore.setType(NANOCORE);
	mCore.allocBonds(radii.size());
	
	//Note: radii.size()==nanoPos.size() now
	//The constants 1200, 45, 250, 50, and 0 are taken from Yu for 0.02 deltaT capable simulation
	for(int i=0;i<radii.size();i++)
	{
		//we need this to modify the types to be all the same
		int nanoOffset=System.readNParticles();
		
		//Janus particle from Yu Zhu
		janus<T>(System,1200,45,250,50,0,nanoPos[i],radii[i],nTessellations[i]);
		
		//prototype for NANOCORE force/potential constants
		//std::vector<T> laradjiSpanglerFCrow(T rc, T rm, T rad, T density, T Umax, T Umin)
		
		//This is a repulsive NANOCORE to all other particles in the system
		//with Umax=30 and Umin=0. The size is slightly smaller than the Janus shell
		auto C=mpd::laradjiSpanglerFCrow(CUTOFF,RMIN,radii[i]-1.3,5.88,30.0,0.0);
		fourVector<int> bond;
		//The last particle in the Janus nanoparticle is the core, 
		//see note at bottom of context
		bond.x=System.readNParticles()-1;
		mCore.addBond(bond);
		for(auto constant:C)
			mCore.addConstant(constant);
		
		//Alternative NP, problem with this one is that the interaction doesn't work with
		//different sizes. You might be able to change that in root/include/MD.h
		//continuumSphere<T>(System, &(nanoPos[0]), nanoPos.size(), radii[i], &C[0], 
		//			System.readNTypes()*System.readNTypes()*nBEADCONST, NANOTYPE);
		
		//Change all the types to the same type, or use a unique type
		for(int j=nanoOffset;j<System.readNParticles();j++)
			System.getPositions()[j].type=NANOTYPE;
		
		// Note: Yu's particle is a tessellated dodecahedron with a center
		// anchor particle. The layout is the shell first and then the anchor.
		// Another note: change the interaction with the general Umin/Umax parameters.
	}
	
	//This adds the core to the janus particle
	System.addMolecule(mCore);
}

//The meat of the generation
int main(int argc, char **argv)
{
	/// Reading in parameters for this simulation
	if(argc!=14)
	{
		//14, brush + NP
		std::cout << "Usage: " << argv[0] << " name seed nBrushes arealDensity UminBrushNano UmaxBrushNano nNanoparticles nanoRadius UminBrushBrush UmaxBrushBrush nMonomers bondLength temperature" << std::endl;
		return 0;
	}
	char *name=argv[1];
	std::stringstream cmdArg;
	
	for(int i=2;i<argc;i++)
		cmdArg << argv[i] << ' ';
	
	int seed,nBrushes,nMonomers,nNanoparticles;
	double arealDensity, UminBrushNano, UmaxBrushNano, UminBrushBrush,UmaxBrushBrush,bondLength,temp,nanoRadius;
	
	cmdArg >> seed >> nBrushes >> arealDensity >> UminBrushNano >> UmaxBrushNano >> nNanoparticles >> nanoRadius >> UminBrushBrush >> UmaxBrushBrush >> nMonomers >> bondLength >> temp;
	
	/// Setup for the initial system state
	Blob<double> System;
	
	System.setGamma(1.0); // Damping for Langevin thermostat, smaller=more inertia larger=less inertia
	System.setNTypes(4);  // Number of particle types
	System.setSeed(seed); // Random number seed, this should be changed for any simulation run
	
	// Keep this here, but wrapping is always on whether true or false
	threeVector<bool> wrap; 
	wrap.x=true;
	wrap.y=true;
	wrap.z=true;
	System.setPeriodic(wrap); // Periodic boundary flag, this doesn't work yet, it is always on
	System.setCutoff(2.0);    // Cutoff distance for non-bonded force, for neighbors
	
	/// Note on units:
	/// 1. For lipid systems, 1 "tau" is approximately ~19 ns
	/// 2. Temperature for lipid systems is ~2.7 and less for gel and fluid for up to ~3.4,
	///    and this would be appropriate for ~310 Kelvin DPPC gel-fluid phase transition
	/// 3. Those temperatures reflect lipids, and the conversion might be different for
	///    other polymers
	System.setInitialTime(0);      // Initial time in tau
	System.setFinalTime(10000);    // Final time in tau
	System.setDeltaT(0.02);        // Time step in tau intervals
	System.setStoreInterval(1000); // Storing the system in tau intervals
	System.setMeasureInterval(100);// Measuring kinetic/potential energy in tau intervals
	//System.setDeltaLXY(0.01);    // "Barostat", Resize the system by a maximum of this value
	System.setInitialTemp(temp);   // Initial temperature
	System.setFinalTemp(temp);     // Final Temperature
	
	// We're just going to zero this for now
	threeVector<double>size=0;
	System.setSize(size);
	
	// offset for bottom of polymers, along the z axis
	threeVector<double> pos=0;
	pos.z=1.0;
	
	/// Set some polymer lengths using poisson distribution
	// I'm not sure what distribution to use, but check out the following paper:
	// http://dx.doi.org/10.3390/polym10080887
	std::vector<int> polyLengths;
	std::mt19937 randM(System.readSeed());
	std::poisson_distribution<int> randDist(nMonomers-3.0);
	for(int i=0;i<nBrushes;i++) polyLengths.emplace_back(randDist(randM)+3);
	// Note: Binomial and Geometric distributions are cool.
	
	/// Constants for brush polymers, add 2 more constants for bending
	std::vector<double> constants;
	constants.push_back(bondLength); // Bond length of monomer-monomer harmonic potential
	constants.push_back(100);        // Bond strength, k from k*(d-a_0)^2, for harmonic potential
	//constants.push_back(1);          // Bend angle in -cos(theta_0) units
	//constants.push_back(100);        // Bend strength, k_be from k_be*(cos(theta)-cos(theta_0)^2, potential
	
	// Simple fixed brush from other simulation, here for reference
	//brush<double>(System, nBrushes, nMonomers, pos, bondLength, 
	//	      arealDensity, constants, 4, 1.0, BB, BM,100,3.0);
	
	// Variable brush with random lattice placed by metropolis algorithm
	variableBrush<double>(System, polyLengths, pos, bondLength, 
		      arealDensity, constants, 1.0, BB, BM, 100.0, 3.0);
	
	// Get the current size because the above function can change it.
	size=System.readSize();
	
	/// We need some boundaries for placing the nanoparticles
	threeVector<double> lowerBound, upperBound;
	// These are wrapped, but I recommend treating them like hard bounds
	lowerBound.x=0.0;
	lowerBound.y=0.0;
	lowerBound.z=size.z;
	// Note: upper bounds are not hard bounds, 
	// when System.addParticle is called the system size is enforced to
	// cause particle positions to wrap near the 0 bound.
	upperBound.x=size.x;
	upperBound.y=size.y;
	// Note: This might need a nanoparticle volume density component.
	// Another note: Because the base monomers are near the 0 Z component, 
	// I would treat this like a heard bound.
	upperBound.z=size.z*2.0;
	
	// Change the Z component of size to accomodate the nanoparticles.
	size.z=upperBound.z;
	System.setSize(size);
	
	//Need some parameters for nanoparticles' radius and tesselations
	std::vector<double> radii;
	std::vector<int> nTessellations;
	for(int i=0;i<nNanoparticles;i++)
	{
		radii.emplace_back(nanoRadius+0.2*(i-0.5*nNanoparticles));// Let's fix this for now.
		nTessellations.emplace_back(2);// Any positive integer including zero.
		// Note: My suggestion is to do one of several things (I like 1:b):
		// 1: Random distributions of nanoparticle radii around a mean with:
		//    a. Fixed tesselations and allow density adjusted Umin/Umax
		//       (may create gaps for polymers into nanoparticle)
		//    b. Integer stepped random radii with tessellation density
		//       (may limit nanoparticle sizes, but for small nanoparticles
		//       this might be considered realistic since "atoms" are treated
		//       classically. I'm not adding nuclear physics to this simulation)
		//    c. Using continuum NPs with different sizes
		//       (easiest to use here, hardest to implement in MD.h and System.h
		//       because we need to solve for overlap in a 7 degree polynomial)
		// 2: Same radius, different tessellations of surface density
	}
	
	// Random rectogonal distribution of nanoparticles placed above brush polymers
	distributeNanoparticles<double>(System, radii, nTessellations, lowerBound, upperBound);
	
	// Floating boundary for bottom monomers to adhear to.
	molecule< double, fourVector<int> > boundaryZ;
	boundaryZ.setType(FLOATING_BASE);
	
	// Add all the particles' indices to this boundary
	for(int i=0;i<System.readNParticles();i++)
	{
		fourVector<int> buf;
		buf.s[0]=i;
		boundaryZ.addBond(buf);
	}
	
	// This is the same relative "strength" as the non-bonded force used
	// in the rest of the simulation. This force acts along the Z boundary
	// rm/rc or RMIN/CUTOFF above Z==0
	std::vector<double> UminFB, UmaxFB;
	
	//default is repulsive
	for(int i=0;i<System.readNTypes();i++)
	{
		UminFB.push_back(0);
		UmaxFB.push_back(200);
	}
	
	//Bottom monomer exceptions
	UminFB[BB]=-4; //Adhesive part, -2 is relatively strong adhesion
	UmaxFB[BB]=100; //Repulsive part
	
	double fbDensity=5.88; // particles per r_min area
	std::vector<double> fbConstants=mpd::laradjiPoursoroushFC(UmaxFB,UminFB,CUTOFF,RMIN,fbDensity);
	for(auto &fbC:fbConstants)
		boundaryZ.addConstant(fbC);//for each substrate constant
	
	if(boundaryZ.readNBond()>0)
		System.addMolecule(boundaryZ);
	
	//Harder boundary
	molecule< double, fourVector<int> > boundaryZ2;
	boundaryZ2.setType(BOUNDARY);
	for(int i=0;i<System.readNParticles();i++)
	{
		fourVector<int> buf;
		buf.s[0]=i;
		boundaryZ2.addBond(buf);
	}
	boundaryZ2.addConstant(2);//along Z, direction
	boundaryZ2.addConstant(0.0);//at 0=z position, z_0
	boundaryZ2.addConstant(4);//new:Doesn't matter old:half width of z, cutoff
	boundaryZ2.addConstant(200);//k, epsilon
	
	if(boundaryZ2.readNBond()>0)
		System.addMolecule(boundaryZ2);
	
	// More height? Go ahead and add more if you want. It can be changed anytime in mpd file.
	// Should only slightly affect the pressure of nanoparticles to the polymer brush.
	size=System.readSize();
	size.z*=2.0;
	System.setSize(size);
	
	/// Initialize non-bonded constants (force and potential constants)
	// Note: calls to Janus will set the number of types to 8, uncomment the next line
	// if you want to do something more complicated
	//System.setNTypes(4);  // Number of particle types
	std::vector<double> Umax(System.readNTypes()*System.readNTypes());
	std::vector<double> Umin(System.readNTypes()*System.readNTypes());
	
	// Could be one loop, but it's good for an example
	// of how to set constant "symmetric" matrix
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
	
	// Brush monomers to brush monomers
	Umin[BM+BM*System.readNTypes()]=UminBrushBrush;
	Umax[BM+BM*System.readNTypes()]=UmaxBrushBrush;
	
	// Brush monomers to Nanoparticles
	Umin[BM+NANOTYPE*System.readNTypes()]=UminBrushNano;
	Umax[BM+NANOTYPE*System.readNTypes()]=UmaxBrushNano;
	
	//Mirror of the previous to maintain symmetry
	Umin[NANOTYPE+BM*System.readNTypes()]=Umin[BM+NANOTYPE*System.readNTypes()];
	Umax[NANOTYPE+BM*System.readNTypes()]=Umax[BM+NANOTYPE*System.readNTypes()];
	
	// Convert two body force constants into a format the simulation uses
	std::vector<double> twoBodyFconst=mpd::laradjiRevaleeFC(Umax,Umin,CUTOFF,RMIN);
	for(auto &v:twoBodyFconst)
		System.addTwoBodyFconst(v);//Put in Blob
	
	// Convert two body potential constants into a format the simulation uses
	std::vector<double> twoBodyUconst=mpd::laradjiRevaleePC(Umax,Umin,CUTOFF,RMIN);
	for(auto &v:twoBodyUconst)
		System.addTwoBodyUconst(v);//Put in Blob
	
	std::cout << "Storing configuration...";
	
	//Write configuration
	Script<double,Blob <double> > output(name,std::ios::out,&System);
	output.write();
	
	//This stores the xyz file to check new configuration
	position<double> *p=System.getPositions();
	int nParticles=System.readNParticles();
	
	//to store an xyz file with the current configuration, not needed, 
	//but makes it easy to visualize the initial configuration.
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

