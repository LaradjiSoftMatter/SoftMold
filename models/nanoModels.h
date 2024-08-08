/**
 * \brief Various model geometries associated with nano-particles.
 * Places nano-particles into your systems. 
 * Many options for connectivity as well. Assumes your system blob has some functions to utilize 
 * with molecules, particles, and temperature.
 */

#include "../include/algorithms/functions.h"
#include "../include/algorithms/cellAlgorithm.h"

#ifndef LATTICE_THRESHOLD
#define LATTICE_THRESHOLD 0.01
#endif

/** \brief Searches for bonds on a lattice.
 * Generic search for bonds. It assumes we are searching among neighboring particles
 * for particles that fit approximately on the lattice (bonded) scaled by latticeLength.
 * If you perform 'searchForBonds< T > object(...); molecule< T, fourVector< int > > m=pairSearch(object)();',
 * then it will create a new molecule m based on our lattice parameters, to manipulate.
 */
template <typename T>
class searchForBonds {
	public:
		//constructor
		searchForBonds(position<T> *particles, int nParticles, position<T> *bonded, int nBonded, T *constants, int nConstants, T latticeLength, int offset)
		{
			m.setType(BOND);
			for(int i=0;i<nConstants;i++)
				m.addConstant(constants[i]);
			
			m.allocBonds(nBonded*nParticles);
			
			p=particles;
			nP=nParticles;
			ll=latticeLength;
			b=bonded;
			nB=nBonded;
			this->offset=offset;
			periodic=false;
		};
		
		//with periodic boundaries
		searchForBonds(position<T> *particles, int nParticles, position<T> *bonded, int nBonded, T *constants, int nConstants, T latticeLength, int offset, threeVector<T> size)
		{
			m.setType(BOND);
			for(int i=0;i<nConstants;i++)
				m.addConstant(constants[i]);
			
			m.allocBonds(nBonded*nParticles);
			
			p=particles;
			nP=nParticles;
			ll=latticeLength;
			b=bonded;
			nB=nBonded;
			this->offset=offset;
			periodic=true;
			s=size;
		}
		
		//destructor
		~searchForBonds(){
			
		};
		
		//output the molecule
		molecule< T, fourVector<int> > &getMolecule()
		{
			return m;
		};
		
		//check if any neihgbor fits the constraints
		void operator() (int &i, int &j)
		{
			for(int k=0;k<nB;k++)
			{
				threeVector<T> d;
				
				d.x=p[i].x-p[j].x+b[k].x*ll;
				d.y=p[i].y-p[j].y+b[k].y*ll;
				d.z=p[i].z-p[j].z+b[k].z*ll;
				if(periodic)
				{
					d.x-=(d.x>s.x/2.0)?s.x:0.0;
					d.x+=(d.x<-s.x/2.0)?s.x:0.0;
					d.y-=(d.y>s.y/2.0)?s.y:0.0;
					d.y+=(d.y<-s.y/2.0)?s.y:0.0;
					d.z-=(d.z>s.z/2.0)?s.z:0.0;
					d.z+=(d.z<-s.z/2.0)?s.z:0.0;
				}
				if(d.x*d.x+d.y*d.y+d.z*d.z<LATTICE_THRESHOLD*LATTICE_THRESHOLD)
				{
					fourVector<int> bond;
					bond.s[0]=i+offset;
					bond.s[1]=j+offset;
					m.addBond(bond);
				}
			}
		};
		
	private:
		molecule< T, fourVector<int> > m;
		
		position<T> *p;
		int nP;
		position<T> *b;
		int nB;
		T ll;
		int offset;
		bool periodic;
		threeVector<T> s;
};

/** \brief Adds a hollow nano-shell to the system.
 * The nanoshell is a continuum based model that utilizes a single point and a specialized potential.
 * The potential is derived from the integration of charge points scattered over a surface
 * with a density of sigma. This particular function simply places a single particle into the system.
 */
template <typename T>
int continuumSphere(Blob<T> &System, threeVector<T> pos, T radius, T *constants, int nConstants, int type)
{
	//Error checking
	if(constants==NULL)
	{
		std::cerr << "Error(continuumSphere): Bad allocation of constants!\n";
		return 0;
	}
	
	//Random number generator
	MTRand randNum(System.readSeed());
	
	//Velocity RMS
	T Vrms=sqrt(3.0*System.readInitialTemp());
	
	position<T> p;
	p.x=pos.x;
	p.y=pos.y;
	p.z=pos.z;
	p.type=type;
	T theta,phi;
	
	threeVector<T> v,a;
	a=0;
	//velocity
	theta=M_PI*randNum.rand53();
	phi=M_PI*2*randNum.rand53();
	v.x=Vrms*cos(phi)*sin(theta);
	v.y=Vrms*sin(phi)*sin(theta);
	v.z=Vrms*cos(theta);
	
	System.addParticle(p,v,a);
	
	molecule< T, fourVector<int> > m;
	m.setType(BEAD);
	//this should have enough to fill
	for(int i=0;i<nConstants;i++)
		m.addConstant(constants[i]);
	//m.addConstant(radius);
	//m.addConstant(type);
	fourVector<int> bond;
	bond.s[0]=System.readNParticles()-1;
	m.addBond(bond);
	
	std::cerr << m.readNBond() << std::endl;
	//std::cin.get();
	System.addMolecule(m);
	return 1;
}

/** \brief Adds a hollow nano-shell to the system.
 * The nanoshell is a continuum based model that utilizes a single point and a specialized potential.
 * The potential is derived from the integration of charge points scattered over a surface
 * with a density of sigma. This particular function adds a field to an index into the system.
 */
template <typename T>
int continuumSphere(Blob<T> &System, int posIndex, T radius, T *constants, int nConstants)
{
	//Error checking
	if(constants==NULL)
	{
		std::cerr << "Error(continuumSphere): Bad allocation of constants!\n";
		return 0;
	}
	
	molecule< T, fourVector<int> > m;
	m.setType(BEAD);
	//this should have enough to fill
	for(int i=0;i<nConstants;i++)
		m.addConstant(constants[i]);
	//m.addConstant(radius);
	//m.addConstant(type);
	fourVector<int> bond;
	bond.s[0]=posIndex;
	m.addBond(bond);
	
	threeVector<T> v=0;
	System.getVelocities()[posIndex]=v;

	std::cerr << m.readNBond() << std::endl;
	//std::cin.get();
	System.addMolecule(m);
	return 1;
}

/** \brief Adds a list of hollow nano-shells to the system.
 * The nanoshell is a continuum based model that utilizes a single point and a specialized potential.
 * The potential is derived from the integration of charge points scattered over a surface
 * with a density of sigma. This particular function simply places a group of particles into the system.
 * The list is a pointer to a group of positions. One could cast "std::vector< threeVector< T> > values"
 * with "&(values[0])" and "values.size()" for nSpheres.
 */
template <typename T>
int continuumSphere(Blob<T> &System, threeVector<T> *pos, int nSpheres, T radius, T *constants, int nConstants, int type)
{
	//Error checking
	if(pos==NULL)
	{
		std::cerr << "Error(continuumSphere): Bad allocation of positions!\n";
		return 0;
	}
	if(constants==NULL)
	{
		std::cerr << "Error(continuumSphere): Bad allocation of constants!\n";
		return 0;
	}
	//Random number generator
	MTRand randNum(System.readSeed());
	
	molecule< T, fourVector<int> > m;
	m.setType(BEAD);
	//this should have enough to fill
	for(int i=0;i<nConstants;i++)
		m.addConstant(constants[i]);
	
	//Velocity RMS
	T Vrms=sqrt(6.0*System.readInitialTemp()/(4.0*M_PI*radius*radius));
	for(int i=0;i<nSpheres;i++)
	{
		position<T> p;
		p.x=pos[i].x;
		p.y=pos[i].y;
		p.z=pos[i].z;
		p.type=type;
		T theta,phi;
		
		threeVector<T> v,a;
		a=0;
		//velocity
		theta=M_PI*randNum.rand53();
		phi=M_PI*2*randNum.rand53();
		v.x=Vrms*cos(phi)*sin(theta);
		v.y=Vrms*sin(phi)*sin(theta);
		v.z=Vrms*cos(theta);
		
		System.addParticle(p,v,a);
		
		//m.addConstant(radius);
		//m.addConstant(type);
		fourVector<int> bond;
		bond.s[0]=System.readNParticles()-1;
		m.addBond(bond);
		
		std::cerr << m.readNBond() << std::endl;
	}
	//std::cin.get();
	System.addMolecule(m);
	return 1;
}

/** \brief Adds a spherical nano-particle to the system. 
 * Spherical nano-particle is generated via a unit cell that is placed over and over until it fills the volume.
 * takes obvious parameters like latticeLength and nano-particle radius. pos is the position of the nano-particle in the system.
 * unitCell is the actual unit cell. An FCC lattice would have 3 Face centered positions and a position at the origin.
 * Relative positions (for FCC) would be (0,0,0), (0.5,0.5,0), (0.5,0,0.5), and (0,0.5,0.5). So nUnitCell for a FCC is 4.
 * bonded is the position of bonded particles. Typically, an FCC would have 6: (0.5,0.5,0), (0.5,0,0.5), (0,0.5,0.5),
 * (-0.5,-0.5,0), (-0.5,0,-0.5), and (0,-0.5,-0.5). 6 would be nBonded. The bond length would be 'sqrt(2.0*0.5*0.5)/2.0' for a
 * perfectly face centered lattice. There is normally an error in the above, but in this particular function, it check if those
 * bonds exist. If they don't exist (within LATTICE_THRESHOLD=0.01), then they are discarded, leaving only certain neighbors.
 */
template <typename T>
T nanoSphere(Blob<T> &System, threeVector<T> pos, T radius, T latticeLength, position<T> *unitCell, int nUnitCell, 
	     position<T> *bonded, int nBonded, T *constants, int nConstants)
{
	//I should really run through everything and replace pointers or something...
	//Error checking
	if(!(nConstants==4 || nConstants==2))
	{
		std::cerr << "Error(nanoSphere): nConstants is not 4 or 2!\n";
		return 0;
	}
	if(constants==NULL)
	{
		std::cerr << "Error(nanoSphere): Bad allocation of constants!\n";
		return 0;
	}
	if(radius<=0)
	{
		std::cerr << "Error(nanoSphere): No radius!\n";
		return 0;
	}
	if(nUnitCell<=0)
	{
		std::cerr << "Error(nanoSphere): nUnitCell<=0!\n";
		return 0;
	}
	if(unitCell==NULL)
	{
		std::cerr << "Error(nanoSphere): Bad allocation of unitCell!\n";
		return 0;
	}
	if(nBonded<=0)
	{
		std::cerr << "Error(nanoSphere): nBonded<=0!\n";
		return 0;
	}
	if(bonded==NULL)
	{
		std::cerr << "Error(nanoSphere): Bad allocation of bonded!\n";
		return 0;
	}
	std::cerr << "nConstants is " << nConstants << " in nanoSphere.\n";
	
	//Build and add nano-particle
	//index cells on a cube for faster searching
	fourVector<int> nLattice;
	nLattice.x=floor(radius*2.0/latticeLength)+1;
	nLattice.y=nLattice.x;
	nLattice.z=nLattice.z;
	nLattice.t=nLattice.x*nLattice.y*nLattice.z;
	
	//Random number generator
	MTRand randNum(System.readSeed());
	
	//Velocity RMS
	T Vrms=sqrt(3.0*System.readInitialTemp());
	
	int offset=System.readNParticles();
	
	//this is a guess
	int nNanoSphereParticles=nUnitCell*((4.0/3.0)*M_PI*(radius+latticeLength)*(radius+latticeLength)*(radius+latticeLength))
				/(latticeLength*latticeLength*latticeLength);
	
	//to avoid allocation slowdown
	System.allocParticle(nNanoSphereParticles+offset);
	
	//fill sphere
	for(threeVector<T> clp=0;clp.x<radius*2.0;clp.x+=latticeLength)
	{
		for(clp.y=0;clp.y<radius*2.0;clp.y+=latticeLength)
		{
			for(clp.z=0;clp.z<radius*2.0;clp.z+=latticeLength)
			{
				//set current lattice index to base particle index
				//this is here for reference
				//threeVector<int> cli;
				//cli.x=floor(clp.x/latticeLength);
				//cli.y=floor(clp.y/latticeLength);
				//cli.z=floor(clp.z/latticeLength);
				
				//set all subsequent particles to a position in the unit cell
				for(int i=0;i<nUnitCell;i++)
				{
					position<T> p=unitCell[i];//recieves type information in this step
					//set at current lattice position
					p.x*=latticeLength;
					p.y*=latticeLength;
					p.z*=latticeLength;
					p.x+=clp.x+pos.x-radius;
					p.y+=clp.y+pos.y-radius;
					p.z+=clp.z+pos.z-radius;
					
					T theta,phi;
					
					threeVector<T> v,a;
					a=0;
					//velocity
					theta=M_PI*randNum.rand53();
					phi=M_PI*2*randNum.rand53();
					v.x=Vrms*cos(phi)*sin(theta);
					v.y=Vrms*sin(phi)*sin(theta);
					v.z=Vrms*cos(theta);
					//save it
					//using point (radius,radius,radius) as center of mass
					threeVector<T> d;
					d.x=p.x-pos.x;
					d.y=p.y-pos.y;
					d.z=p.z-pos.z;
					
					//is the current lattice point out of range?
					if(d.x*d.x+d.y*d.y+d.z*d.z<radius*radius)
						System.addParticle(p,v,a);
				}
			}
		}
	}
	
	//set up a pair list
	Cell<double> nearbyPairs(System.getPositions()+offset, System.readNParticles()-offset, latticeLength*1.1, System.readSize());
	nearbyPairs.build();
	
	searchForBonds<double> newSearch(System.getPositions()+offset, System.readNParticles()-offset, bonded, 
					 nBonded, constants, nConstants, latticeLength, offset, System.readSize());
	nearbyPairs.twoWayCompareIndex(newSearch);
	
	molecule< T, fourVector<int> > m;
	m=newSearch.getMolecule();
	std::cerr << m.readNBond() << std::endl;
	//std::cin.get();
	System.addMolecule(m);
	//Compare the simplicity of the lines above to the commented block below.
	/*
	//set up the molecule (this will be extended to use more than one.
	molecule< T, fourVector<int> > m;
	
	m.setType(BOND);
	for(int i=0;i<nConstants;i++)
		m.addConstant(constants[i]);
	position<T> *p=System.getPositions();
	
	m.allocBonds(nBonded*(System.readNParticles()-offset));
	
	//search for nearby bonds using bonded as base
	for(threeVector<int> cli=0;cli.x<nLattice.x;cli.x++)
	{
		for(cli.y=0;cli.y<nLattice.y;cli.y++)
		{
			for(cli.z=0;cli.z<nLattice.z;cli.z++)
			{
				//if it exists
				if(cube[cli.x][cli.y][cli.z]>=0)
				{
					//check all the particles in the unit cell
					for(int i=0;i<nUnitCell;i++)
					{
						//for all chances to bond with a nearby particle
						for(int j=0;j<nBonded;j++)
						{
							//here is our first base index
							int baseIndexA=cube[cli.x][cli.y][cli.z]+i;
							threeVector<int> nearbyIndex;
							nearbyIndex.x=floor((p[baseIndexA].x/latticeLength)+bonded[j].x);
							nearbyIndex.y=floor((p[baseIndexA].y/latticeLength)+bonded[j].y);
							nearbyIndex.z=floor((p[baseIndexA].z/latticeLength)+bonded[j].z);
							
							//check if lattice point exists
							if(nearbyIndex.x>=0 && nearbyIndex.x<nLattice.x &&
								nearbyIndex.y>=0 && nearbyIndex.y<nLattice.y &&
								nearbyIndex.z>=0 && nearbyIndex.z<nLattice.z)
							//Make sure it isn't repeated
							//if(cli.
							{
								//make sure it isn't masked
								if(cube[nearbyIndex.x][nearbyIndex.y][nearbyIndex.z]>=0)
								{
									//go through the cell to locate one that matches the bond specifications
									for(int k=0;k<nUnitCell;k++)
									{
										int baseIndexB=cube[nearbyIndex.x][nearbyIndex.y][nearbyIndex.z]+k;
										//Make sure it isn't already bonded
										if(baseIndexB>baseIndexA)
										{
											threeVector<T> d;
											d.x=p[baseIndexB].x-p[baseIndexA].x-bonded[j].x*latticeLength;
											d.y=p[baseIndexB].y-p[baseIndexA].y-bonded[j].y*latticeLength;
											d.z=p[baseIndexB].z-p[baseIndexA].z-bonded[j].z*latticeLength;
											
											if(d.x*d.x+d.y*d.y+d.z*d.z<LATTICE_THRESHOLD*LATTICE_THRESHOLD)
											{
											//	std::cerr << d.x << '\t' << d.y << '\t' << d.z << '\n';
											//	std::cerr << baseIndexA << '\t' << baseIndexB << '\n';
											//	std::cin.get();
												
												fourVector<int> bond;
												bond.s[0]=baseIndexA;
												bond.s[1]=baseIndexB;
												m.addBond(bond);
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	
	std::cerr << m.readNBond() << '\n';
	
	System.addMolecule(m);
	*/
	return radius;
}


/** \brief Adds a spherical nano-particle to the system with a bifrucated surface.
 * Spherical nano-particle is generated via a unit cell that is placed over and over until it fills the volume.
 * takes obvious parameters like latticeLength and nano-particle radius. pos is the position of the nano-particle in the system.
 * unitCell is the actual unit cell. An FCC lattice would have 3 Face centered positions and a position at the origin.
 * Relative positions (for FCC) would be (0,0,0), (0.5,0.5,0), (0.5,0,0.5), and (0,0.5,0.5). So nUnitCell for a FCC is 4.
 * bonded is the position of bonded particles. Typically, an FCC would have 6: (0.5,0.5,0), (0.5,0,0.5), (0,0.5,0.5),
 * (-0.5,-0.5,0), (-0.5,0,-0.5), and (0,-0.5,-0.5). 6 would be nBonded. The bond length would be 'sqrt(2.0*0.5*0.5)/2.0' for a
 * perfectly face centered lattice. There is normally an error in the above, but in this particular function, it check if those
 * bonds exist. If they don't exist (within LATTICE_THRESHOLD=0.01), then they are discarded, leaving only certain neighbors.
 * planeHeight is the distance from the lowest z value in the nanoparticle. alterType is the alternate type above the planeHeight.
 */
template <typename T>
T janusSphere(Blob<T> &System, threeVector<T> pos, T radius, T latticeLength, position<T> *unitCell, int nUnitCell, 
	     position<T> *bonded, int nBonded, T *constants, int nConstants, T planeHeight, int alterType)//std::vector< threeVector<double> > &planes)
{
	//I should really run through everything and replace pointers or something...
	//Error checking
	if(!(nConstants==4 || nConstants==2))
	{
		std::cerr << "Error(nanoSphere): nConstants is not 4 or 2!\n";
		return 0;
	}
	if(constants==NULL)
	{
		std::cerr << "Error(nanoSphere): Bad allocation of constants!\n";
		return 0;
	}
	if(radius<=0)
	{
		std::cerr << "Error(nanoSphere): No radius!\n";
		return 0;
	}
	if(nUnitCell<=0)
	{
		std::cerr << "Error(nanoSphere): nUnitCell<=0!\n";
		return 0;
	}
	if(unitCell==NULL)
	{
		std::cerr << "Error(nanoSphere): Bad allocation of unitCell!\n";
		return 0;
	}
	if(nBonded<=0)
	{
		std::cerr << "Error(nanoSphere): nBonded<=0!\n";
		return 0;
	}
	if(bonded==NULL)
	{
		std::cerr << "Error(nanoSphere): Bad allocation of bonded!\n";
		return 0;
	}
	std::cerr << "nConstants is " << nConstants << " in nanoSphere.\n";
	
	//Build and add nano-particle
	//index cells on a cube for faster searching
	fourVector<int> nLattice;
	nLattice.x=floor(radius*2.0/latticeLength)+1;
	nLattice.y=nLattice.x;
	nLattice.z=nLattice.z;
	nLattice.t=nLattice.x*nLattice.y*nLattice.z;
	
	//Random number generator
	MTRand randNum(System.readSeed());
	
	//Velocity RMS
	T Vrms=sqrt(3.0*System.readInitialTemp());
	
	int offset=System.readNParticles();
	
	//this is a guess
	int nNanoSphereParticles=nUnitCell*((4.0/3.0)*M_PI*(radius+latticeLength)*(radius+latticeLength)*(radius+latticeLength))
				/(latticeLength*latticeLength*latticeLength);
	
	//to avoid allocation slowdown
	System.allocParticle(nNanoSphereParticles+offset);
	
	//fill sphere
	for(threeVector<T> clp=0;clp.x<radius*2.0;clp.x+=latticeLength)
	{
		for(clp.y=0;clp.y<radius*2.0;clp.y+=latticeLength)
		{
			for(clp.z=0;clp.z<radius*2.0;clp.z+=latticeLength)
			{
				//set current lattice index to base particle index
				threeVector<int> cli;
				cli.x=floor(clp.x/latticeLength);
				cli.y=floor(clp.y/latticeLength);
				cli.z=floor(clp.z/latticeLength);
				
				//set all subsequent particles to a position in the unit cell
				for(int i=0;i<nUnitCell;i++)
				{
					position<T> p=unitCell[i];//recieves type information in this step
					//set at current lattice position
					p.x*=latticeLength;
					p.y*=latticeLength;
					p.z*=latticeLength;
					p.x+=clp.x+pos.x-radius;
					p.y+=clp.y+pos.y-radius;
					p.z+=clp.z+pos.z-radius;
					
					T theta,phi;
					
					threeVector<T> v,a;
					a=0;
					//velocity
					theta=M_PI*randNum.rand53();
					phi=M_PI*2*randNum.rand53();
					v.x=Vrms*cos(phi)*sin(theta);
					v.y=Vrms*sin(phi)*sin(theta);
					v.z=Vrms*cos(theta);
					//save it
					//using point (radius,radius,radius) as center of mass
					threeVector<T> d;
					d.x=p.x-pos.x;
					d.y=p.y-pos.y;
					d.z=p.z-pos.z;
					
					if(p.z>planeHeight+pos.z-radius)
						p.type=alterType;
					
					//is the current lattice point out of range?
					if(d.x*d.x+d.y*d.y+d.z*d.z<radius*radius)
						System.addParticle(p,v,a);
				}
			}
		}
	}
	
	//set up a pair list
	Cell<double> nearbyPairs(System.getPositions()+offset, System.readNParticles()-offset, latticeLength*1.1, System.readSize());
	nearbyPairs.build();
	
	searchForBonds<double> newSearch(System.getPositions()+offset, System.readNParticles()-offset, bonded, 
					 nBonded, constants, nConstants, latticeLength, offset, System.readSize());
	nearbyPairs.twoWayCompareIndex(newSearch);
	
	molecule< T, fourVector<int> > m;
	m=newSearch.getMolecule();
	std::cerr << m.readNBond() << std::endl;
	//std::cin.get();
	System.addMolecule(m);
	
	return radius;
}

template <typename T>
T rigidSphere(Blob<T> &System, threeVector<T> pos, T radius, T latticeLength, position<T> *unitCell, int nUnitCell, 
	      T *constants, int nConstants)
{
	//I should really run through everything and replace pointers or something...
	//Error checking
	if(nConstants!=1)
	{
		std::cerr << "Error(nanoSphere): nConstants is not 1!\n";
		return 0;
	}
	if(constants==NULL)
	{
		std::cerr << "Error(nanoSphere): Bad allocation of constants!\n";
		return 0;
	}
	if(radius<=0)
	{
		std::cerr << "Error(nanoSphere): No radius!\n";
		return 0;
	}
	if(nUnitCell<=0)
	{
		std::cerr << "Error(nanoSphere): nUnitCell<=0!\n";
		return 0;
	}
	if(unitCell==NULL)
	{
		std::cerr << "Error(nanoSphere): Bad allocation of unitCell!\n";
		return 0;
	}
	std::cerr << "nConstants is " << nConstants << " in nanoSphere.\n";
	
	//Build and add nano-particle
	//index cells on a cube for faster searching
	fourVector<int> nLattice;
	nLattice.x=floor(radius*2.0/latticeLength)+1;
	nLattice.y=nLattice.x;
	nLattice.z=nLattice.z;
	nLattice.t=nLattice.x*nLattice.y*nLattice.z;
	
	//Random number generator
	MTRand randNum(System.readSeed());
	
	//Velocity RMS
	T Vrms=sqrt(3.0*System.readInitialTemp());
	
	int offset=System.readNParticles();
	
	//this is a guess
	int nNanoSphereParticles=nUnitCell*((4.0/3.0)*M_PI*(radius+latticeLength)*(radius+latticeLength)*(radius+latticeLength))
				/(latticeLength*latticeLength*latticeLength);
	
	//to avoid allocation slowdown
	System.allocParticle(nNanoSphereParticles+offset);
	
	//fill sphere
	for(threeVector<T> clp=0;clp.x<radius*2.0;clp.x+=latticeLength)
	{
		for(clp.y=0;clp.y<radius*2.0;clp.y+=latticeLength)
		{
			for(clp.z=0;clp.z<radius*2.0;clp.z+=latticeLength)
			{
				//set current lattice index to base particle index
				threeVector<int> cli;
				cli.x=floor(clp.x/latticeLength);
				cli.y=floor(clp.y/latticeLength);
				cli.z=floor(clp.z/latticeLength);
				
				//set all subsequent particles to a position in the unit cell
				for(int i=0;i<nUnitCell;i++)
				{
					position<T> p=unitCell[i];//recieves type information in this step
					//set at current lattice position
					p.x*=latticeLength;
					p.y*=latticeLength;
					p.z*=latticeLength;
					p.x+=clp.x+pos.x-radius;
					p.y+=clp.y+pos.y-radius;
					p.z+=clp.z+pos.z-radius;
					
					T theta,phi;
					
					threeVector<T> v,a;
					a=0;
					//velocity
					theta=M_PI*randNum.rand53();
					phi=M_PI*2*randNum.rand53();
					v.x=Vrms*cos(phi)*sin(theta);
					v.y=Vrms*sin(phi)*sin(theta);
					v.z=Vrms*cos(theta);
					//save it
					//using point (radius,radius,radius) as center of mass
					threeVector<T> d;
					d.x=p.x-pos.x;
					d.y=p.y-pos.y;
					d.z=p.z-pos.z;
					
					T s=d.x*d.x+d.y*d.y+d.z*d.z;
					//is the current lattice point out of range?
					if(s<radius*radius && s>(radius-2.0)*(radius-2.0))
						System.addParticle(p,v,a);
				}
			}
		}
	}
	
	molecule< T, fourVector<int> > m;
	
	m.setType(SOLID);
	m.addConstant(constants[0]);
		
	m.allocBonds(System.readNParticles()-offset);
	
	for(int i=offset;i<System.readNParticles();i++)
	{
		fourVector<int> bond;
		bond.s[0]=i;
		m.addBond(bond);
	}
	
	std::cerr << m.readNBond() << std::endl;
	//std::cin.get();
	System.addMolecule(m);
	return radius;
}

/** \brief Adds a generic nano-particle to the system. 
 * Generic nano-particle is generated via a unit cell that is placed over and over until it fills the volume.
 * Takes only parameters like latticeLength. pos is the position of the nano-particle in the system. exclude template
 * class handles elimination of out of bound particles. Bond search might leave some particles out of bonding if
 * excluded volume is too complicated. Center of particle might not be center of mass, just center of lattice. nCells
 * controls the number of cells along each direction of a rectangular polygon initially centered at zero coordinate.
 * unitCell is the actual unit cell. An FCC lattice would have 3 Face centered positions and a position at the origin.
 * Relative positions (for FCC) would be (0,0,0), (0.5,0.5,0), (0.5,0,0.5), and (0,0.5,0.5). So nUnitCell for a FCC is 4.
 * bonded is the position of bonded particles. Typically, an FCC would have 6: (0.5,0.5,0), (0.5,0,0.5), (0,0.5,0.5),
 * (-0.5,-0.5,0), (-0.5,0,-0.5), and (0,-0.5,-0.5). 6 would be nBonded. The bond length would be 'sqrt(2.0*0.5*0.5)/2.0' for a
 * perfectly face centered lattice. There is normally an error in the above, but in this particular function, it check if those
 * bonds exist. If they don't exist (within LATTICE_THRESHOLD=0.01), then they are discarded, leaving only certain neighbors.
 */
template <typename T, class exclude>
T nanoLatticeExclude(Blob<T> &System, threeVector<T> pos, T latticeLength, threeVector<int> nCells, position<T> *unitCell, 
		     int nUnitCell, position<T> *bonded, int nBonded, T *constants, int nConstants)
{
	//I should really run through everything and replace pointers or something...
	//Error checking
	if(!(nConstants==4 || nConstants==2))
	{
		std::cerr << "Error(nanoSphere): nConstants is not 4 or 2!\n";
		return 0;
	}
	if(constants==NULL)
	{
		std::cerr << "Error(nanoSphere): Bad allocation of constants!\n";
		return 0;
	}
	if(nUnitCell<=0)
	{
		std::cerr << "Error(nanoSphere): nUnitCell<=0!\n";
		return 0;
	}
	if(unitCell==NULL)
	{
		std::cerr << "Error(nanoSphere): Bad allocation of unitCell!\n";
		return 0;
	}
	if(nBonded<=0)
	{
		std::cerr << "Error(nanoSphere): nBonded<=0!\n";
		return 0;
	}
	if(bonded==NULL)
	{
		std::cerr << "Error(nanoSphere): Bad allocation of bonded!\n";
		return 0;
	}
	std::cerr << "nConstants is " << nConstants << " in nanoSphere.\n";
	
	//Random number generator
	MTRand randNum(System.readSeed());
	
	//Velocity RMS
	T Vrms=sqrt(3.0*System.readInitialTemp());
	
	int offset=System.readNParticles();
	
	//to avoid allocation slowdown
	System.allocParticle(nCells.x*nCells.y*nCells.z*nUnitCell+offset);
	
	//fill sphere
	for(threeVector<T> clp=0;clp.x<static_cast<T>(nCells.x)*latticeLength;clp.x+=latticeLength)
	{
		for(clp.y=0;clp.y<static_cast<T>(nCells.y)*latticeLength;clp.y+=latticeLength)
		{
			for(clp.z=0;clp.z<static_cast<T>(nCells.z)*latticeLength;clp.z+=latticeLength)
			{
				//set current lattice index to base particle index
				threeVector<int> cli;
				cli.x=floor(clp.x/latticeLength);
				cli.y=floor(clp.y/latticeLength);
				cli.z=floor(clp.z/latticeLength);
				
				//set all subsequent particles to a position in the unit cell
				for(int i=0;i<nUnitCell;i++)
				{
					position<T> p=unitCell[i];//recieves type information in this step
					//set at current lattice position
					p.x*=latticeLength;
					p.y*=latticeLength;
					p.z*=latticeLength;
					p.x+=clp.x-static_cast<T>(nCells.x)*latticeLength/2.0;
					p.y+=clp.y-static_cast<T>(nCells.y)*latticeLength/2.0;
					p.z+=clp.z-static_cast<T>(nCells.z)*latticeLength/2.0;
					
					T theta,phi;
					
					threeVector<T> v,a;
					a=0;
					//velocity
					theta=M_PI*randNum.rand53();
					phi=M_PI*2*randNum.rand53();
					v.x=Vrms*cos(phi)*sin(theta);
					v.y=Vrms*sin(phi)*sin(theta);
					v.z=Vrms*cos(theta);
					//save it
					//is the current lattice point out of range?
					if(!exclude(p))
					{
						p.x+=pos.x;
						p.y+=pos.y;
						p.z+=pos.z;
						System.addParticle(p,v,a);
					}
				}
			}
		}
	}
	
	//set up a pair list
	Cell<double> nearbyPairs(System.getPositions()+offset, System.readNParticles()-offset, latticeLength*1.1, System.readSize());
	nearbyPairs.build();
	
	searchForBonds<double> newSearch(System.getPositions()+offset, System.readNParticles()-offset, bonded, 
					 nBonded, constants, nConstants, latticeLength, offset, System.readSize());
	nearbyPairs.twoWayCompareIndex(newSearch);
	
	molecule< T, fourVector<int> > m;
	m=newSearch.getMolecule();
	std::cerr << m.readNBond() << std::endl;
	//std::cin.get();
	System.addMolecule(m);
	return 1;
}

/*
//searches for midpoints in a walk around a polyhedra, derived from someone elses code for generating tesselated polyhedra
template <typename T>
int search_midpoint(position<T> *vertices, threeVector<int> *edge, int &edge_walk, int &nVertices, int index_start, int index_end) 
{ 
	//Note these values for the edges
	//edge[i].x=START;
	//edge[i].y=MIDDLE;
	//edge[i].z=END;
	
	for (int i=0; i<edge_walk; i++) 
		if ((edge[i].x == index_start && edge[i].z == index_end) || 
			(edge[i].x == index_end && edge[i].z == index_start)) 
		{
			int res = edge[i].y;
			

			edge[i].x=edge[edge_walk-1].x;
			edge[i].z=edge[edge_walk-1].z;
			edge[i].y=edge[edge_walk-1].y;
			edge_walk--;
			return res; 
		}


	edge[edge_walk].x = index_start;
	edge[edge_walk].z = index_end; 
	edge[edge_walk].y = nVertices; 
	

	vertices[nVertices].x = (vertices[index_start].x + vertices[index_end].x) / 2.0;
	vertices[nVertices].y = (vertices[index_start].y + vertices[index_end].y) / 2.0;
	vertices[nVertices].z = (vertices[index_start].z + vertices[index_end].z) / 2.0;
	

	T length = sqrt (vertices[nVertices].x * vertices[nVertices].x +
				vertices[nVertices].y * vertices[nVertices].y +
				vertices[nVertices].z * vertices[nVertices].z);
	length = 1/length;
	vertices[nVertices].x *= length;
	vertices[nVertices].y *= length;
	vertices[nVertices].z *= length;
	
	nVertices++;
	edge_walk++;
	
	return edge[edge_walk-1].y;
}
*/

template <typename T>
T shellNano(Blob<T> &System, threeVector<T> pos, T radius, int nTess, int type)
{
	//Error checking
	if(radius<=0)
	{
		std::cerr << "Error(sphericalCyto): No radius!\n";
		return 0;
	}
	if(nTess<0)
	{
		std::cerr << "Error(sphericalCyto): No tesselation!\n";
		return 0;
	}
	
	//tesselate an icosahedron
	int nVertices = 12, nEdges=30, nFaces = 20;
	
	//these are normalized constants
	T t = (1+(T)sqrt(5.0))/2;
	T tau = t/(T)sqrt(1.0+t*t);
	T one = 1/(T)sqrt(1.0+t*t);
	
	position <T> *anc = (position <T> *)malloc(nVertices*sizeof(position<T>));
	if(anc==NULL)
	{
		std::cerr << "Not enough memory!\n";
		return 0;
	}
	
	//vertices of an icosahedron
	anc[0].x=tau;
	anc[0].y=one;
	anc[0].z=0.0;
	anc[1].x=-tau;
	anc[1].y=one;
	anc[1].z=0.0;
	anc[2].x=-tau;
	anc[2].y=-one;
	anc[2].z=0.0;
	anc[3].x=tau;
	anc[3].y=-one;
	anc[3].z=0.0;
	anc[4].x=one;
	anc[4].y=0.0;
	anc[4].z=tau;
	anc[5].x=one;
	anc[5].y=0.0;
	anc[5].z=-tau;
	anc[6].x=-one;
	anc[6].y=0.0;
	anc[6].z=-tau;
	anc[7].x=-one;
	anc[7].y=0.0;
	anc[7].z=tau;
	anc[8].x=0.0;
	anc[8].y=tau;
	anc[8].z=one;
	anc[9].x=0.0;
	anc[9].y=-tau;
	anc[9].z=one;
	anc[10].x=0.0;
	anc[10].y=-tau;
	anc[10].z=-one;
	anc[11].x=0.0;
	anc[11].y=tau;
	anc[11].z=-one;
	
	threeVector<int> *faces = (threeVector<int>*)malloc(nFaces*sizeof(threeVector<int>));
	if(faces==NULL)
	{
		std::cerr << "Not enough memory!\n";
		return 0;
	}
	
	faces[0].s[0]=4;
	faces[0].s[1]=8;
	faces[0].s[2]=7;
	faces[1].s[0]=4;
	faces[1].s[1]=7;
	faces[1].s[2]=9;
	faces[2].s[0]=5;
	faces[2].s[1]=6;
	faces[2].s[2]=11;	
	faces[3].s[0]=5;
	faces[3].s[1]=10;
	faces[3].s[2]=6;
	faces[4].s[0]=0;
	faces[4].s[1]=4;
	faces[4].s[2]=3;
	faces[5].s[0]=0;
	faces[5].s[1]=3;
	faces[5].s[2]=5;
	faces[6].s[0]=2;
	faces[6].s[1]=7;
	faces[6].s[2]=1;
	faces[7].s[0]=2;
	faces[7].s[1]=1;
	faces[7].s[2]=6;
	faces[8].s[0]=8;
	faces[8].s[1]=0;
	faces[8].s[2]=11;
	faces[9].s[0]=8;
	faces[9].s[1]=11;
	faces[9].s[2]=1;
	faces[10].s[0]=9;
	faces[10].s[1]=10;
	faces[10].s[2]=3;
	faces[11].s[0]=9;
	faces[11].s[1]=2;
	faces[11].s[2]=10;
	faces[12].s[0]=8;
	faces[12].s[1]=4;
	faces[12].s[2]=0;
	faces[13].s[0]=11;
	faces[13].s[1]=0;
	faces[13].s[2]=5;
	faces[14].s[0]=4;
	faces[14].s[1]=9;
	faces[14].s[2]=3;
	faces[15].s[0]=5;
	faces[15].s[1]=3;
	faces[15].s[2]=10;
	faces[16].s[0]=7;
	faces[16].s[1]=8;
	faces[16].s[2]=1;
	faces[17].s[0]=6;
	faces[17].s[1]=1;
	faces[17].s[2]=11;
	faces[18].s[0]=7;
	faces[18].s[1]=2;
	faces[18].s[2]=9;
	faces[19].s[0]=6;
	faces[19].s[1]=10;
	faces[19].s[2]=2;
	
	for (int i=0; i<nTess; i++)
	{
		int n_a_new = nVertices+2*nEdges; 
		int n_faces_new = 4*nFaces; 
		
		int edge_walk = 0; 
		nEdges = 2*nVertices + 3*nFaces; 
		threeVector<int> *edge = new threeVector<int>[nEdges]; 
		
		for(int j=0;j<nEdges;j++)
		{
			edge[j].s[0]=-1;
			edge[j].s[1]=-1;
			edge[j].s[2]=-1;
		}
		
		threeVector<int> *faces_old=new threeVector<int>[nFaces]; 
		for(int j=0;j<nFaces;j++)
			faces_old[j]=faces[j];
		anc=(position<T>*)realloc ((void*)anc, n_a_new*sizeof(position<T>));
		if(anc==NULL)
		{
			std::cerr << "Not enough memory!\n";
			return 0;
		}
		
		faces=(threeVector<int>*)realloc ((void*)faces, n_faces_new*sizeof(threeVector<int>));
		if(faces==NULL)
		{
			std::cerr << "Not enough memory!\n";
			return 0;
		}
		
		n_faces_new=0; 
		
		for (int j=0; j<nFaces; j++) 
		{ 
			int xa = faces_old[j].s[0]; 
			int xb = faces_old[j].s[1]; 
			int xc = faces_old[j].s[2]; 
			
			int ab_midpoint = search_midpoint (anc, edge, edge_walk, nVertices, xb, xa); 
			int bc_midpoint = search_midpoint (anc, edge, edge_walk, nVertices, xc, xb); 
			int ca_midpoint = search_midpoint (anc, edge, edge_walk, nVertices, xa, xc); 
			
			faces[n_faces_new].s[0] = xa; 
			faces[n_faces_new].s[1] = ab_midpoint; 
			faces[n_faces_new].s[2] = ca_midpoint; 
			n_faces_new++;
			
			faces[n_faces_new].s[0] = ca_midpoint; 
			faces[n_faces_new].s[1] = ab_midpoint; 
			faces[n_faces_new].s[2] = bc_midpoint; 
			n_faces_new++;
			
			faces[n_faces_new].s[0] = ca_midpoint; 
			faces[n_faces_new].s[1] = bc_midpoint; 
			faces[n_faces_new].s[2] = xc;
			n_faces_new++;
			
			faces[n_faces_new].s[0] = ab_midpoint; 
			faces[n_faces_new].s[1] = xb; 
			faces[n_faces_new].s[2] = bc_midpoint; 
			n_faces_new++; 
		}
		
		nFaces = n_faces_new;
		
		delete faces_old;
		delete edge;
	}
	
	std::cerr << "nFaces: " << nFaces << " nVertices: " << nVertices << '\n';
	
	//expand to radius-1.0 (just below radius of a liposome), assumes it is normalized
	
	for(int i=0;i<nVertices;i++)
	{
		anc[i].x*=(radius);
		anc[i].y*=(radius);
		anc[i].z*=(radius);
	} 
	
	/*
	//Adding new molecules to system
	fourVector<int> bond;
	molecule<T,fourVector<int> > chainMol;
	
	chainMol.setType(CHAIN);
	
	//constants for CHAIN
	chainMol.addConstant(constants[0]);
	chainMol.addConstant(constants[1]);
	chainMol.addConstant(constants[2]);
	chainMol.addConstant(constants[3]);
	
	//this is used for the anchors as well
	int startOffset=System.readNParticles();
	
	//bond for a chain type
	bond.s[START]=startOffset+nVertices;
	bond.s[NCHAINS]=(nFaces+nVertices)-2;
	bond.s[CHAINLENGTH]=nMonomers;
	chainMol.addBond(bond);
	System.addMolecule(chainMol);
	*/
	//This is important to reduce excessive reallocations
	System.allocParticle(System.readNParticles()+nVertices);
	
	//Velocities
	MTRand *randNum=new MTRand(System.readSeed());
	T Vrms=sqrt(3.0*System.readInitialTemp());
	
	//molecule<T,fourVector<int> > anchorMol;
	
	//put vertices in system
	for(int i=0;i<nVertices;i++)
	{
		position<T> p;
		threeVector<T> v,a;
		a.x=0;
		a.y=0;
		a.z=0;
		p.x=anc[i].x+pos.x;
		p.y=anc[i].y+pos.y;
		p.z=anc[i].z+pos.z;
		p.type=type;
		//velocity
		T theta=M_PI*randNum->rand53();
		T phi=M_PI*2*randNum->rand53();
		v.x=Vrms*cos(phi)*sin(theta);
		v.y=Vrms*sin(phi)*sin(theta);
		v.z=Vrms*cos(theta);
		
		System.addParticle(p,v,a);
	}
	/*
	anchorMol.setType(BOND);
	
	//constants for BOND
	anchorMol.addConstant(constants[0]);
	anchorMol.addConstant(constants[1]);
	
	//Just allocate more than you expect (12 of them have 5 each, but midpoints have 6)
	anchorMol.allocBonds(nVertices*6);
	
	
	//place monomers between vertices
	for(int i=0;i<nFaces;i++)
	{
		//indices to points on the face
		int xa=faces[i].s[0];
		int xb=faces[i].s[1];
		int xc=faces[i].s[2];
		
		//checking for edge's existance on previous faces
		bool ab=true;
		bool bc=true;
		bool ca=true;
		
		//null vector
		threeVector<T> a;
		a.x=0;
		a.y=0;
		a.z=0;
		
		//check if edge pairs exist on other faces
		for(int j=0;j<i;j++)
		{
			int ya=faces[j].s[0];
			int yb=faces[j].s[1];
			int yc=faces[j].s[2];

			if((xa==ya || xa==yb || xa==yc) && (xb==ya || xb==yb || xb==yc))
				ab=false;
			if((xb==ya || xb==yb || xb==yc) && (xc==ya || xc==yb || xc==yc))
				bc=false;
			if((xc==ya || xc==yb || xc==yc) && (xa==ya || xa==yb || xa==yc))
				ca=false;
		}
		
		//place monomers between vertices if edge isn't a duplicate
		if(ab)
		{
			//bond for a bond type
			bond.s[0]=startOffset+xa;
			bond.s[1]=System.readNParticles()-1;
			anchorMol.addBond(bond);
			
			//bond for a bond type
			bond.s[0]=startOffset+xb;
			bond.s[1]=System.readNParticles()-nMonomers;
			anchorMol.addBond(bond);
		}
		
		if(bc)
		{
			//bond for a bond type
			bond.s[0]=startOffset+xb;
			bond.s[1]=System.readNParticles()-1;
			anchorMol.addBond(bond);
			
			//bond for a bond type
			bond.s[0]=startOffset+xc;
			bond.s[1]=System.readNParticles()-nMonomers;
			anchorMol.addBond(bond);
		}
		
		if(ca)
		{
			//bond for a bond type
			bond.s[0]=startOffset+xc;
			bond.s[1]=System.readNParticles()-1;
			anchorMol.addBond(bond);
			
			//bond for a bond type
			bond.s[0]=startOffset+xa;
			bond.s[1]=System.readNParticles()-nMonomers;
			anchorMol.addBond(bond);
		}
	}
	
	//add the vertex bonds to the system
	System.addMolecule(anchorMol);
	*/
	free(faces);
	free(anc);
	return radius;
}
