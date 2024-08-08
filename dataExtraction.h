/**
 * \brief Data extraction while the program runs.
 * This is the data extraction header. It is currently used for molecular dynamics
 * (see MD.cpp), but it is probably general enough for other usage, like DPD. This
 * is a user modifiable header that can extract quantities like kinetic or potential
 * energy, center of mass, tension, etc... Some objects are defined twice, once here
 * and once in MD.cpp; that shouldn't be a problem unless you run out of memory. The
 * basic usage of this is as follows:
 * 	1. Build your system and this object.
 * 	2. Initialize some values.
 *	3. Run MD until measure interval is complete.
 * 	4. Execute the compute member function of this object.
 * 	5. repeat step 2 unless system is done.
 * 	6. Destroy objects and exit.
 * 
 * What you are modifying here is three parts:
 * 	1. The constructor member is run during build phase (step 1 above).
 * 	2. The initialize member is run during initialization (step 2 above).
 * 	2. The compute member is executed every measure interval (step 4 above).\n
 * 	3. The destructor member is run prior to exiting (step 6 above).
 * 
 * Also, note that you are passing the blob object into the object's constructor; so
 * at an interval, you should have all the information about the system. 
 * 
 * This is also a template. T is the numeric type (float, double, int, etc...) and
 * U is the system blob type (see systemMD.h and systemDPD.h for an example of the blob).
 * 
 * If you plan on modifying the object, then you should probably put quickly computed 
 * variables in the compute member, and persistent variables, which require information
 * from the previous steps, as private members of this object. 
 */

template <typename T, typename U>
class dataExtraction {
	public:
		/**
		 * The constructor, where blob is a pointer (don't forget the reference operator (&)!)
		 * and name is a pointer to a '\0' terminated char string. Execute this after you build
		 * the system (after reading in an mpd file, for example). 
		 */
		dataExtraction(U *blob, char *name);
		
		/**
		 * Alternative constructor, where blob is a pointer (don't forget the reference operator (&)!)
		 * and name is a pointer to a '\0' terminated char string. Execute this after you build
		 * the system (after reading in an mpd file, for example). Also accepts the absolute position
		 * of a particle list, so that you can track diffusion.
		 */
		dataExtraction(U *blob, char *name, position<T> *absolutePosition);
		
		/**
		 * The destructor, where old memory goes to die. This should be executed after it falls
		 * out of scope, like when a function (e.g. main()) exits.
		 */
		~dataExtraction();
		
		/**
		 * An initialize member. Call this after you build it. It automatically calls compute too.
		 * Mainly, it just does some extra setup prior to your molecular dynamics loop.
		 */
		void initialize();
		
		/**
		 * This does some computations faster. Data extraction is still put off until compute is called.
		 */
		void computeFast();
		
		/**
		 * A compute member. Call this every measure interval to update the data extraction. When
		 * you modify it, don't forget to close files when you are done using them. Typical extraction
		 * is as follows:
		 * 	1. Set values you want.
		 * 	2. Compute against data in the system blob.
		 * 	3. Open a file.
		 * 	4. Write computed data to the file.
		 * 	5. Close the file.
		 * 
		 * Although, you could swap steps 2 and 3.
		 */
		void compute();
		
		/**
		 * Starts calculation of diffusion.
		 */
		void startDiffusion();
		
		
		#ifdef ANCHOR_DATA
			/**
			* An extra member to view the current state of a parameter. I use this to cut the simulation
			* short if there are too many anchors broken. I recommend adding whatever you want.
			* Not active unless ANCHOR type is defined!
			*/
			int readNAnchors()
			{
				return nAnchor;
			};
			
			/**
			* An extra member to view the current state of a parameter. I use this to cut the simulation
			* short if there are too many anchors broken. I recommend adding whatever you want.
			* Not active unless ANCHOR type is defined!
			*/
			int readNBrokenAnchors()
			{
				return nBroken;
			};
		#endif
		
	private:
		//These variables are typically set during construction:
			//The system blob.
			U *System;
			
			//The name of the system.
			char *name;
			
			//Non-bonded pair interactions. Persistant object, but you have to execute
			//the build member prior to use. 
			//See 'include/algorithms/cellOpt.h'.
			CellOpt<T, Potential<T>, Force <T> > pairInteractions;
			
			//An object that computes kinetic energy. Doesn't really need to be persistent.
			//See 'include/algorithms/dataCollection.h'.
			Kinetic<T> kinetic;
			
			//Tracks the kinetic energy density distribution as a histagram
			std::vector<int> kEnergyDensity;
			T kEnergyDensityPartition;//actual length of energy ranges
			
			//An object that computes inner volumes. This one is quite complex. 
			//See 'include/algorithms/volume.h'.
			//VolumeExtraction<T> volumize;
			
			//An object for locating neighbors with cell lists. (Standard neighbor list for system)
			//See 'include/algorithms/cell.h'.
			Cell<T> neighbors;
			
			//An object for getting nearby bleb particles. (Different indices)
			//See 'include/algorithms/cell.h'.
			Cell<T> blebNeighbors;
			
			//A stack for recursion.
			int *stack;
			
			//Flags for our stack.
			int *flag;
			
			//Bleb particle indices.
			int *blebParticles;
			
			//Flip flop object? (nComponents, inner ratio, outer ratio, ...)
			
			
		//Useful references for past (persistent) data
			//#ifdef SOLVENT
				//Number of particles exchanged between inside and outside.
				//Need to know what it was before to compare.
				int nExchanged;
			//#endif
				
				//For flip flop rate calculations?
				//bool *flipFlop;
				
			#ifdef ANCHOR_DATA
				//A list of anchor indices. Not necessarily persistent, just difficult to compute
				//in large systems.
				int *anchorIndex;
				
				//Tracks broken anchors.
				bool *brokenAnchor;
				
				//Total number of anchors.
				int nAnchor;
				
				//Total number of broken anchors
				int nBroken;
				
				//Indices of clusters.
				int *clusterIndices;
				
				//Initial anchor distance.
				T initialAnchorDistance;
				
				//Average anchor distance.
				T anchorDistance;
				
				//A vector containing all anchor connections, shares indices with anchorIndex above.
				std::vector< std::vector <int> > anchorConnections;
				
				//A vector containing all anchor and cytoskeleton particles, We are using this in
				// place of a neighbor list under the assumption that there are far fewer cytoskeleton
				// particles than there are neighbors when we search for exclusions.
				std::vector<int> cytoAnchorList;
				
				//Track mean squared displacement of cytoskeleton anchor, for execluded volume.
				std::vector< std::vector<T> > msCytoAnchor;
				
				//Track number of mean square displacement datapoints.
				int nMSCytoAnchorSteps;
				
				//Old anchor positions
				std::vector< position<T> > oldAnchorPositions;
				
				//list of cytoskeleton molecules
				std::vector<int> cytoList;
			#endif
			
			#ifdef FLAT_MEMBRANE
				//Height map for flat membrane fluctuations.
				T *heightMap;
				
				//Height map number of elements per cell.
				T *heightMapElements;
				
				//Number of heightMap cells projected onto the Z plane.
				twoVector<int> heightMapCells;
				
				//Size of each heightMap Cell.
				twoVector<T> heightMapCellSize;
			#endif
			
			#ifdef NANOPARTICLE
				int nanoParticleOffset;
				int nNanoParticleElements;
			#endif
			
			//For diffusion, it is the absolute position of a particle (no pbc)
			position<T> *aP;
			//For diffusion, stores the start positions and signifies the start of diffusion calculations
			position<T> *aPStart;
			
};

template <typename T, typename U>
dataExtraction<T,U>::dataExtraction(U *blob, char *name, position<T> *absolutePosition)
{
	System=blob;
	this->name=name;
	
	//initialize data structures here, some of these are created twice per execution (once in main, once here)
	kinetic.initialize((*System).getVelocities(), (*System).readNParticles());
	//potential info
	pairInteractions.initialize((*System).getPositions(), (*System).getAccelerations(), 
		(*System).getTwoBodyFconst(), (*System).getTwoBodyUconst(), (*System).readNParticles(), 
		(*System).readNTypes(), (*System).readSize(), (*System).readPeriodic(), (*System).readCutoff());
	
	//For volume
	//int *excludeType=new int[(*System).readNTypes()];
	//for(int i=1;i<(*System).readNTypes();i++)
	//	excludeType[i-1]=i;
	
	//volumize.initialize((*System).getPositions(), (*System).readNParticles(), (*System).readSize(), 
	//	(*System).readCutoff(), excludeType, (*System).readNTypes()-1, (*System).readSeed());
	//delete excludeType;
	
	//neighbor list
	neighbors.initialize((*System).getPositions(), (*System).readNParticles(), (*System).readCutoff(), (*System).readSize());
	
	#ifdef ANCHOR_DATA
		nAnchor=0;
		
		//count number of anchors
		for(int i=0;i<(*System).readNParticles();i++)
		{
			if((*System).getPositions()[i].type==ANCHOR)
				nAnchor++;
			if((*System).getPositions()[i].type==ANCHOR || (*System).getPositions()[i].type==CYTO)
				cytoAnchorList.push_back(i);
		}
		
		if(nAnchor!=0)
		{
			anchorIndex=new int[nAnchor];
			brokenAnchor=new bool[nAnchor];
		}
		else
		{
			brokenAnchor=NULL;
			anchorIndex=NULL;
		}
	
		nAnchor=0;
		
		//Fill anchor indices
		for(int i=0;i<(*System).readNParticles();i++)
			if((*System).getPositions()[i].type==ANCHOR)
				anchorIndex[nAnchor++]=i;
		
		//For cluster finding during blebbing
		stack=new int[(*System).readNParticles()];
		flag=new int[(*System).readNParticles()];
		blebParticles=new int[(*System).readNParticles()];
		
		blebNeighbors.initialize((*System).getPositions(), (*System).readNParticles(), (*System).readCutoff(), \
			(*System).readSize(), blebParticles, (*System).readNParticles());
		
		for(int i=0;i<(*System).readNMolecules();i++)
		{
			bool containsCyto=false;
			//pick a structure by type
			switch((*System).getMolecule()[i].readType())
			{
				case BOND:
				{
					for(int l=0;l<(*System).getMolecule()[i].readNBond();l++)
					{
						//These are the first and second particles of the bond
						int firstParticle=(*System).getMolecule()[i].getBonds()[l].s[0];
						int secondParticle=(*System).getMolecule()[i].getBonds()[l].s[1];
						if((*System).getPositions()[firstParticle].type==CYTO ||
						   (*System).getPositions()[firstParticle].type==ANCHOR ||
						   (*System).getPositions()[firstParticle].type==MONOMER)
							containsCyto=true;
							
						if((*System).getPositions()[secondParticle].type==CYTO ||
						   (*System).getPositions()[secondParticle].type==ANCHOR ||
						   (*System).getPositions()[secondParticle].type==MONOMER)
							containsCyto=true;
					}
					break;
				}
				case BEND:
				{
					for(int l=0;l<(*System).getMolecule()[i].readNBond();l++)
					{
						//These are the first and second particles of the bond
						int firstParticle=(*System).getMolecule()[i].getBonds()[l].s[0];
						int secondParticle=(*System).getMolecule()[i].getBonds()[l].s[1];
						int thirdParticle=(*System).getMolecule()[i].getBonds()[l].s[2];
						if((*System).getPositions()[firstParticle].type==CYTO ||
						   (*System).getPositions()[firstParticle].type==ANCHOR ||
						   (*System).getPositions()[firstParticle].type==MONOMER)
							containsCyto=true;
							
						if((*System).getPositions()[secondParticle].type==CYTO ||
						   (*System).getPositions()[secondParticle].type==ANCHOR ||
						   (*System).getPositions()[secondParticle].type==MONOMER)
							containsCyto=true;
							
						if((*System).getPositions()[thirdParticle].type==CYTO ||
						   (*System).getPositions()[thirdParticle].type==ANCHOR ||
						   (*System).getPositions()[thirdParticle].type==MONOMER)
							containsCyto=true;
					}
					break;
				}
				case CHAIN:
				{
					//Go through all bond descriptions
					for(int l=0;l<(*System).getMolecule()[i].readNBond();l++)
					{
						fourVector<int> *bond=(*System).getMolecule()[i].getBonds();
						//bond info
						int start=bond[l].s[START];
						int nChains=bond[l].s[NCHAINS];
						int length=bond[l].s[CHAINLENGTH];
						//go through all chain lengths
						for(int k=start; k<start+length*nChains; k++)
						{
							if((*System).getPositions()[k].type==CYTO ||
							   (*System).getPositions()[k].type==ANCHOR ||
							   (*System).getPositions()[k].type==MONOMER)
								containsCyto=true;
						}
					}
					break;
				}
				default:
				{
					//does nothing
					break;
				}
			}
			if(containsCyto)
				cytoList.push_back(i);
		}
	#endif
		
	#ifdef FLAT_MEMBRANE
		heightMapCells.x=256;
		heightMapCells.y=256;
		heightMapCellSize.x=(*System).readSize().x/T(heightMapCells.x);
		heightMapCellSize.y=(*System).readSize().y/T(heightMapCells.y);
		heightMap=new T[heightMapCells.x*heightMapCells.y];
		heightMapElements=new T[heightMapCells.x*heightMapCells.y];
	#endif
	aP=absolutePosition;
	aPStart=NULL;
	kEnergyDensityPartition=0.0001;
}

template <typename T, typename U>
dataExtraction<T,U>::dataExtraction(U *blob, char *name)
{
	System=blob;
	this->name=name;
	
	//initialize data structures here, some of these are created twice per execution (once in main, once here)
	kinetic.initialize((*System).getVelocities(), (*System).readNParticles());
	//potential info
	pairInteractions.initialize((*System).getPositions(), (*System).getAccelerations(), 
		(*System).getTwoBodyFconst(), (*System).getTwoBodyUconst(), (*System).readNParticles(), 
		(*System).readNTypes(), (*System).readSize(), (*System).readPeriodic(), (*System).readCutoff());
	
	//For volume
	//int *excludeType=new int[(*System).readNTypes()];
	//for(int i=1;i<(*System).readNTypes();i++)
	//	excludeType[i-1]=i;
	
	//volumize.initialize((*System).getPositions(), (*System).readNParticles(), (*System).readSize(), 
	//	(*System).readCutoff(), excludeType, (*System).readNTypes()-1, (*System).readSeed());
	//delete excludeType;
	
	//neighbor list
	neighbors.initialize((*System).getPositions(), (*System).readNParticles(), (*System).readCutoff(), (*System).readSize());
	
	#ifdef ANCHOR_DATA
		nAnchor=0;
		
		//count number of anchors
		for(int i=0;i<(*System).readNParticles();i++)
		{
			if((*System).getPositions()[i].type==ANCHOR)
				nAnchor++;
			if((*System).getPositions()[i].type==ANCHOR || (*System).getPositions()[i].type==CYTO)
				cytoAnchorList.push_back(i);
		}
		
		if(nAnchor!=0)
		{
			anchorIndex=new int[nAnchor];
			brokenAnchor=new bool[nAnchor];
		}
		else
		{
			brokenAnchor=NULL;
			anchorIndex=NULL;
		}
	
		nAnchor=0;
		
		//Fill anchor indices
		for(int i=0;i<(*System).readNParticles();i++)
			if((*System).getPositions()[i].type==ANCHOR)
				anchorIndex[nAnchor++]=i;
		
		//For cluster finding during blebbing
		stack=new int[(*System).readNParticles()];
		flag=new int[(*System).readNParticles()];
		blebParticles=new int[(*System).readNParticles()];
		
		blebNeighbors.initialize((*System).getPositions(), (*System).readNParticles(), (*System).readCutoff(), \
			(*System).readSize(), blebParticles, (*System).readNParticles());
		
				for(int i=0;i<(*System).readNMolecules();i++)
		{
			bool containsCyto=false;
			//pick a structure by type
			switch((*System).getMolecule()[i].readType())
			{
				case BOND:
				{
					for(int l=0;l<(*System).getMolecule()[i].readNBond();l++)
					{
						//These are the first and second particles of the bond
						int firstParticle=(*System).getMolecule()[i].getBonds()[l].s[0];
						int secondParticle=(*System).getMolecule()[i].getBonds()[l].s[1];
						if((*System).getPositions()[firstParticle].type==CYTO ||
						   (*System).getPositions()[firstParticle].type==ANCHOR ||
						   (*System).getPositions()[firstParticle].type==MONOMER)
							containsCyto=true;
							
						if((*System).getPositions()[secondParticle].type==CYTO ||
						   (*System).getPositions()[secondParticle].type==ANCHOR ||
						   (*System).getPositions()[secondParticle].type==MONOMER)
							containsCyto=true;
					}
					break;
				}
				case BEND:
				{
					for(int l=0;l<(*System).getMolecule()[i].readNBond();l++)
					{
						//These are the first and second particles of the bond
						int firstParticle=(*System).getMolecule()[i].getBonds()[l].s[0];
						int secondParticle=(*System).getMolecule()[i].getBonds()[l].s[1];
						int thirdParticle=(*System).getMolecule()[i].getBonds()[l].s[2];
						if((*System).getPositions()[firstParticle].type==CYTO ||
						   (*System).getPositions()[firstParticle].type==ANCHOR ||
						   (*System).getPositions()[firstParticle].type==MONOMER)
							containsCyto=true;
							
						if((*System).getPositions()[secondParticle].type==CYTO ||
						   (*System).getPositions()[secondParticle].type==ANCHOR ||
						   (*System).getPositions()[secondParticle].type==MONOMER)
							containsCyto=true;
							
						if((*System).getPositions()[thirdParticle].type==CYTO ||
						   (*System).getPositions()[thirdParticle].type==ANCHOR ||
						   (*System).getPositions()[thirdParticle].type==MONOMER)
							containsCyto=true;
					}
					break;
				}
				case CHAIN:
				{
					//Go through all bond descriptions
					for(int l=0;l<(*System).getMolecule()[i].readNBond();l++)
					{
						fourVector<int> *bond=(*System).getMolecule()[i].getBonds();
						//bond info
						int start=bond[l].s[START];
						int nChains=bond[l].s[NCHAINS];
						int length=bond[l].s[CHAINLENGTH];
						//go through all chain lengths
						for(int k=start; k<start+length*nChains; k++)
						{
							if((*System).getPositions()[k].type==CYTO ||
							   (*System).getPositions()[k].type==ANCHOR ||
							   (*System).getPositions()[k].type==MONOMER)
								containsCyto=true;
						}
					}
					break;
				}
				default:
				{
					//does nothing
					break;
				}
			}
			if(containsCyto)
				cytoList.push_back(i);
		}
	#endif
		
	#ifdef FLAT_MEMBRANE
		heightMapCells.x=256;
		heightMapCells.y=256;
		heightMapCellSize.x=(*System).readSize().x/float(heightMapCells.x);
		heightMapCellSize.y=(*System).readSize().y/float(heightMapCells.y);
		heightMap=new T[heightMapCells.x*heightMapCells.y];
		heightMapElements=new T[heightMapCells.x*heightMapCells.y];
	#endif
	aP=NULL;
	aPStart=NULL;
	kEnergyDensityPartition=0.0001;
}

template <typename T, typename U>
dataExtraction<T,U>::~dataExtraction()
{
	std::fstream dataFile;
	std::string buf("kEnergyDensity_");
	buf+=name;
	buf+=".dat";
	dataFile.open(buf.c_str(), std::ios::out);
	long int partitionSum=0;
	for(int i=0;i<kEnergyDensity.size();i++)
		partitionSum+=kEnergyDensity[i];
	for(int i=0;i<kEnergyDensity.size();i++)
	{
		dataFile << static_cast<float>(i)*kEnergyDensityPartition << '\t';
		dataFile << static_cast<float>(kEnergyDensity[i])/static_cast<float>(partitionSum) << std::endl;
	}
	dataFile.close();
	#ifdef ANCHOR_DATA
		if(anchorIndex!=NULL)
			delete anchorIndex;
		if(brokenAnchor!=NULL)
			delete brokenAnchor;
		if(stack!=NULL)
			delete stack;
		if(flag!=NULL)
			delete flag;
		if(blebParticles!=NULL)
			delete blebParticles;
	#endif
	
	#ifdef FLAT_MEMBRANE
		if(heightMap!=NULL)
			delete heightMap;
		if(heightMapElements!=NULL)
			delete heightMapElements;
	#endif
	if(aPStart!=NULL)
		delete aPStart;
}

template <typename T, typename U>
void dataExtraction<T,U>::initialize()
{
	//do initial data collection here (usually before system is run)
	nExchanged=0;
	
	//useful references
	position<T> *p=(*System).getPositions();
	threeVector<T> *v=(*System).getVelocities();
	threeVector<T> *a=(*System).getAccelerations();
	threeVector<T> s=(*System).readSize();
	T cutoffSqr=(*System).readCutoff();
	cutoffSqr*=cutoffSqr;
	
	#ifdef ANCHOR_DATA
	if(nAnchor>0)
	{
		int nConnections=0;
		
		//let's just assume no anchors are broken when it starts
		for(int i=0;i<nAnchor;i++)
			brokenAnchor[i]=false;
		
		//quick molecule access
		molecule<T, fourVector<int> > *m=(*System).getMolecule();
		
		initialAnchorDistance=0;
		
		//Quick overview of what will be done:
		// 1. Find a list of BOND types that contains the anchor.
		// 2. Find a  all CHAINs that contain the anchor's bonded particle (likely a cytoskeleton chain).
		// 3. Find a BOND type that contains the CHAIN type's other end.
		// 4. Check if it is bonded to another anchor.
		// 5. Determine the distance between that anchor and the current anchor.
		//Note: this assumes each chain that connects anchors is unique!
		
		//check initial distance between anchors. 
		//This is deeply nested, but should be fairly quick because there are far fewer anchors and 
		// chains than particles in the system.
		for(int anchor=0;anchor<nAnchor;anchor++)
		{
			std::vector<int> connected;
			//1. Find a BOND type that contains the anchor.
			for(int mol=0;mol<(*System).readNMolecules();mol++)
			{
				if(m[mol].readType()==BOND)
				{
					//reset every time we start a new search
					m[mol].resetFind();
					for(twoVector<int> buf=m[mol].findBond(anchorIndex[anchor]);buf.x!=-1;buf=m[mol].findBond(anchorIndex[anchor]))
					{
						//We found a bond
						if(buf.x!=-1)
						{
							//we will push back the current cytoskeleton chain start particle
							if(buf.y==0)
								connected.push_back(m[mol].getBonds()[buf.x].s[1]);
							else
								connected.push_back(m[mol].getBonds()[buf.x].s[0]);
						}
					}
				}
			}
			//2. Find a CHAIN type that contains the anchor's bonded particle
			for(int mol=0;mol<(*System).readNMolecules();mol++)
			{
				if(m[mol].readType()==CHAIN)
				{
					//do it for each connected chain particle
					for(int i=0;i<connected.size();i++)
					{
						//reset every time we start a new search
						m[mol].resetFind();
						twoVector<int> buf=m[mol].findBond(connected[i]);
						//We found a chain
						if(buf.x!=-1)
						{
							int start=m[mol].getBonds()[buf.x].s[START];
							int length=m[mol].getBonds()[buf.x].s[CHAINLENGTH];
							//we will adjust this to the the other end of the chain.
							if((connected[i]-start)%length==0)
								connected[i]+=(length);//connected[i] is at the beginning, grab the one at the end.
							else
								connected[i]-=(length);//connected[i] is at the end, grab the one at the beginning.
						}
						else
						{
							#ifdef WARNINGS_ENABLED
								//std::cout << "Warning (dataExtraction): Can't locate chain!\n";
							#endif
						}
					}
				}
			}
			//3. Find a BOND type that contains the CHAIN type's other end.
			for(int mol=0;mol<(*System).readNMolecules();mol++)
			{
				if(m[mol].readType()==BOND)
				{
					//do it for each connected chain particle
					for(int i=0;i<connected.size();i++)
					{
						//reset every time we start a new search
						m[mol].resetFind();
						twoVector<int> buf=m[mol].findBond(connected[i]);
						//We found a bond
						if(buf.x!=-1)
						{
							//we will push back the current cytoskeleton chain start particle
							if(buf.y==0)
								connected[i]=m[mol].getBonds()[buf.x].s[1];
							else
								connected[i]=m[mol].getBonds()[buf.x].s[0];
						}
						else
						{
							#ifdef WARNINGS_ENABLED
								//std::cout << "Warning (dataExtraction): Can't locate other anchor!\n";
							#endif
						}
					}
				}
			}
			//4. Check if it is bonded to another anchor.
			for(int i=0;i<connected.size();i++)
				if(p[connected[i]].type!=ANCHOR)//is it an anchor?
					connected.erase(connected.begin()+i);//it isn't an anchor, remove it
					
			//We don't have to repeat the above for every compute() call, only the 5th step is repeated in compute().
			anchorConnections.push_back(connected);
			
			//5. Determine the distance between that anchor and the current anchor.
			for(int i=0;i<anchorConnections[anchor].size();i++)
			{
				threeVector<T> d;
				//calculate bond length using minimum image
				d.x=p[anchorIndex[anchor]].x-p[anchorConnections[anchor][i]].x;
				d.y=p[anchorIndex[anchor]].y-p[anchorConnections[anchor][i]].y;
				d.z=p[anchorIndex[anchor]].z-p[anchorConnections[anchor][i]].z;
				
				if(d.x>=s.x/2.0) d.x-=s.x;
				if(d.x<=-s.x/2.0) d.x+=s.x;
				if(d.y>=s.y/2.0) d.y-=s.y;
				if(d.y<=-s.y/2.0) d.y+=s.y;
				if(d.z>=s.z/2.0) d.z-=s.z;
				if(d.z<=-s.z/2.0) d.z+=s.z;
				
				anchorDistance+=sqrt(d.x*d.x+d.y*d.y+d.z*d.z);
				nConnections++;
			}
		}
		if(nConnections>0)
			anchorDistance/=float(nConnections);
		else
			anchorDistance=1.0;
		initialAnchorDistance=anchorDistance;
		
		if(initialAnchorDistance<10.0)//There is a thickness associated with the bilayer, this 2x that thickness
			initialAnchorDistance=10.0;
		
		std::cout << "Using anchor to anchor distance for bleb search cutoff: " << initialAnchorDistance << '\n';
		std::cout << "There are " << nAnchor << " anchors.\n";
		//Initialize mean square displacement
		if(aP!=NULL)
		{
			for(int i=0;i<nAnchor;i++)
				msCytoAnchor.push_back(std::vector<T>());
			
			for(int i=0;i<nAnchor;i++)
				oldAnchorPositions.push_back(aP[cytoAnchorList[i]]);
			
			nMSCytoAnchorSteps=0;
		}
	}
	
	
	
	#endif
	#ifdef NANOPARTICLE
		nanoParticleOffset=-1;
		nNanoParticleElements=0;
		for(int i=0;i<(*System).readNParticles();i++)
		{
			if(p[i].type==NANOPARTICLE && nanoParticleOffset==-1)
			{
				nanoParticleOffset=i;
			}
			if(p[i].type!=NANOPARTICLE && nanoParticleOffset>=0)
			{
				nNanoParticleElements=i-nanoParticleOffset;
				break;
			}
		}
		if(nNanoParticleElements==0)
			nNanoParticleElements=(*System).readNParticles()-nanoParticleOffset;
	#endif
}

template <typename T, typename U>
void dataExtraction<T,U>::startDiffusion()
{
	if(aP!=NULL && aPStart==NULL)
	{
		aPStart=new position<T>[(*System).readNParticles()];
		for(int i=0;i<(*System).readNParticles();i++)
		{
			aPStart[i]=aP[i];
			//Alternative version (Doesn't work the same way!):
			//aPStart[i]=(*System).getPositions()[i];
		}
	}
	
}

template <typename T, typename U>
void dataExtraction<T,U>::computeFast()
{
	#ifdef ANCHOR_DATA
	/*
	//track anchor deviation
	if(aP!=NULL)
	{
		for(int i=0;i<nAnchor;i++)
		{
			//Compute the MSD from their start
			threeVector<T> d;
			d.x=aP[cytoAnchorList[i]].x-oldAnchorPositions[i].x;
			d.y=aP[cytoAnchorList[i]].y-oldAnchorPositions[i].y;
			d.z=aP[cytoAnchorList[i]].z-oldAnchorPositions[i].z;
			msCytoAnchor[i].push_back(sqrt(d.x*d.x+d.y*d.y+d.z*d.z));
		}
	}
	*/
	#endif
}

template <typename T, typename U>
void dataExtraction<T,U>::compute()
{
	//useful references
	position<T> *p=(*System).getPositions();
	threeVector<T> *v=(*System).getVelocities();
	threeVector<T> *a=(*System).getAccelerations();
	threeVector<T> s=(*System).readSize();
	T cutoffSqr=(*System).readCutoff();
	cutoffSqr*=cutoffSqr;
	
	//Cuda works for this
	pairInteractions.resize((*System).readSize());//Can't forget this!
	pairInteractions.build();
	T potential=pairInteractions.computePotential();
	
	//Cuda does not work on this loop yet
	//Better version of molecule interactions
	T lBond=0;
	int nBond=0;
	T costhetaBend=0;
	int nBend=0;
	twoVector<T> lBend;
	lBend.x=0;
	lBend.y=0;
	T beadPotential=0;
	int nBeads=0;
	
	//Go through all molecule structures
	for(int k=0;k<(*System).readNMolecules();k++)
	{
		//pick a structure by type
		switch((*System).getMolecule()[k].readType())
		{
			case BOND:
			{
				potential+=(*System).doBondPotential(k);
				//Go through all in bond list
				T lBondPrivate=0;
				//#pragma omp parallel for reduction(+:lBondPrivate)
				for(int l=0;l<(*System).getMolecule()[k].readNBond();l++)
				{
					//These are the first and second particles of the bond
					int firstParticle=(*System).getMolecule()[k].getBonds()[l].s[0];
					int secondParticle=(*System).getMolecule()[k].getBonds()[l].s[1];
					
					//calculate bond length using minimum image
					T dx=p[firstParticle].x-p[secondParticle].x;
					T dy=p[firstParticle].y-p[secondParticle].y;
					T dz=p[firstParticle].z-p[secondParticle].z;
					
					if(dx>=s.x/2.0) dx-=s.x;
					if(dx<=-s.x/2.0) dx+=s.x;
					if(dy>=s.y/2.0) dy-=s.y;
					if(dy<=-s.y/2.0) dy+=s.y;
					if(dz>=s.z/2.0) dz-=s.z;
					if(dz<=-s.z/2.0) dz+=s.z;
					
					lBondPrivate+=sqrt(dx*dx+dy*dy+dz*dz);
				}
				lBond+=lBondPrivate;
				nBond+=(*System).getMolecule()[k].readNBond();
				break;
			}
			case BEND:
			{
				for(int l=0;l<(*System).getMolecule()[k].readNBond();l++)
				{
					//These are the first and second particles of the bond
					int firstParticle=(*System).getMolecule()[k].getBonds()[l].s[0];
					int secondParticle=(*System).getMolecule()[k].getBonds()[l].s[1];
					int thirdParticle=(*System).getMolecule()[k].getBonds()[l].s[2];
					
					twoVector<T> dr;
					threeVector<T> da, db;
					
					//calculate bond length using minimum image
					da.x=p[firstParticle].x-p[secondParticle].x;
					da.y=p[firstParticle].y-p[secondParticle].y;
					da.z=p[firstParticle].z-p[secondParticle].z;
					
					db.x=p[secondParticle].x-p[thirdParticle].x;
					db.y=p[secondParticle].y-p[thirdParticle].y;
					db.z=p[secondParticle].z-p[thirdParticle].z;
					
					if(da.x>=s.x/2.0) da.x-=s.x;
					if(da.x<=-s.x/2.0) da.x+=s.x;
					if(da.y>=s.y/2.0) da.y-=s.y;
					if(da.y<=-s.y/2.0) da.y+=s.y;
					if(da.z>=s.z/2.0) da.z-=s.z;
					if(da.z<=-s.z/2.0) da.z+=s.z;
					
					if(db.x>=s.x/2.0) db.x-=s.x;
					if(db.x<=-s.x/2.0) db.x+=s.x;
					if(db.y>=s.y/2.0) db.y-=s.y;
					if(db.y<=-s.y/2.0) db.y+=s.y;
					if(db.z>=s.z/2.0) db.z-=s.z;
					if(db.z<=-s.z/2.0) db.z+=s.z;
					
					dr.s[0]=sqrt(da.x*da.x+da.y*da.y+da.z*da.z);
					dr.s[1]=sqrt(db.x*db.x+db.y*db.y+db.z*db.z);
					
					lBend.s[0]+=dr.s[0];
					lBend.s[1]+=dr.s[1];
					
					costhetaBend+=(da.x*db.x+da.y*db.y+da.z*db.z)/(dr.s[0]*dr.s[1]);
				}
				nBend+=(*System).getMolecule()[k].readNBond();
				potential+=(*System).doBendPotential(k);
				break;
			}
			case CHAIN:
			{
				potential+=(*System).doChainPotential(k);
				break;
			}
			case BEAD:
			{
				T tempBeadPotential=(*System).doBeadPotential(k);
				beadPotential+=tempBeadPotential;
				potential+=tempBeadPotential;
				nBeads++;
				break;
			}
			case BOUNDARY:
			{
				potential+=(*System).doBoundaryPotential(k);
				break;
			}
			case FLOATING_BASE:
			{
				potential+=(*System).doFloatingBasePotential(k);
				break;
			}
			case ZTORQUE:
			{
				potential+=(*System).doZTorquePotential(k);
				break;
			}
			case ZPOWERPOTENTIAL:
			{
				potential+=(*System).doZPowerPotential(k);
				break;
			}
			case NANOCORE:
			{
				T tempBeadPotential=(*System).doNanoCorePotential(k);
				beadPotential+=tempBeadPotential;
				potential+=tempBeadPotential;
				nBeads++;
				break;
			}
			case BALL:
			{
				potential+=(*System).doBallPotential(k);
				break;
			}
			default:
			{
				//does nothing
				break;
			}
		}
	}
	
	std::fstream dataFile;
	std::string buf("potential_");
	buf+=name;
	buf+=".dat";
	dataFile.open(buf.c_str(), std::ios::app | std::ios::out);
	dataFile << (*System).readInitialTime() << '\t' << potential << std::endl;
	dataFile.close();
	
	//Just copy this one and use it as a template
	buf.clear();
	buf="size_";
	buf+=name;
	buf+=".dat";
	dataFile.open(buf.c_str(), std::ios::app | std::ios::out);
	dataFile << (*System).readInitialTime() << '\t' << s.x << '\t' << s.y << '\t' << s.z << std::endl;
	dataFile.close();
	
	buf.clear();
	buf="lBond_";
	buf+=name;
	buf+=".dat";
	dataFile.open(buf.c_str(), std::ios::app | std::ios::out);
	dataFile << (*System).readInitialTime() << '\t' << (lBond/(T)nBond) << std::endl;
	dataFile.close();
	
	if(nBend>1)
	{
		costhetaBend/=static_cast<T>(nBend);
		lBend.s[0]/=static_cast<T>(nBend);
		lBend.s[1]/=static_cast<T>(nBend);
	}
	buf.clear();
	buf="bend_";
	buf+=name;
	buf+=".dat";
	dataFile.open(buf.c_str(), std::ios::app | std::ios::out);
	dataFile << (*System).readInitialTime() << '\t' << costhetaBend << '\t' << lBend.s[0] << '\t' << lBend.s[1] << std::endl;
	dataFile.close();
	
	if(nBeads>0)
	{
		buf.clear();
		buf="beadPotential_";
		buf+=name;
		buf+=".dat";
		dataFile.open(buf.c_str(), std::ios::app | std::ios::out);
		dataFile << (*System).readInitialTime() << '\t' << beadPotential << std::endl;
		dataFile.close();
	}
	
	buf.clear();
	buf="temp_";
	buf+=name;
	buf+=".dat";
	dataFile.open(buf.c_str(), std::ios::app | std::ios::out);
	dataFile << (*System).readInitialTime() << '\t' << (*System).readInitialTemp() << std::endl;
	dataFile.close();
	
	kinetic.output((*System).readInitialTime(),name);
	
	//do this before the next 2 outputs
	/*volumize.build();
	
	nExchanged+=volumize.nExchanged();
	if((*System).readRemoveSolvent()>0)
	{
		buf.clear();
		buf="nExchange_";
		buf+=name;
		buf+=".dat";
		dataFile.open(buf.c_str(), std::ios::app | std::ios::out);
		dataFile << (*System).readInitialTime() << '\t' << nExchanged << '\n';
		dataFile.close();
		
		nExchanged=0;
		
		T vInner=0;
		for(int k=1;k<volumize.readNVolumes();k++)
			if(k!=volumize.readOuter())
				vInner=volumize.readVolumeSize(k);
			buf.clear();
		buf="vInner_";
		buf+=name;
		buf+=".dat";
		dataFile.open(buf.c_str(), std::ios::app | std::ios::out);
		dataFile << (*System).readInitialTime() << '\t' << vInner << '\n';
		dataFile.close();
		
		buf.clear();
		buf="vExclude_";
		buf+=name;
		buf+=".dat";
		dataFile.open(buf.c_str(), std::ios::app | std::ios::out);
		dataFile << (*System).readInitialTime() << '\t' << volumize.readVolumeSize(0) << '\n';
		dataFile.close();
	}
	*/
	
	#ifdef ANCHOR_DATA
	//track anchor breaks, and determine changes in elasticity and bleb size
	if(nAnchor>0)
	{
		int nConnections=0;
		
		//neighbor list for anchor breakage
		neighbors.build();
		
		nBroken=0;
		for(int anchor=0;anchor<nAnchor;anchor++)
		{
			bool broken=true;
			threeVector<T> minImg=0;
			
			//check to see if there is a head within cutoff distance
			for(int j=neighbors.query(anchorIndex[anchor],minImg);j!=-1;j=neighbors.query(anchorIndex[anchor],minImg))
			{
				if(p[j].type==HEAD)
				{
					threeVector<T> d;
					
					//calculate distance from anchor
					d.x=p[j].x-p[anchorIndex[anchor]].x;
					d.y=p[j].y-p[anchorIndex[anchor]].y;
					d.z=p[j].z-p[anchorIndex[anchor]].z;
					
					d.x+=((d.x>s.x/2.0)?-s.x:0);
					d.x+=((d.x<-s.x/2.0)?s.x:0);
					
					d.y+=((d.y>s.y/2.0)?-s.y:0);
					d.y+=((d.y<-s.y/2.0)?s.y:0);
					
					d.z+=((d.z>s.z/2.0)?-s.z:0);
					d.z+=((d.z<-s.z/2.0)?s.z:0);
					
					//is it in range of an anchor?
					if(d.x*d.x+d.y*d.y+d.z*d.z<cutoffSqr)
						broken=false;
				}
			}
			
			//Output info about anchor
			if(brokenAnchor[anchor]!=broken)
			{
				if(broken)
					std::cout << "Anchor " << anchor << " broke!\n";
				else
					std::cout << "Anchor " << anchor << " reattached!\n";
			}
			
			//set the current state
			brokenAnchor[anchor]=broken;
			
			//I'm not going to depend on the fact that true==1 in most(all???) cases. (e.g. nBroken+=broken)
			nBroken+=(broken)?1:0;
			
			//5. Determine the distance between that anchor and the current anchor.
			for(int j=0;j<anchorConnections[anchor].size();j++)
			{
				threeVector<T> d;
				//calculate bond length using minimum image
				d.x=p[anchorIndex[anchor]].x-p[anchorConnections[anchor][j]].x;
				d.y=p[anchorIndex[anchor]].y-p[anchorConnections[anchor][j]].y;
				d.z=p[anchorIndex[anchor]].z-p[anchorConnections[anchor][j]].z;
				
				if(d.x>=s.x/2.0) d.x-=s.x;
				if(d.x<=-s.x/2.0) d.x+=s.x;
				if(d.y>=s.y/2.0) d.y-=s.y;
				if(d.y<=-s.y/2.0) d.y+=s.y;
				if(d.z>=s.z/2.0) d.z-=s.z;
				if(d.z<=-s.z/2.0) d.z+=s.z;
				
				anchorDistance+=sqrt(d.x*d.x+d.y*d.y+d.z*d.z);
				nConnections++;
			}
		}
		if(nConnections>0)
			anchorDistance/=float(nConnections);
		else
			anchorDistance=1.0;
		float anchorDistanceSqr=initialAnchorDistance*initialAnchorDistance;
		
		buf.clear();
		buf="brokenAnchors_";
		buf+=name;
		buf+=".dat";
		dataFile.open(buf.c_str(), std::ios::app | std::ios::out);
		dataFile << (*System).readInitialTime() << '\t' << nBroken << std::endl;
		dataFile.close();
		
		//this will help us determine any changes in tension
		buf.clear();
		buf="anchorDistance_";
		buf+=name;
		buf+=".dat";
		dataFile.open(buf.c_str(), std::ios::app | std::ios::out);
		dataFile << (*System).readInitialTime() << '\t' << anchorDistance << std::endl;
		dataFile.close();
		
		//determine potential of cytoskeleton and anchors for tension
		T cytoPotential=0;
		for(int i=0;i<(*System).readNParticles();i++)
			if(p[i].type==CYTO || p[i].type==ANCHOR || p[i].type==MONOMER)
				cytoPotential+=pairInteractions.computePotential(i);
		
		for(int k=0;k<cytoList.size();k++)
		{
			//pick a structure by type
			switch((*System).getMolecule()[cytoList[k]].readType())
			{
				case BOND:
				{
					cytoPotential+=(*System).doBondPotential(cytoList[k]);
					break;
				}
				case BEND:
				{
					cytoPotential+=(*System).doBendPotential(cytoList[k]);
					break;
				}
				case CHAIN:
				{
					cytoPotential+=(*System).doChainPotential(cytoList[k]);
					break;
				}
				case BEAD:
				{
					cytoPotential+=(*System).doBeadPotential(cytoList[k]);
					break;
				}
				case BOUNDARY:
				{
					cytoPotential+=(*System).doBoundaryPotential(cytoList[k]);
					break;
				}	
				case NANOCORE:
				{
					cytoPotential+=(*System).doNanoCorePotential(k);
					break;
				}
				default:
				{
					//does nothing
					break;
				}
			}
		}
			
		buf.clear();
		buf="cytoPotential_";
		buf+=name;
		buf+=".dat";
		dataFile.open(buf.c_str(), std::ios::app | std::ios::out);
		dataFile << (*System).readInitialTime() << '\t' << cytoPotential << std::endl;
		dataFile.close();
		
		//determine potential of cytoskeleton and anchors for tension
		T anchorPotential=0;
		for(int i=0;i<(*System).readNParticles();i++)
			if(p[i].type==ANCHOR)
				anchorPotential+=pairInteractions.computePotential(i);
			
		buf.clear();
		buf="anchorPotential_";
		buf+=name;
		buf+=".dat";
		dataFile.open(buf.c_str(), std::ios::app | std::ios::out);
		dataFile << (*System).readInitialTime() << '\t' << (anchorPotential/T(nAnchor)) << std::endl;
		dataFile.close();
		
		//make sure to rebuild our neighbor search!
		int nBlebParticles=0;
		//Now lets check for blebbing!
		//We will assume the blebs do not contain cytoskeleton or anchors
		for(int i=0;i<(*System).readNParticles();i++)
		{
			if(p[i].type==TAIL)//is it already a cytoskeleton or anchor particle?
			{
				bool range=false;//nope, lets check what's in range.
				//Lets utilize our neighbor list
				for(int j=0;j<cytoAnchorList.size();j++)
				{
					threeVector<T> d;
					
					//calculate distance between neighboring particle
					d.x=p[i].x-p[cytoAnchorList[j]].x;
					d.y=p[i].y-p[cytoAnchorList[j]].y;
					d.z=p[i].z-p[cytoAnchorList[j]].z;
					
					d.x-=(d.x>s.x/2.0)?s.x:0.0;
					d.x+=(d.x<-s.x/2.0)?s.x:0.0;
					
					d.y-=(d.y>s.y/2.0)?s.y:0.0;
					d.y+=(d.y<-s.y/2.0)?s.y:0.0;
					
					d.z-=(d.z>s.z/2.0)?s.z:0.0;
					d.z+=(d.z<-s.z/2.0)?s.z:0.0;
					
					//is it in range of an anchor or cytoskeleton particle?
					if(d.x*d.x+d.y*d.y+d.z*d.z<anchorDistanceSqr)
					{
						range=true;//yes it is in range
						break;
					}
				}
				if(!range)
					blebParticles[nBlebParticles++]=i;//push it on our blebParticles list
			}
		}

		blebNeighbors.changeN(nBlebParticles,blebParticles);
		blebNeighbors.build();
		
		for(int i=0;i<nBlebParticles;i++)
			flag[i]=-1;
		
		std::vector< std::vector<int> > blebs;

		//we've got some blebs! Now let's count them
		for(int i=0;i<nBlebParticles;i++)
		{
			if(flag[i]==-1)
			{
				std::vector<int> bleb;
				stack[0]=i;
				flag[i]=blebs.size();
				for(int j=0;j>-1;j)
				{
					int current=stack[j--];
					bleb.push_back(current);
					blebNeighbors.resetState();
					threeVector<T> minImg;
					
					for(int k=blebNeighbors.queryIndex(current,minImg);k!=-1;k=blebNeighbors.queryIndex(current,minImg))
					{
						if(flag[k]==-1)
						{
							int iIndex=blebParticles[current];
							int kIndex=blebParticles[k];
							threeVector<float> d;
							d.x=p[iIndex].x-p[kIndex].x+minImg.x;
							d.y=p[iIndex].y-p[kIndex].y+minImg.y;
							d.z=p[iIndex].z-p[kIndex].z+minImg.z;
							
							//is it within the range?
							if(d.x*d.x+d.y*d.y+d.z*d.z<cutoffSqr*0.75)
							{
								//exclude it and push it on the stack
								flag[k]=blebs.size();
								stack[++j]=k;
							}
						}
					}
				}
				if(bleb.size()>300)
					blebs.push_back(bleb);
			}
		}
		/*
		for(int i=0;i<blebs.size();i++)
		{
			for(int j=0;j<
		}
		*/
		buf.clear();
		buf="nBlebs_";
		buf+=name;
		buf+=".dat";
		dataFile.open(buf.c_str(), std::ios::app | std::ios::out);
		dataFile << (*System).readInitialTime() << '\t' << blebs.size() << std::endl;
		dataFile.close();
		
		buf.clear();
		buf="blebSize_";
		buf+=name;
		buf+=".dat";
		dataFile.open(buf.c_str(), std::ios::app | std::ios::out);
		if(blebs.size()>0)
			dataFile << (*System).readInitialTime() << '\t' << (nBlebParticles/blebs.size()) << std::endl;
		else
			dataFile << (*System).readInitialTime() << '\t' << 0 << std::endl;
		dataFile.close();
	}
	/*//Need a better measure of anchor volume!
	if(aP!=NULL)
	{
		buf.clear();
		buf="exVol_";
		buf+=name;
		buf+=".dat";
		dataFile.open(buf.c_str(), std::ios::app | std::ios::out);
		//dataFile << (*System).readInitialTime();
		
		std::vector<T> dispHist;
		T stepSize=0.1;
		T range=100.0;
		int nSteps=range/stepSize;
		
		for(int i=0;i<nSteps;i++)
			dispHist.push_back(0.0);
		
		for(int i=0;i<msCytoAnchor.size();i++)
		{
			for(int j=0;j<msCytoAnchor[i].size();j++)
				dispHist[static_cast<int>(msCytoAnchor[i][j]/stepSize)]++;
		*/	/*
			T mean=0;
			for(int j=0;j<msCytoAnchor[i].size();j++)
				mean+=msCytoAnchor[i][j];
			
			if(msCytoAnchor[i].size()!=0)
				mean/=msCytoAnchor[i].size();
			
			T meanSquare=0;
			
			for(int j=0;j<msCytoAnchor[i].size();j++)
				meanSquare+=((msCytoAnchor[i][j]-mean)*(msCytoAnchor[i][j]-mean));
			
			if(msCytoAnchor[i].size()!=0)
				meanSquare/=msCytoAnchor[i].size();
			
			dataFile << '\t' << meanSquare;
			*//*
			msCytoAnchor[i].clear();
			oldAnchorPositions[i]=aP[cytoAnchorList[i]];
		}
		
		for(int i=0;i<nSteps;i++)
			dataFile << static_cast<float>(i)*stepSize+(stepSize/2.0) << '\t' << dispHist[i] << '\n';
		
		dataFile << '\n';
		dataFile.close();
	}
	*/
	#endif
	
	#ifdef FLAT_MEMBRANE
		//initialize heightMap and heightMapElements to 0
		for(int i=0;i<heightMapCells.x*heightMapCells.y;i++)
		{
			heightMap[i]=0;
			heightMapElements[i]=0;
		}
		
		//sum heights from Z plane
		for(int i=0;i<(*System).readNParticles();i++)
		{
			if(p[i].type==HEAD || p[i].type==TAIL)
			{
				//get heightMap index
				int x=p[i].x/heightMapCellSize.x;
				int y=p[i].y/heightMapCellSize.y;
				
				int index=x+y*heightMapCells.x;
				
				heightMap[index]+=p[i].z;
				heightMapElements[index]++;
			}
		}
		
		//Calculate average heightMap
		for(int i=0;i<heightMapCells.x*heightMapCells.y;i++)
			if(heightMapElements[i]>0)
				heightMap[i]/=heightMapElements[i];
		
	#else
		threeVector<T> minPos, maxPos;
		minPos=(*System).readSize();
		maxPos.x=0;
		maxPos.y=0;
		maxPos.z=0;
		
		for(int i=0;i<(*System).readNParticles();i++)
		{
			for(int dim=0;dim<3;dim++)
			{
				if(p[i].s[dim]<minPos.s[dim])
					minPos.s[dim]=p[i].s[dim];
				if(p[i].s[dim]>maxPos.s[dim])
					maxPos.s[dim]=p[i].s[dim];
			}
		}
		
		threeVector<T> flick;
		flick.x=maxPos.x-minPos.x;
		flick.y=maxPos.y-minPos.y;
		flick.z=maxPos.z-minPos.z;
		
		//flick.x+=(flick.x>size.x/2.0)?size.x
		
		buf.clear();
		buf="flicker_";
		buf+=name;
		buf+=".dat";
		dataFile.open(buf.c_str(), std::ios::app | std::ios::out);
		dataFile << (*System).readInitialTime() << '\t' << flick.x << '\t' << flick.y << '\t' << flick.z << std::endl;
		dataFile.close();
		
	#endif
	
	#ifdef NANOPARTICLE
		position<T> nanoCOM=com< position<T> >(p+nanoParticleOffset,nNanoParticleElements);
		buf.clear();
		buf="nanoCOM_";
		buf+=name;
		buf+=".dat";
		dataFile.open(buf.c_str(), std::ios::app | std::ios::out);
		dataFile << (*System).readInitialTime() << '\t' << nanoCOM.x << '\t' << nanoCOM.y << '\t' << nanoCOM.z << std::endl;
		dataFile.close();
		
		position<T> systemCOM=com< position<T> >(p,(*System).readNParticles());
		buf.clear();
		buf="systemCOM_";
		buf+=name;
		buf+=".dat";
		dataFile.open(buf.c_str(), std::ios::app | std::ios::out);
		dataFile << (*System).readInitialTime() << '\t' << systemCOM.x << '\t' << systemCOM.y << '\t' << systemCOM.z << std::endl;
		dataFile.close();
	#endif
	
	if(aP!=NULL)
	{
		for(int i=0;i<(*System).readNParticles();i++)
		{
			int partition=(0.5*(v[i].x*v[i].x+v[i].y*v[i].y+v[i].z*v[i].z)/kEnergyDensityPartition);
			while(kEnergyDensity.size()<=partition)
				kEnergyDensity.push_back(0);
			kEnergyDensity[partition]++;
		}
	}
		
	//If you are interested in lateral diffusion, don't worry about whether you need 2 or 3 of
	// the dimensions in d.x, d.y, and d.z. Just change your dimension, d, in <r^2>=d*D*t.
	//If you see a contribution from another dimension, it means you haven't set up your system correctly.
	if(aPStart!=NULL && aP!=NULL)
	{
		
		buf.clear();
		buf="meanSquareDisplacement_";
		buf+=name;
		buf+=".dat";
		dataFile.open(buf.c_str(), std::ios::app | std::ios::out);
		dataFile << (*System).readInitialTime();
		
		for(int molIndex=0;molIndex<(*System).readNMolecules();molIndex++)
		{
			T msD=0;
			int nParticles=0;
			molecule<T, fourVector<int> > *m=(*System).getMolecule();
			
			//pick a structure by type
			switch(m[molIndex].readType())
			{
				case BOND:
				{
					for(int l=0;l<m[molIndex].readNBond();l++)
					{
						//These are the first and second particles of the bond
						int firstParticle=m[molIndex].getBonds()[l].s[0];
						int secondParticle=m[molIndex].getBonds()[l].s[1];
						
						//Compute the MSD from their start
						threeVector<T> d;
						d.x=aP[firstParticle].x-aPStart[firstParticle].x;
						d.y=aP[firstParticle].y-aPStart[firstParticle].y;
						d.z=aP[firstParticle].z-aPStart[firstParticle].z;
						msD+=(d.x*d.x+d.y*d.y+d.z*d.z);
						
						d.x=aP[secondParticle].x-aPStart[secondParticle].x;
						d.y=aP[secondParticle].y-aPStart[secondParticle].y;
						d.z=aP[secondParticle].z-aPStart[secondParticle].z;
						msD+=(d.x*d.x+d.y*d.y+d.z*d.z);
					}
					nParticles+=m[molIndex].readNBond()*2;
					break;
				}
				case BEND:
				{
					for(int l=0;l<m[molIndex].readNBond();l++)
					{
						//These are the first, second, and third particles of the bend
						int firstParticle=m[molIndex].getBonds()[l].s[0];
						int secondParticle=m[molIndex].getBonds()[l].s[1];
						int thirdParticle=m[molIndex].getBonds()[l].s[2];
						
						threeVector<T> d;
						d.x=aP[firstParticle].x-aPStart[firstParticle].x;
						d.y=aP[firstParticle].y-aPStart[firstParticle].y;
						d.z=aP[firstParticle].z-aPStart[firstParticle].z;
						msD+=(d.x*d.x+d.y*d.y+d.z*d.z);
						
						d.x=aP[secondParticle].x-aPStart[secondParticle].x;
						d.y=aP[secondParticle].y-aPStart[secondParticle].y;
						d.z=aP[secondParticle].z-aPStart[secondParticle].z;
						msD+=(d.x*d.x+d.y*d.y+d.z*d.z);
						
						d.x=aP[thirdParticle].x-aPStart[thirdParticle].x;
						d.y=aP[thirdParticle].y-aPStart[thirdParticle].y;
						d.z=aP[thirdParticle].z-aPStart[thirdParticle].z;
						msD+=(d.x*d.x+d.y*d.y+d.z*d.z);
						
					}
					nParticles+=m[molIndex].readNBond()*3;
					break;
				}
				case CHAIN:
				{
					
					for(int l=0;l<m[molIndex].readNBond();l++)
					{
						int start=m[molIndex].getBonds()[l].s[START];
						int length=m[molIndex].getBonds()[l].s[CHAINLENGTH];
						int nChains=m[molIndex].getBonds()[l].s[NCHAINS];
						for(int j=start;j<start+length*nChains;j++)
						{
							threeVector<T> d;
							d.x=aP[j].x-aPStart[j].x;
							d.y=aP[j].y-aPStart[j].y;
							d.z=aP[j].z-aPStart[j].z;
							msD+=(d.x*d.x+d.y*d.y+d.z*d.z);
							
						}
						nParticles+=length*nChains;
					}
					break;
				}
				case BEAD:
				{
					for(int l=0;l<m[molIndex].readNBond();l++)
					{
						//These are the first and second particles of the bond
						int firstParticle=m[molIndex].getBonds()[l].s[0];
						
						threeVector<T> d;
						d.x=aP[firstParticle].x-aPStart[firstParticle].x;
						d.y=aP[firstParticle].y-aPStart[firstParticle].y;
						d.z=aP[firstParticle].z-aPStart[firstParticle].z;
						msD+=(d.x*d.x+d.y*d.y+d.z*d.z);
					}
					nParticles+=m[molIndex].readNBond();
					break;
				}
				case FLOATING_BASE:
				{
					break;
				}
				case NANOCORE:
				{
					for(int l=0;l<m[molIndex].readNBond();l++)
					{
						//These are the first and second particles of the bond
						int firstParticle=m[molIndex].getBonds()[l].s[0];
						
						threeVector<T> d;
						d.x=aP[firstParticle].x-aPStart[firstParticle].x;
						d.y=aP[firstParticle].y-aPStart[firstParticle].y;
						d.z=aP[firstParticle].z-aPStart[firstParticle].z;
						msD+=(d.x*d.x+d.y*d.y+d.z*d.z);
					}
					nParticles+=m[molIndex].readNBond();
					break;
				}
				default:
				{
					//does nothing
					break;
				}
			}
			if(nParticles!=0)
				dataFile << '\t' << (msD/static_cast<T>(nParticles));
		}
		dataFile << '\n';
		dataFile.close();
		
		#ifdef SOLVENT_TRACER
			T msD=0;
			int nParticles=0;
			for(int i=0;i<(*System).readNParticles();i++)
			{
				if(aP[i].type==SOLVENT_TRACER)
				{
					threeVector<T> d;
					d.x=aP[i].x-aPStart[i].x;
					d.y=aP[i].y-aPStart[i].y;
					d.z=aP[i].z-aPStart[i].z;
					msD+=(d.x*d.x+d.y*d.y+d.z*d.z);
					nParticles++;
				}
			}
			
			if(nParticles!=0)
			{
				buf.clear();
				buf="MSD_1_";
				buf+=name;
				buf+=".dat";
				dataFile.open(buf.c_str(), std::ios::app | std::ios::out);
				dataFile << (*System).readInitialTime() << '\t' << (msD/static_cast<T>(nParticles)) << '\n';
				dataFile.close();
			}
		#endif
	}
}
