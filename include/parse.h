/**
 * \brief Parsing object. Hopefully this will replace some of this script and memory handling mess.
 * It should be inherited into another class, likely a class that handles memory.
 */
template <typename T>
class Parse {
	public:
		Parse();
		~Parse();
		void initialize();
		void reset();
		bool insertCommand();
		
	private:
		std::map<std::string,int> commandMap;
		bool *required;
		int nCommands;
		std::string *commands;
};

template <typename T>
Blob<T>::Blob()
{
	//set most of these to a default value
	p=NULL;
	v=NULL;
	a=NULL;
	m=NULL;
	twoBodyFconst=NULL;
	twoBodyUconst=NULL;
	twoBodyMFconst=NULL;
	twoBodyMUconst=NULL;
	threeBodyMFconst=NULL;
	threeBodyMUconst=NULL;
	seed=0;
	nTypes=0;
	nMolecules=0;
	nParticles=0;
	periodic.x=true;
	periodic.y=true;
	periodic.z=true;
	cutoff=0;
	size.x=0;
	size.y=0;
	size.z=0;
	initialTime=0;
	finalTime=0;
	deltaT=0;
	storeInterval=0;
	measureInterval=0;
	deltaLXY=0;
	removeSolvent=0;
	tempStepInterval=0;
	//Keep this comment here:
	//COMMANDINPUTHANDLER
	solventGamma=0;
	
	//list of commands this simulation utilizes
	nCommands=NENUMCOMMANDID;
	commands=new std::string[nCommands];
	required=new bool[nCommands];
	nAllocated=new int[nCommands];
	
	
	//these are the mappings of the command ID's strings
	commands[GAMMA]="gamma";
	commands[INITIALTEMP]="initialTemp";
	commands[FINALTEMP]="finalTemp";
	commands[SEED]="seed";
	commands[NTYPES]="nTypes";
	commands[NMOLECULES]="nMolecules";
	commands[NPARTICLES]="nParticles";
	commands[PERIODIC]="periodic";
	commands[CUTOFF]="cutoff";
	commands[SIZE]="size";
	commands[INITIALTIME]="initialTime";
	commands[FINALTIME]="finalTime";
	commands[DELTAT]="deltaT";
	commands[STOREINTERVAL]="storeInterval";
	commands[MEASUREINTERVAL]="measureInterval";
	//these are really required to be after the nWhatevers for output
	commands[TWOBODYFCONST]="twoBodyFconst";
	commands[TWOBODYUCONST]="twoBodyUconst";
	commands[TWOBODYMFCONST]="twoBodyMFconst";
	commands[TWOBODYMUCONST]="twoBodyMUconst";
	commands[THREEBODYMFCONST]="threeBodyMFconst";
	commands[THREEBODYMUCONST]="threeBodyMUconst";
	commands[POSITIONS]="positions";
	commands[VELOCITIES]="velocities";
	commands[MOLECULE]="molecule";
	//No more nWhatevers, I should regroup these
	commands[BANANA]="banana";
	commands[DELTALXY]="deltaLXY";
	commands[REMOVESOLVENT]="removeSolvent";
	commands[TEMPSTEPINTERVAL]="tempStepInterval";
	//Add more after this
	//Keep this comment here:
	//COMMANDINPUTHANDLER
	commands[SOLVENTGAMMA]="solventGamma";
	
	//map the commands to their indexes to make them easy to find
	for(int i=0;i<nCommands;i++)
	{
		commandMap[commands[i]]=i;
		required[i]=true;
		nAllocated[i]=0;
	}
	
	required[NMOLECULES]=false;
	required[MOLECULE]=false;
	required[BANANA]=false;
	required[TWOBODYFCONST]=false;
	required[TWOBODYUCONST]=false;
	required[TWOBODYMFCONST]=false;
	required[TWOBODYMUCONST]=false;
	required[THREEBODYMFCONST]=false;
	required[THREEBODYMUCONST]=false;
	required[DELTALXY]=false;
	required[REMOVESOLVENT]=false;
	required[TEMPSTEPINTERVAL]=false;
	//Keep this comment here:
	//COMMANDINPUTHANDLER
	required[SOLVENTGAMMA]=false;
}

//Mostly memory handling, deletes etc...
template <typename T>
Blob<T>::~Blob()
{
	if(nAllocated[POSITIONS]>0)
		delete p;
	if(nAllocated[VELOCITIES]>0)
		delete v;
	if(nAllocated[POSITIONS]>0)
		delete a;
	if(nAllocated[TWOBODYFCONST]>0)
		delete twoBodyFconst;
	if(nAllocated[TWOBODYUCONST]>0)
		delete twoBodyUconst;
	if(nAllocated[TWOBODYMFCONST]>0)
		delete twoBodyMFconst;
	if(nAllocated[TWOBODYMUCONST]>0)
		delete twoBodyMUconst;
	if(nAllocated[THREEBODYMFCONST]>0)
		delete threeBodyMFconst;
	if(nAllocated[THREEBODYMUCONST]>0)
		delete threeBodyMUconst;
	if(nAllocated[MOLECULE]>0)
		delete[] m;
	
	delete[] commands;//the mysterious dynamic deallocator of arrays [], actually, this calls each destructor properly
	delete required;
	delete nAllocated;
}

template <typename T>
void Blob<T>::errorChecking()
{
	for(int i=0;i<nParticles;i++)
	{
		if(p[i].x>size.x || p[i].x<0)
		{
			std::cout << "X position of particle " << i << " is out of bounds.\n";
			throw 0;
		}
		if(p[i].y>size.y || p[i].y<0)
		{
			std::cout << "Y position of particle " << i << " is out of bounds.\n";
			throw 0;
		}
		if(p[i].z>size.z || p[i].z<0)
		{
			std::cout << "Z position of particle " << i << " is out of bounds.\n";
			throw 0;
		}
	}
	
	for(int i=0;i<nMolecules;i++)
	{
		int nBonded=m[i].readNBond();
		fourVector<int> *bond=m[i].getBonds();
		int type=m[i].readType();
		
		for(int j=0;j<nBonded;j++)
		{
			switch(type)
			{
				case BOND:
					if(bond[j].s[0]>nParticles || bond[j].s[1]>nParticles || bond[j].s[0]<0 || bond[j].s[1]<0)
					{
						std::cout << "BOND Molecule " << i << ", bond " << j << " is out of bounds!\n";
						throw 0;
					}
					break;
				case BEND:
					if(bond[j].s[0]>nParticles || bond[j].s[1]>nParticles || 
						bond[j].s[2]>nParticles || bond[j].s[0]<0 ||
						bond[j].s[1]<0 || bond[j].s[2]<0)
					{
						std::cout << "BEND Molecule " << i << ", bond " << j << " is out of bounds!\n";
						throw 0;
					}
					break;
				case CHAIN:
					if(bond[j].s[START]+bond[j].s[NCHAINS]*bond[j].s[LENGTH]>nParticles ||
						bond[j].s[START]+bond[j].s[NCHAINS]*bond[j].s[LENGTH]<0)
					{
						std::cout << "CHAIN Molecule " << i << ", bond " << j << " is out of bounds!\n";
						throw 0;
					}
					break;
				default:
					//does nothing
					break;
			}
		}
	}
}

template <typename T>
void Blob<T>::reset()
{
	wP=-1;
	it=0;
}

template <typename T>
bool Blob<T>::input(std::string in)
{
	if(wP==-1)
	{
		it=0;
		//this checks for end
		if(in==END_INPUT)
		{
			for(int i=0;i<nCommands;i++)
			{
				if((nAllocated[i]<=0) && required[i])
				{
					#ifdef WARNINGS_ENABLED
						std::cout << "Warning (Blob): Script missing " << commands[i] << "!\n";
						//std::cout << "Press any key to continue loading or ctrl-z to stop!";
						//std::cin.get();
					#endif
					///Maybe in the future I should add warning and error levels for crap like this
					//throw 0;
				}
				//go on as if nothing happend
			}
			this->errorChecking();
			return true;//stream is complete
		}
		
		std::map<std::string,int>::iterator found;
		found=commandMap.find(in);
		if(found!=commandMap.end())
		{
			wP=found->second;//save the command index
		}
		else
		{
			std::cout << in << " is not a recognized command!\n";
			std::cout << "Locate command before " << in << "!\n";
			std::cout << "You probably have too many positions, velocities, or molecules!\n";
			std::cout << "Also, check nMolecules and nParticles.\n";
			throw 0;
		}
		//return false;//just recieved command, wait until next input to get parameter
	}
	
	//generate a string stream to extract the data from the string safely
	std::stringstream reformatString(in);
	
	//get the values for the command
	switch(wP)
	{
		case GAMMA://gamma
			if(it>0)
			{
				reformatString >> gamma;
				if(reformatString.fail())
				{
					std::cout << "Value after gamma is not a numerical type!\n";
					std::cout << "Format:\ngamma [float]\n";
					throw 0;
				}
				nAllocated[wP]=1;
				wP=-1;//clear command index
			}
			break;
			
		case INITIALTEMP://initialTemp
			if(it>0)
			{
				reformatString >> initialTemp;
				if(reformatString.fail())
				{
					std::cout << "Value after initialTemp is not a numerical type!\n";
					std::cout << "Format:\ninitialTemp [float]\n";
					throw 0;
				}
				nAllocated[wP]=1;
				wP=-1;
			}
			break;
			
		case FINALTEMP://finalTemp
			if(it>0)
			{
				reformatString >> finalTemp;
				if(reformatString.fail())
				{
					std::cout << "Value after finalTemp is not a numerical type!\n";
					std::cout << "Format:\nfinalTemp [float]\n";
					throw 0;
				}
				nAllocated[wP]=1;
				wP=-1;
			}
			break;
			
		case SEED://seed
			if(it>0)
			{
				reformatString >> seed;
				if(reformatString.fail())
				{
					std::cout << "Value after seed is not a numerical type!\n";
					std::cout << "Format:\nseed [integer]\n";
					throw 0;
				}
				nAllocated[wP]=1;
				wP=-1;
			}
			break;
			
		case NTYPES://nTypes
			if(it>0)
			{
				reformatString >> nTypes;
				if(reformatString.fail())
				{
					std::cout << "Value after nTypes is not a numerical type!\n";
					std::cout << "Format:\nnTypes [integer]\n";
					throw 0;
				}
				nAllocated[wP]=1;
				wP=-1;
			}
			break;
			
		case NMOLECULES://nMolecules
			if(it>0)
			{
				reformatString >> nMolecules;
				if(reformatString.fail())
				{
					std::cout << "Value after nMolecules is not a numerical type!\n";
					std::cout << "Format:\nnMolecules [integer]\n";
					throw 0;
				}
				nAllocated[wP]=1;
				wP=-1;
			}
			break;
			
		case NPARTICLES://nParticles
			if(it>0)
			{
				reformatString >> nParticles;
				if(reformatString.fail())
				{
					std::cout << "Value after nParticles is not a numerical type!\n";
					std::cout << "Format:\nnParticles [integer]\n";
					throw 0;
				}
				nAllocated[wP]=1;
				wP=-1;
			}
			break;
			
		case PERIODIC://periodic
			if(it>0)
			{
				//this doesn't work anyway
				//reformatString >> nMolecules;
				//if(reformatString.fail())
				//{
				//	std::cout << "Value after nMolecules is not a numerical type!\n";
				//	throw 0;
				//}
				periodic.x=true;
				periodic.y=true;
				periodic.z=true;
				nAllocated[wP]=1;
				wP=-1;
			}
			break;
			
		case CUTOFF://cutoff
			if(it>0)
			{
				reformatString >> cutoff;
				if(reformatString.fail())
				{
					std::cout << "Value after cutoff is not a numerical type!\n";
					std::cout << "Format:\ncutoff [float]\n";
					throw 0;
				}
				cutoffSquared=cutoff*cutoff;
				nAllocated[wP]=1;
				wP=-1;
			}
			break;
			
		case SIZE://size
			if(it>0)
			{
				//this uses the it iterative variable
				//also, this demonstrates the union concept, size.s[0]==size.x
				reformatString >> size.s[it-1];
				if(reformatString.fail())
				{
					std::cout << "Value after size is not a numerical type!\n";
					std::cout << "Format:\nsize [float] [float] [float]\n";
					throw 0;
				}
				if(it==3)
				{
					nAllocated[wP]=1;
					wP=-1;
				}
			}
			break;
			
		case INITIALTIME://initialTime
			if(it>0)
			{
				reformatString >> initialTime;
				if(reformatString.fail())
				{
					std::cout << "Value after initialTime is not a numerical type!\n";
					std::cout << "Format:\ninitialTime [float]\n";
					throw 0;
				}
				nAllocated[wP]=1;
				wP=-1;
			}
			break;
			
		case FINALTIME://finalTime
			if(it>0)
			{
				reformatString >> finalTime;
				if(reformatString.fail())
				{
					std::cout << "Value after finalTime is not a numerical type!\n";
					std::cout << "Format:\nfinalTime [float]\n";
					throw 0;
				}
				nAllocated[wP]=1;
				wP=-1;
			}
			break;
			
		case DELTAT://deltaT
			if(it>0)
			{
				reformatString >> deltaT;
				if(reformatString.fail())
				{
					std::cout << "Value after deltaT is not a numerical type!\n";
					std::cout << "Format:\ndeltaT [float]\n";
					throw 0;
				}
				nAllocated[wP]=1;
				wP=-1;
			}
			break;
			
		case STOREINTERVAL://storeInterval
			if(it>0)
			{
				reformatString >> storeInterval;
				if(reformatString.fail())
				{
					std::cout << "Value after storeInterval is not a numerical type!\n";
					std::cout << "Format:\nstoreInterval [float]\n";
					throw 0;
				}
				nAllocated[wP]=1;
				wP=-1;
			}
			break;
			
		case MEASUREINTERVAL://measureInterval
			if(it>0)
			{
				reformatString >> measureInterval;
				if(reformatString.fail())
				{
					std::cout << "Value after measureInterval is not a numerical type!\n";
					std::cout << "Format:\nmeasureInterval [float]\n";
					throw 0;
				}
				nAllocated[wP]=1;
				wP=-1;
			}
			break;
			
		case TWOBODYFCONST://twoBodyFconst
			if(it>0)//more advanced "it" usage
			{
				reformatString >> twoBodyFconst[it-1];
				if(reformatString.fail())//you broke it...
				{
					std::cout << "Value after twoBodyFconst is not a numerical type!\n";
					std::cout << "Check that nTypes is correct!\n";
					std::cout << "Format:\ntwoBodyFconst\n [float]...\n";
					throw 0;
				}
				if(it==nTWOBODYFCONST*nTypes*nTypes)//last number
				{
					wP=-1;
				}
			}
			else
			{
				if(nAllocated[NTYPES]<1)//check for the prescence of types
				{
					std::cout << "nTypes was not present before twoBodyFconst!\n";
					throw 0;
				}
				if(nAllocated[wP]==0)//check for any allocation
				{
					twoBodyFconst=new T[nTWOBODYFCONST*nTypes*nTypes];
					nAllocated[wP]=nTWOBODYFCONST*nTypes*nTypes;
				}
				if(nAllocated[wP]==-1)//check if allocation has changed
				{
					delete twoBodyFconst;
					twoBodyFconst=new T[nTWOBODYFCONST*nTypes*nTypes];
					nAllocated[wP]=nTWOBODYFCONST*nTypes*nTypes;
				}
			}
			break;
			
		case TWOBODYUCONST://twoBodyUconst
			if(it>0)
			{
				reformatString >> twoBodyUconst[it-1];
				if(reformatString.fail())
				{
					std::cout << "Value after twoBodyUconst is not a numerical type!\n";
					std::cout << "Check that nTypes is correct!\n";
					std::cout << "Format:\ntwoBodyUconst\n [float]...\n";
					throw 0;
				}
				if(it==nTWOBODYUCONST*nTypes*nTypes)
				{
					wP=-1;
				}
			}
			else
			{
				if(nAllocated[NTYPES]<1)
				{
					std::cout << "nTypes was not present before twoBodyUconst!\n";
					throw 0;
				}
				if(nAllocated[wP]==0)//check for any allocation
				{
					twoBodyUconst=new T[nTWOBODYUCONST*nTypes*nTypes];
					nAllocated[wP]=nTWOBODYUCONST*nTypes*nTypes;
				}
				if(nAllocated[wP]==-1)//check if allocation has changed
				{
					delete twoBodyUconst;
					twoBodyUconst=new T[nTWOBODYUCONST*nTypes*nTypes];
					nAllocated[wP]=nTWOBODYUCONST*nTypes*nTypes;
				}
			}
			break;
			
		case TWOBODYMFCONST://twoBodyMFconst
			if(it>0)
			{
				reformatString >> twoBodyMFconst[it-1];
				if(reformatString.fail())
				{
					std::cout << "Value after twoBodyMFconst is not a numerical type!\n";
					std::cout << "Check that nTypes is correct!\n";
					std::cout << "Format:\ntwoBodyMFconst\n [float]...\n";
					throw 0;
				}
				if(it==nTWOBODYMFCONST*nTypes*nTypes)
				{
					wP=-1;
				}
			}
			else
			{
				if(nAllocated[NTYPES]<1)
				{
					std::cout << "nTypes was not present before twoBodyMFconst!\n";
					throw 0;
				}
				if(nAllocated[wP]==0)//check for any allocation
				{
					twoBodyMFconst=new T[nTWOBODYMFCONST*nTypes*nTypes];
					nAllocated[wP]=nTWOBODYMFCONST*nTypes*nTypes;
				}
				if(nAllocated[wP]==-1)//check if allocation has changed
				{
					delete twoBodyMFconst;
					twoBodyMFconst=new T[nTWOBODYMFCONST*nTypes*nTypes];
					nAllocated[wP]=nTWOBODYMFCONST*nTypes*nTypes;
				}
			}
			break;
			
		case TWOBODYMUCONST://twoBodyMUconst
			if(it>0)
			{
				reformatString >> twoBodyMUconst[it-1];
				if(reformatString.fail())
				{
					std::cout << "Value after twoBodyMUconst is not a numerical type!\n";
					std::cout << "Check that nTypes is correct!\n";
					std::cout << "Format:\ntwoBodyMUconst\n [float]...\n";
					throw 0;
				}
				if(it==nTWOBODYMUCONST*nTypes*nTypes)
				{
					wP=-1;
				}
			}
			else
			{
				if(nAllocated[NTYPES]<1)
				{
					std::cout << "nTypes was not present before twoBodyMUconst!\n";
					throw 0;
				}
				if(nAllocated[wP]==0)//check for any allocation
				{
					twoBodyMUconst=new T[nTWOBODYMUCONST*nTypes*nTypes];
					nAllocated[wP]=nTWOBODYMUCONST*nTypes*nTypes;
				}
				if(nAllocated[wP]==-1)//check if allocation has changed
				{
					delete twoBodyMUconst;
					twoBodyMUconst=new T[nTWOBODYMUCONST*nTypes*nTypes];
					nAllocated[wP]=nTWOBODYMUCONST*nTypes*nTypes;
				}
			}
			break;
			
		case THREEBODYMFCONST://threeBodyMFconst
			if(it>0)
			{
				reformatString >> threeBodyMFconst[it-1];
				if(reformatString.fail())
				{
					std::cout << "Value after threeBodyMFconst is not a numerical type!\n";
					std::cout << "Check that nTypes is correct!\n";
					std::cout << "Format:\nthreeBodyMFconst\n [float]...\n";
					throw 0;
				}
				if(it==nTHREEBODYMFCONST*nTypes*nTypes*nTypes)
				{
					wP=-1;
				}
			}
			else
			{
				if(nAllocated[NTYPES]<1)
				{
					std::cout << "nTypes was not present before threeBodyMFconst!\n";
					throw 0;
				}
				if(nAllocated[wP]==0)//check for any allocation
				{
					threeBodyMFconst=new T[nTHREEBODYMFCONST*nTypes*nTypes*nTypes];
					nAllocated[wP]=nTHREEBODYMFCONST*nTypes*nTypes*nTypes;
				}
				if(nAllocated[wP]==-1)//check if allocation has changed
				{
					delete threeBodyMFconst;
					threeBodyMFconst=new T[nTHREEBODYMFCONST*nTypes*nTypes*nTypes];
					nAllocated[wP]=nTHREEBODYMFCONST*nTypes*nTypes*nTypes;
				}
			}
			break;
			
		case THREEBODYMUCONST://threeBodyMUconst
			if(it>0)
			{
				//Extract
				reformatString >> threeBodyMUconst[it-1];
				//Test
				if(reformatString.fail())
				{
					std::cout << "Value after threeBodyMUconst is not a numerical type!\n";
					std::cout << "Check that nTypes is correct!\n";
					std::cout << "Format:\ntwoBodyMUconst\n [float]...\n";
					throw 0;
				}
				//Exit upon end of stream condition
				if(it==nTHREEBODYMUCONST*nTypes*nTypes*nTypes)
				{
					wP=-1;
				}
			}
			else//if iterator is 0 or less
			{
				//Dependencies
				if(nAllocated[NTYPES]<1)
				{
					std::cout << "nTypes was not present before threeBodyMUconst!\n";
					throw 0;
				}
				//Allocation
				if(nAllocated[wP]==0)//check for any allocation
				{
					threeBodyMUconst=new T[nTHREEBODYMUCONST*nTypes*nTypes*nTypes];
					nAllocated[wP]=nTHREEBODYMUCONST*nTypes*nTypes*nTypes;
				}
				if(nAllocated[wP]==-1)//check if allocation has changed
				{
					delete threeBodyMUconst;
					threeBodyMUconst=new T[nTHREEBODYMUCONST*nTypes*nTypes*nTypes];
					nAllocated[wP]=nTHREEBODYMUCONST*nTypes*nTypes*nTypes;
				}
			}
			break;
			
		case POSITIONS://positions
			if(it>0)
			{
				int particle=int((it-1)/4);
				if((it-1)%4==0)//these are the only values
					//type may change to a character string, but that is messy...
					reformatString >> p[particle].type;
				else
					reformatString >> p[particle].s[(it-1)%4-1];
				
				if(reformatString.fail())
				{
					std::cout << "Value after positions is not a numerical type!\n";
					std::cout << "This likely means that nParticles is larger than the number of positions in file.\n";
					std::cout << "Format:\npositions\n [integer] [float] [float] [float]\n...\n";
					throw 0;
				}
				if(it==nParticles*4)
				{
					wP=-1;
				}
			}
			else
			{
				if(nAllocated[NPARTICLES]<1)
				{
					std::cout << "nParticles was not present before positions!\n";
					throw 0;
				}
				if(nAllocated[wP]==0)//check for any allocation
				{
					p=new position<T>[nParticles];
					a=new threeVector<T>[nParticles];
					nAllocated[wP]=nParticles;
				}
				if(nAllocated[wP]==-1)//check if allocation has changed
				{
					delete p,a;
					p=new position<T>[nParticles];
					a=new threeVector<T>[nParticles];
					nAllocated[wP]=nParticles;
				}
			}
			break;
			
		case VELOCITIES://velocities
			if(it>0)
			{
				int particle=int((it-1)/3);
				reformatString >> v[particle].s[(it-1)%3];
				if(reformatString.fail())
				{
					std::cout << "Value after velocities is not a numerical type!\n";
					std::cout << "This likely means that nParticles is larger than the number of velocities in file.\n";
					std::cout << "Format:\nvelocities\n [float] [float] [float]\n...\n";
					throw 0;
				}
				if(it==nParticles*3)
				{
					wP=-1;
				}
			}
			else
			{
				if(nAllocated[NPARTICLES]<1)
				{
					std::cout << "nParticles was not present before velocities!\n";
					throw 0;
				}
				if(nAllocated[wP]==0)//check for any allocation
				{
					v=new threeVector<T>[nParticles];
					nAllocated[wP]=nParticles;
				}
				if(nAllocated[wP]==-1)//check if allocation has changed
				{
					delete v;
					v=new threeVector<T>[nParticles];
					nAllocated[wP]=nParticles;
				}
			}
			break;
			
		case MOLECULE://molecules
			if(it>0)
			{
				if(mol!=current)
				{
					int buf;
					reformatString >> buf;
					m[mol].setType(buf);
					molState=0;
					bondState=0;
					bondNumber=0;
					current=mol;
				}
				else
				{
					int nMoleculeParticles, nMoleculeConstants;
					switch(m[mol].readType())
					{
						//I seriously need to implement the next 2
						case TORSION:
							//m[mol].allocConstant(0);
							nMoleculeConstants=0;
							nMoleculeParticles=0;
							break;
						case DIHEDRAL:
							//m[mol].allocConstant(0);
							nMoleculeConstants=0;
							nMoleculeParticles=0;
							break;
						case BOND:
							nMoleculeConstants=nBONDCONST;
							nMoleculeParticles=2;//2 particles per bond
							break;
						case BEND:
							nMoleculeConstants=nBONDCONST;
							nMoleculeParticles=3;//3 particles per bond
							break;
						case CHAIN:
							nMoleculeConstants=nCHAINCONST;
							nMoleculeParticles=3;//3 index relations
							break;
						default://default formatting
							//m[mol].allocConstant(0);
							nMoleculeConstants=0;
							nMoleculeParticles=0;
							#ifdef WARNINGS_ENABLED
								std::cout << "Warning (Blob): default formatting while reading molecule!\n";
							#endif
							break;
					}
					if(nMoleculeParticles==0)
					{
						int buf;
						reformatString >> buf;
						bondBuf.s[bondState++]=buf;
						if(bondState>=4)//4 elements in the default
							m[mol++].addBond(bondBuf);
					}
					else
					{
						if(molState==0)//nBonds
						{
							int buf;
							reformatString >> buf;
							m[mol].allocBonds(buf);
						}
						else
						{
							if(molState<nMoleculeConstants+1)//constants for bond
							{
								double buf;
								reformatString >> buf;
								m[mol].addConstant(buf);
							}
							else //if(molState>m[mol].readNConstantsAlloc()+1)//bonds themselves
							{
								int buf;
								reformatString >> buf;
								bondBuf.s[bondState++]=buf;
								if(bondState>=nMoleculeParticles)//read a bond component
								{
									m[mol].addBond(bondBuf);
									bondState=0;
									bondNumber++;
									if(bondNumber>=m[mol].readNBondAlloc())
										mol++;
								}
							}
						}
					}
					molState++;
				}
				
				if(reformatString.fail())
				{
					std::cout << "Value after molecules is not a numerical type!\n";
					std::cout << "Read: " << in << '\n';
					throw 0;
				}
				
				if(mol==nMolecules)
				{
					wP=-1;
				}
			}
			else
			{
				if(nAllocated[NMOLECULES]<1)
				{
					std::cout << "nMolecules was not present before molecules!\n";
					throw 0;
				}
				if(nAllocated[MOLECULE]==0)
				{
					m=new molecule<T,fourVector<int> >[nMolecules];
					nAllocated[MOLECULE]=nMolecules;
				}
				if(nAllocated[MOLECULE]==-1)
				{
					delete[] m;
					m=new molecule<T,fourVector<int> >[nMolecules];
					nAllocated[MOLECULE]=nMolecules;
				}
				current=-1;
				mol=0;//set number of mols to 0
			}
			break;
			
		case BANANA://banana!
			//this is an example of a flag that can be set, notice the distinct lack of an iterator
			std::cout << "Monkeys must have messed with your script because I found a banana!\n";
			nAllocated[wP]=1;
			wP=-1;
			break;
			
		case DELTALXY://deltaLXY
			if(it>0)
			{
				//this uses the it iterative variable
				//also, this demonstrates the union concept, size.s[0]==size.x
				reformatString >> deltaLXY;
				if(reformatString.fail())
				{
					std::cout << "Value after deltaLXY is not a numerical type!\n";
					std::cout << "Format:\n\tdeltaLXY [float]\n";
					throw 0;
				}
				nAllocated[wP]=1;
				wP=-1;
			}
			break;
		
		case REMOVESOLVENT:
			if(it>0)
			{
				reformatString >> removeSolvent;
				if(reformatString.fail())
				{
					std::cout << "Value after removeSolvent is not a numerical type!\n";
					std::cout << "Format:\n\tremoveSolvent [float]\n";
					std::cout << "Float is the fraction of the total solvent, not just the inner volume of solvent\n";
					throw 0;
				}
				nAllocated[wP]=1;
				wP=-1;
			}
			break;
			
		case TEMPSTEPINTERVAL:
			if(it>0)
			{
				reformatString >> tempStepInterval;
				if(reformatString.fail())
				{
					std::cout << "Value after tempStepInterval is not a numerical type!\n";
					std::cout << "Format:\ntempStepInterval [float]\n";
					throw 0;
				}
				nAllocated[wP]=1;
				wP=-1;
			}
			break;
		
		//Keep this comment here:
		//COMMANDINPUTHANDLER
		case SOLVENTGAMMA:
			if(it>0)
			{
				reformatString >> solventGamma;
				if(reformatString.fail())
				{
					std::cout << "Value after solventGamma is not a numerical type!\n";
					std::cout << "Format:\nsolventGamma [float]\n";
					throw 0;
				}
				nAllocated[wP]=1;
				wP=-1;
			}
			break;
			
		default://other commands
			std::cout << "This command doesn't yet exist in the input() handler.\n";
			throw 0;
			break;
	}
	it++;//increment iterator
	return false;
}

template <typename T>
std::string Blob<T>::output()
{
	int particle;
	if(wP==-1)
		wP=0;
	if(wP==nCommands)
		return END_INPUT;
	
	std::stringstream out("");
	if(it==0)
	{
		///This will enable the command to output with or without it being required
		if(nAllocated[wP]>0)
		{
			if(wP!=0)
				out << '\n';
			out << commands[wP];
			it++;
		}
		else
		{
			it=0;
			wP++;
		}
	}
	else
	{
		switch(wP)
		{
			case GAMMA://gamma
			out << std::setprecision(15) << gamma;
			wP++;
			it=0;
			break;
			
			case INITIALTEMP://initialTemp
			out << std::setprecision(15) << initialTemp;
			wP++;
			it=0;
			break;
			
			case FINALTEMP://finalTemp
			out << std::setprecision(15) << finalTemp;
			wP++;
			it=0;
			break;
			
			case SEED://seed
			out << std::setprecision(15) << seed;
			wP++;
			it=0;
			break;
			
			case NTYPES://nTypes
			out << std::setprecision(15) << nTypes;
			wP++;
			it=0;
			break;
			
			case NMOLECULES://nMolecules
			out << std::setprecision(15) << nMolecules;
			wP++;
			it=0;
			if(nMolecules>0)
				required[23]=true;
			break;
			
			case NPARTICLES://nParticles
			out << std::setprecision(15) << nParticles;
			wP++;
			it=0;
			break;
			
			case PERIODIC://periodic
			out << std::setprecision(15) << periodic.x;
			wP++;
			it=0;
			break;
			
			case CUTOFF://cutoff
			out << std::setprecision(15) << cutoff;
			wP++;
			it=0;
			break;
			
			case SIZE://size
			out << std::setprecision(15) << size.s[it-1];
			it++;
			if(it-1==3)
			{
				wP++;
				it=0;
			}
			break;
			
			case INITIALTIME://initialTime
			out << std::setprecision(15) << initialTime;
			wP++;
			it=0;
			break;
			
			case FINALTIME://finalTime
			out << std::setprecision(15) << finalTime;
			wP++;
			it=0;
			break;
			
			case DELTAT://deltaT
			out << std::setprecision(15) << deltaT;
			wP++;
			it=0;
			break;
			
			case STOREINTERVAL://storeInterval
			out << std::setprecision(15) << storeInterval;
			wP++;
			it=0;
			break;
			
			case MEASUREINTERVAL://measureInterval
			out << std::setprecision(15) << measureInterval;
			wP++;
			it=0;
			break;
			
			case TWOBODYFCONST://twoBodyFconst
			if(it==1)
				out << std::setprecision(15) << "\n ";
			out << std::setprecision(15) << twoBodyFconst[it-1];
			
			if((it)%nTWOBODYFCONST==0)
				out << std::setprecision(15) << "\n";
			
			it++;
			if(it-1==nTWOBODYFCONST*nTypes*nTypes)
			{
				wP++;
				it=0;
			}
			break;
			
			case TWOBODYUCONST://twoBodyUconst
			if(it==1)
				out << std::setprecision(15) << "\n ";
			out << std::setprecision(15) << twoBodyUconst[it-1];
			if((it)%nTWOBODYUCONST==0)
				out << std::setprecision(15) << "\n";
			it++;
			if(it-1==nTWOBODYUCONST*nTypes*nTypes)
			{
				wP++;
				it=0;
			}
			break;
			
			case TWOBODYMFCONST://twoBodyMFconst
			if(it==1)
				out << std::setprecision(15) << "\n "; 
			out << std::setprecision(15) << twoBodyMFconst[it-1];
			if((it)%nTWOBODYMFCONST==0)
				out << std::setprecision(15) << "\n";
			it++;
			if(it-1==nTWOBODYMFCONST*nTypes*nTypes)
			{
				wP++;
				it=0;
			}
			break;
			
			case TWOBODYMUCONST://twoBodyMUconst
			if(it==1)
				out << std::setprecision(15) << "\n "; 
			out << std::setprecision(15) << twoBodyMUconst[it-1];
			if((it)%nTWOBODYMUCONST==0)
				out << std::setprecision(15) << "\n";
			it++;
			if((it-1)==nTWOBODYMUCONST*nTypes*nTypes)
			{
				wP++;
				it=0;
			}
			break;
			
			case THREEBODYMFCONST://threeBodyMFconst
			if(it==1)
				out << std::setprecision(15) << "\n ";
			out << std::setprecision(15) << threeBodyMFconst[it-1];
			if((it)%nTHREEBODYMFCONST==0)
				out << std::setprecision(15) << "\n";
			it++;
			if((it-1)==nTHREEBODYMFCONST*nTypes*nTypes*nTypes)
			{
				wP++;
				it=0;
			}
			break;
			
			case THREEBODYMUCONST://threeBodyMUconst
			if(it==1)
				out << std::setprecision(15) << "\n "; 
			out << std::setprecision(15) << threeBodyMUconst[it-1];
			if((it)%nTHREEBODYMUCONST==0)
				out << std::setprecision(15) << "\n";
			it++;
			if((it-1)==nTHREEBODYMUCONST*nTypes*nTypes*nTypes)
			{
				wP++;
				it=0;
			}
			break;
			
			case POSITIONS://positions
			{
				int particle=(it-1)/4;
				if(it==1)
					out << std::setprecision(15) << "\n ";
				if((it-1)%4==0)
					out << std::setprecision(15) << p[particle].type;
				else
					out << std::setprecision(15) << p[particle].s[(it-1)%4-1];
				if((it-1)%4==3)
					out << std::setprecision(15) << '\n';
				it++;
				if((it-1)==nParticles*4)
				{
					wP++;
					it=0;
				}
				
			}
			break;
			
			case VELOCITIES://velocities
			{
				int particle=(it-1)/3;
				if(it==1)
					out << std::setprecision(15) << "\n ";
				out << std::setprecision(15) << v[particle].s[(it-1)%3];
				if((it-1)%3==2)
					out << std::setprecision(15) << '\n';
				it++;
				if((it-1)==nParticles*3)
				{
					wP++;
					it=0;
				}
				
				
			}
			break;
			
			case MOLECULE://molecules
			if(it==1)
			{
				out << std::setprecision(15) << '\n';
				mol=0;
				current=-1;
			}
			if(it>0)
			{
				if(mol!=current)
				{
					out << std::setprecision(15) << m[mol].readType() << '\t';
					molState=0;
					bondState=0;
					bondNumber=0;
					current=mol;
				}
				else
				{
					int nMoleculeParticles, nMoleculeConstants;
					switch(m[mol].readType())
					{
						//I seriously need to implement the next 2
						case TORSION:
							nMoleculeParticles=0;//4;
							nMoleculeConstants=0;
							break;
						case DIHEDRAL:
							nMoleculeParticles=0;//4;?
							nMoleculeConstants=0;
							break;
						case BOND:
							nMoleculeParticles=2;//2 particles per bond
							nMoleculeConstants=nBONDCONST;
							break;
						case BEND:
							nMoleculeParticles=3;//3 particles per bond
							nMoleculeConstants=nBENDCONST;
							break;
						case CHAIN:
							nMoleculeParticles=3;//3 index relations
							nMoleculeConstants=nCHAINCONST;
							break;
						default://default formatting
							//m[mol].allocConstant(0);
							nMoleculeParticles=0;
							#ifdef WARNINGS_ENABLED
							std::cout << "Warning (Blob): default formatting while writing molecule!\n";
							#endif
							break;
					}
					if(nMoleculeParticles==0)//Default
					{
						fourVector<int> *bond=m[mol].getBonds();
						out << std::setprecision(15) << bond[0].s[bondState++] << ' ';
						if(bondState>=4)//4 elements in the default
							mol++;
					}
					else
					{
						if(molState==0)//nBonds
						{
							out << std::setprecision(15) << m[mol].readNBond() << '\n';
						}
						else
						{
							if(molState<nMoleculeConstants+1)//constants for bond
							{
								T *constant=m[mol].getConstants();
								out << std::setprecision(15) << constant[molState-1] << ' ';
								if(molState==nMoleculeConstants)
									out << std::setprecision(15) << '\n';
							}
							else //if(molState>m[mol].readNConstantsAlloc()+1)//the bonds themselves
							{
								fourVector<int> *bond=m[mol].getBonds();
								out << std::setprecision(15) << bond[bondNumber].s[bondState++] << ' ';//particle index of bond
								if(bondState>=nMoleculeParticles)
								{
									out << std::setprecision(15) << '\n';
									bondState=0;
									bondNumber++;
									if(bondNumber>=m[mol].readNBond())
										mol++;
								}
							}
						}
					}
					molState++;
				}
			}
			else
			{
				if(nAllocated[NMOLECULES]<1)
				{
					std::cout << "nMolecules was not present before molecules!\n";
					throw 0;
				}
				if(nAllocated[MOLECULE]==0)
				{
					m=new molecule<T,fourVector<int> >[nMolecules];
					nAllocated[MOLECULE]=nMolecules;
				}
				if(nAllocated[MOLECULE]==-1)
				{
					delete[] m;
					m=new molecule<T,fourVector<int> >[nMolecules];
					nAllocated[MOLECULE]=nMolecules;
				}
				current=-1;
				mol=0;//set number of mols to 0
			}
			it++;
			if(mol==nMolecules)
			{
				wP++;
				it=0;
			}
			break;
			
			case BANANA://banana!
			it=0;
			wP++;
			out << std::setprecision(15) << "\n";
			break;
			
			case DELTALXY://deltaLXY
			out << std::setprecision(15) << deltaLXY;
			wP++;
			it=0;
			break;
			
			case REMOVESOLVENT://remove solvent flag
			out << std::setprecision(15) << removeSolvent;
			wP++;
			it=0;
			break;
			
			case TEMPSTEPINTERVAL://deltaLXY
			out << std::setprecision(15) << tempStepInterval;
			wP++;
			it=0;
			break;
			
			
			//Keep this comment here:
			//COMMANDINPUTHANDLER
			case SOLVENTGAMMA://deltaLXY
			out << std::setprecision(15) << solventGamma;
			wP++;
			it=0;
			break;
			
			
			default://other commands
			std::cout << "This command doesn't yet exist in the output() handler.\n";
			throw 0;
			//it=0;
			//wP++;
			//return "\n";
			break;
		}
	}
	return out.str();
}


///Nonbonded forces


///These are just conveniences, so you don't have to worry about typing the same crap over and over
template <typename T>
inline void Blob<T>::doForce(int i, int j, threeVector<T> &minImg)
{
	threeVector<T> d;
	
	d.x=p[i].x-p[j].x+minImg.x;
	d.y=p[i].y-p[j].y+minImg.y;
	d.z=p[i].z-p[j].z+minImg.z;
	
	int cindex=nTWOBODYFCONST*((p[i].type*nTypes)+p[j].type);
	
	this->nonBondedF(d,a[i],a[j],&twoBodyFconst[cindex]);
}

template <typename T>
T Blob<T>::doPotential(int i, int j, threeVector<T> &minImg)
{
	threeVector<T> d;
	
	d.x=p[i].x-p[j].x+minImg.x;
	d.y=p[i].y-p[j].y+minImg.y;
	d.z=p[i].z-p[j].z+minImg.z;
	
	int cindex=nTWOBODYUCONST*((p[i].type*nTypes)+p[j].type);
	
	return this->nonBondedP(d,&twoBodyUconst[cindex]);
}

//These are molecule related things
template <typename T>
void Blob<T>::doChainForce(int i)
{
	//this is just going to check it without an error condition, it ignores it
	if(m[i].readType()==CHAIN)
	{
		fourVector<int> *bond=m[i].getBonds();
		int first, second, third;
		threeVector<T> da,db;
		//Go through all bond descriptions
		for(int j=0;j<m[i].readNBond();j++)
		{
			//bond info
			int start=bond[j].s[START];
			int nChains=bond[j].s[NCHAINS];
			int length=bond[j].s[LENGTH];
			//go through all chain lengths
			for(int k=start; k<start+length*nChains; k+=length)
			{
				//go through all particles in chain
				for(int l=k;l<k+length-3;l++)
				{
					int first=l;
					int second=l+1;
					int third=l+2;
					
					da.x=p[first].x-p[second].x;
					if(da.x>size.x/2.0) da.x-=size.x;
					if(da.x<-size.x/2.0) da.x+=size.x;
					da.y=p[first].y-p[second].y;
					if(da.y>size.y/2.0) da.y-=size.y;
					if(da.y<-size.y/2.0) da.y+=size.y;
					da.z=p[first].z-p[second].z;
					if(da.z>size.z/2.0) da.z-=size.z;
					if(da.z<-size.z/2.0) da.z+=size.z;
					
					db.x=p[second].x-p[third].x;
					if(db.x>size.x/2.0) db.x-=size.x;
					if(db.x<-size.x/2.0) db.x+=size.x;
					db.y=p[second].y-p[third].y;
					if(db.y>size.y/2.0) db.y-=size.y;
					if(db.y<-size.y/2.0) db.y+=size.y;
					db.z=p[second].z-p[third].z;
					if(db.z>size.z/2.0) db.z-=size.z;
					if(db.z<-size.z/2.0) db.z+=size.z;
					
					this->harmonicF(da, a[first], a[second], &(m[i].getConstants()[CHAINBOND]));
					this->bendF(da, db, a[first], a[second], a[third], &(m[i].getConstants()[CHAINBEND]));
				}
				//last three are a special case, requires bond between last two particles in chain
				first=k+length-3;
				second=k+length-2;
				third=k+length-1;
				
				da.x=p[first].x-p[second].x;
				if(da.x>size.x/2.0) da.x-=size.x;
				if(da.x<-size.x/2.0) da.x+=size.x;
				da.y=p[first].y-p[second].y;
				if(da.y>size.y/2.0) da.y-=size.y;
				if(da.y<-size.y/2.0) da.y+=size.y;
				da.z=p[first].z-p[second].z;
				if(da.z>size.z/2.0) da.z-=size.z;
				if(da.z<-size.z/2.0) da.z+=size.z;
				
				db.x=p[second].x-p[third].x;
				if(db.x>size.x/2.0) db.x-=size.x;
				if(db.x<-size.x/2.0) db.x+=size.x;
				db.y=p[second].y-p[third].y;
				if(db.y>size.y/2.0) db.y-=size.y;
				if(db.y<-size.y/2.0) db.y+=size.y;
				db.z=p[second].z-p[third].z;
				if(db.z>size.z/2.0) db.z-=size.z;
				if(db.z<-size.z/2.0) db.z+=size.z;
				
				this->harmonicF(da, a[first], a[second], &(m[i].getConstants()[CHAINBOND]));
				this->harmonicF(db, a[second], a[third], &(m[i].getConstants()[CHAINBOND]));
				this->bendF(da, db, a[first], a[second], a[third], &(m[i].getConstants()[CHAINBEND]));
			}
		}
	}
}

template <typename T>
void Blob<T>::doDihedralForce(int i)
{
	//not implemented
}

template <typename T>
void Blob<T>::doTorsionForce(int i)
{
	//not implemented
}

template <typename T>
void Blob<T>::doBondForce(int i)
{
	//this is just going to check it without an error condition, it ignores it
	if(m[i].readType()==BOND)
	{
		fourVector<int> *bond=m[i].getBonds();
		threeVector<T> d;
		for(int j=0;j<m[i].readNBond();j++)
		{
			int first=bond[j].s[0];
			int second=bond[j].s[1];
			
			d.x=p[first].x-p[second].x;
			if(d.x>size.x/2.0) d.x-=size.x;
			if(d.x<-size.x/2.0) d.x+=size.x;
			d.y=p[first].y-p[second].y;
			if(d.y>size.y/2.0) d.y-=size.y;
			if(d.y<-size.y/2.0) d.y+=size.y;
			d.z=p[first].z-p[second].z;
			if(d.z>size.z/2.0) d.z-=size.z;
			if(d.z<-size.z/2.0) d.z+=size.z;
			
			harmonicF(d, a[first], a[second], m[i].getConstants());
		}
	}
}

template <typename T>
void Blob<T>::doBendForce(int i)
{
	//this is just going to check it without an error condition, it ignores it
	if(m[i].readType()==BEND)
	{
		fourVector<int> *bond=m[i].getBonds();
		threeVector<T> da,db;
		for(int j=0;j<m[i].readNBond();j++)
		{
			int first=bond[j].s[0];
			int second=bond[j].s[1];
			int third=bond[j].s[2];
			
			da.x=p[first].x-p[second].x;
			if(da.x>size.x/2.0) da.x-=size.x;
			if(da.x<-size.x/2.0) da.x+=size.x;
			da.y=p[first].y-p[second].y;
			if(da.y>size.y/2.0) da.y-=size.y;
			if(da.y<-size.y/2.0) da.y+=size.y;
			da.z=p[first].z-p[second].z;
			if(da.z>size.z/2.0) da.z-=size.z;
			if(da.z<-size.z/2.0) da.z+=size.z;
			
			db.x=p[second].x-p[third].x;
			if(db.x>size.x/2.0) db.x-=size.x;
			if(db.x<-size.x/2.0) db.x+=size.x;
			db.y=p[second].y-p[third].y;
			if(db.y>size.y/2.0) db.y-=size.y;
			if(db.y<-size.y/2.0) db.y+=size.y;
			db.z=p[second].z-p[third].z;
			if(db.z>size.z/2.0) db.z-=size.z;
			if(db.z<-size.z/2.0) db.z+=size.z;
			
			bendF(da, db, a[first], a[second], a[third], m[i].getConstants());
		}
	}
}

template <typename T>
T Blob<T>::doChainPotential(int i)
{
	T potential=0;
	//this is just going to check it without an error condition, it ignores it
	if(m[i].readType()==CHAIN)
	{
		fourVector<int> *bond=m[i].getBonds();
		int first, second, third;
		threeVector<T> da,db;
		for(int j=0;j<m[i].readNBond();j++)
		{
			int start=bond[j].s[START];
			int nChains=bond[j].s[NCHAINS];
			int length=bond[j].s[LENGTH];
			for(int k=start; k<start+length*nChains; k+=length)
			{
				for(int l=k;l<k+length-3;l++)
				{
					int first=l;
					int second=l+1;
					int third=l+2;
					
					da.x=p[first].x-p[second].x;
					if(da.x>size.x/2.0) da.x-=size.x;
					if(da.x<-size.x/2.0) da.x+=size.x;
					da.y=p[first].y-p[second].y;
					if(da.y>size.y/2.0) da.y-=size.y;
					if(da.y<-size.y/2.0) da.y+=size.y;
					da.z=p[first].z-p[second].z;
					if(da.z>size.z/2.0) da.z-=size.z;
					if(da.z<-size.z/2.0) da.z+=size.z;
					
					db.x=p[second].x-p[third].x;
					if(db.x>size.x/2.0) db.x-=size.x;
					if(db.x<-size.x/2.0) db.x+=size.x;
					db.y=p[second].y-p[third].y;
					if(db.y>size.y/2.0) db.y-=size.y;
					if(db.y<-size.y/2.0) db.y+=size.y;
					db.z=p[second].z-p[third].z;
					if(db.z>size.z/2.0) db.z-=size.z;
					if(db.z<-size.z/2.0) db.z+=size.z;
					
					potential+=harmonicP(da, &(m[i].getConstants()[CHAINBOND]));
					potential+=bendP(da, db, &(m[i].getConstants()[CHAINBEND]));
				}
				
				first=k+length-3;
				second=k+length-2;
				third=k+length-1;
				
				da.x=p[first].x-p[second].x;
				if(da.x>size.x/2.0) da.x-=size.x;
				if(da.x<-size.x/2.0) da.x+=size.x;
				da.y=p[first].y-p[second].y;
				if(da.y>size.y/2.0) da.y-=size.y;
				if(da.y<-size.y/2.0) da.y+=size.y;
				da.z=p[first].z-p[second].z;
				if(da.z>size.z/2.0) da.z-=size.z;
				if(da.z<-size.z/2.0) da.z+=size.z;
				
				db.x=p[second].x-p[third].x;
				if(db.x>size.x/2.0) db.x-=size.x;
				if(db.x<-size.x/2.0) db.x+=size.x;
				db.y=p[second].y-p[third].y;
				if(db.y>size.y/2.0) db.y-=size.y;
				if(db.y<-size.y/2.0) db.y+=size.y;
				db.z=p[second].z-p[third].z;
				if(db.z>size.z/2.0) db.z-=size.z;
				if(db.z<-size.z/2.0) db.z+=size.z;
				
				potential+=harmonicP(da, &(m[i].getConstants()[CHAINBOND]));
				potential+=harmonicP(db, &(m[i].getConstants()[CHAINBOND]));
				potential+=bendP(da, db, &(m[i].getConstants()[CHAINBEND]));
			}
		}
	}
	return potential;
}

template <typename T>
T Blob<T>::doDihedralPotential(int i)
{
	return 0;//not yet implemented
}

template <typename T>
T Blob<T>::doTorsionPotential(int i)
{
	return 0;//not yet implemented
}

template <typename T>
T Blob<T>::doBondPotential(int i)
{
	T potential=0;
	//this is just going to check it without an error condition, it ignores it
	if(m[i].readType()==BOND)
	{
		fourVector<int> *bond=m[i].getBonds();
		threeVector<T> d;
		for(int j=0;j<m[i].readNBond();j++)
		{
			int first=bond[j].s[0];
			int second=bond[j].s[1];
			
			d.x=p[first].x-p[second].x;
			if(d.x>size.x/2.0) d.x-=size.x;
			if(d.x<-size.x/2.0) d.x+=size.x;
			d.y=p[first].y-p[second].y;
			if(d.y>size.y/2.0) d.y-=size.y;
			if(d.y<-size.y/2.0) d.y+=size.y;
			d.z=p[first].z-p[second].z;
			if(d.z>size.z/2.0) d.z-=size.z;
			if(d.z<-size.z/2.0) d.z+=size.z;
			
			potential+=harmonicP(d, m[i].getConstants());
		}
	}
	return potential;
}

template <typename T>
T Blob<T>::doBendPotential(int i)
{
	T potential=0;
	//this is just going to check it without an error condition, it ignores it
	if(m[i].readType()==BEND)
	{
		fourVector<int> *bond=m[i].getBonds();
		threeVector<T> da,db;
		for(int j=0;j<m[i].readNBond();j++)
		{
			int first=bond[j].s[0];
			int second=bond[j].s[1];
			int third=bond[j].s[2];
			
			da.x=p[first].x-p[second].x;
			if(da.x>size.x/2.0) da.x-=size.x;
			if(da.x<-size.x/2.0) da.x+=size.x;
			da.y=p[first].y-p[second].y;
			if(da.y>size.y/2.0) da.y-=size.y;
			if(da.y<-size.y/2.0) da.y+=size.y;
			da.z=p[first].z-p[second].z;
			if(da.z>size.z/2.0) da.z-=size.z;
			if(da.z<-size.z/2.0) da.z+=size.z;
			
			db.x=p[second].x-p[third].x;
			if(db.x>size.x/2.0) db.x-=size.x;
			if(db.x<-size.x/2.0) db.x+=size.x;
			db.y=p[second].y-p[third].y;
			if(db.y>size.y/2.0) db.y-=size.y;
			if(db.y<-size.y/2.0) db.y+=size.y;
			db.z=p[second].z-p[third].z;
			if(db.z>size.z/2.0) db.z-=size.z;
			if(db.z<-size.z/2.0) db.z+=size.z;
			
			potential+=bendP(da, db, m[i].getConstants());
		}
	}
	return potential;
}


//These are for the change in potential with a given change in size.
template <typename T>
T Blob<T>::doChainDPotential(int i, threeVector<T> aSize)
{
	T dPotential=0;
	//this is just going to check it without an error condition, it ignores it
	if(m[i].readType()==CHAIN)
	{
		fourVector<int> *bond=m[i].getBonds();
		int first, second, third;
		threeVector<T> da,db;
		for(int j=0;j<m[i].readNBond();j++)
		{
			int start=bond[j].s[START];
			int nChains=bond[j].s[NCHAINS];
			int length=bond[j].s[LENGTH];
			for(int k=start; k<start+length*nChains; k+=length)
			{
				for(int l=k;l<k+length-3;l++)
				{
					int first=l;
					int second=l+1;
					int third=l+2;
					
					da.x=p[first].x-p[second].x;
					if(da.x>size.x/2.0) da.x-=size.x;
					if(da.x<-size.x/2.0) da.x+=size.x;
					da.y=p[first].y-p[second].y;
					if(da.y>size.y/2.0) da.y-=size.y;
					if(da.y<-size.y/2.0) da.y+=size.y;
					da.z=p[first].z-p[second].z;
					if(da.z>size.z/2.0) da.z-=size.z;
					if(da.z<-size.z/2.0) da.z+=size.z;
					
					db.x=p[second].x-p[third].x;
					if(db.x>size.x/2.0) db.x-=size.x;
					if(db.x<-size.x/2.0) db.x+=size.x;
					db.y=p[second].y-p[third].y;
					if(db.y>size.y/2.0) db.y-=size.y;
					if(db.y<-size.y/2.0) db.y+=size.y;
					db.z=p[second].z-p[third].z;
					if(db.z>size.z/2.0) db.z-=size.z;
					if(db.z<-size.z/2.0) db.z+=size.z;
					
					double oldPotential=harmonicP(da, &(m[i].getConstants()[CHAINBOND]));
					oldPotential+=bendP(da, db, &(m[i].getConstants()[CHAINBEND]));
						
					da.x*=aSize.x;
					da.y*=aSize.y;
					da.z*=aSize.z;
					
					db.x*=aSize.x;
					db.y*=aSize.y;
					db.z*=aSize.z;
					
					double newPotential=harmonicP(da, &(m[i].getConstants()[CHAINBOND]));
					newPotential+=bendP(da, db, &(m[i].getConstants()[CHAINBEND]));
					dPotential+=(newPotential-oldPotential);
				}
				
				first=k+length-3;
				second=k+length-2;
				third=k+length-1;
				
				da.x=p[first].x-p[second].x;
				if(da.x>size.x/2.0) da.x-=size.x;
				if(da.x<-size.x/2.0) da.x+=size.x;
				da.y=p[first].y-p[second].y;
				if(da.y>size.y/2.0) da.y-=size.y;
				if(da.y<-size.y/2.0) da.y+=size.y;
				da.z=p[first].z-p[second].z;
				if(da.z>size.z/2.0) da.z-=size.z;
				if(da.z<-size.z/2.0) da.z+=size.z;
				
				db.x=p[second].x-p[third].x;
				if(db.x>size.x/2.0) db.x-=size.x;
				if(db.x<-size.x/2.0) db.x+=size.x;
				db.y=p[second].y-p[third].y;
				if(db.y>size.y/2.0) db.y-=size.y;
				if(db.y<-size.y/2.0) db.y+=size.y;
				db.z=p[second].z-p[third].z;
				if(db.z>size.z/2.0) db.z-=size.z;
				if(db.z<-size.z/2.0) db.z+=size.z;
				
				double oldPotential=harmonicP(da, &(m[i].getConstants()[CHAINBOND]));
				oldPotential+=harmonicP(db, &(m[i].getConstants()[CHAINBOND]));
				oldPotential+=bendP(da, db, &(m[i].getConstants()[CHAINBEND]));
					
				da.x*=aSize.x;
				da.y*=aSize.y;
				da.z*=aSize.z;
				
				db.x*=aSize.x;
				db.y*=aSize.y;
				db.z*=aSize.z;
				
				double newPotential=harmonicP(da, &(m[i].getConstants()[CHAINBOND]));
				newPotential+=harmonicP(db, &(m[i].getConstants()[CHAINBOND]));				
				newPotential+=bendP(da, db, &(m[i].getConstants()[CHAINBEND]));
				dPotential+=(newPotential-oldPotential);
			}
		}
	}
	return dPotential;
}

template <typename T>
T Blob<T>::doDihedralDPotential(int i, threeVector<T> aSize)
{
	return 0;//not yet implemented
}

template <typename T>
T Blob<T>::doTorsionDPotential(int i, threeVector<T> aSize)
{
	return 0;//not yet implemented
}

template <typename T>
T Blob<T>::doBondDPotential(int i, threeVector<T> aSize)
{
	T dPotential=0;
	//this is just going to check it without an error condition, it ignores it
	if(m[i].readType()==BOND)
	{
		fourVector<int> *bond=m[i].getBonds();
		threeVector<T> d;
		for(int j=0;j<m[i].readNBond();j++)
		{
			int first=bond[j].s[0];
			int second=bond[j].s[1];
			
			d.x=p[first].x-p[second].x;
			if(d.x>size.x/2.0) d.x-=size.x;
			if(d.x<-size.x/2.0) d.x+=size.x;
			d.y=p[first].y-p[second].y;
			if(d.y>size.y/2.0) d.y-=size.y;
			if(d.y<-size.y/2.0) d.y+=size.y;
			d.z=p[first].z-p[second].z;
			if(d.z>size.z/2.0) d.z-=size.z;
			if(d.z<-size.z/2.0) d.z+=size.z;
			
			double oldPotential=harmonicP(d, m[i].getConstants());
			
			d.x*=aSize.x;
			d.y*=aSize.y;
			d.z*=aSize.z;
			
			double newPotential=harmonicP(d,m[i].getConstants());
			
			dPotential+=(newPotential-oldPotential);
		}
	}
	return dPotential;
}

template <typename T>
T Blob<T>::doBendDPotential(int i, threeVector<T> aSize)
{
	T dPotential=0;
	//this is just going to check it without an error condition, it ignores it
	if(m[i].readType()==BEND)
	{
		fourVector<int> *bond=m[i].getBonds();
		threeVector<T> da,db;
		for(int j=0;j<m[i].readNBond();j++)
		{
			int first=bond[j].s[0];
			int second=bond[j].s[1];
			int third=bond[j].s[2];
			
			da.x=p[first].x-p[second].x;
			if(da.x>size.x/2.0) da.x-=size.x;
			if(da.x<-size.x/2.0) da.x+=size.x;
			da.y=p[first].y-p[second].y;
			if(da.y>size.y/2.0) da.y-=size.y;
			if(da.y<-size.y/2.0) da.y+=size.y;
			da.z=p[first].z-p[second].z;
			if(da.z>size.z/2.0) da.z-=size.z;
			if(da.z<-size.z/2.0) da.z+=size.z;
			
			db.x=p[second].x-p[third].x;
			if(db.x>size.x/2.0) db.x-=size.x;
			if(db.x<-size.x/2.0) db.x+=size.x;
			db.y=p[second].y-p[third].y;
			if(db.y>size.y/2.0) db.y-=size.y;
			if(db.y<-size.y/2.0) db.y+=size.y;
			db.z=p[second].z-p[third].z;
			if(db.z>size.z/2.0) db.z-=size.z;
			if(db.z<-size.z/2.0) db.z+=size.z;
			
			double oldPotential=bendP(da, db, m[i].getConstants());
			
			da.x*=aSize.x;
			da.y*=aSize.y;
			da.z*=aSize.z;
			
			db.x*=aSize.x;
			db.y*=aSize.y;
			db.z*=aSize.z;
			
			double newPotential=bendP(da, db, m[i].getConstants());
			dPotential+=(newPotential-oldPotential);
		}
	}
	return dPotential;
}

//These are pretty straight forward:
// threeVector d, da, and db are difference vectors (a vector normal to the conservative force between point 1 and 2);
// a1 is the acceleration/force along the difference vector(s);
// a2 is the acceleration/force against the difference vector(s);
// when a3 is present (three-body), a1+a3 are opposite of a2, and any other combination;
// there might eventually be a a2 less one, but I don't care right now.
//They allow you to do your own calculation without the need for a specific particle.
//pair potentials and forces as described by revalee, et. al., forces derived by spangler
//floating point operations (flops) only count required operations
template <typename T>
inline void Blob<T>::nonBondedF(threeVector<T> d, threeVector<T> &a1, threeVector<T> &a2, T *constants)
{
	//3 flops from the difference vector d
	T dr=d.x*d.x+d.y*d.y+d.z*d.z;//5 flops
	//is it in range?
	if(dr<cutoffSquared)
	{
		dr=sqrt(dr);
		
		int offset=int(dr/(constants[0]))*3;
		T magnitude=constants[offset]-dr;//1 flop
		magnitude=((constants[offset+1]-constants[offset+2]*magnitude)*magnitude)/dr;//4 flops
		
		d.x*=magnitude;//1 flop
		d.y*=magnitude;//1 flop
		d.z*=magnitude;//1 flop
		
		a1.x+=d.x;//1 flop
		a1.y+=d.y;//1 flop
		a1.z+=d.z;//1 flop
		#pragma omp atomic
		a2.x-=d.x;//1 flop
		#pragma omp atomic
		a2.y-=d.y;//1 flop
		#pragma omp atomic
		a2.z-=d.z;//1 flop
	}
	//22 "useful" flops total, compare to a lennard jones 23
}

template <typename T>
T Blob<T>::nonBondedP(threeVector<T> d, T *constants)
{
	T potential=0;
	T dr=d.x*d.x+d.y*d.y+d.z*d.z;
	//is it in range?
	if(dr<cutoffSquared)
	{
		dr=sqrt(dr);
		
		if(dr<=constants[0])
		{
			potential=constants[0]-dr;
			potential=constants[1]*potential*potential+constants[2];
		}
		else
		{
			potential=constants[3]-dr;
			potential=potential*potential*(constants[4]-potential*constants[5]);
		}
	}
	return potential;
}

template <typename T>
void Blob<T>::harmonicF(threeVector<T> d, threeVector<T> &a1, threeVector<T> &a2, T *constants)
{
	T dr=sqrt(d.x*d.x+d.y*d.y+d.z*d.z);
	
	T magnitude=dr-constants[0];//1 flop
	magnitude=-magnitude*constants[1]/dr;//2 flops
	
	d.x*=magnitude;
	d.y*=magnitude;
	d.z*=magnitude;
	
	a1.x+=d.x;
	a1.y+=d.y;
	a1.z+=d.z;
	
	a2.x-=d.x;
	a2.y-=d.y;
	a2.z-=d.z;
}

template <typename T>
T Blob<T>::harmonicP(threeVector<T> d, T *constants)
{
	T dr=sqrt(d.x*d.x+d.y*d.y+d.z*d.z);
	
	T potential=dr-constants[0];
	potential=0.5*constants[1]*potential*potential;
	
	return potential;
}

template <typename T>
void Blob<T>::bendF(threeVector<T> da, threeVector<T> db, threeVector<T> &a1, 
		   threeVector<T> &a2,threeVector<T> &a3,T *constants)
{
	T dra=sqrt(da.x*da.x+da.y*da.y+da.z*da.z);
	T drb=sqrt(db.x*db.x+db.y*db.y+db.z*db.z);
	
	da.x/=dra;
	da.y/=dra;
	da.z/=dra;
	
	db.x/=drb;
	db.y/=drb;
	db.z/=drb;
	
	T costheta=(da.x*db.x)+(da.y*db.y)+(da.z*db.z);
	
	T magnitude=constants[0]-costheta;
	magnitude*=constants[1];
	
	T bufa=magnitude*(db.x-(da.x*costheta))/dra;
	T bufb=magnitude*(da.x-(db.x*costheta))/drb;
	
	a1.x+=bufa;
	a2.x+=(bufb-bufa);
	a3.x-=bufb;
	
	bufa=magnitude*(db.y-(da.y*costheta))/dra;
	bufb=magnitude*(da.y-(db.y*costheta))/drb;
	
	a1.y+=bufa;
	a2.y+=(bufb-bufa);
	a3.y-=bufb;
	
	bufa=magnitude*(db.z-(da.z*costheta))/dra;
	bufb=magnitude*(da.z-(db.z*costheta))/drb;
	
	a1.z+=bufa;
	a2.z+=(bufb-bufa);
	a3.z-=bufb;
}

template <typename T>
T Blob<T>::bendP(threeVector<T> da, threeVector<T> db, T *constants)
{
	T dra=sqrt(da.x*da.x+da.y*da.y+da.z*da.z);
	T drb=sqrt(db.x*db.x+db.y*db.y+db.z*db.z);
	
	da.x/=dra;
	da.y/=dra;
	da.z/=dra;
	
	db.x/=drb;
	db.y/=drb;
	db.z/=drb;
	
	T costheta=(da.x*db.x)+(da.y*db.y)+(da.z*db.z);
	
	T potential=constants[0]-costheta;
	potential=constants[1]*potential*potential*0.5;
	
	return potential;
}

template <typename T>
void Blob<T>::addMolecule(molecule<T, fourVector<int> > &value)
{
	if(nAllocated[MOLECULE]==0)
	{
		m=new molecule<T, fourVector<int> >[1];
		nAllocated[MOLECULE]=1;
	}
	else if(nAllocated[MOLECULE]<nMolecules+1)
	{
		//++nAllocated[...] isn't an error, its a preincrement
		molecule<T, fourVector<int> > *buf=new molecule<T, fourVector<int> >[++nAllocated[MOLECULE]];
		for(int i=0;i<nMolecules;i++)
			buf[i]=m[i];
		delete[] m;
		m=buf;
	}
	m[nMolecules++]=value;
	nAllocated[NMOLECULES]=1;
}

template <typename T>
void Blob<T>::setMolecule(molecule<T,fourVector<int> > &value, int index)
{
	#ifdef ERRORS_ENABLED
	if(nAllocated[MOLECULE]==0)
	{
		std::cout << "Error (Blob): molecule doesn't appear to be allocated in setMolecule!\n";
		throw 0;
	}
	if(index>nMolecules-1 || index<0)
	{
		std::cout << "Error (Blob): molecule " << index << " is out of bounds in setMolecule!\n";
		throw 0;
	}
	#endif
	m[index]=value;
}

template <typename T>
void Blob<T>::delMolecule(int index)
{
	#ifdef ERRORS_ENABLED
	if(nAllocated[MOLECULE]==0)
	{
		std::cout << "Error (Blob): molecule doesn't appear to be allocated in delMolecule!\n";
		throw 0;
	}
	if(index>nMolecules-1 || index<0)
	{
		std::cout << "Error (Blob): molecule " << index << " is out of bounds in delMolecule!\n";
		throw 0;
	}
	#endif
	//due to the behavior of molecule, delete is as simple as shifting all of the molecules back one
	for(int i=index;i<nMolecules-1;i++)
		m[i]=m[i+1];
	nMolecules--;
}

template <typename T>
molecule<T,fourVector<int> > * Blob<T>::getMolecule()
{
	#ifdef WARNINGS_ENABLED
	if(nAllocated[MOLECULE]==0)
		std::cout << "Warning (Blob): molecule doesn't appear to be present in getMolecule!\n";
	#endif
		return m;
}

template <typename T>
void Blob<T>::addParticle(position<T> pos, threeVector<T> vel, threeVector<T> acc)
{
	#ifdef WARNINGS_ENABLED
	if(nAllocated[POSITIONS]!=nAllocated[VELOCITIES])
		std::cout << "Warning (Blob): number of velocity and position elements allocated are not equal in addParticle!\n";
	#endif
	if(nAllocated[POSITIONS]==0)
	{
		p=new position<T>[1];
		a=new threeVector<T>[1];
		nAllocated[POSITIONS]=1;
		v=new threeVector<T>[1];
		nAllocated[VELOCITIES]=1;
	} 
	else if(nAllocated[POSITIONS]<nParticles+1)
	{
		position<T> *bufP=new position<T>[++nAllocated[POSITIONS]];
		threeVector<T> *bufA=new threeVector<T>[nAllocated[POSITIONS]];
		threeVector<T> *bufV=new threeVector<T>[++nAllocated[VELOCITIES]];
		for(int i=0;i<nParticles;i++)
		{
			bufP[i]=p[i];
			bufV[i]=v[i];
			bufA[i]=a[i];
		}
		delete p, a, v;
		p=bufP;
		v=bufV;
		a=bufA;
	}
	p[nParticles]=pos;
	v[nParticles]=vel;
	a[nParticles++]=acc;
	nAllocated[NPARTICLES]=1;
}

template <typename T>
void Blob<T>::setParticle(int index, position<T> pos, threeVector<T> vel, threeVector<T> acc)
{
	#ifdef ERRORS_ENABLED
	if(nAllocated[POSITIONS]==0)
	{
		std::cout << "Error (Blob): position elements do not appear to be allocated in setParticle!\n";
		if(nAllocated[VELOCITIES]==0)
			std::cout << "Error (Blob): velocity elements do not appear to be allocated in setParticle!\n";
		throw 0;
	}
	if(nAllocated[VELOCITIES]==0)
	{
		std::cout << "Error (Blob): velocity elements do not appear to be allocated in setParticle!\n";
		if(nAllocated[POSITIONS]==0)
			std::cout << "Error (Blob): position elements do not appear to be allocated in setParticle!\n";
		throw 0;
	}
	if(index>=nParticles || index<0)
	{
		std::cout << "Error (Blob): particle " << index << " is out of bounds in setParticle!\n";
		throw 0;
	}
	#endif
	p[index]=pos;
	v[index]=vel;
	a[index]=acc;
}

template <typename T>
void Blob<T>::delParticle(int index)
{
	#ifdef ERRORS_ENABLED
	if(nAllocated[POSITIONS]==0)
	{
		std::cout << "Error (Blob): position elements do not appear to be allocated in delParticle!\n";
		if(nAllocated[VELOCITIES]==0)
			std::cout << "Error (Blob): velocity elements do not appear to be allocated in delParticle!\n";
		throw 0;
	}
	if(nAllocated[VELOCITIES]==0)
	{
		std::cout << "Error (Blob): velocity elements do not appear to be allocated in delParticle!\n";
		if(nAllocated[POSITIONS]==0)
			std::cout << "Error (Blob): position elements do not appear to be allocated in delParticle!\n";
		throw 0;
	}
	if(index>=nParticles || index<0)
	{
		std::cout << "Error (Blob): particle " << index << " is out of bounds in delParticle!\n";
		throw 0;
	}
	#endif
	for(int i=index;i<nParticles-1;i++)
	{
		p[i]=p[i+1];
		v[i]=v[i+1];
		a[i]=a[i+1];
	}
	nParticles--;
}

template <typename T>
void Blob<T>::allocParticle(int n)
{
	if(nAllocated[POSITIONS]<n || nAllocated[VELOCITIES]<n)
	{
		position<T> *pBuf=new position<T>[n];
		threeVector<T> *vBuf=new threeVector<T>[n];
		threeVector<T> *aBuf=new threeVector<T>[n];
		for(int i=0;i<nParticles;i++)
		{
			pBuf[i]=p[i];
			vBuf[i]=v[i];
			aBuf[i]=a[i];
		}
		if(p!=NULL)
			delete p;
		if(v!=NULL)
			delete v;
		if(a!=NULL)
			delete a;
		p=pBuf;
		v=vBuf;
		a=aBuf;
		nAllocated[POSITIONS]=n;
		nAllocated[VELOCITIES]=n;
	}
}

template <typename T>
position<T> * Blob<T>::getPositions()
{
	#ifdef WARNINGS_ENABLED
	if(nAllocated[POSITIONS]==0)
		std::cout << "Warning (Blob): positions do not appear to be allocated in getPositions!\n";
	#endif
		return p;
}

template <typename T>
threeVector<T> * Blob<T>::getVelocities()
{
	#ifdef WARNINGS_ENABLED
	if(nAllocated[VELOCITIES]==0)
		std::cout << "Warning (Blob): velocities do not appear to be allocated in getVelocities!\n";
	#endif
		return v;
}

template <typename T>
threeVector<T> * Blob<T>::getAccelerations()
{
	#ifdef WARNINGS_ENABLED
	if(nAllocated[POSITIONS]==0)
		std::cout << "Warning (Blob): accelerations do not appear to be allocated in getAccelerations!\n";
	#endif
		return a;
}

template <typename T>
void Blob<T>::addTwoBodyFconst(T value)
{
	T *buf=new T[nAllocated[TWOBODYFCONST]+1];
	for(int i=0;i<nAllocated[TWOBODYFCONST];i++)
		buf[i]=twoBodyFconst[i];
	if(nAllocated[TWOBODYFCONST]!=0)
		delete twoBodyFconst;
	twoBodyFconst=buf;
	
	twoBodyFconst[nAllocated[TWOBODYFCONST]++]=value;
}

template <typename T>
void Blob<T>::setTwoBodyFconst(int index, T value)
{
	#ifdef ERRORS_ENABLED
	if(nAllocated[TWOBODYFCONST]==0)
	{
		std::cout << "Error (Blob): twoBodyFconst doesn't appear to be allocated in setTwoBodyFconst!\n";
		throw 0;
	}
	if(index>=nAllocated[TWOBODYFCONST] || index<0)
	{
		std::cout << "Error (Blob): index " << index << " is out of bounds in setTwoBodyFconst!\n";
		throw 0;
	}
	#endif
	twoBodyFconst[index]=value;
}

template <typename T>
void Blob<T>::delTwoBodyFconst(int index)
{
	#ifdef ERRORS_ENABLED
	if(nAllocated[TWOBODYFCONST]==0)
	{
		std::cout << "Error (Blob): twoBodyFconst doesn't appear to be allocated in delTwoBodyFconst!\n";
		throw 0;
	}
	if(index>=nAllocated[TWOBODYFCONST] || index<0)
	{
		std::cout << "Error (Blob): index " << index << " is out of bounds in delTwoBodyFconst!\n";
		throw 0;
	}
	#endif
	T *buf=new T[nAllocated[TWOBODYFCONST]-1];
	for(int i=0;i<index;i++)
		buf[i]=twoBodyFconst[i];
	for(int i=index+1;i<nAllocated[TWOBODYFCONST];i++)
		buf[i-1]=twoBodyFconst[i];
	if(twoBodyFconst!=NULL)
		delete twoBodyFconst;
	twoBodyFconst=buf;
	nAllocated[TWOBODYFCONST]--;
}

template <typename T>
T * Blob<T>::getTwoBodyFconst()
{
	#ifdef WARNINGS_ENABLED
	if(nAllocated[TWOBODYFCONST]==0)
		std::cout << "Warning (Blob): twoBodyFconst does not appear to be allocated in getTwoBodyFconst!\n";
	#endif
		return twoBodyFconst;
}

template <typename T>
void Blob<T>::addTwoBodyUconst(T value)
{
	T *buf=new T[nAllocated[TWOBODYUCONST]+1];
	for(int i=0;i<nAllocated[TWOBODYUCONST];i++)
		buf[i]=twoBodyUconst[i];
	if(nAllocated[TWOBODYUCONST]>0)
		delete twoBodyUconst;
	twoBodyUconst=buf;
	
	twoBodyUconst[nAllocated[TWOBODYUCONST]++]=value;
}

template <typename T>
void Blob<T>::setTwoBodyUconst(int index, T value)
{
	#ifdef ERRORS_ENABLED
	if(nAllocated[TWOBODYUCONST]==0)
	{
		std::cout << "Error (Blob): twoBodyUconst doesn't appear to be allocated in setTwoBodyUconst!\n";
		throw 0;
	}
	if(index>=nAllocated[TWOBODYUCONST] || index<0)
	{
		std::cout << "Error (Blob): index " << index << " is out of bounds in setTwoBodyUconst!\n";
		throw 0;
	}
	#endif
	twoBodyUconst[index]=value;
}

template <typename T>
void Blob<T>::delTwoBodyUconst(int index)
{
	#ifdef ERRORS_ENABLED
	if(nAllocated[TWOBODYUCONST]==0)
	{
		std::cout << "Error (Blob): twoBodyUconst doesn't appear to be allocated in delTwoBodyUconst!\n";
		throw 0;
	}
	if(index>=nAllocated[TWOBODYUCONST] || index<0)
	{
		std::cout << "Error (Blob): index " << index << " is out of bounds in delTwoBodyUconst!\n";
		throw 0;
	}
	#endif
	T *buf=new T[nAllocated[TWOBODYUCONST]-1];
	for(int i=0;i<index;i++)
		buf[i]=twoBodyUconst[i];
	for(int i=index+1;i<nAllocated[TWOBODYUCONST];i++)
		buf[i-1]=twoBodyUconst[i];
	if(twoBodyUconst!=NULL)
		delete twoBodyUconst;
	twoBodyUconst=buf;
	nAllocated[TWOBODYUCONST]--;
}

template <typename T>
T * Blob<T>::getTwoBodyUconst()
{
	#ifdef WARNINGS_ENABLED
	if(nAllocated[TWOBODYUCONST]==0)
		std::cout << "Warning (Blob): twoBodyUconst does not appear to be allocated in getTwoBodyUconst!\n";
	#endif
		return twoBodyUconst;
}

template <typename T>
void Blob<T>::addTwoBodyMFconst(T value)
{
	T *buf=new T[nAllocated[TWOBODYMFCONST]+1];
	for(int i=0;i<nAllocated[TWOBODYMFCONST];i++)
		buf[i]=twoBodyMFconst[i];
	if(nAllocated[TWOBODYMFCONST]>0)
		delete twoBodyMFconst;
	twoBodyMFconst=buf;
	twoBodyMFconst[nAllocated[TWOBODYMFCONST]++]=value;
}

template <typename T>
void Blob<T>::setTwoBodyMFconst(int index, T value)
{
	#ifdef ERRORS_ENABLED
	if(nAllocated[TWOBODYMFCONST]==0)
	{
		std::cout << "Error (Blob): twoBodyMFconst doesn't appear to be allocated in setTwoBodyMFconst!\n";
		throw 0;
	}
	if(index>=nAllocated[TWOBODYMFCONST] || index<0)
	{
		std::cout << "Error (Blob): index " << index << " is out of bounds in setTwoBodyMFconst!\n";
		throw 0;
	}
	#endif
	twoBodyMFconst[index]=value;
}

template <typename T>
void Blob<T>::delTwoBodyMFconst(int index)
{
	#ifdef ERRORS_ENABLED
	if(nAllocated[TWOBODYMFCONST]==0)
	{
		std::cout << "Error (Blob): twoBodyMFconst doesn't appear to be allocated in delTwoBodyMFconst!\n";
		throw 0;
	}
	if(index>=nAllocated[TWOBODYMFCONST] || index<0)
	{
		std::cout << "Error (Blob): index " << index << " is out of bounds in delTwoBodyMFconst!\n";
		throw 0;
	}
	#endif
	T *buf=new T[nAllocated[TWOBODYMFCONST]-1];
	for(int i=0;i<index;i++)
		buf[i]=twoBodyMFconst[i];
	for(int i=index+1;i<nAllocated[TWOBODYMFCONST];i++)
		buf[i-1]=twoBodyMFconst[i];
	if(twoBodyMFconst!=NULL)
		delete twoBodyMFconst;
	twoBodyMFconst=buf;
	nAllocated[TWOBODYMFCONST]--;
}

template <typename T>
T * Blob<T>::getTwoBodyMFconst()
{
	#ifdef WARNINGS_ENABLED
	if(nAllocated[TWOBODYMFCONST]==0)
		std::cout << "Warning (Blob): twoBodyMFconst does not appear to be allocated in getTwoBodyMFconst!\n";
	#endif
		return twoBodyMFconst;
}

template <typename T>
void Blob<T>::addTwoBodyMUconst(T value)
{
	T *buf=new T[nAllocated[TWOBODYMUCONST]+1];
	for(int i=0;i<nAllocated[TWOBODYMUCONST];i++)
		buf[i]=twoBodyMUconst[i];
	if(nAllocated[TWOBODYMUCONST]>0)
		delete twoBodyMUconst;
	twoBodyMUconst=buf;
	twoBodyMUconst[nAllocated[TWOBODYMUCONST]++]=value;
}

template <typename T>
void Blob<T>::setTwoBodyMUconst(int index, T value)
{
	#ifdef ERRORS_ENABLED
	if(nAllocated[TWOBODYMUCONST]==0)
	{
		std::cout << "Error (Blob): twoBodyMUconst doesn't appear to be allocated in setTwoBodyMUconst!\n";
		throw 0;
	}
	if(index>=nAllocated[TWOBODYMUCONST] || index<0)
	{
		std::cout << "Error (Blob): index " << index << " is out of bounds in setTwoBodyMUconst!\n";
		throw 0;
	}
	#endif
	twoBodyMUconst[index]=value;
}

template <typename T>
void Blob<T>::delTwoBodyMUconst(int index)
{
	#ifdef ERRORS_ENABLED
	if(nAllocated[TWOBODYMUCONST]==0)
	{
		std::cout << "Error (Blob): twoBodyMUconst doesn't appear to be allocated in delTwoBodyMUconst!\n";
		throw 0;
	}
	if(index>=nAllocated[TWOBODYMUCONST] || index<0)
	{
		std::cout << "Error (Blob): index " << index << " is out of bounds in delTwoBodyMUconst!\n";
		throw 0;
	}
	#endif
	T *buf=new T[nAllocated[TWOBODYMUCONST]-1];
	for(int i=0;i<index;i++)
		buf[i]=twoBodyMUconst[i];
	for(int i=index+1;i<nAllocated[TWOBODYMUCONST];i++)
		buf[i-1]=twoBodyMUconst[i];
	if(twoBodyMUconst!=NULL)
		delete twoBodyMUconst;
	twoBodyMUconst=buf;
	nAllocated[TWOBODYMUCONST]--;
}

template <typename T>
T * Blob<T>::getTwoBodyMUconst()
{
	#ifdef WARNINGS_ENABLED
	if(nAllocated[TWOBODYMUCONST]==0)
		std::cout << "Warning (Blob): twoBodyMUconst does not appear to be allocated in getTwoBodyMUconst!\n";
	#endif
		return twoBodyMUconst;
}

template <typename T>
void Blob<T>::addThreeBodyMFconst(T value)
{
	T *buf=new T[nAllocated[THREEBODYMFCONST]+1];
	for(int i=0;i<nAllocated[THREEBODYMFCONST];i++)
		buf[i]=threeBodyMFconst[i];
	if(nAllocated[THREEBODYMFCONST]>0)
		delete threeBodyMFconst;
	threeBodyMFconst=buf;
	threeBodyMFconst[nAllocated[THREEBODYMFCONST]++]=value;
}

template <typename T>
void Blob<T>::setThreeBodyMFconst(int index, T value)
{
	#ifdef ERRORS_ENABLED
	if(nAllocated[THREEBODYMFCONST]==0)
	{
		std::cout << "Error (Blob): threeBodyMFconst doesn't appear to be allocated in setThreeBodyMFconst!\n";
		throw 0;
	}
	if(index>=nAllocated[THREEBODYMFCONST] || index<0)
	{
		std::cout << "Error (Blob): index " << index << " is out of bounds in setThreeBodyMFconst!\n";
		throw 0;
	}
	#endif
	threeBodyMFconst[index]=value;
}

template <typename T>
void Blob<T>::delThreeBodyMFconst(int index)
{
	#ifdef ERRORS_ENABLED
	if(nAllocated[THREEBODYMFCONST]==0)
	{
		std::cout << "Error (Blob): threeBodyMFconst doesn't appear to be allocated in delTwoBodyMFconst!\n";
		throw 0;
	}
	if(index>=nAllocated[THREEBODYMFCONST] || index<0)
	{
		std::cout << "Error (Blob): index " << index << " is out of bounds in delTwoBodyMFconst!\n";
		throw 0;
	}
	#endif
	T *buf=new T[nAllocated[THREEBODYMFCONST]-1];
	for(int i=0;i<index;i++)
		buf[i]=threeBodyMFconst[i];
	for(int i=index+1;i<nAllocated[THREEBODYMFCONST];i++)
		buf[i-1]=threeBodyMFconst[i];
	if(threeBodyMFconst!=NULL)
		delete threeBodyMFconst;
	threeBodyMFconst=buf;
	nAllocated[THREEBODYMFCONST]--;
}

template <typename T>
T * Blob<T>::getThreeBodyMFconst()
{
	#ifdef WARNINGS_ENABLED
	if(nAllocated[THREEBODYMFCONST]==0)
		std::cout << "Warning (Blob): threeBodyMFconst does not appear to be allocated in getThreeBodyMFconst!\n";
	#endif
		return threeBodyMFconst;
}

template <typename T>
void Blob<T>::addThreeBodyMUconst(T value)
{
	T *buf=new T[nAllocated[THREEBODYMUCONST]+1];
	for(int i=0;i<nAllocated[THREEBODYMUCONST];i++)
		buf[i]=threeBodyMUconst[i];
	if(nAllocated[THREEBODYMUCONST]>0)
		delete threeBodyMUconst;
	threeBodyMUconst=buf;
	threeBodyMUconst[nAllocated[THREEBODYMUCONST]++]=value;
}

template <typename T>
void Blob<T>::setThreeBodyMUconst(int index, T value)
{
	#ifdef ERRORS_ENABLED
	if(nAllocated[THREEBODYMUCONST]==0)
	{
		std::cout << "Error (Blob): threeBodyMUconst doesn't appear to be allocated in setThreeBodyMUconst!\n";
		throw 0;
	}
	if(index>=nAllocated[THREEBODYMUCONST] || index<0)
	{
		std::cout << "Error (Blob): index " << index << " is out of bounds in setThreeBodyMUconst!\n";
		throw 0;
	}
	#endif
	threeBodyMUconst[index]=value;
}

template <typename T>
void Blob<T>::delThreeBodyMUconst(int index)
{
	#ifdef ERRORS_ENABLED
	if(nAllocated[THREEBODYMUCONST]==0)
	{
		std::cout << "Error (Blob): threeBodyMUconst doesn't appear to be allocated in delTwoBodyMUconst!\n";
		throw 0;
	}
	if(index>=nAllocated[THREEBODYMUCONST] || index<0)
	{
		std::cout << "Error (Blob): index " << index << " is out of bounds in delTwoBodyMUconst!\n";
		throw 0;
	}
	#endif
	T *buf=new T[nAllocated[THREEBODYMUCONST]-1];
	for(int i=0;i<index;i++)
		buf[i]=threeBodyMUconst[i];
	for(int i=index+1;i<nAllocated[THREEBODYMUCONST];i++)
		buf[i-1]=threeBodyMUconst[i];
	if(threeBodyMFconst!=NULL)
		delete threeBodyMUconst;
	threeBodyMUconst=buf;
	nAllocated[THREEBODYMUCONST]--;
}

template <typename T>
T * Blob<T>::getThreeBodyMUconst()
{
	#ifdef WARNINGS_ENABLED
	if(nAllocated[THREEBODYMUCONST]==0)
		std::cout << "Warning (Blob): threeBodyMUconst does not appear to be allocated in getThreeBodyMUconst!\n";
	#endif
		return threeBodyMUconst;
}

template <typename T>
int Blob<T>::readNParticles()
{
	#ifdef WARNINGS_ENABLED
	if(0==nAllocated[NPARTICLES])
		std::cout << "Warning (Blob): nParticles doesn't appear to be present in readNParticles!\n";
	#endif
		return nParticles;
}

template <typename T>
int Blob<T>::readNMolecules()
{
	//#ifdef WARNINGS_ENABLED
	//if(0==nAllocated[NMOLECULES])
	//	std::cout << "Warning (Blob): nMolecules doesn't appear to be present in readNMolecules!\n";
	//#endif
		return nMolecules;
}

template <typename T>
T Blob<T>::readGamma()
{
	#ifdef WARNINGS_ENABLED
	if(0==nAllocated[GAMMA])
		std::cout << "Warning (Blob): gamma doesn't appear to be present in readGamma!\n";
	#endif
		return gamma;
}

template <typename T>
T Blob<T>::readInitialTemp()
{
	#ifdef WARNINGS_ENABLED
	if(0==nAllocated[INITIALTEMP])
		std::cout << "Warning (Blob): initialTemp doesn't appear to be present in readInitialTemp!\n";
	#endif
		return initialTemp;
}

template <typename T>
T Blob<T>::readFinalTemp()
{
	#ifdef WARNINGS_ENABLED
	if(0==nAllocated[FINALTEMP])
		std::cout << "Warning (Blob): finalTemp doesn't appear to be present in readFinalTemp!\n";
	#endif
		return finalTemp;
}

template <typename T>
int Blob<T>::readSeed()
{
	#ifdef WARNINGS_ENABLED
	if(0==nAllocated[SEED])
		std::cout << "Warning (Blob): seed doesn't appear to be present in readSeed!\n";
	#endif
		return seed;
}

template <typename T>
int Blob<T>::readNTypes()
{
	#ifdef WARNINGS_ENABLED
	if(0==nAllocated[NTYPES])
		std::cout << "Warning (Blob): nTypes doesn't appear to be present in readNTypes!\n";
	#endif
		return nTypes;
}

template <typename T>
threeVector<bool> Blob<T>::readPeriodic()
{
	#ifdef WARNINGS_ENABLED
	if(0==nAllocated[PERIODIC])
		std::cout << "Warning (Blob): periodic doesn't appear to be present in readPeriodic!\n";
	#endif
		return periodic;
}

template <typename T>
T Blob<T>::readCutoff()
{
	#ifdef WARNINGS_ENABLED
	if(0==nAllocated[CUTOFF])
		std::cout << "Warning (Blob): cutoff doesn't appear to be present in readCutoff!\n";
	#endif
		return cutoff;
}

template <typename T>
threeVector<T> Blob<T>::readSize()
{
	#ifdef WARNINGS_ENABLED
	if(0==nAllocated[SIZE])
		std::cout << "Warning (Blob): size doesn't appear to be present in readSize!\n";
	#endif
		return size;
}

template <typename T>
T Blob<T>::readInitialTime()
{
	#ifdef WARNINGS_ENABLED
	if(0==nAllocated[INITIALTIME])
		std::cout << "Warning (Blob): initialTime doesn't appear to be present in readInitialTime!\n";
	#endif
		return initialTime;
}

template <typename T>
T Blob<T>::readFinalTime()
{
	#ifdef WARNINGS_ENABLED
	if(0==nAllocated[FINALTIME])
		std::cout << "Warning (Blob): finalTime doesn't appear to be present in readFinalTime!\n";
	#endif
		return finalTime;
}

template <typename T>
T Blob<T>::readDeltaT()
{
	#ifdef WARNINGS_ENABLED
	if(0==nAllocated[DELTAT])
		std::cout << "Warning (Blob): deltaT doesn't appear to be present in readDeltaT!\n";
	#endif
		return deltaT;
}

template <typename T>
T Blob<T>::readStoreInterval()
{
	#ifdef WARNINGS_ENABLED
	if(0==nAllocated[STOREINTERVAL])
		std::cout << "Warning (Blob): storeInterval doesn't appear to be present in readStoreInterval!\n";
	#endif
		return storeInterval;
}

template <typename T>
T Blob<T>::readMeasureInterval()
{
	#ifdef WARNINGS_ENABLED
	if(0==nAllocated[MEASUREINTERVAL])
		std::cout << "Warning (Blob): measureInterval doesn't appear to be present in readMeasureInterval!\n";
	#endif
		return measureInterval;
}

template <typename T>
T Blob<T>::readDeltaLXY()
{
	//#ifdef WARNINGS_ENABLED
	//if(0==nAllocated[DELTALXY])
	//	std::cout << "Warning (Blob): deltaLXY doesn't appear to be present in readDeltaLXY!\n";
	//#endif
		return deltaLXY;
}

template <typename T>
T Blob<T>::readTempStepInterval()
{
	return tempStepInterval;
}

//Keep this comment here:
//COMMANDINPUTHANDLER
template <typename T>
T Blob<T>::readSolventGamma()
{
	return solventGamma;
}

template <typename T>
T Blob<T>::readRemoveSolvent()
{
	return removeSolvent;
}

template <typename T>
void Blob<T>::setGamma(T value)
{
	gamma=value;
	nAllocated[GAMMA]=1;
}

template <typename T>
void Blob<T>::setInitialTemp(T value)
{
	initialTemp=value;
	nAllocated[INITIALTEMP]=1;
}

template <typename T>
void Blob<T>::setFinalTemp(T value)
{
	finalTemp=value;
	nAllocated[FINALTEMP]=1;
}

template <typename T>
void Blob<T>::setSeed(int value)
{
	seed=value;
	nAllocated[SEED]=1;
}

template <typename T>
void Blob<T>::setNTypes(int value)
{
	nTypes=value;
	nAllocated[NTYPES]=1;
}

template <typename T>
void Blob<T>::setPeriodic(threeVector<bool> value)
{
	periodic=value;
	nAllocated[PERIODIC]=1;
}

template <typename T>
void Blob<T>::setCutoff(T value)
{
	cutoff=value;
	nAllocated[CUTOFF]=1;
}

template <typename T>
void Blob<T>::setSize(threeVector<T> value)
{
	size=value;
	nAllocated[SIZE]=1;
}

template <typename T>
void Blob<T>::setInitialTime(T value)
{
	initialTime=value;
	nAllocated[INITIALTIME]=1;
}

template <typename T>
void Blob<T>::setFinalTime(T value)
{
	finalTime=value;
	nAllocated[FINALTIME]=1;
}

template <typename T>
void Blob<T>::setDeltaT(T value)
{
	deltaT=value;
	nAllocated[DELTAT]=1;
}

template <typename T>
void Blob<T>::setStoreInterval(T value)
{
	storeInterval=value;
	nAllocated[STOREINTERVAL]=1;
}

template <typename T>
void Blob<T>::setMeasureInterval(T value)
{
	measureInterval=value;
	nAllocated[MEASUREINTERVAL]=1;
}

template <typename T>
void Blob<T>::setDeltaLXY(T value)
{
	deltaLXY=value;
	nAllocated[DELTALXY]=1;
}

template <typename T>
void Blob<T>::setRemoveSolvent(T value)
{
	removeSolvent=value;
	nAllocated[REMOVESOLVENT]=1;
}

template <typename T>
void Blob<T>::setTempStepInterval(T value)
{
	tempStepInterval=value;
	nAllocated[TEMPSTEPINTERVAL]=1;
}

//Keep this comment here:
//COMMANDINPUTHANDLER
template <typename T>
void Blob<T>::setSolventGamma(T value)
{
	solventGamma=value;
	nAllocated[SOLVENTGAMMA]=1;
}

//These are ultra optimized for use in CellOpt
//Non-bonded pair force function as described by Revalee, et. al., derived by Spangler.
//It is optimized for utilization in cellOpt.
//Note that there is a slow down using this function. CellOpt is capable of being 4 times faster without this function.
template <typename T>
inline void Force(T cutoffSquared, int nT, position<T> &p1, position<T> &p2, threeVector<T> &a1, threeVector<T> &a2, T *fC)
{
	T dx,dy,dz,dr;

	dx=p1.x-p2.x;//1 flop
	dy=p1.y-p2.y;//1 flop
	dz=p1.z-p2.z;//1 flop
	
	dr=dx*dx+dy*dy+dz*dz;//5 flops
	
	//is it in range?
	#ifndef SOLVENT_FLAG
		if(dr<cutoffSquared)
	#else
		int cindex=nTWOBODYFCONST*((p1.type*nT)+p2.type);
		//also, don't do solvent solvent interactions
		if(dr<cutoffSquared && cindex!=nTWOBODYFCONST*((SOLVENT_FLAG*nT)+SOLVENT_FLAG))
	#endif
	{
		dr=sqrt(dr);
		#ifndef SOLVENT_FLAG
			int cindex=nTWOBODYFCONST*((p1.type*nT)+p2.type);
		#endif
		//floor() is worse than int()
		cindex+=(int(dr/(fC[cindex]))*3);
		T magnitude=fC[cindex]-dr;//1 flop
		magnitude=((fC[cindex+1]-fC[cindex+2]*magnitude)*magnitude)/dr;//4 flops
		
		dx*=magnitude;//1 flop
		dy*=magnitude;//1 flop
		dz*=magnitude;//1 flop
		
		a1.x+=dx;//1 flop
		a1.y+=dy;//1 flop
		a1.z+=dz;//1 flop
		a2.x-=dx;//1 flop
		a2.y-=dy;//1 flop
		a2.z-=dz;//1 flop
		
	}
	//22 "useful" flops total, compare to a lennard jones 23
}

///Non-bonded pair potential function as described by Revalee, et. al..
///It is optimized for use in cellOpt.
template <typename T>
inline T Potential(T cutoffSquared, int nT, position<T> &p1, position<T> &p2, T *uC)
{
	T dx,dy,dz,dr,potential=0;
	dx=p1.x-p2.x;
	dy=p1.y-p2.y;
	dz=p1.z-p2.z;
	
	dr=dx*dx+dy*dy+dz*dz;
	//is it in range?
	#ifndef SOLVENT_FLAG
		if(dr<cutoffSquared)
	#else
		int cindex=nTWOBODYUCONST*((p1.type*nT)+p2.type);
		//also, don't do solvent solvent interactions
		if(dr<cutoffSquared && cindex!=nTWOBODYUCONST*((SOLVENT_FLAG*nT)+SOLVENT_FLAG))
	#endif
	{
		dr=sqrt(dr);
		#ifndef SOLVENT_FLAG
			int cindex=nTWOBODYUCONST*((p1.type*nT)+p2.type);
		#endif
		
		if(dr<=uC[cindex])
		{
			potential=uC[cindex]-dr;
			potential=uC[cindex+1]*potential*potential+uC[cindex+2];
		}
		else
		{
			potential=uC[cindex+3]-dr;
			potential=potential*potential*(uC[cindex+4]-potential*uC[cindex+5]);
		}
	}
	return potential;
}

///Counter non-bonded pair force. To subtract bonds if it isn't easy to do otherwise. Optimized for cellOpt.
template <typename T>
inline void CounterForce(T cutoffSquared, int nT, position<T> &p1, position<T> &p2, threeVector<T> &a1, threeVector<T> &a2, T *fC)
{
	T dx,dy,dz,dr;

	dx=p1.x-p2.x;//1 flop
	dy=p1.y-p2.y;//1 flop
	dz=p1.z-p2.z;//1 flop
	
	dr=dx*dx+dy*dy+dz*dz;//5 flops
	//is it in range?
	if(dr<cutoffSquared)
	{
		dr=sqrt(dr);
		
		int cindex=nTWOBODYFCONST*((p1.type*nT)+p2.type);
		cindex+=(int(dr/(fC[cindex]))*3);
		T magnitude=fC[cindex]-dr;//1 flop
		magnitude=((fC[cindex+1]-fC[cindex+2]*magnitude)*magnitude)/dr;//4 flops
		
		dx*=magnitude;//1 flop
		dy*=magnitude;//1 flop
		dz*=magnitude;//1 flop
		
		a1.x-=dx;//1 flop
		a1.y-=dy;//1 flop
		a1.z-=dz;//1 flop
		a2.x+=dx;//1 flop
		a2.y+=dy;//1 flop
		a2.z+=dz;//1 flop
		
	}
	//22 "useful" flops total, compare to a lennard jones 23
}

///Counter non-bonded pair potential. To subtract bonds if it isn't easy to do otherwise. Optimized for cellOpt.
template <typename T>
inline T CounterPotential(T cutoffSquared, int nT, position<T> &p1, position<T> &p2, T *uC)
{
	T dx,dy,dz,dr,potential=0;
	dx=p1.x-p2.x;
	dy=p1.y-p2.y;
	dz=p1.z-p2.z;
	
	dr=dx*dx+dy*dy+dz*dz;
	//is it in range?
	if(dr<cutoffSquared)
	{
		dr=sqrt(dr);
		int cindex=nTWOBODYUCONST*((p1.type*nT)+p2.type);
		
		if(dr<=uC[cindex])
		{
			potential=uC[cindex]-dr;
			potential=uC[cindex+1]*potential*potential+uC[cindex+2];
		}
		else
		{
			potential=uC[cindex+3]-dr;
			potential=potential*potential*(uC[cindex+4]-potential*uC[cindex+5]);
		}
	}
	return 0-potential;
}
