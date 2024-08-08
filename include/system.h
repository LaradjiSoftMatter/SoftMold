//this is for some lipid specific data collection
//#include "lipids/lipids.h"
#include "fileFormats/scriptFormat.h"
#include "../analysis/cellInclude/cell.h"
#include <array>
/**
 * \brief Molecular dynamics system blob.
 * This is an example system blob. Works with Revalee's, et al., model
 * This blob is system specific, ergo the filename name systemMD.h. It also contains some
 * of the script handling required by scriptFormat.h. Flags can be set here as well.
 * This object maintains all the variables and functions related to the simulation.
 * Most of this is user defined, but in some instances it must have certain public members.
 */
template <typename T>
class Blob {
	public:
		//the Script template shouldn't see the constructors or destructors
		///The constructor is generic, it just sets the values to 0, NULL, or a default value.
		Blob();
		///Destructor. Deletes memory that was allocated.
		~Blob();
		
		//the required functions for script/parameters, etc...
		///Resets the read write stream for scriptFormat
		void reset();
		///outputs one element of the stream based on wP and it, returns "end" when done
		std::string output();
		///inputs one element of the stream based on wP and it, returns false if the input ended
		bool input(std::string in);
		
		///This just checks if anything in the simulation has gone out of some error bound
		void errorChecking();
		
		//System specific variables
		///Add a molecule residue to the system.
		molecule<T, fourVector<int> > * addMolecule(molecule<T, fourVector<int> > &value);
		///Set or change a molecule residue in the system.
		molecule<T, fourVector<int> > * setMolecule(molecule<T,fourVector<int> > &value, int index);
		///Delete a molecule residue in the system.
		molecule<T, fourVector<int> > * delMolecule(int index);
		///Return the pointer to the molecule residues in the system.
		molecule<T, fourVector<int> > * getMolecule();
		
		//if you add a particle you include the complete state of the particle as well
		///Add a particle to the system. Includes the whole initial state for equations of motion.
		void addParticle(position<T> pos, threeVector<T> vel, threeVector<T> acc);
		///Set or change a particle in the system.
		void setParticle(int index, position<T> pos, threeVector<T> vel, threeVector<T> acc);
		///Delete a particle in the system.
		void delParticle(int index);
		///Allocate some particles in the system. This is here to reduce excessive allocations.
		void allocParticle(int n);
		///Get the pointer to the positions of particles in the system.
		position<T> * getPositions();
		///Get the pointer to the velocities of the particles in the system.
		threeVector<T> * getVelocities();
		///Get the pointer to the accelerations of the particles in the system.
		threeVector<T> * getAccelerations();
		
		///Add a non-bonded force constant to the system.
		T * addTwoBodyFconst(T value);
		///Set or change a non-bonded force constant in the system.
		void setTwoBodyFconst(int index, T value);
		///Delete a non-bonded force constant in the system.
		T * delTwoBodyFconst(int index);
		///Get a pointer to the non-bonded force constants in the system.
		T * getTwoBodyFconst();
		
		///Add a non-bonded potential constant to the system.
		T * addTwoBodyUconst(T value);
		///Set or change a non-bonded potential constant in the system.
		void setTwoBodyUconst(int index, T value);
		///Delete a non-bonded potential constant in the system.
		T * delTwoBodyUconst(int index);
		///Get a pointer to the non-bonded potential constants in the system.
		T * getTwoBodyUconst();
		
		///Read the number of particles in the system.
		int readNParticles();
		///Read the number of molecular residues in the system.
		int readNMolecules();
		///Read gamma for the Langevin thermostat.
		T readGamma();
		///Read the initial temperature of the system, in reduced units.
		T readInitialTemp();
		///Read the final temperature of the system, in reduced units.
		T readFinalTemp();
		///Read the random number seed of the system.
		int readSeed();
		///Read the number of particle types in the system.
		int readNTypes();
		///Read the periodic boundaries flag in the system.
		threeVector<bool> readPeriodic();
		///Read the force and potential cutoff of the non-bonded interactions.
		T readCutoff();
		///Read the size of the system, in reduced units.
		threeVector<T> readSize();
		///Read the initial time of the system, in tau.
		T readInitialTime();
		///Read the ending time of the system, in tau.
		T readFinalTime();
		///Read the timestep of the system, in tau.
		T readDeltaT();
		///Read the store interval of the system, in tau.
		T readStoreInterval();
		///Read the measure interval of the system, in tau.
		T readMeasureInterval();
		///Read the maximum step size of the constant pressure length change, in reduced units.
		T readDeltaLXY();
		///Read the amount of solvent to remove, as a fraction of total solvent.
		T readRemoveSolvent();
		///Read temperature step interval.
		T readTempStepInterval();
		///Read solvent gamma.
		twoVector<T> readSolventGamma();
		///Read gamma type.
		T readGammaType(int type);
		///Get gammaType pointer
		T * getGammaType();
		///Read constant tension parameter.
		T readTension();
		
		///Set or change gamma.
		void setGamma(T value);
		///Set or change the initial temperature, in reduced units.
		void setInitialTemp(T value);
		///Set or change the final temperature, in reduced units.
		void setFinalTemp(T value);
		///Set or change the random number seed.
		void setSeed(int value);
		///Set or change the number of particle types.
		void setNTypes(int value);
		///Set or change the periodic boundaries flag.
		void setPeriodic(threeVector<bool> value);
		///Set or change the non-bonded cutoff.
		void setCutoff(T value);
		///Set or change the system size, in reduced units.
		void setSize(threeVector<T> value);
		///Set or change the initial time, in tau.
		void setInitialTime(T value);
		///Set or change the final time, in tau.
		void setFinalTime(T value);
		///Set or change the timestep, in tau.
		void setDeltaT(T value);
		///Set or change the store interval, in tau.
		void setStoreInterval(T value);
		///Set or change the measure interval, in tau.
		void setMeasureInterval(T value);
		///Set or change the maximum step size of the constant pressure length change, in reduced units.
		void setDeltaLXY(T value);
		///Set or change the fraction of solvent removal, as a fraction of total solvent.
		void setRemoveSolvent(T value);
		///Set temperature step interval.
		void setTempStepInterval(T value);
		///Set solvent gamma.
		void setSolventGamma(twoVector<T> value);
		///Set gamma by type.
		void setGammaType(int type, T value);
		///Set constant tension parameter
		void setTension(T value);
		
		//Functions for the simulation
		///Non-bonded force function between particles i and j. Utilizing the minimum image for boundaries.
		void doForce(int i, int j, threeVector<T> &minImg);
		///Non-bonded potential function between particles i and j. Utilizing the minimum image for boundaries.
		T doPotential(int i, int j, threeVector<T> &minImg);
		
		//Types of molecules have been moved to dataTypes
		//After years and years of trying to make a molecules class... This is what you get:
		///Do CHAIN force bonded interactions on residue i. CHAIN is a combination of two and three body bonded forces.
		void doChainForce(int i);
		///Do DIHEADRAL force bonded interactions on residue i.
		void doDihedralForce(int i);
		///Do TORSION force bonded interactions on residue i.
		void doTorsionForce(int i);
		///Do BOND force bonded interactions on residue i.
		void doBondForce(int i);
		///Do BEND force bonded interactions on residue i.
		void doBendForce(int i);
		///Do BEAD force bonded interactions on residue i.
		void doBeadForce(int i);
		///Do NANOCORE force bonded interactions on residue i. Same as BEAD, but lacks BEAD-BEAD
		void doNanoCoreForce(int i);
		///Do BOUNDARY force on some particles in residue i
		void doBoundaryForce(int i);
		void doOffsetBoundaryForce(int i);
		///Do RIGIDBEND force on some particles in residue i
		void doRigidBendForce(int i);
		///Do PULLBEAD force on some particles in residue i
		void doPullBeadForce(int i);
		void doFloatingBaseForce(int i);
		void doZTorqueForce(int i);
		void doZPowerForce(int i);
		void doBallForce(int i);
		
		///Extract neighbor info
		fourVector<T> doBeadNeighbors(int i);
		std::vector<threeVector<int> > doBeadChain(int i);
		
		///Return the CHAIN potential for residue i. 
		T doChainPotential(int i);
		///Return the DIHEADRAL potential for residue i.
		T doDihedralPotential(int i);
		///Return the Torsion Potential for residue i.
		T doTorsionPotential(int i);
		///Return the Bond potential for residue i.
		T doBondPotential(int i);
		///Return the Bend potential for residue i.
		T doBendPotential(int i);
		///Return the Bead potential for residue i.
		T doBeadPotential(int i);
		T doBeadBeadPotential(int i);
		///Return the nanocore potential for residue i. Same as BEAD, but lacks BEAD-BEAD
		T doNanoCorePotential(int i);
		///Return the BOUNDARY potential on some particles in residue i.
		T doBoundaryPotential(int i);
		T doOffsetBoundaryPotential(int i);
		///Return the RIGIDBEND potential on some particles in residue i, this is zero.
		T doRigidBendPotential(int i);
		///Return the PULLBEAD potential on some particles in residue i.
		T doPullBeadPotential(int i);
		T doFloatingBasePotential(int i);
		T doZTorquePotential(int i);
		T doZPowerPotential(int i);
		T doBallPotential(int i);
		
		///Return the change in CHAIN potential for residue i when the size changes to a new state.
		T doChainDPotential(int i, threeVector<T> aSize);
		///Return the change in DIHEADRAL potential for residue i when the size changes to a new state.
		T doDihedralDPotential(int i, threeVector<T> aSize);
		///Return the change in Torsion potential for residue i when the size changes to a new state.
		T doTorsionDPotential(int i, threeVector<T> aSize);
		///Return the change in Bond potential for residue i when the size changes to a new state.
		T doBondDPotential(int i, threeVector<T> aSize);
		///Return the change in Bend potential for residue i when the size changes to a new state.
		T doBendDPotential(int i, threeVector<T> aSize);
		///Return the change in Bead potential for residue i when the size changes to a new state.
		T doBeadDPotential(int i, threeVector<T> aSize);
		///Return the change in Rigid potential for residue i when the size changes to a new state, this is zero.
		T doRigidBendDPotential(int i, threeVector<T> aSize);
		///Return the change in Bead potential for residue i when the size changes to a new state.
		T doPullBeadDPotential(int i, threeVector<T> aSize);
		T doNanoCoreDPotential(int i, threeVector<T> aSize);
		T doBallDPotential(int i, threeVector<T> aSize);
		
	private:
		std::vector< position<T> > p;
		std::vector< threeVector<T> > v;
		std::vector< threeVector<T> > a;
		std::vector< molecule<T, fourVector<int> > > m;
		
		int nParticles;
		int nMolecules;
		
		//These are supposedly internal
		std::vector<T> twoBodyFconst;
		std::vector<T> twoBodyUconst;
		
		//constants loaded at run time
		T gamma;
		T initialTemp;
		T finalTemp;
		
		//basic properties of the particle system
		int seed;              //random number seed
		int nTypes;            //number of particle types
		
		threeVector<bool> periodic;//periodic boundaries
		T cutoff;              //cutoff for algorithms
		threeVector<T> size;   //component sizes of system
		T initialTime;         //initial time
		T finalTime;           //length of run
		T deltaT;              //time step size, depends on algorithm in process
		T storeInterval;       //period for storing data for continuation, float?
		T measureInterval;     //measure and gather data at time steps, float?
		T deltaLXY;            //Change in the length of the x-y plane
		T removeSolvent;//percentage of solvent to move to outer
		T tempStepInterval;
		//Keep this comment here:
		//COMMANDINPUTHANDLER
		twoVector<T> solventGamma;
		std::vector<T> gammaType;
		T tension;
		
		T cutoffSquared;
		
		#ifdef OMP_MOLECULES
			threeVector<T> **omp_a;
		#endif
			
		//this is a script engine (finite state machine), here are the command IDs (states)
		enum commandId {
		GAMMA,INITIALTEMP,FINALTEMP,SEED,NTYPES,NMOLECULES,NPARTICLES,PERIODIC,CUTOFF,SIZE,INITIALTIME,FINALTIME,
		DELTAT,STOREINTERVAL,MEASUREINTERVAL,TWOBODYFCONST,TWOBODYUCONST,
		POSITIONS,VELOCITIES,MOLECULE,BANANA,DELTALXY,REMOVESOLVENT,TEMPSTEPINTERVAL
			//Keep this comment here:
			//COMMANDINPUTHANDLER
			,SOLVENTGAMMA
			,GAMMATYPE
			,TENSION
			,NENUMCOMMANDID
		};
		
		//map of the states and their values (hash table)
		std::map< std::string, int > commandMap;
		std::vector<bool> commandRequired;
		std::vector<bool> commandPresent;
		std::vector<std::string> commands;
		
		int nP;//number of parameters
		int wP;//which parameter is currently being used
		int it;//which iteration of current parameter
		int molState,mol,bondState,current,bondNumber;//these ones control molecule iterations
		fourVector<int> bondBuf;//global container for a bond
};

template <typename T>
Blob<T>::Blob()
{
	//set most of these to a default value
	seed=0;
	nTypes=0;
	nParticles=0;
	nMolecules=0;
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
	solventGamma.x=0;
	solventGamma.y=0;
	gamma=0;
	tension=0;
	
	#ifdef OMP_MOLECULES
		omp_a=new threeVector<T> *[omp_get_max_threads()];
		#pragma omp parallel
		{
			//std::cerr << omp_get_thread_num() << '\n';
			omp_a[omp_get_thread_num()]=NULL;
		}
		//std::cin.get();
	#endif
	
	//list of commands this simulation utilizes
	for(int i=0;i<NENUMCOMMANDID;i++)
		commands.push_back(std::string());
	
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
	commands[GAMMATYPE]="gammaType";
	commands[TENSION]="tension";
	
	//map the commands to their indexes to make them easy to find
	for(int i=0;i<commands.size();i++)
	{
		commandMap[commands[i]]=i;
		#ifdef NO_REQUIRED_COMMANDS
			commandRequired.push_back(false);
		#else
			commandRequired.push_back(true);
		#endif
		commandPresent.push_back(false);
	}
	
	commandRequired[MEASUREINTERVAL]=false;
	commandRequired[NMOLECULES]=false;
	commandRequired[MOLECULE]=false;
	commandRequired[BANANA]=false;
	commandRequired[TWOBODYFCONST]=false;
	commandRequired[TWOBODYUCONST]=false;
	commandRequired[DELTALXY]=false;
	commandRequired[REMOVESOLVENT]=false;
	commandRequired[TEMPSTEPINTERVAL]=false;
	//Keep this comment here:
	//COMMANDINPUTHANDLER
	commandRequired[SOLVENTGAMMA]=false;
	commandRequired[GAMMATYPE]=false;
	commandRequired[GAMMA]=false;
	commandRequired[TENSION]=false;
}

//Mostly memory handling, deletes etc...
template <typename T>
Blob<T>::~Blob()
{
	#ifdef OMP_MOLECULES
		#pragma omp parallel
		{
			if(omp_a[omp_get_thread_num()]!=NULL)
				delete omp_a[omp_get_thread_num()];
		}
		delete omp_a;
	#endif
}

template <typename T>
void Blob<T>::errorChecking()
{
	if(p.size()!=a.size() || p.size()!=v.size() || a.size()!=v.size() ||
	nParticles!=p.size() || nParticles!=a.size() || nParticles!=v.size())
	{
		std::cerr << "Mismatch of parameters:\n";
		std::cerr << "\tNumber of particles assumed: " << nParticles << std::endl;
		std::cerr << "\tNumber of particle positions: " << p.size() << std::endl;
		std::cerr << "\tNumber of particle velocities: " << v.size() << std::endl;
		std::cerr << "\tNumber of particle accelerations: " << a.size() << std::endl;
		throw 0;
	}
	
	for(int i=0;i<p.size();i++)
	{
		if(p[i].x>size.x || p[i].x<0)
		{
			std::cerr << "X position of particle " << i << " is out of bounds.\n";
			throw 0;
		}
		if(p[i].y>size.y || p[i].y<0)
		{
			std::cerr << "Y position of particle " << i << " is out of bounds.\n";
			throw 0;
		}
		if(p[i].z>size.z || p[i].z<0)
		{
			std::cerr << "Z position of particle " << i << " is out of bounds.\n";
			throw 0;
		}
	}
	
	//for locating bonded particles
	std::vector< std::vector<int> > bonded;
	std::vector< std::vector<int> > bended;
	
	//set up the reserve
	for(int i=0;i<nParticles;i++)
	{
		bonded.push_back(std::vector<int>());
		bended.push_back(std::vector<int>());
	}
	
	//check all molecules
	for(int i=0;i<m.size();i++)
	{
		int nBonded=m[i].readNBond();
		fourVector<int> *bond=m[i].getBonds();
		int type=m[i].readType();
		
		for(int j=0;j<nBonded;j++)
		{
			switch(type)
			{
				case BOND:
				{
					//check the bond index ranges
					if(bond[j].s[0]>nParticles || bond[j].s[1]>nParticles || bond[j].s[0]<0 || bond[j].s[1]<0)
					{
						std::cerr << "BOND Molecule " << i << ", bond " << j << " is out of bounds!\n";
						throw 0;
					}
					
					//collect bonds to check for duplication.
					bonded[bond[j].s[0]].push_back(bond[j].s[1]);
					bonded[bond[j].s[1]].push_back(bond[j].s[0]);
					
				}
				break;
				
				case BEND:
				{
					//check the bend index ranges
					if(bond[j].s[0]>nParticles || bond[j].s[1]>nParticles || 
						bond[j].s[2]>nParticles || bond[j].s[0]<0 ||
						bond[j].s[1]<0 || bond[j].s[2]<0)
					{
						std::cerr << "BEND Molecule " << i << ", bond " << j << " is out of bounds!\n";
						throw 0;
					}
					
					//collect bends to check for duplication.
					bended[bond[j].s[0]].push_back(bond[j].s[1]);
					bended[bond[j].s[0]].push_back(bond[j].s[2]);
					
					bended[bond[j].s[2]].push_back(bond[j].s[1]);
					bended[bond[j].s[2]].push_back(bond[j].s[0]);
				}
				break;
				
				case CHAIN:
				{
					if(bond[j].s[START]+bond[j].s[NCHAINS]*bond[j].s[CHAINLENGTH]>nParticles ||
						bond[j].s[START]+bond[j].s[NCHAINS]*bond[j].s[CHAINLENGTH]<0)
					{
						std::cerr << "CHAIN Molecule " << i << ", bond " << j << " is out of bounds!\n";
						throw 0;
					}
				}
				break;
				
				default:
					//does nothing
				break;
			}
		}
	}
	
	#ifdef WARNINGS_ENABLED
	for(int i=0;i<bonded.size();i++)
	{
		for(int j=0;j<bonded[i].size();j++)
		{
			int nDuplicateBonds=0;
			for(int k=j+1;k<bonded[i].size();k++)
				if(bonded[i][j]==bonded[i][k])
					nDuplicateBonds++;
			if(nDuplicateBonds>0)
			{
				std::cerr << "Warning(Blob): Particle " << i << " has " << nDuplicateBonds << " bonded duplicates of ";
				std::cerr << "particle " << bonded[i][j] << ".\n";
			}
		}
	}
	
	for(int i=0;i<bended.size();i++)
	{
		for(int j=0;j<bended[i].size();j+=2)
		{
			int nDuplicateBends=0;
			for(int k=j+2;k<bended[i].size();k+=2)
				if(bended[i][j]==bended[i][k] && bended[i][j+1]==bended[i][k+1])
					nDuplicateBends++;
			if(nDuplicateBends>0)
			{
				std::cerr << "Warning(Blob): Particle " << i << " has " << nDuplicateBends << " bend duplicates of ";
				std::cerr << "particles " << bonded[i][j] << " and " << bonded[i][j+1] << ".\n";
			}
		}
	}
	#endif
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
			for(int i=0;i<NENUMCOMMANDID;i++)
			{
				if(!commandPresent[i] && commandRequired[i])
				{
					#ifdef WARNINGS_ENABLED
						std::cerr << "Warning (Blob): Script missing " << commands[i] << "!\n";
						//std::cerr << "Press any key to continue loading or ctrl-z to stop!";
						//std::cin.get();
					#endif
					//Maybe in the future I should add warning and error levels for crap like this
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
			std::cerr << in << " is not a recognized command!\n";
			std::cerr << "Locate command before " << in << "!\n";
			std::cerr << "You probably have too many positions, velocities, or molecules!\n";
			std::cerr << "Also, check nMolecules and nParticles.\n";
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
					std::cerr << "Value after gamma is not a numerical type!\n";
					std::cerr << "Format:\ngamma [float]\n";
					throw 0;
				}
				commandPresent[wP]=true;
				wP=-1;//clear command index
			}
			break;
			
		case INITIALTEMP://initialTemp
			if(it>0)
			{
				reformatString >> initialTemp;
				if(reformatString.fail())
				{
					std::cerr << "Value after initialTemp is not a numerical type!\n";
					std::cerr << "Format:\ninitialTemp [float]\n";
					throw 0;
				}
				commandPresent[wP]=true;
				wP=-1;
			}
			break;
			
		case FINALTEMP://finalTemp
			if(it>0)
			{
				reformatString >> finalTemp;
				if(reformatString.fail())
				{
					std::cerr << "Value after finalTemp is not a numerical type!\n";
					std::cerr << "Format:\nfinalTemp [float]\n";
					throw 0;
				}
				commandPresent[wP]=true;
				wP=-1;
			}
			break;
			
		case SEED://seed
			if(it>0)
			{
				reformatString >> seed;
				if(reformatString.fail())
				{
					std::cerr << "Value after seed is not a numerical type!\n";
					std::cerr << "Format:\nseed [integer]\n";
					throw 0;
				}
				commandPresent[wP]=true;
				wP=-1;
			}
			break;
			
		case NTYPES://nTypes
			if(it>0)
			{
				reformatString >> nTypes;
				if(reformatString.fail())
				{
					std::cerr << "Value after nTypes is not a numerical type!\n";
					std::cerr << "Format:\nnTypes [integer]\n";
					throw 0;
				}
				
				commandPresent[wP]=true;
				wP=-1;
			}
			break;
			
		case NMOLECULES://nMolecules
			if(it>0)
			{
				reformatString >> nMolecules;
				if(reformatString.fail())
				{
					std::cerr << "Value after nMolecules is not a numerical type!\n";
					std::cerr << "Format:\nnMolecules [integer]\n";
					throw 0;
				}
				commandPresent[wP]=true;
				wP=-1;
			}
			break;
			
		case NPARTICLES://nParticles
			if(it>0)
			{
				reformatString >> nParticles;
				if(reformatString.fail())
				{
					std::cerr << "Value after nParticles is not a numerical type!\n";
					std::cerr << "Format:\nnParticles [integer]\n";
					throw 0;
				}
				commandPresent[wP]=true;
				wP=-1;
			}
			break;
			
		case PERIODIC://periodic
			if(it>0)
			{
				periodic.x=true;
				periodic.y=true;
				periodic.z=true;
				commandPresent[wP]=true;
				wP=-1;
			}
			break;
			
		case CUTOFF://cutoff
			if(it>0)
			{
				reformatString >> cutoff;
				if(reformatString.fail())
				{
					std::cerr << "Value after cutoff is not a numerical type!\n";
					std::cerr << "Format:\ncutoff [float]\n";
					throw 0;
				}
				cutoffSquared=cutoff*cutoff;
				commandPresent[wP]=true;
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
					std::cerr << "Value after size is not a numerical type!\n";
					std::cerr << "Format:\nsize [float] [float] [float]\n";
					throw 0;
				}
				if(it==3)
				{
					commandPresent[wP]=true;
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
					std::cerr << "Value after initialTime is not a numerical type!\n";
					std::cerr << "Format:\ninitialTime [float]\n";
					throw 0;
				}
				commandPresent[wP]=true;
				wP=-1;
			}
			break;
			
		case FINALTIME://finalTime
			if(it>0)
			{
				reformatString >> finalTime;
				if(reformatString.fail())
				{
					std::cerr << "Value after finalTime is not a numerical type!\n";
					std::cerr << "Format:\nfinalTime [float]\n";
					throw 0;
				}
				commandPresent[wP]=true;
				wP=-1;
			}
			break;
			
		case DELTAT://deltaT
			if(it>0)
			{
				reformatString >> deltaT;
				if(reformatString.fail())
				{
					std::cerr << "Value after deltaT is not a numerical type!\n";
					std::cerr << "Format:\ndeltaT [float]\n";
					throw 0;
				}
				commandPresent[wP]=true;
				wP=-1;
			}
			break;
			
		case STOREINTERVAL://storeInterval
			if(it>0)
			{
				reformatString >> storeInterval;
				if(reformatString.fail())
				{
					std::cerr << "Value after storeInterval is not a numerical type!\n";
					std::cerr << "Format:\nstoreInterval [float]\n";
					throw 0;
				}
				commandPresent[wP]=true;
				wP=-1;
			}
			break;
			
		case MEASUREINTERVAL://measureInterval
			if(it>0)
			{
				reformatString >> measureInterval;
				if(reformatString.fail())
				{
					std::cerr << "Value after measureInterval is not a numerical type!\n";
					std::cerr << "Format:\nmeasureInterval [float]\n";
					throw 0;
				}
				commandPresent[wP]=true;
				wP=-1;
			}
			break;
			
		case TWOBODYFCONST://twoBodyFconst
			if(it>0)//more advanced "it" usage
			{
				T buf;
				reformatString >> buf;
				
				if(reformatString.fail())//you broke it...
				{
					std::cerr << "Value after twoBodyFconst is not a numerical type!\n";
					std::cerr << "Check that nTypes is correct!\n";
					std::cerr << "Format:\ntwoBodyFconst\n [float]...\n";
					throw 0;
				}
				
				if(twoBodyFconst.size()<=it-1)
					twoBodyFconst.push_back(buf);
				else
					twoBodyFconst[it-1]=buf;
				
				if(it==nTWOBODYFCONST*nTypes*nTypes)//last number
				{
					commandPresent[wP]=true;
					wP=-1;
				}
			}
			else
			{
				if(!commandPresent[NTYPES])
				{
					std::cerr << "nTypes was not present before twoBodyFconst!\n";
					throw 0;
				}
			}
			break;
			
		case TWOBODYUCONST://twoBodyUconst
			if(it>0)
			{
				T buf;
				reformatString >> buf;
				
				if(reformatString.fail())
				{
					std::cerr << "Value after twoBodyUconst is not a numerical type!\n";
					std::cerr << "Check that nTypes is correct!\n";
					std::cerr << "Format:\ntwoBodyUconst\n [float]...\n";
					throw 0;
				}
				
				if(twoBodyUconst.size()<=it-1)
					twoBodyUconst.push_back(buf);
				else
					twoBodyUconst[it-1]=buf;
				
				if(it==nTWOBODYUCONST*nTypes*nTypes)
				{
					commandPresent[wP]=true;
					wP=-1;
				}
			}
			else
			{
				if(!commandPresent[NTYPES])
				{
					std::cerr << "nTypes was not present before twoBodyUconst!\n";
					throw 0;
				}
			}
			break;
			
		case POSITIONS://positions
			if(it>0)
			{
				int particle=int((it-1)/4);
				if(p.size()<=particle)
				{
					position<T> buf;
					p.push_back(buf);
				}
				
				if((it-1)%4==0)//these are the only values
					reformatString >> p[particle].type;
				else
					reformatString >> p[particle].s[(it-1)%4-1];
				
				if(reformatString.fail())
				{
					std::cerr << "Value after positions is not a numerical type!\n";
					std::cerr << "This likely means that nParticles is larger than the number of positions in file.\n";
					std::cerr << "Format:\npositions\n [integer] [float] [float] [float]\n...\n";
					throw 0;
				}
				
				if(it==nParticles*4)
				{
					commandPresent[wP]=true;
					wP=-1;
				}
			}
			else
			{
				if(!commandPresent[NPARTICLES])
				{
					std::cerr << "nParticles was not present before positions!\n";
					throw 0;
				}
			}
			break;
			
		case VELOCITIES://velocities
			if(it>0)
			{
				int particle=int((it-1)/3);
				if(v.size()<=particle)
				{
					threeVector<T> buf;
					v.push_back(buf);
				}
				if(a.size()<=particle)
				{
					threeVector<T> buf;
					a.push_back(buf);
				}
				
				reformatString >> v[particle].s[(it-1)%3];
				
				if(reformatString.fail())
				{
					std::cerr << "Value after velocities is not a numerical type!\n";
					std::cerr << "This likely means that nParticles is larger than the number of velocities in file.\n";
					std::cerr << "Format:\nvelocities\n [float] [float] [float]\n...\n";
					throw 0;
				}
				if(it==nParticles*3)
				{
					commandPresent[wP]=true;
					wP=-1;
				}
			}
			else
			{
				if(!commandPresent[NPARTICLES])
				{
					std::cerr << "nParticles was not present before velocities!\n";
					throw 0;
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
					if(m.size()<=mol)
					{
						m.push_back(molecule<T, fourVector<int> >());
					}
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
							nMoleculeConstants=nBENDCONST;
							nMoleculeParticles=3;//3 particles per bond
							break;
						case CHAIN:
							nMoleculeConstants=nCHAINCONST;
							nMoleculeParticles=3;//3 index relations
							break;
						case BEAD:
							nMoleculeParticles=1;//1 index
							if(nTypes<=0)
							{
								std::cerr << "Error (Blob): nTypes undefined for BEAD molecule!" << std::endl;
								throw 1;
							}
							nMoleculeConstants=nBEADCONST*nTypes*nTypes;//radius, type
							break;
						case SOLID:
							nMoleculeParticles=1;//1 index
							nMoleculeConstants=1;//r0
							break;
						case BOUNDARY:
							nMoleculeParticles=1;
							nMoleculeConstants=nBOUNDARYCONST;
							break;
						case OFFSET_BOUNDARY:
							nMoleculeParticles=1;
							nMoleculeConstants=nOFFSET_BOUNDARYCONST;
							break;
						case RIGIDBEND:
							nMoleculeConstants=nRIGIDBENDCONST;
							nMoleculeParticles=2;//2 particles per bond
							break;
						case PULLBEAD:
							nMoleculeConstants=nPULLBEADCONST;
							nMoleculeParticles=1;//1 particles per bond
							break;
						case FLOATING_BASE:
							nMoleculeConstants=nFLOATING_BASECONST*nTypes;
							nMoleculeParticles=1;
							break;
						case ZTORQUE:
							nMoleculeConstants=nZTORQUECONST;
							nMoleculeParticles=3;//3 index relations
							break;
						case ZPOWERPOTENTIAL:
							nMoleculeConstants=nZPOWERPOTENTIALCONST;
							nMoleculeParticles=2;//2 index relations
							break;
						case NANOCORE:
							nMoleculeConstants=nBEADCONST;
							nMoleculeParticles=1;//2 index relations
							//Correction for this type after nBonds read
							if(molState>0) nMoleculeConstants*=m[mol].readNBondAlloc();
							break;
						case BALL:
							nMoleculeConstants=nBONDCONST;
							nMoleculeParticles=2;//2 particles per bond
							break;
						default://default formatting
							//m[mol].allocConstant(0);
							nMoleculeConstants=0;
							nMoleculeParticles=0;
							#ifdef WARNINGS_ENABLED
								std::cerr << "Warning (Blob): default formatting while reading molecule!\n";
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
								T buf;
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
					std::cerr << "Value after molecules is not a numerical type!\n";
					std::cerr << "Read: " << in << '\n';
					throw 0;
				}
				
				if(mol==nMolecules)
				{
					commandPresent[MOLECULE]=true;
					wP=-1;
				}
			}
			else
			{
				if(!commandPresent[NMOLECULES])
				{
					std::cerr << "nMolecules was not present before molecules!\n";
					throw 0;
				}
				current=-1;
				mol=0;//set number of mols to 0
				
				//while(m.size()<nMolecules)
				//	m.push_back(molecule<T, fourVector<int> >());
			}
			break;
			
		case BANANA://banana!
			//this is an example of a flag that can be set, notice the distinct lack of an iterator
			std::cerr << "Monkeys must have messed with your script because I found a banana!\n";
			commandPresent[wP]=true;
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
					std::cerr << "Value after deltaLXY is not a numerical type!\n";
					std::cerr << "Format:\n\tdeltaLXY [float]\n";
					throw 0;
				}
				commandPresent[wP]=true;
				wP=-1;
			}
			break;
		
		case REMOVESOLVENT:
			if(it>0)
			{
				reformatString >> removeSolvent;
				if(reformatString.fail())
				{
					std::cerr << "Value after removeSolvent is not a numerical type!\n";
					std::cerr << "Format:\n\tremoveSolvent [float]\n";
					std::cerr << "Float is the fraction of the total solvent, not just the inner volume of solvent\n";
					throw 0;
				}
				commandPresent[wP]=true;
				wP=-1;
			}
			break;
			
		case TEMPSTEPINTERVAL:
			if(it>0)
			{
				reformatString >> tempStepInterval;
				if(reformatString.fail())
				{
					std::cerr << "Value after tempStepInterval is not a numerical type!\n";
					std::cerr << "Format:\ntempStepInterval [float]\n";
					throw 0;
				}
				commandPresent[wP]=true;
				wP=-1;
			}
			break;
		
		//Keep this comment here:
		//COMMANDINPUTHANDLER
		case SOLVENTGAMMA:
			if(it>0)
			{
				reformatString >> solventGamma.s[it-1];
				if(reformatString.fail())
				{
					std::cerr << "Value after solventGamma is not a numerical type!\n";
					std::cerr << "Format:\nsolventGamma [inner float] [outer float]\n";
					throw 0;
				}
				if(it==2)
				{
					commandPresent[wP]=true;
					wP=-1;
				}
			}
			break;
		case GAMMATYPE://gamma
			if(it>0)
			{
				reformatString >> gammaType[it-1];
				if(reformatString.fail())
				{
					std::cerr << "Value after gammaType is not a numerical type!\n";
					std::cerr << "Format:\ngammaType [float] [float] ...\n";
					throw 0;
				}
				if(it==nTypes)
				{
					commandPresent[wP]=true;
					wP=-1;
				}
			}
			else
			{
				if(!commandPresent[NTYPES])
				{
					std::cerr << "nTypes was not present before gammaType!\n";
					throw 0;
				}
			}
			break;
		case TENSION://tension
			if(it>0)
			{
				reformatString >> tension;
				if(reformatString.fail())
				{
					std::cerr << "Value after tension is not a numerical type!\n";
					std::cerr << "Format:\ntension [float]\n";
					throw 0;
				}
				commandPresent[wP]=true;
				wP=-1;
			}
			break;
		
		default://other commands
			std::cerr << "This command doesn't yet exist in the input() handler.\n";
			throw 0;
			//break;
	}
	it++;//increment iterator
	return false;
}

template <typename T>
std::string Blob<T>::output()
{
	if(wP==-1)
		wP=0;
	if(wP==NENUMCOMMANDID)
		return END_INPUT;
	
	std::stringstream out("");
	if(it==0)
	{
		//This will enable the command to output with or without it being required
		if(commandPresent[wP])
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
				commandRequired[MOLECULE]=true;
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
			if(sizeof(T)==8)
				out << std::setprecision(15) << size.s[it-1];
			else
				out << std::setprecision(8) << size.s[it-1];
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
						case BEAD:
							nMoleculeParticles=1;//1 index
							if(nTypes<=0)
							{
								std::cerr << "Error (Blob): nTypes undefined for BEAD molecule!" << std::endl;
								throw 1;
							}
							nMoleculeConstants=nBEADCONST*nTypes*nTypes;//radius, type
							break;
						case SOLID:
							nMoleculeParticles=1;//1 index
							nMoleculeConstants=1;//r0
							break;
						case BOUNDARY:
							nMoleculeParticles=1;
							nMoleculeConstants=nBOUNDARYCONST;
							break;
						case OFFSET_BOUNDARY:
							nMoleculeParticles=1;
							nMoleculeConstants=nOFFSET_BOUNDARYCONST;
							break;
						case RIGIDBEND:
							nMoleculeConstants=nRIGIDBENDCONST;
							nMoleculeParticles=2;//2 particles per bond
							break;
						case PULLBEAD:
							nMoleculeConstants=nPULLBEADCONST;
							nMoleculeParticles=1;//1 particles per bond
							break;
						case FLOATING_BASE:
							nMoleculeConstants=nFLOATING_BASECONST*nTypes;
							nMoleculeParticles=1;
							break;
						case ZTORQUE:
							nMoleculeConstants=nZTORQUECONST;
							nMoleculeParticles=3;//3 index relations
							break;
						case ZPOWERPOTENTIAL:
							nMoleculeConstants=nZPOWERPOTENTIALCONST;
							nMoleculeParticles=2;//2 index relations
							break;
						case NANOCORE:
							nMoleculeConstants=nBEADCONST*m[mol].readNBond();
							nMoleculeParticles=1;//2 index relations
							break;
						case BALL:
							nMoleculeParticles=2;//2 particles per bond
							nMoleculeConstants=nBONDCONST;
							break;
						default://default formatting
							//m[mol].allocConstant(0);
							nMoleculeParticles=0;
							#ifdef WARNINGS_ENABLED
							std::cerr << "Warning (Blob): default formatting while writing molecule!\n";
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
								if((molState)%nBEADCONST==0 && m[mol].readType()==BEAD)
									out << std::endl;
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
				/*//I have no idea why this was here before...
				if(commandPresent[NMOLECULES]<1)
				{
					std::cerr << "nMolecules was not present before molecules!\n";
					throw 0;
				}
				if(commandPresent[MOLECULE]==0)
				{
					m=new molecule<T,fourVector<int> >[nMolecules];
					commandPresent[MOLECULE]=nMolecules;
				}
				if(commandPresent[MOLECULE]==-1)
				{
					delete[] m;
					m=new molecule<T,fourVector<int> >[nMolecules];
					commandPresent[MOLECULE]=nMolecules;
				}
				*/
				current=-1;
				mol=0;//set number of mols to 0
			}
			it++;
			if(mol==m.size())
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
			case SOLVENTGAMMA://gamma for inner and outer solvent
			out << std::setprecision(15) << solventGamma.s[it-1];
			it++;
			if(it-1==2)
			{
				wP++;
				it=0;
			}
			break;
			
			case GAMMATYPE://size
			out << std::setprecision(15) << gammaType[it-1];
			it++;
			if(it-1==nTypes)
			{
				wP++;
				it=0;
			}
			break;
			
			case TENSION://tension
			out << std::setprecision(15) << tension;
			wP++;
			it=0;
			break;
			
			default://other commands
			std::cerr << "This command doesn't yet exist in the output() handler.\n";
			throw 0;
			//it=0;
			//wP++;
			//return "\n";
			//break;
		}
	}
	return out.str();
}


//These are just conveniences, so you don't have to worry about typing the same crap over and over
template <typename T>
inline void Blob<T>::doForce(int i, int j, threeVector<T> &minImg)
{
	threeVector<T> d;
	
	d.x=p[i].x-p[j].x+minImg.x;
	d.y=p[i].y-p[j].y+minImg.y;
	d.z=p[i].z-p[j].z+minImg.z;
	
	int cindex=nTWOBODYFCONST*((p[i].type*nTypes)+p[j].type);
	
	nonBondedF(d,a[i],a[j],&twoBodyFconst[cindex]);
}

template <typename T>
T Blob<T>::doPotential(int i, int j, threeVector<T> &minImg)
{
	threeVector<T> d;
	
	d.x=p[i].x-p[j].x+minImg.x;
	d.y=p[i].y-p[j].y+minImg.y;
	d.z=p[i].z-p[j].z+minImg.z;
	
	int cindex=nTWOBODYUCONST*((p[i].type*nTypes)+p[j].type);
	
	return nonBondedP(d,&twoBodyUconst[cindex]);
}

//These are molecule related things
template <typename T>
void Blob<T>::doChainForce(int i)
{
	//this is just going to check it without an error condition, it ignores it
	if(m[i].readType()==CHAIN)
	{
		
		//Go through all bond descriptions
		for(int j=0;j<m[i].readNBond();j++)
		{
			fourVector<int> *bond=m[i].getBonds();
			//bond info
			int start=bond[j].s[START];
			int nChains=bond[j].s[NCHAINS];
			int length=bond[j].s[CHAINLENGTH];
			//go through all chain lengths
			#pragma omp parallel for
			for(int k=start; k<start+length*nChains; k+=length)
			{
				threeVector<T> da,db;
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
					
					harmonicF(da, a[first], a[second], &(m[i].getConstants()[CHAINBOND]));
					//CounterForce(cutoffSquared, nTypes, da, a[first], a[second], &twoBodyFconst[0], p[first].type, p[second].type);
					bendF(da, db, a[first], a[second], a[third], &(m[i].getConstants()[CHAINBEND]));
				}
				//last three are a special case, requires bond between last two particles in chain
				int first=k+length-3;
				int second=k+length-2;
				int third=k+length-1;
				
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
				
				harmonicF(da, a[first], a[second], &(m[i].getConstants()[CHAINBOND]));
				harmonicF(db, a[second], a[third], &(m[i].getConstants()[CHAINBOND]));
				//CounterForce(cutoffSquared, nTypes, da, a[first], a[second], &twoBodyFconst[0],p[first].type,p[second].type);
				//CounterForce(cutoffSquared, nTypes, db, a[second], a[third], &twoBodyFconst[0],p[second].type,p[third].type);
				bendF(da, db, a[first], a[second], a[third], &(m[i].getConstants()[CHAINBEND]));
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
	#ifdef OMP_MOLECULES
		#pragma omp parallel
		{
			for(int j=0;j<nParticles;j++)
			{
				omp_a[omp_get_thread_num()][j]=0;
			}
		}
	#endif
	//this is just going to check it without an error condition, it ignores it
	if(m[i].readType()==BOND)
	{
		fourVector<int> *bond=m[i].getBonds();
		
		#ifdef OMP_MOLECULES
			#pragma omp parallel for
		#endif
		for(int j=0;j<m[i].readNBond();j++)
		{
			int first=bond[j].s[0];
			int second=bond[j].s[1];
			threeVector<T> d;
			
			d.x=p[first].x-p[second].x;
			if(d.x>size.x/2.0) d.x-=size.x;
			if(d.x<-size.x/2.0) d.x+=size.x;
			d.y=p[first].y-p[second].y;
			if(d.y>size.y/2.0) d.y-=size.y;
			if(d.y<-size.y/2.0) d.y+=size.y;
			d.z=p[first].z-p[second].z;
			if(d.z>size.z/2.0) d.z-=size.z;
			if(d.z<-size.z/2.0) d.z+=size.z;
			#ifdef OMP_MOLECULES
				harmonicF(d, omp_a[omp_get_thread_num()][first], omp_a[omp_get_thread_num()][second], m[i].getConstants());
			#else
				harmonicF(d, a[first], a[second], m[i].getConstants());
			#endif
		}
	}
	#ifdef OMP_MOLECULES
		#pragma omp parallel for
		for(int j=0;j<nParticles;j++)
		{
			for(int k=0;k<omp_get_max_threads();k++)
			{
				a[j]+=omp_a[k][j];
				//std::cerr << omp_a[j][i].x << '\t' << omp_a[j][i].y << '\t' << omp_a[j][i].z << '\n';
			}
			//std::cin.get();
		}
	#endif
}

template <typename T>
void Blob<T>::doBallForce(int i)
{
	//this is just going to check it without an error condition, it ignores it
	if(m[i].readType()==BALL)
	{
		fourVector<int> *bond=m[i].getBonds();
		
		int first=bond[0].s[0];
		T ax=0, ay=0, az=0;
		#pragma omp parallel for reduction(+:ax) reduction(+:ay) reduction(+:az)
		for(int j=0;j<m[i].readNBond();j++)
		{
			threeVector<T> aFirst=0;
			int second=bond[j].s[1];
			threeVector<T> d;
			
			d.x=p[first].x-p[second].x;
			if(d.x>size.x/2.0) d.x-=size.x;
			if(d.x<-size.x/2.0) d.x+=size.x;
			d.y=p[first].y-p[second].y;
			if(d.y>size.y/2.0) d.y-=size.y;
			if(d.y<-size.y/2.0) d.y+=size.y;
			d.z=p[first].z-p[second].z;
			if(d.z>size.z/2.0) d.z-=size.z;
			if(d.z<-size.z/2.0) d.z+=size.z;
			harmonicHalfF(d, aFirst, a[second], m[i].getConstants());
			ax+=aFirst.x;
			ay+=aFirst.y;
			az+=aFirst.z;
		}
		a[first].x+=ax;
		a[first].y+=ay;
		a[first].z+=az;
	}
}



template <typename T>
void Blob<T>::doBendForce(int i)
{
	#ifdef OMP_MOLECULES
		#pragma omp parallel
		{
			for(int j=0;j<nParticles;j++)
			{
				omp_a[omp_get_thread_num()][j]=0;
			}
		}
	#endif
	//this is just going to check it without an error condition, it ignores it
	if(m[i].readType()==BEND)
	{
		fourVector<int> *bond=m[i].getBonds();
		#ifdef OMP_MOLECULES
			#pragma omp parallel for
		#endif
		for(int j=0;j<m[i].readNBond();j++)
		{
			int first=bond[j].s[0];
			int second=bond[j].s[1];
			int third=bond[j].s[2];
			threeVector<T> da,db;
			
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
			
			#ifdef OMP_MOLECULES
				bendF(da, db, omp_a[omp_get_thread_num()][first], omp_a[omp_get_thread_num()][second], omp_a[omp_get_thread_num()][third], m[i].getConstants());
			#else
				bendF(da, db, a[first], a[second], a[third], m[i].getConstants());
			#endif
		}
	}
	#ifdef OMP_MOLECULES
		#pragma omp parallel for
		for(int j=0;j<nParticles;j++)
		{
			for(int k=0;k<omp_get_max_threads();k++)
			{
				a[j]+=omp_a[k][j];
				//std::cerr << omp_a[j][i].x << '\t' << omp_a[j][i].y << '\t' << omp_a[j][i].z << '\n';
			}
			//std::cin.get();
		}
	#endif
}

//the range of this might change
template <typename T>
void Blob<T>::doBeadForce(int i)
{
	//this is just going to check it without an error condition, it ignores it
	if(m[i].readType()==BEAD)
	{
		//x=index, y=bead-bead type, z=bead-particle type, t=radius
		std::vector<threeVector<int> > beadInfo;
		std::vector<T> beadRadius;
		
		for(int j=i;j<m.size();j++)
		{
			if(m[j].readType()==BEAD)
			{
				for(int k=0;k<m[j].readNBond();k++)
				{
					//x is index, y is type 1, z is type 2
					threeVector<int> info;
					info.x=m[j].getBonds()[k].x;
					info.y=0;//m[i].getConstants()[1];
					info.z=p[info.x].type;//we will do this for now, but we may want multiple types later
					//nBEADCONST*(info.z*nTypes)+info.z);
					beadRadius.push_back(m[i].getConstants()[BEADRADIUS]);
					beadInfo.push_back(info);
					//std::cout << info.x << std::endl;
				}
			}
		}
		T *C=m[i].getConstants();	
		
		//for all other beads
		for(int j=0;j<m[i].readNBond();j++)
		{
			for(int k=j+1;k<beadInfo.size();k++)
			{
				threeVector<T> d;
				int first=beadInfo[j].x;
				int second=beadInfo[k].x;
				
				d.x=p[first].x-p[second].x;
				if(d.x>size.x/2.0) d.x-=size.x;
				if(d.x<-size.x/2.0) d.x+=size.x;
				d.y=p[first].y-p[second].y;
				if(d.y>size.y/2.0) d.y-=size.y;
				if(d.y<-size.y/2.0) d.y+=size.y;
				d.z=p[first].z-p[second].z;
				if(d.z>size.z/2.0) d.z-=size.z;
				if(d.z<-size.z/2.0) d.z+=size.z;
				
				T cutoffSquared2=beadRadius[j]+beadRadius[k]+2.0;
				cutoffSquared2*=cutoffSquared2;
				
				int cindex=nBEADCONST*((beadInfo[j].z*nTypes)+beadInfo[k].z);
				beadBeadForce(d, a[first], a[second], &C[cindex], cutoffSquared2);
			}
		}
		
		T cutoffRad=m[i].getConstants()[0];//radius+rc
		T cutoffSquaredRad=cutoffRad*cutoffRad;
		
		//p, a, C, C, nParticles,
		//nTypes,size,periodic,beadRadius[0]);
		if(m[i].readNBond()>20)
		{
		threeVector<int> nCells;
		nCells.x=size.x/(cutoffRad*2.0);
		nCells.y=size.y/(cutoffRad*2.0);
		nCells.z=size.z/(cutoffRad*2.0);
		if(nCells.x<3) std::cerr << "Warning: nCells.x is " << nCells.x << std::endl;
		if(nCells.y<3) std::cerr << "Warning: nCells.y is " << nCells.y << std::endl;
		if(nCells.z<3) std::cerr << "Warning: nCells.z is " << nCells.z << std::endl;
		threeVector<T> cSize;
		cSize.x=size.x/nCells.x;
		cSize.y=size.y/nCells.y;
		cSize.z=size.z/nCells.z;
		auto iMap=createIMap(&p[0],nParticles,nCells,cSize);
		
		//for all other particles
		#pragma omp parallel for
		for(int j=0;j<m[i].readNBond();j++)
		{
			int first=beadInfo[j].x;
			auto iCell=getCell(p[first],cSize);
			int cHash=hash(iCell,nCells);
			auto neigh=neighIndices(cHash, nCells);
			T ax=0, ay=0, az=0;
			//#pragma parallel for reduction(+:ax,ay,az)
			//for(int k=0;k<this->nParticles;k++)
			for(auto nHash:neigh)
			{
			auto iMapIt=iMap.find(nHash);
			if(iMapIt!=iMap.end())
			for(auto &k:iMapIt->second)
			{
				int second=k;
				if(first!=second)
				{
					threeVector<T> d;
					threeVector<T> aBead=0;
					d.x=p[first].x-p[second].x;
					if(d.x>size.x/2.0) d.x-=size.x;
					if(d.x<-size.x/2.0) d.x+=size.x;
					d.y=p[first].y-p[second].y;
					if(d.y>size.y/2.0) d.y-=size.y;
					if(d.y<-size.y/2.0) d.y+=size.y;
					d.z=p[first].z-p[second].z;
					if(d.z>size.z/2.0) d.z-=size.z;
					if(d.z<-size.z/2.0) d.z+=size.z;
					
					int cindex=nBEADCONST*((beadInfo[j].z*nTypes)+p[k].type);
					beadForce(d, aBead, a[second], &(C[cindex]), cutoffSquaredRad);
					ax+=aBead.x;
					ay+=aBead.y;
					az+=aBead.z;
				}
			}
			}
			a[first].x+=ax;
			a[first].y+=ay;
			a[first].z+=az;
		}
		}
		else
		{
		std::array<int,20> pExclude;
		//std::vector<int> pExclude;
		for(int j=0;j<m[i].readNBond();j++)
			//pExclude.push_back(beadInfo[j].x);
			pExclude[j]=beadInfo[j].x;
		for(int j=m[i].readNBond();j<20;j++)
			pExclude[j]=-1;
		
		//for all other particles
		for(int j=0;j<m[i].readNBond();j++)
		{
			int first=beadInfo[j].x;
			T ax=0, ay=0, az=0;
			#pragma omp parallel for reduction(+:ax,ay,az)
			for(int k=0;k<this->nParticles;k++)
			{
				int second=k;
				//if(first!=second)
				if(std::find(pExclude.begin(),pExclude.end(),k)==pExclude.end())
				{
					threeVector<T> d;
					threeVector<T> aBead=0;
					d.x=p[first].x-p[second].x;
					if(d.x>size.x/2.0) d.x-=size.x;
					if(d.x<-size.x/2.0) d.x+=size.x;
					d.y=p[first].y-p[second].y;
					if(d.y>size.y/2.0) d.y-=size.y;
					if(d.y<-size.y/2.0) d.y+=size.y;
					d.z=p[first].z-p[second].z;
					if(d.z>size.z/2.0) d.z-=size.z;
					if(d.z<-size.z/2.0) d.z+=size.z;
					
					int cindex=nBEADCONST*((beadInfo[j].z*nTypes)+p[k].type);
					beadForce(d, aBead, a[second], &(C[cindex]), cutoffSquaredRad);
					ax+=aBead.x;
					ay+=aBead.y;
					az+=aBead.z;
				}
			}
			a[first].x+=ax;
			a[first].y+=ay;
			a[first].z+=az;
		}
		}
	}
}

//the range of this might change
template <typename T>
void Blob<T>::doNanoCoreForce(int i)
{
	//this is just going to check it without an error condition, it ignores it
	if(m[i].readType()==NANOCORE)
	{	
		if(m[i].readNBond()<50)
		{
		for(int j=0;j<m[i].readNBond();j++)
		{
			int pIndex=m[i].getBonds()[j].x;
			T *C=&m[i].getConstants()[j*nBEADCONST];
			T beadRadius=C[BEADRADIUS];
			
			//for all other beads
			// This is the part we are skipping compared to BEAD type
			
			T cutoffRad=C[0];//radius+rc
			T cutoffSquaredRad=cutoffRad*cutoffRad;
			
			//for all other particles
			T ax=0, ay=0, az=0;
			#pragma omp parallel for reduction(+:ax,ay,az)
			for(int k=0;k<this->nParticles;k++)
			{
				if(k!=pIndex)
				{
				threeVector<T> d;
				threeVector<T> aBead=0;
				d.x=p[pIndex].x-p[k].x;
				if(d.x>size.x/2.0) d.x-=size.x;
				if(d.x<-size.x/2.0) d.x+=size.x;
				d.y=p[pIndex].y-p[k].y;
				if(d.y>size.y/2.0) d.y-=size.y;
				if(d.y<-size.y/2.0) d.y+=size.y;
				d.z=p[pIndex].z-p[k].z;
				if(d.z>size.z/2.0) d.z-=size.z;
				if(d.z<-size.z/2.0) d.z+=size.z;
				beadForce(d, aBead, a[k], C, cutoffSquaredRad);
				ax+=aBead.x;
				ay+=aBead.y;
				az+=aBead.z;
				}
			}
			a[pIndex].x+=ax;
			a[pIndex].y+=ay;
			a[pIndex].z+=az;
		}
		}
		else
		{
		int cutoffRad=0;
		for(int j=0;j<m[i].readNBond();j++)
		{
			T *C=&m[i].getConstants()[j*nBEADCONST];
			T rad=C[0];
			if(rad>cutoffRad) cutoffRad=rad;
		}
		
		threeVector<int> nCells;
		nCells.x=size.x/(cutoffRad*2.0);
		nCells.y=size.y/(cutoffRad*2.0);
		nCells.z=size.z/(cutoffRad*2.0);
		if(nCells.x<3) std::cerr << "Warning: nCells.x is " << nCells.x << std::endl;
		if(nCells.y<3) std::cerr << "Warning: nCells.y is " << nCells.y << std::endl;
		if(nCells.z<3) std::cerr << "Warning: nCells.z is " << nCells.z << std::endl;
		threeVector<T> cSize;
		cSize.x=size.x/nCells.x;
		cSize.y=size.y/nCells.y;
		cSize.z=size.z/nCells.z;
		auto iMap=createIMap(&p[0],nParticles,nCells,cSize);
		
		//for all other particles
		#pragma omp parallel for
		for(int j=0;j<m[i].readNBond();j++)
		{
			int first=m[i].getBonds()[j].x;
			T *C=&m[i].getConstants()[j*nBEADCONST];
			T cutoffSquaredRad=C[0]*C[0];
			auto iCell=getCell(p[first],cSize);
			int cHash=hash(iCell,nCells);
			auto neigh=neighIndices(cHash, nCells);
			T ax=0, ay=0, az=0;
			//#pragma parallel for reduction(+:ax,ay,az)
			//for(int k=0;k<this->nParticles;k++)
			for(auto nHash:neigh)
			{
			auto iMapIt=iMap.find(nHash);
			if(iMapIt!=iMap.end())
			for(auto &k:iMapIt->second)
			{
				if(first!=k)
				{
				threeVector<T> d;
				threeVector<T> aBead=0;
				d.x=p[first].x-p[k].x;
				if(d.x>size.x/2.0) d.x-=size.x;
				if(d.x<-size.x/2.0) d.x+=size.x;
				d.y=p[first].y-p[k].y;
				if(d.y>size.y/2.0) d.y-=size.y;
				if(d.y<-size.y/2.0) d.y+=size.y;
				d.z=p[first].z-p[k].z;
				if(d.z>size.z/2.0) d.z-=size.z;
				if(d.z<-size.z/2.0) d.z+=size.z;
				beadForce(d, aBead, a[k], C, cutoffSquaredRad);
				ax+=aBead.x;
				ay+=aBead.y;
				az+=aBead.z;
				}
			}
			}
			a[first].x+=ax;
			a[first].y+=ay;
			a[first].z+=az;
		}
		}
	}
}

template <typename T>
void Blob<T>::doBoundaryForce(int i)
{
	//this is just going to check it without an error condition, it ignores it
	if(m[i].readType()==BOUNDARY)
	{
		fourVector<int> *bond=m[i].getBonds();
		T *C=m[i].getConstants();
		#pragma omp parallel for
		for(int j=0;j<m[i].readNBond();j++)
		{
			boundaryF(p[bond[j].s[0]],a[bond[j].s[0]],&C[0],size);
		}
	}
}
		
template <typename T>
void Blob<T>::doOffsetBoundaryForce(int i)
{
	if(m[i].readType()==OFFSET_BOUNDARY)
	{
		fourVector<int> *bond=m[i].getBonds();
		T *C=m[i].getConstants();
		#pragma omp parallel for
		for(int j=0;j<m[i].readNBond();j++)
		{
			offsetBoundaryF(p[bond[j].s[0]],a[bond[j].s[0]],&C[0],size);
		}
	}
}

template <typename T>
void Blob<T>::doRigidBendForce(int i)
{
	//this is just going to check it without an error condition, it ignores it
	if(m[i].readType()==RIGIDBEND)
	{
		T *mC=m[i].getConstants();
		fourVector<int> *bond=m[i].getBonds();
		for(int j=0;j<m[i].readNBond();j++)
		{
			int first=bond[j].s[0];
			int second=bond[j].s[1];
			
			threeVector<T> d;
			d.x=p[first].x-p[second].x;
			if(d.x>size.x/2.0) d.x-=size.x;
			if(d.x<-size.x/2.0) d.x+=size.x;
			d.y=p[first].y-p[second].y;
			if(d.y>size.y/2.0) d.y-=size.y;
			if(d.y<-size.y/2.0) d.y+=size.y;
			d.z=p[first].z-p[second].z;
			if(d.z>size.z/2.0) d.z-=size.z;
			if(d.z<-size.z/2.0) d.z+=size.z;
			
			threeVector<T> force=kmaxTorqueF(d, mC);
			
			a[first].x+=force.x;
			a[first].y+=force.y;
			a[first].z+=force.z;
			
			a[second].x-=force.x;
			a[second].y-=force.y;
			a[second].z-=force.z;
		}
	}
}

template <typename T>
void Blob<T>::doFloatingBaseForce(int i)
{
	//this is just going to check it without an error condition, it ignores it
	if(m[i].readType()==FLOATING_BASE)
	{
		T *mC=m[i].getConstants();
		fourVector<int> *bond=m[i].getBonds();
		#pragma omp parallel for
		for(int j=0;j<m[i].readNBond();j++)
		{
			int toPull=bond[j].s[0];
			threeVector<T> d;
			d.z=p[toPull].z;
			//if(d.z>size.z/2.0) d.z-=size.z;
			//if(d.z<-size.z/2.0) d.z+=size.z;
			int cindex=nFLOATING_BASECONST;
			floatingBaseForce(d, a[toPull], mC+cindex*p[toPull].type,size);
		}
	}
}

template <typename T>
T Blob<T>::doFloatingBasePotential(int i)
{
	T potential=0;
	//this is just going to check it without an error condition, it ignores it
	if(m[i].readType()==FLOATING_BASE)
	{
		T *mC=m[i].getConstants();
		fourVector<int> *bond=m[i].getBonds();
		for(int j=0;j<m[i].readNBond();j++)
		{
			int toPull=bond[j].s[0];
			threeVector<T> d;
			d.z=p[toPull].z-mC[0];
			//if(d.z>size.z/2.0) d.z-=size.z;
			//if(d.z<-size.z/2.0) d.z+=size.z;
			int cindex=nFLOATING_BASECONST;
			
			potential+=floatingBasePotential(d, mC+cindex*p[toPull].type,size);
		}
	}
	return potential;
}

template <typename T>
void Blob<T>::doPullBeadForce(int i)
{
	//this is just going to check it without an error condition, it ignores it
	if(m[i].readType()==PULLBEAD)
	{
		T *mC=m[i].getConstants();
		fourVector<int> *bond=m[i].getBonds();
		for(int j=0;j<m[i].readNBond();j++)
		{
			int toPull=bond[j].s[0];
			threeVector<T> d;
			d.x=p[toPull].x-mC[0];
			if(d.x>size.x/2.0) d.x-=size.x;
			if(d.x<-size.x/2.0) d.x+=size.x;
			d.y=p[toPull].y-mC[1];
			if(d.y>size.y/2.0) d.y-=size.y;
			if(d.y<-size.y/2.0) d.y+=size.y;
			d.z=p[toPull].z-mC[2];
			if(d.z>size.z/2.0) d.z-=size.z;
			if(d.z<-size.z/2.0) d.z+=size.z;
			
			threeVector<T> force=0;
			harmonicFZ(d, force, a[toPull], &mC[3]);
			
			force/=static_cast<T>(nParticles-1);
			for(int k=0;k<nParticles;k++)
			{
				//if(p[k].type!=4)
				{
					a[k].x+=force.x;
					a[k].y+=force.y;
					a[k].z+=force.z;
				}
			}
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
		#pragma omp parallel for reduction(+:potential)
		for(int j=0;j<m[i].readNBond();j++)
		{
			threeVector<T> da,db;
			int start=bond[j].s[START];
			int nChains=bond[j].s[NCHAINS];
			int length=bond[j].s[CHAINLENGTH];
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
					//potential+=CounterPotential(cutoffSquared,nTypes,da,&twoBodyUconst[0],p[first].type,p[second].type);
					potential+=bendP(da, db, &(m[i].getConstants()[CHAINBEND]));
				}
				
				int first=k+length-3;
				int second=k+length-2;
				int third=k+length-1;
				
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
				//potential+=CounterPotential(cutoffSquared,nTypes,da,&twoBodyUconst[0],p[first].type,p[second].type);
				//potential+=CounterPotential(cutoffSquared,nTypes,db,&twoBodyUconst[0],p[second].type,p[third].type);
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
T Blob<T>::doZTorquePotential(int i)
{
	T potential=0;
	//this is just going to check it without an error condition, it ignores it
	if(m[i].readType()==ZTORQUE)
	{
		fourVector<int> *bond=m[i].getBonds();
		#pragma omp parallel for reduction(+:potential)
		for(int j=0;j<m[i].readNBond();j++)
		{
			int start=bond[j].s[START];
			int nChains=bond[j].s[NCHAINS];
			int length=bond[j].s[CHAINLENGTH];
			for(int k=start; k<start+length*nChains; k+=length)
			{
				for(int l=k;l<k+length-2;l++)
				{
					int first=l;
					int second=l+1;
					threeVector<T> d;
					
					d.x=p[first].x-p[second].x;
					if(d.x>size.x/2.0) d.x-=size.x;
					if(d.x<-size.x/2.0) d.x+=size.x;
					d.y=p[first].y-p[second].y;
					if(d.y>size.y/2.0) d.y-=size.y;
					if(d.y<-size.y/2.0) d.y+=size.y;
					d.z=p[first].z-p[second].z;
					if(d.z>size.z/2.0) d.z-=size.z;
					if(d.z<-size.z/2.0) d.z+=size.z;
			
					T zHeight=p[first].z+d.z/2.0;
					while(zHeight>=size.z)zHeight-=size.z;
					while(zHeight<0)zHeight+=size.z;
					potential+=zTorquePotential(d, zHeight, m[i].getConstants());
				}
			}
		}
	}
	return potential;
}

template <typename T>
void Blob<T>::doZTorqueForce(int i)
{
	//this is just going to check it without an error condition, it ignores it
	if(m[i].readType()==ZTORQUE)
	{
		fourVector<int> *bond=m[i].getBonds();
		for(int j=0;j<m[i].readNBond();j++)
		{
			int start=bond[j].s[START];
			int nChains=bond[j].s[NCHAINS];
			int length=bond[j].s[CHAINLENGTH];
			for(int k=start; k<start+length*nChains; k+=length)
			{
				for(int l=k;l<k+length-2;l++)
				{
					int first=l;
					int second=l+1;
					threeVector<T> d;
					
					d.x=p[first].x-p[second].x;
					if(d.x>size.x/2.0) d.x-=size.x;
					if(d.x<-size.x/2.0) d.x+=size.x;
					d.y=p[first].y-p[second].y;
					if(d.y>size.y/2.0) d.y-=size.y;
					if(d.y<-size.y/2.0) d.y+=size.y;
					d.z=p[first].z-p[second].z;
					if(d.z>size.z/2.0) d.z-=size.z;
					if(d.z<-size.z/2.0) d.z+=size.z;
			
					T zHeight=p[first].z+d.z/2.0;
					while(zHeight>=size.z)zHeight-=size.z;
					while(zHeight<0)zHeight+=size.z;
					zTorqueForce(d, zHeight,a[first],a[second], m[i].getConstants());
				}
			}
		}
	}
}

template <typename T>
T Blob<T>::doZPowerPotential(int i)
{
	T potential=0;
	//this is just going to check it without an error condition, it ignores it
	if(m[i].readType()==ZPOWERPOTENTIAL)
	{
		fourVector<int> *bond=m[i].getBonds();
		#pragma omp parallel for reduction(+:potential)
		for(int j=0;j<m[i].readNBond();j++)
		{
			int start=bond[j].s[START];
			int nZBond=bond[j].s[NZBOND];
			//int length=bond[j].s[CHAINLENGTH];
			for(int k=start; k<start+nZBond; k++)
			{
				threeVector<T> d;
				d.x=p[k].x;
				d.y=p[k].y;
				d.z=p[k].z;
				potential+=zPowerPotential(d,m[i].getConstants());
			}
		}
	}
	return potential;
}

template <typename T>
void Blob<T>::doZPowerForce(int i)
{
	//this is just going to check it without an error condition, it ignores it
	if(m[i].readType()==ZPOWERPOTENTIAL)
	{
		fourVector<int> *bond=m[i].getBonds();
		for(int j=0;j<m[i].readNBond();j++)
		{
			int start=bond[j].s[START];
			int nZBond=bond[j].s[NZBOND];
			//int length=bond[j].s[CHAINLENGTH];
			#pragma omp parallel for
			for(int k=start; k<start+nZBond; k++)
			{
				threeVector<T> d;
				d.x=p[k].x;
				d.y=p[k].y;
				d.z=p[k].z;
				zPowerForce(d,a[k],m[i].getConstants());
			}
		}
	}
}

template <typename T>
T Blob<T>::doBondPotential(int i)
{
	T potential=0;
	//this is just going to check it without an error condition, it ignores it
	if(m[i].readType()==BOND)
	{
		fourVector<int> *bond=m[i].getBonds();
		
		#pragma omp parallel for reduction(+:potential)
		for(int j=0;j<m[i].readNBond();j++)
		{
			threeVector<T> d;
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
T Blob<T>::doBallPotential(int i)
{
	T potential=0;
	//this is just going to check it without an error condition, it ignores it
	if(m[i].readType()==BALL)
	{
		fourVector<int> *bond=m[i].getBonds();
		
		#pragma omp parallel for reduction(+:potential)
		for(int j=0;j<m[i].readNBond();j++)
		{
			threeVector<T> d;
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
			
			potential+=harmonicHalfP(d, m[i].getConstants());
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
		
		#pragma omp parallel for reduction(+:potential)
		for(int j=0;j<m[i].readNBond();j++)
		{
			threeVector<T> da,db;
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

//It is a non-conservative force
template <typename T>
T Blob<T>::doRigidBendPotential(int i)
{
	T potential=0;
	//this is just going to check it without an error condition, it ignores it
	//if(m[i].readType()==RIGIDBEND)
	//{
	//}
	return potential;
}

template <typename T>
T Blob<T>::doPullBeadPotential(int i)
{
	T potential=0;
	//this is just going to check it without an error condition, it ignores it
	if(m[i].readType()==PULLBEAD)
	{
		T *mC=m[i].getConstants();
		fourVector<int> *bond=m[i].getBonds();
		for(int j=0;j<m[i].readNBond();j++)
		{
			int toPull=bond[j].s[0];
			threeVector<T> d;
			d.x=p[toPull].x-mC[0];
			if(d.x>size.x/2.0) d.x-=size.x;
			if(d.x<-size.x/2.0) d.x+=size.x;
			d.y=p[toPull].y-mC[1];
			if(d.y>size.y/2.0) d.y-=size.y;
			if(d.y<-size.y/2.0) d.y+=size.y;
			d.z=p[toPull].z-mC[2];
			if(d.z>size.z/2.0) d.z-=size.z;
			if(d.z<-size.z/2.0) d.z+=size.z;
			
			potential+=harmonicPZ(d, &mC[3]);
			
		}
	}
	return potential;
}

//Yank one molecule from another by center of mass
/*template <typename T>
T Blob<T>::doYankPotential(int i)
{
	T potential=0;
	//this is just going to check it without an error condition, it ignores it
	if(m[i].readType()==YANK)
	{
		fourVector<int> *bond=m[i].getBonds();
		
		#pragma omp parallel for reduction(+:potential)
		for(int j=0;j<m[i].readNBond();j++)
		{
			threeVector<T> comA=0,comB=0;
			for(int k=0;k<2;k++)
			switch(m[bond.s[k]].readType())
			{
				case BOND:
				{
					fourVector<int> *altBond=m[bond.s[k]].getBonds();
					for(int l=0;l<m[bond.s[k]].readNBond();l++)
					{
						
					}
					break;
				}
				case BEND:
				{
					
					break;
				}
				case CHAIN:
				{
					
					break;
				}
				case BEAD:
				{
					
				}
				//Holy crap this doesn't work right now
				case SOLID:
				{
					
					break;
				}
				case BOUNDARY:
				{
					
					break;
				}
				case RIGIDBEND:
				{
					
					break;
				}
				default:
				{
					//does nothing
					break;
				}
			}
		}
	}
	return potential;
}
*/

//the range of this might change
template <typename T>
T Blob<T>::doBeadPotential(int i)
{
	T potential=0;
	//this is just going to check it without an error condition, it ignores it
	if(m[i].readType()==BEAD)
	{
		//x=index, y=bead-bead type, z=bead-particle type, t=radius
		std::vector<threeVector<int> > beadInfo;
		std::vector<T> beadRadius;
		
		for(int j=i;j<m.size();j++)
		{
			if(m[j].readType()==BEAD)
			{
				for(int k=0;k<m[i].readNBond();k++)
				{
					threeVector<int> info;
					info.x=m[j].getBonds()[k].x;
					info.y=0;//m[i].getConstants()[1];
					info.z=p[info.x].type;
					//this is actually repeated across all particles
					beadRadius.push_back(m[i].getConstants()[BEADRADIUS]);
					beadInfo.push_back(info);
				}
			}
		}
		T *C=m[i].getConstants();	
		
		//for all other beads
		for(int j=0;j<m[i].readNBond();j++)
		{
			#pragma omp parallel for reduction(+:potential)
			for(int k=j+1;k<beadInfo.size();k++)
			{
				threeVector<T> d;
				int first=beadInfo[j].x;
				int second=beadInfo[k].x;
				
				d.x=p[first].x-p[second].x;
				if(d.x>size.x/2.0) d.x-=size.x;
				if(d.x<-size.x/2.0) d.x+=size.x;
				d.y=p[first].y-p[second].y;
				if(d.y>size.y/2.0) d.y-=size.y;
				if(d.y<-size.y/2.0) d.y+=size.y;
				d.z=p[first].z-p[second].z;
				if(d.z>size.z/2.0) d.z-=size.z;
				if(d.z<-size.z/2.0) d.z+=size.z;
				
				T cutoffSquared2=beadRadius[j]+beadRadius[k]+2.0;
				cutoffSquared2*=cutoffSquared2;
				
				int cindex=nBEADCONST*((beadInfo[j].z*nTypes)+beadInfo[k].z);
				
				potential+=beadBeadPotential(d,&(C[cindex]), cutoffSquared2);
			}
		}
		
		T cutoffSquaredRad=m[i].getConstants()[0];//radius+rc
		cutoffSquaredRad*=cutoffSquaredRad;
		
		//for all other particles
		for(int j=0;j<m[i].readNBond();j++)
		{
			#pragma omp parallel for reduction(+:potential)
			for(int k=0;k<this->nParticles;k++)
			{
				threeVector<T> d;
				int first=beadInfo[j].x;
				int second=k;
				if(first!=second)
				{
					d.x=p[first].x-p[second].x;
					if(d.x>size.x/2.0) d.x-=size.x;
					if(d.x<-size.x/2.0) d.x+=size.x;
					d.y=p[first].y-p[second].y;
					if(d.y>size.y/2.0) d.y-=size.y;
					if(d.y<-size.y/2.0) d.y+=size.y;
					d.z=p[first].z-p[second].z;
					if(d.z>size.z/2.0) d.z-=size.z;
					if(d.z<-size.z/2.0) d.z+=size.z;
					
					int cindex=nBEADCONST*((beadInfo[j].z*nTypes)+p[k].type);
					
					potential+=beadPotential(d,&(C[cindex]), cutoffSquaredRad);
				}
			}
		}
	}
	return potential;
}

//the range of this might change
template <typename T>
T Blob<T>::doNanoCorePotential(int i)
{
	T potential=0;
	//this is just going to check it without an error condition, it ignores it
	if(m[i].readType()==NANOCORE)
	{
		for(int j=0;j<m[i].readNBond();j++)
		{
			int pIndex=m[i].getBonds()[j].x;
			//T beadRadius=m[i].getConstants()[BEADRADIUS];
			T *C=&m[i].getConstants()[j*nBEADCONST];
			T beadRadius=C[BEADRADIUS];	
			
			//for all other beads
			// This is the part we are skipping compared to BEAD type
			
			T cutoffRad=C[0];//radius+rc
			T cutoffSquaredRad=cutoffRad*cutoffRad;
			
			//for all other particles
			#pragma omp parallel for reduction(+:potential)
			for(int k=0;k<this->nParticles;k++)
			{
				if(k!=pIndex)
				{
				threeVector<T> d;
				threeVector<T> aBead=0;
				d.x=p[pIndex].x-p[k].x;
				if(d.x>size.x/2.0) d.x-=size.x;
				if(d.x<-size.x/2.0) d.x+=size.x;
				d.y=p[pIndex].y-p[k].y;
				if(d.y>size.y/2.0) d.y-=size.y;
				if(d.y<-size.y/2.0) d.y+=size.y;
				d.z=p[pIndex].z-p[k].z;
				if(d.z>size.z/2.0) d.z-=size.z;
				if(d.z<-size.z/2.0) d.z+=size.z;
				
				potential+=beadPotential(d, C, cutoffSquaredRad);
				}
			}
		}
	}
	return potential;
}

//the range of this might change
template <typename T>
T Blob<T>::doBeadBeadPotential(int i)
{
	T potential=0;
	//this is just going to check it without an error condition, it ignores it
	if(m[i].readType()==BEAD)
	{
		//x=index, y=bead-bead type, z=bead-particle type, t=radius
		std::vector<threeVector<int> > beadInfo;
		std::vector<T> beadRadius;
		
		for(int j=i;j<m.size();j++)
		{
			if(m[j].readType()==BEAD)
			{
				for(int k=0;k<m[i].readNBond();k++)
				{
					threeVector<int> info;
					info.x=m[j].getBonds()[k].x;
					info.y=0;//m[i].getConstants()[1];
					info.z=p[info.x].type;
					//this is actually repeated across all particles
					beadRadius.push_back(m[i].getConstants()[BEADRADIUS]);
					//std::cout << info.t << std::endl;
					beadInfo.push_back(info);
				}
			}
		}
		T *C=m[i].getConstants();	
		
		//for all other beads
		for(int j=0;j<m[i].readNBond();j++)
		{
			#pragma omp parallel for reduction(+:potential)
			for(int k=j+1;k<beadInfo.size();k++)
			{
				threeVector<T> d;
				int first=beadInfo[j].x;
				int second=beadInfo[k].x;
				
				d.x=p[first].x-p[second].x;
				if(d.x>size.x/2.0) d.x-=size.x;
				if(d.x<-size.x/2.0) d.x+=size.x;
				d.y=p[first].y-p[second].y;
				if(d.y>size.y/2.0) d.y-=size.y;
				if(d.y<-size.y/2.0) d.y+=size.y;
				d.z=p[first].z-p[second].z;
				if(d.z>size.z/2.0) d.z-=size.z;
				if(d.z<-size.z/2.0) d.z+=size.z;
				
				T cutoffSquared2=beadRadius[j]+beadRadius[k]+2.0;
				cutoffSquared2*=cutoffSquared2;
				
				int cindex=nBEADCONST*((beadInfo[j].z*nTypes)+beadInfo[k].z);
				
				potential+=beadBeadPotential(d,&(C[cindex]), cutoffSquared2);
			}
		}
	}
	return potential;
}

template <typename T>
T Blob<T>::doBoundaryPotential(int i)
{
	T potential=0;
	//this is just going to check it without an error condition, it ignores it
	if(m[i].readType()==BOUNDARY)
	{
		fourVector<int> *bond=m[i].getBonds();
		T *C=m[i].getConstants();
		#pragma omp parallel for reduction(+:potential)
		for(int j=0;j<m[i].readNBond();j++)
			potential+=boundaryP(p[bond[j].s[0]],&C[0],size);
	}
	return potential;
}
	
template <typename T>
T Blob<T>::doOffsetBoundaryPotential(int i)
{
	T potential=0;
	//this is just going to check it without an error condition, it ignores it
	if(m[i].readType()==OFFSET_BOUNDARY)
	{
		fourVector<int> *bond=m[i].getBonds();
		T *C=m[i].getConstants();
		#pragma omp parallel for reduction(+:potential)
		for(int j=0;j<m[i].readNBond();j++)
			potential+=offsetBoundaryP(p[bond[j].s[0]],&C[0],size);
	}
	return potential;
}

template <typename T>
fourVector<T> Blob<T>::doBeadNeighbors(int i)
{
	fourVector<T> neighbors;
	neighbors.x=0;//number of bead-bead neighbors
	neighbors.y=0;//average distance between bead bead neighbors
	neighbors.z=0;//number of bead-particle neibhros
	neighbors.t=0;//average distance between bead particle neibhors
	
	//this is just going to check it without an error condition, it ignores it
	if(m[i].readType()==BEAD)
	{
		//x=index, y=bead-bead type, z=bead-particle type, t=radius
		std::vector<threeVector<int> > beadInfo;
		std::vector<T> beadRadius;
		
		for(int j=i;j<m.size();j++)
		{
			if(m[j].readType()==BEAD)
			{
				for(int k=0;k<m[i].readNBond();k++)
				{
					threeVector<int> info;
					info.x=m[j].getBonds()[k].x;
					info.y=0;//m[i].getConstants()[1];
					info.z=p[info.x].type;
					//this is actually repeated across all particles
					beadRadius.push_back(m[i].getConstants()[BEADRADIUS]);
					//std::cout << info.t << std::endl;
					beadInfo.push_back(info);
				}
			}
		}
		T *C=m[i].getConstants();	
		
		//for all other beads
		for(int j=0;j<m[i].readNBond();j++)
		{
			for(int k=j+1;k<beadInfo.size();k++)
			{
				threeVector<T> d;
				int first=beadInfo[j].x;
				int second=beadInfo[k].x;
				
				d.x=p[first].x-p[second].x;
				if(d.x>size.x/2.0) d.x-=size.x;
				if(d.x<-size.x/2.0) d.x+=size.x;
				d.y=p[first].y-p[second].y;
				if(d.y>size.y/2.0) d.y-=size.y;
				if(d.y<-size.y/2.0) d.y+=size.y;
				d.z=p[first].z-p[second].z;
				if(d.z>size.z/2.0) d.z-=size.z;
				if(d.z<-size.z/2.0) d.z+=size.z;
				
				T cutoffSquared2=beadRadius[j]+beadRadius[k]+2.0;
				cutoffSquared2*=cutoffSquared2;
				
				T beadBeadDistance=d.x*d.x+d.y*d.y+d.z*d.z;
				if(beadBeadDistance<cutoffSquared2)
				{
					neighbors.x++;
					neighbors.y+=sqrt(beadBeadDistance);
				}
			}
		}
		
		T cutoffSquaredRad=m[i].getConstants()[0];//radius+rc
		cutoffSquaredRad*=cutoffSquaredRad;
		
		//for all other particles
		for(int j=0;j<m[i].readNBond();j++)
		{
			//#pragma omp parallel for reduction(+:potential)
			for(int k=0;k<this->nParticles;k++)
			{
				threeVector<T> d;
				int first=beadInfo[j].x;
				int second=k;
				if(first!=second)
				{
					d.x=p[first].x-p[second].x;
					if(d.x>size.x/2.0) d.x-=size.x;
					if(d.x<-size.x/2.0) d.x+=size.x;
					d.y=p[first].y-p[second].y;
					if(d.y>size.y/2.0) d.y-=size.y;
					if(d.y<-size.y/2.0) d.y+=size.y;
					d.z=p[first].z-p[second].z;
					if(d.z>size.z/2.0) d.z-=size.z;
					if(d.z<-size.z/2.0) d.z+=size.z;
					
					T beadBeadDistance=d.x*d.x+d.y*d.y+d.z*d.z;
					if(beadBeadDistance<cutoffSquaredRad)
					{
						neighbors.z++;
						neighbors.t+=sqrt(beadBeadDistance);
					}
					
					//potential+=beadPotential(d,&(C[cindex]), cutoffSquaredRad);
				}
			}
		}
	}
	if(neighbors.x!=0)
		neighbors.y/=neighbors.x;
	if(neighbors.z!=0)
		neighbors.t/=neighbors.z;
	return neighbors;
	//return potential;
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
		
		#pragma omp parallel for reduction(+:dPotential)
		for(int j=0;j<m[i].readNBond();j++)
		{
			threeVector<T> da,db;
			int start=bond[j].s[START];
			int nChains=bond[j].s[NCHAINS];
			int length=bond[j].s[CHAINLENGTH];
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
					
					T oldPotential=harmonicP(da, &(m[i].getConstants()[CHAINBOND]));
					oldPotential+=bendP(da, db, &(m[i].getConstants()[CHAINBEND]));
						
					da.x*=aSize.x;
					da.y*=aSize.y;
					da.z*=aSize.z;
					
					db.x*=aSize.x;
					db.y*=aSize.y;
					db.z*=aSize.z;
					
					T newPotential=harmonicP(da, &(m[i].getConstants()[CHAINBOND]));
					newPotential+=bendP(da, db, &(m[i].getConstants()[CHAINBEND]));
					dPotential+=(oldPotential-newPotential);
				}
				
				int first=k+length-3;
				int second=k+length-2;
				int third=k+length-1;
				
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
				
				T oldPotential=harmonicP(da, &(m[i].getConstants()[CHAINBOND]));
				oldPotential+=harmonicP(db, &(m[i].getConstants()[CHAINBOND]));
				oldPotential+=bendP(da, db, &(m[i].getConstants()[CHAINBEND]));
					
				da.x*=aSize.x;
				da.y*=aSize.y;
				da.z*=aSize.z;
				
				db.x*=aSize.x;
				db.y*=aSize.y;
				db.z*=aSize.z;
				
				T newPotential=harmonicP(da, &(m[i].getConstants()[CHAINBOND]));
				newPotential+=harmonicP(db, &(m[i].getConstants()[CHAINBOND]));				
				newPotential+=bendP(da, db, &(m[i].getConstants()[CHAINBEND]));
				dPotential+=(oldPotential-newPotential);
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
		#pragma omp parallel for reduction(+:dPotential)
		for(int j=0;j<m[i].readNBond();j++)
		{
			int first=bond[j].s[0];
			int second=bond[j].s[1];
			threeVector<T> d;
			
			d.x=p[first].x-p[second].x;
			if(d.x>size.x/2.0) d.x-=size.x;
			if(d.x<-size.x/2.0) d.x+=size.x;
			d.y=p[first].y-p[second].y;
			if(d.y>size.y/2.0) d.y-=size.y;
			if(d.y<-size.y/2.0) d.y+=size.y;
			d.z=p[first].z-p[second].z;
			if(d.z>size.z/2.0) d.z-=size.z;
			if(d.z<-size.z/2.0) d.z+=size.z;
			
			T oldPotential=harmonicP(d, m[i].getConstants());
			
			d.x*=aSize.x;
			d.y*=aSize.y;
			d.z*=aSize.z;
			
			T newPotential=harmonicP(d,m[i].getConstants());
			
			dPotential+=(oldPotential-newPotential);
		}
	}
	return dPotential;
}

template <typename T>
T Blob<T>::doBallDPotential(int i, threeVector<T> aSize)
{
	T dPotential=0;
	//this is just going to check it without an error condition, it ignores it
	if(m[i].readType()==BALL)
	{
		fourVector<int> *bond=m[i].getBonds();
		#pragma omp parallel for reduction(+:dPotential)
		for(int j=0;j<m[i].readNBond();j++)
		{
			int first=bond[j].s[0];
			int second=bond[j].s[1];
			threeVector<T> d;
			
			d.x=p[first].x-p[second].x;
			if(d.x>size.x/2.0) d.x-=size.x;
			if(d.x<-size.x/2.0) d.x+=size.x;
			d.y=p[first].y-p[second].y;
			if(d.y>size.y/2.0) d.y-=size.y;
			if(d.y<-size.y/2.0) d.y+=size.y;
			d.z=p[first].z-p[second].z;
			if(d.z>size.z/2.0) d.z-=size.z;
			if(d.z<-size.z/2.0) d.z+=size.z;
			
			T oldPotential=harmonicHalfP(d, m[i].getConstants());
			
			d.x*=aSize.x;
			d.y*=aSize.y;
			d.z*=aSize.z;
			
			T newPotential=harmonicHalfP(d,m[i].getConstants());
			
			dPotential+=(oldPotential-newPotential);
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
		#pragma omp parallel for reduction(+:dPotential)
		for(int j=0;j<m[i].readNBond();j++)
		{
			int first=bond[j].s[0];
			int second=bond[j].s[1];
			int third=bond[j].s[2];
			threeVector<T> da,db;
			
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
			
			T oldPotential=bendP(da, db, m[i].getConstants());
			
			da.x*=aSize.x;
			da.y*=aSize.y;
			da.z*=aSize.z;
			
			db.x*=aSize.x;
			db.y*=aSize.y;
			db.z*=aSize.z;
			
			T newPotential=bendP(da, db, m[i].getConstants());
			dPotential+=(oldPotential-newPotential);
		}
	}
	return dPotential;
}

//It is a non-conservative force
template <typename T>
T Blob<T>::doRigidBendDPotential(int i, threeVector<T> aSize)
{
	T dPotential=0;
	//this is just going to check it without an error condition, it ignores it
	//if(m[i].readType()==RIGIDBEND)
	//{
	//	
	//}
	return dPotential;
}

template <typename T>
T Blob<T>::doPullBeadDPotential(int i, threeVector<T> aSize)
{
	T dPotential=0;
	//this is just going to check it without an error condition, it ignores it
	if(m[i].readType()==PULLBEAD)
	{	/*
		T *mC=m[i].getConstants();
		fourVector<int> *bond=m[i].getBonds();
		for(int j=0;j<m[i].readNBond();j++)
		{
			int toPull=bond[j].s[0];
			threeVector<T> d;
			d.x=p[toPull].x-mC[0];
			if(d.x>size.x/2.0) d.x-=size.x;
			if(d.x<-size.x/2.0) d.x+=size.x;
			d.y=p[toPull].y-mC[1];
			if(d.y>size.y/2.0) d.y-=size.y;
			if(d.y<-size.y/2.0) d.y+=size.y;
			d.z=p[toPull].z-mC[2];
			if(d.z>size.z/2.0) d.z-=size.z;
			if(d.z<-size.z/2.0) d.z+=size.z;
			
			T oldPotential=harmonicPZ(d, &mC[3]);
			
			d.x*=aSize.x;
			d.y*=aSize.y;
			d.z*=aSize.z;
			
			T newPotential=harmonicPZ(d,m[i].getConstants());
			
			dPotential+=(oldPotential-newPotential);
			
		}
		*/
	}
	return dPotential;
}


//the range of this might change
template <typename T>
T Blob<T>::doBeadDPotential(int i, threeVector<T> aSize)
{
	T dPotential=0;
	//this is just going to check it without an error condition, it ignores it
	if(m[i].readType()==BEAD)
	{
		//x=index, y=bead-bead type, z=bead-particle type, t=radius
		std::vector<threeVector<int> > beadInfo;
		std::vector<T> beadRadius;
		
		for(int j=i;j<m.size();j++)
		{
			if(m[j].readType()==BEAD)
			{
				for(int k=0;k<m[i].readNBond();k++)
				{
					threeVector<int> info;
					info.x=m[j].getBonds()[k].x;
					info.y=0;//m[i].getConstants()[1];
					info.z=p[info.x].type;
					//this is actually repeated across all particles
					beadRadius.push_back(m[i].getConstants()[BEADRADIUS]);
					beadInfo.push_back(info);
				}
			}
		}
		T *C=m[i].getConstants();	
		//for all other beads
		for(int j=0;j<m[i].readNBond();j++)
		{
			//#pragma omp parallel for reduction(+:dPotential)
			for(int k=j+1;k<beadInfo.size();k++)
			{
				threeVector<T> d;
				int first=beadInfo[j].x;
				int second=beadInfo[k].x;
				
				//T radiusOffset=beadInfo[j].t+beadInfo[k].t;
				
				d.x=p[first].x-p[second].x;
				if(d.x>size.x/2.0) d.x-=size.x;
				if(d.x<-size.x/2.0) d.x+=size.x;
				d.y=p[first].y-p[second].y;
				if(d.y>size.y/2.0) d.y-=size.y;
				if(d.y<-size.y/2.0) d.y+=size.y;
				d.z=p[first].z-p[second].z;
				if(d.z>size.z/2.0) d.z-=size.z;
				if(d.z<-size.z/2.0) d.z+=size.z;
				
				T cutoffSquared2=beadRadius[j]+beadRadius[k]+2.0;
				cutoffSquared2*=cutoffSquared2;
				
				int cindex=nBEADCONST*((beadInfo[j].z*nTypes)+beadInfo[k].z);
				
				T oldPotential=beadBeadPotential(d,&(C[cindex]), cutoffSquared2);
				
				//int cindex=nTWOBODYUCONST*((beadInfo[j].y*nTypes)+beadInfo[k].y);
				//T oldPotential=nonBondedP(d,&twoBodyUconst[cindex],cutoffSquared);
				
				d.x*=aSize.x;
				d.y*=aSize.y;
				d.z*=aSize.z;
				
				T newPotential=beadBeadPotential(d,&(C[cindex]), cutoffSquared2);
				//T newPotential=nonBondedP(d,&twoBodyUconst[cindex],cutoffSquared);
				
				dPotential+=(oldPotential-newPotential);
			}
		}
		
		T cutoffSquaredRad=m[i].getConstants()[0];//radius+rc
		cutoffSquaredRad*=cutoffSquaredRad;
		
		//for all other particles
		for(int j=0;j<m[i].readNBond();j++)
		{
			//#pragma omp parallel for reduction(+:dPotential)
			for(int k=0;k<this->nParticles;k++)
			{
				threeVector<T> d;
				int first=beadInfo[j].x;
				int second=k;
				
				//T radiusOffset=beadInfo[j].t+beadInfo[k].t;
				if(first!=second)
				{
					d.x=p[first].x-p[second].x;
					if(d.x>size.x/2.0) d.x-=size.x;
					if(d.x<-size.x/2.0) d.x+=size.x;
					d.y=p[first].y-p[second].y;
					if(d.y>size.y/2.0) d.y-=size.y;
					if(d.y<-size.y/2.0) d.y+=size.y;
					d.z=p[first].z-p[second].z;
					if(d.z>size.z/2.0) d.z-=size.z;
					if(d.z<-size.z/2.0) d.z+=size.z;
					
					//T radiusOffset=beadInfo[j].t;
					//threeVector<T> offsetVec=unitVector(d)*radiusOffset;
					//d-=offsetVec;
					
					//if(magnitude(d)<2.0)
					//	std::cout << j*2+1 << '\t' << p[second].x << '\t' << p[second].y << '\t' << p[second].z << std::endl;
					
					//int cindex=nTWOBODYUCONST*((beadInfo[j].z*nTypes)+p[second].type);
					int cindex=nBEADCONST*((beadInfo[j].z*nTypes)+p[k].type);
					
					//potential+=beadPotential(d,&(m[i].getConstants()[cindex]));
					
					T oldPotential=beadPotential(d,&(C[cindex]), cutoffSquaredRad);
					
					d.x*=aSize.x;
					d.y*=aSize.y;
					d.z*=aSize.z;
					
					T newPotential=beadPotential(d,&(C[cindex]), cutoffSquaredRad);
					
					dPotential+=(oldPotential-newPotential);
				}
			}
		}
	}
	return dPotential;
}

template <typename T>
T Blob<T>::doNanoCoreDPotential(int i, threeVector<T> aSize)
{
	T dPotential=0;
	//this is just going to check it without an error condition, it ignores it
	if(m[i].readType()==NANOCORE)
	{
		for(int j=0;j<m[i].readNBond();j++)
		{
			int pIndex=m[i].getBonds()[j].x;
			//T beadRadius=m[i].getConstants()[BEADRADIUS];
			T *C=&m[i].getConstants()[j*nBEADCONST];
			T beadRadius=C[BEADRADIUS];
			
			//for all other beads
			// This is the part we are skipping compared to BEAD type
			
			T cutoffRad=C[0];//radius+rc
			T cutoffSquaredRad=cutoffRad*cutoffRad;
			
			//for all other particles
			#pragma omp parallel for reduction(+:dPotential)
			for(int k=0;k<this->nParticles;k++)
			{
				if(k!=pIndex)
				{
				threeVector<T> d;
				threeVector<T> aBead=0;
				d.x=p[pIndex].x-p[k].x;
				if(d.x>size.x/2.0) d.x-=size.x;
				if(d.x<-size.x/2.0) d.x+=size.x;
				d.y=p[pIndex].y-p[k].y;
				if(d.y>size.y/2.0) d.y-=size.y;
				if(d.y<-size.y/2.0) d.y+=size.y;
				d.z=p[pIndex].z-p[k].z;
				if(d.z>size.z/2.0) d.z-=size.z;
				if(d.z<-size.z/2.0) d.z+=size.z;
				
				T oldPotential=beadPotential(d, C, cutoffSquaredRad);
				d.x*=aSize.x;
				d.y*=aSize.y;
				d.z*=aSize.z;
				T newPotential=beadPotential(d, C, cutoffSquaredRad);
				dPotential+=(oldPotential-newPotential);
				}
			}
		}
	}
	return dPotential;
}

template <typename T>
molecule<T, fourVector<int> > * Blob<T>::addMolecule(molecule<T, fourVector<int> > &value)
{
	m.push_back(value);
	nMolecules++;
	commandPresent[NMOLECULES]=true;
	commandPresent[MOLECULE]=true;
	return &m[0];
}

template <typename T>
molecule<T, fourVector<int> > * Blob<T>::setMolecule(molecule<T,fourVector<int> > &value, int index)
{
	#ifdef ERRORS_ENABLED
	if(!commandPresent[MOLECULE])
	{
		std::cerr << "Error (Blob): molecule doesn't appear to be allocated in setMolecule!\n";
		throw 0;
	}
	if(index>nMolecules-1 || index<0)
	{
		std::cerr << "Error (Blob): molecule " << index << " is out of bounds in setMolecule!\n";
		throw 0;
	}
	#endif
	m[index]=value;
	return &m[0];
}

template <typename T>
molecule<T, fourVector<int> > * Blob<T>::delMolecule(int index)
{
	#ifdef ERRORS_ENABLED
	if(!commandPresent[MOLECULE])
	{
		std::cerr << "Error (Blob): molecule doesn't appear to be allocated in delMolecule!\n";
		throw 0;
	}
	if(index>=m.size() || index<0)
	{
		std::cerr << "Error (Blob): molecule " << index << " is out of bounds in delMolecule!\n";
		throw 0;
	}
	#endif
	//due to the behavior of molecule, delete is as simple as shifting all of the molecules back one
	//for(int i=index;i<nMolecules-1;i++)
	//	m[i]=m[i+1];
	
	//copy the end to the middle then pop_back
	if(index!=m.size()-1)
		m[index]=m[m.size()-1];
	m.pop_back();
	nMolecules--;
	if(m.size()==0)
	{
		commandPresent[MOLECULE]=false;
		commandPresent[NMOLECULES]=false;
		commandRequired[MOLECULE]=false;
		commandRequired[NMOLECULES]=false;
		return NULL;
	}
	return &m[0];
}

template <typename T>
molecule<T,fourVector<int> > * Blob<T>::getMolecule()
{
	#ifdef WARNINGS_ENABLED
	if(!commandPresent[MOLECULE])
		std::cerr << "Warning (Blob): molecule doesn't appear to be present in getMolecule!\n";
	#endif
	if(m.size()==0)
		return NULL;
	return &m[0];
}

template <typename T>
void Blob<T>::addParticle(position<T> pos, threeVector<T> vel, threeVector<T> acc)
{
	#ifdef WARNINGS_ENABLED
	if(p.size()!=v.size() || p.size()!=a.size())
		std::cerr << "Warning (Blob): number of acceleration, velocity, and position elements allocated are not equal in addParticle!\n";
	#endif
	
	//Make sure it fits within our system
	if(size.x>0 && size.y>0 && size.z>0)
	{
		while(pos.x>=size.x)
			pos.x-=size.x;
		while(pos.y>=size.y)
			pos.y-=size.y;
		while(pos.z>=size.z)
			pos.z-=size.z;
		
		while(pos.x<0)
			pos.x+=size.x;
		while(pos.y<0)
			pos.y+=size.y;
		while(pos.z<0)
			pos.z+=size.z;
	}
	
	p.push_back(pos);
	v.push_back(vel);
	a.push_back(acc);
	nParticles++;
	commandPresent[NPARTICLES]=true;
	commandPresent[POSITIONS]=true;
	commandPresent[VELOCITIES]=true;
	commandPresent[NPARTICLES]=true;
}

template <typename T>
void Blob<T>::setParticle(int index, position<T> pos, threeVector<T> vel, threeVector<T> acc)
{
	//Make sure it fits within our system
	if(size.x>0 && size.y>0 && size.z>0)
	{
		while(pos.x>=size.x)
			pos.x-=size.x;
		while(pos.y>=size.y)
			pos.y-=size.y;
		while(pos.z>=size.z)
			pos.z-=size.z;
		
		while(pos.x<0)
			pos.x+=size.x;
		while(pos.y<0)
			pos.y+=size.y;
		while(pos.z<0)
			pos.z+=size.z;
	}
	
	#ifdef WARNINGS_ENABLED
	if(!commandPresent[POSITIONS])
		std::cerr << "Warning (Blob): position elements do not appear to be present in setParticle!\n";
	if(!commandPresent[VELOCITIES])
		std::cerr << "Warning (Blob): velocity elements do not appear to be present in setParticle!\n";
	if(!commandPresent[NPARTICLES])
		std::cerr << "Warning (Blob): nParticles do not appear to be present in setParticle!\n";
	#endif
	
	#ifdef ERRORS_ENABLED
	if(p.size()==0)
	{
		std::cerr << "Error (Blob): positions elements do not appear to be allocated in setParticle!\n";
		throw 0;
	}
	if(v.size()==0)
	{
		std::cerr << "Error (Blob): velocity elements do not appear to be allocated in setParticle!\n";
		throw 0;
	}
	if(a.size()==0)
	{
		std::cerr << "Error (Blob): acceleration elements do not appear to be allocated in setParticle!\n";
		throw 0;
	}
	if(index>=nParticles || index<0)
	{
		std::cerr << "Error (Blob): particle " << index << " is out of bounds (0," << nParticles-1 << ") in setParticle!\n";
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
	if(p.size()==0)
	{
		std::cerr << "Error (Blob): position elements do not appear to be allocated in delParticle!\n";
		if(v.size()==0)
			std::cerr << "Error (Blob): velocity elements do not appear to be allocated in delParticle!\n";
		if(a.size()==0)
			std::cerr << "Error (Blob): acceleration elements do not appear to be allocated in delParticle!\n";
		throw 0;
	}
	if(v.size()==0)
	{
		if(a.size()==0)
			std::cerr << "Error (Blob): acceleration elements do not appear to be allocated in delParticle!\n";
		std::cerr << "Error (Blob): velocity elements do not appear to be allocated in delParticle!\n";
		throw 0;
	}
	if(a.size()==0)
	{
		std::cerr << "Error (Blob): acceleration elements do not appear to be allocated in delParticle!\n";
		throw 0;
	}
	if(index>=nParticles || index<0)
	{
		std::cerr << "Error (Blob): particle " << index << " is out of bounds in delParticle!\n";
		throw 0;
	}
	#endif
	
	//Should probably run through all the molecule indices to remove instances of the particular element
	
	p.erase(index);
	v.erase(index);
	a.erase(index);
	nParticles--;
	if(p.size()==0)
		commandPresent[POSITIONS]=false;
	if(v.size()==0)
		commandPresent[VELOCITIES]=false;
	if(nParticles==0)
		commandPresent[NPARTICLES]=false;
}

template <typename T>
void Blob<T>::allocParticle(int n)
{
	if(p.size()<n || v.size()<n || a.size()<n)
	{
		p.reserve(n);
		v.reserve(n);
		a.reserve(n);
	}
}

template <typename T>
position<T> * Blob<T>::getPositions()
{
	#ifdef WARNINGS_ENABLED
	if(!commandPresent[POSITIONS])
		std::cerr << "Warning (Blob): positions do not appear to be present in getPositions!\n";
	#endif
	if(p.size()==0)
		return NULL;
	return &p[0];
}

template <typename T>
threeVector<T> * Blob<T>::getVelocities()
{
	#ifdef WARNINGS_ENABLED
	if(!commandPresent[VELOCITIES])
		std::cerr << "Warning (Blob): velocities do not appear to be present in getVelocities!\n";
	#endif
	if(v.size()==0)
		return NULL;
	return &v[0];
}

template <typename T>
threeVector<T> * Blob<T>::getAccelerations()
{
	#ifdef WARNINGS_ENABLED
	if(!commandPresent[POSITIONS])
		std::cerr << "Warning (Blob): accelerations do not appear to be present in getAccelerations!\n";
	#endif
	if(a.size()==0)
		return NULL;
	return &a[0];
}

template <typename T>
T * Blob<T>::addTwoBodyFconst(T value)
{
	//std::cout << twoBodyFconst.size() << std::endl;
	twoBodyFconst.push_back(value);
	commandPresent[TWOBODYFCONST]=true;
	return &twoBodyFconst[0];
}

template <typename T>
void Blob<T>::setTwoBodyFconst(int index, T value)
{
	#ifdef ERRORS_ENABLED
	if(!commandPresent[TWOBODYFCONST])
	{
		std::cerr << "Error (Blob): twoBodyFconst doesn't appear to be present in setTwoBodyFconst!\n";
		throw 0;
	}
	if(index>=twoBodyFconst.size() || index<0)
	{
		std::cerr << "Error (Blob): index " << index << " is out of bounds (0," << twoBodyFconst.size()-1 << ") in setTwoBodyFconst!\n";
		throw 0;
	}
	#endif
	
	twoBodyFconst[index]=value;
}

template <typename T>
T * Blob<T>::delTwoBodyFconst(int index)
{
	#ifdef ERRORS_ENABLED
	if(!commandPresent[TWOBODYFCONST])
	{
		std::cerr << "Error (Blob): twoBodyFconst doesn't appear to be present in delTwoBodyFconst!\n";
		throw 0;
	}
	if(index>=twoBodyFconst.size() || index<0)
	{
		std::cerr << "Error (Blob): index " << index << " is out of bounds (0," << twoBodyFconst.size()-1 << ") in delTwoBodyFconst!\n";
		throw 0;
	}
	#endif
	
	twoBodyFconst.erase(index);
	if(twoBodyFconst.size()==0)
	{
		commandPresent[TWOBODYFCONST]=false;
		return NULL;
	}
	return &twoBodyFconst[0];
}

template <typename T>
T * Blob<T>::getTwoBodyFconst()
{
	#ifdef WARNINGS_ENABLED
	if(!commandPresent[TWOBODYFCONST])
		std::cerr << "Warning (Blob): twoBodyFconst does not appear to be present in getTwoBodyFconst!\n";
	#endif
	
	if(twoBodyFconst.size()==0)
		return NULL;
	return &twoBodyFconst[0];
}

template <typename T>
T * Blob<T>::addTwoBodyUconst(T value)
{
	twoBodyUconst.push_back(value);
	commandPresent[TWOBODYUCONST]=true;
	return &twoBodyUconst[0];
}

template <typename T>
void Blob<T>::setTwoBodyUconst(int index, T value)
{
	#ifdef ERRORS_ENABLED
	if(!commandPresent[TWOBODYUCONST])
	{
		std::cerr << "Error (Blob): twoBodyUconst doesn't appear to be allocated in setTwoBodyUconst!\n";
		throw 0;
	}
	if(index>=twoBodyUconst.size() || index<0)
	{
		std::cerr << "Error (Blob): index " << index << " is out of bounds in setTwoBodyUconst!\n";
		throw 0;
	}
	#endif
	twoBodyUconst[index]=value;
}

template <typename T>
T * Blob<T>::delTwoBodyUconst(int index)
{
	#ifdef ERRORS_ENABLED
	if(!commandPresent[TWOBODYUCONST])
	{
		std::cerr << "Error (Blob): twoBodyUconst doesn't appear to be allocated in delTwoBodyUconst!\n";
		throw 0;
	}
	if(index>=twoBodyUconst.size() || index<0)
	{
		std::cerr << "Error (Blob): index " << index << " is out of bounds in delTwoBodyUconst!\n";
		throw 0;
	}
	#endif
	twoBodyUconst.erase(index);
	if(twoBodyUconst.size()==0)
	{
		commandPresent[TWOBODYUCONST]=false;
		return NULL;
	}
	return &twoBodyUconst[0];
}

template <typename T>
T * Blob<T>::getTwoBodyUconst()
{
	#ifdef WARNINGS_ENABLED
	if(!commandPresent[TWOBODYUCONST])
		std::cerr << "Warning (Blob): twoBodyUconst does not appear to be allocated in getTwoBodyUconst!\n";
	#endif
	
	if(twoBodyUconst.size()==0)
		return NULL;
	return &twoBodyUconst[0];
}

template <typename T>
int Blob<T>::readNParticles()
{
	if(!commandPresent[NPARTICLES])
		return 0;
	return nParticles;
}

template <typename T>
int Blob<T>::readNMolecules()
{
	if(!commandPresent[NMOLECULES])
		return 0;
	return nMolecules;
}

template <typename T>
T Blob<T>::readGamma()
{
	#ifdef WARNINGS_ENABLED
	if(!commandPresent[GAMMA])
		std::cerr << "Warning (Blob): gamma doesn't appear to be present in readGamma!\n";
	#endif
	if(!commandPresent[GAMMA])
		return 0;
	return gamma;
}

template <typename T>
T Blob<T>::readInitialTemp()
{
	#ifdef WARNINGS_ENABLED
	if(!commandPresent[INITIALTEMP])
		std::cerr << "Warning (Blob): initialTemp doesn't appear to be present in readInitialTemp!\n";
	#endif
	if(!commandPresent[INITIALTEMP])
		return 0;
	return initialTemp;
}

template <typename T>
T Blob<T>::readFinalTemp()
{
	#ifdef WARNINGS_ENABLED
	if(!commandPresent[FINALTEMP])
		std::cerr << "Warning (Blob): finalTemp doesn't appear to be present in readFinalTemp!\n";
	#endif
	if(!commandPresent[FINALTEMP])
		return 0;
	return finalTemp;
}

template <typename T>
int Blob<T>::readSeed()
{
	#ifdef WARNINGS_ENABLED
	if(!commandPresent[SEED])
		std::cerr << "Warning (Blob): seed doesn't appear to be present in readSeed!\n";
	#endif
	if(!commandPresent[SEED])
		return 0;
	return seed;
}

template <typename T>
int Blob<T>::readNTypes()
{
	if(!commandPresent[NTYPES])
		return 0;
	return nTypes;
}

template <typename T>
threeVector<bool> Blob<T>::readPeriodic()
{
	#ifdef WARNINGS_ENABLED
	if(!commandPresent[PERIODIC])
		std::cerr << "Warning (Blob): periodic doesn't appear to be present in readPeriodic!\n";
	#endif
	return periodic;
}

template <typename T>
T Blob<T>::readCutoff()
{
	#ifdef WARNINGS_ENABLED
	if(!commandPresent[CUTOFF])
		std::cerr << "Warning (Blob): cutoff doesn't appear to be present in readCutoff!\n";
	#endif
	if(!commandPresent[CUTOFF])
		return 0;
	return cutoff;
}

template <typename T>
threeVector<T> Blob<T>::readSize()
{
	#ifdef WARNINGS_ENABLED
	if(!commandPresent[SIZE])
		std::cerr << "Warning (Blob): size doesn't appear to be present in readSize!\n";
	#endif
	return size;
}

template <typename T>
T Blob<T>::readInitialTime()
{
	#ifdef WARNINGS_ENABLED
	if(!commandPresent[INITIALTIME])
		std::cerr << "Warning (Blob): initialTime doesn't appear to be present in readInitialTime!\n";
	#endif
	if(!commandPresent[INITIALTIME])
		return 0;
	return initialTime;
}

template <typename T>
T Blob<T>::readFinalTime()
{
	#ifdef WARNINGS_ENABLED
	if(!commandPresent[FINALTIME])
		std::cerr << "Warning (Blob): finalTime doesn't appear to be present in readFinalTime!\n";
	#endif
	if(!commandPresent[FINALTIME])
		return 0;
	return finalTime;
}

template <typename T>
T Blob<T>::readDeltaT()
{
	#ifdef WARNINGS_ENABLED
	if(!commandPresent[DELTAT])
		std::cerr << "Warning (Blob): deltaT doesn't appear to be present in readDeltaT!\n";
	#endif
	if(!commandPresent[DELTAT])
		return 0;
	return deltaT;
}

template <typename T>
T Blob<T>::readStoreInterval()
{
	if(!commandPresent[STOREINTERVAL])
		return 0;
	return storeInterval;
}

template <typename T>
T Blob<T>::readMeasureInterval()
{
	if(!commandPresent[MEASUREINTERVAL])
		return 0;
	return measureInterval;
}

template <typename T>
T Blob<T>::readDeltaLXY()
{
	if(!commandPresent[DELTALXY])
		return 0;
	return deltaLXY;
}

template <typename T>
T Blob<T>::readTempStepInterval()
{
	if(!commandPresent[TEMPSTEPINTERVAL])
		return 0;
	return tempStepInterval;
}

//Keep this comment here:
//COMMANDINPUTHANDLER
template <typename T>
twoVector<T> Blob<T>::readSolventGamma()
{
	return solventGamma;
}

//This member is problematic, it has an answer for everything, use getGammaType if you need to check for presense
template <typename T>
T Blob<T>::readGammaType(int type)
{
	//check if it exists, otherwise fallback on regular gamma value
	if(!commandPresent[GAMMATYPE])
	{
		if(commandPresent[GAMMA])
		{
			return gamma;
		}
		#ifdef ERRORS_ENABLED
			else
			{
				std::cerr << "Error (Blob): type " << type << " doesn't appear to be present in gammaType!\n";
				throw 0;
			}
		#else
			else
			{
				return 1.0;//Horrifying, I know!
			}
		#endif
	}
	#ifdef ERRORS_ENABLED
		if(type<0)
		{
			std::cerr << "Error (Blob): type " << type << " doesn't appear to be present in gammaType!\n";
			throw 0;
		}
	#endif
	return gammaType[type];
}

//get the pointer to gamma type
template <typename T>
T * Blob<T>::getGammaType()
{
	if(!commandPresent[GAMMATYPE])
		return NULL;
	return &gammaType[0];
}

template <typename T>
T Blob<T>::readRemoveSolvent()
{
	if(!commandPresent[REMOVESOLVENT])
		return 0;
	return removeSolvent;
}

template <typename T>
T Blob<T>::readTension()
{
	if(!commandPresent[TENSION])
		return 0;
	return tension;
}

template <typename T>
void Blob<T>::setGamma(T value)
{
	gamma=value;
	commandPresent[GAMMA]=true;
}

template <typename T>
void Blob<T>::setInitialTemp(T value)
{
	initialTemp=value;
	commandPresent[INITIALTEMP]=true;
}

template <typename T>
void Blob<T>::setFinalTemp(T value)
{
	finalTemp=value;
	commandPresent[FINALTEMP]=true;
}

template <typename T>
void Blob<T>::setSeed(int value)
{
	seed=value;
	commandPresent[SEED]=true;
}

template <typename T>
void Blob<T>::setNTypes(int value)
{
	nTypes=value;
	commandPresent[NTYPES]=true;
}

template <typename T>
void Blob<T>::setPeriodic(threeVector<bool> value)
{
	periodic=value;
	commandPresent[PERIODIC]=true;
}

template <typename T>
void Blob<T>::setCutoff(T value)
{
	cutoff=value;
	commandPresent[CUTOFF]=true;
}

template <typename T>
void Blob<T>::setSize(threeVector<T> value)
{
	size=value;
	commandPresent[SIZE]=true;
}

template <typename T>
void Blob<T>::setInitialTime(T value)
{
	initialTime=value;
	commandPresent[INITIALTIME]=true;
}

template <typename T>
void Blob<T>::setFinalTime(T value)
{
	finalTime=value;
	commandPresent[FINALTIME]=true;
}

template <typename T>
void Blob<T>::setDeltaT(T value)
{
	deltaT=value;
		deltaT=value;
	if(value<=0)
		commandPresent[DELTAT]=false;
	else
		commandPresent[DELTAT]=true;
}

template <typename T>
void Blob<T>::setStoreInterval(T value)
{
	storeInterval=value;
	commandPresent[STOREINTERVAL]=true;
}

template <typename T>
void Blob<T>::setMeasureInterval(T value)
{
	measureInterval=value;
	commandPresent[MEASUREINTERVAL]=true;
}

template <typename T>
void Blob<T>::setDeltaLXY(T value)
{
	deltaLXY=value;
	if(deltaLXY==0)
		commandPresent[DELTALXY]=false;
	else
		commandPresent[DELTALXY]=true;
}

template <typename T>
void Blob<T>::setRemoveSolvent(T value)
{
	removeSolvent=value;
	commandPresent[REMOVESOLVENT]=true;
}

template <typename T>
void Blob<T>::setTempStepInterval(T value)
{
	tempStepInterval=value;
	commandPresent[TEMPSTEPINTERVAL]=true;
}

//Keep this comment here:
//COMMANDINPUTHANDLER
template <typename T>
void Blob<T>::setSolventGamma(twoVector<T> value)
{
	solventGamma=value;
	commandPresent[SOLVENTGAMMA]=true;
}

template <typename T>
void Blob<T>::setGammaType(int type, T value)
{
	if(!commandPresent[GAMMATYPE])
	{
		if(commandPresent[NTYPES])
		{
			gammaType=std::vector<double>(nTypes);
		}
		#ifdef ERRORS_ENABLED
		else
		{
			std::cerr << "Error (Blob): nTypes not set!\n";
		}
		#endif
	}
	gammaType[type]=value;
	commandPresent[GAMMATYPE]=true;
}

template <typename T>
void Blob<T>::setTension(T value)
{
	tension=value;
	if(tension==0)
		commandPresent[TENSION]=false;
	else
		commandPresent[TENSION]=true;
}
