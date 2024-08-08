#include "../algorithms/functions.h"
//mpd specific formating


template <typename S>
class mpdFormat {
	public:
		mpdFormat(S &systemBlob, std::string file, std::ios_base::openmode mode);
		mpdFormat(S &systemBlob);
		~mpdFormat();
	
		void open(std::string file, std::ios_base::openmode mode);
		void close();
		void store();
		void safeStore();
		bool load();
		
	private:
		S *sB;
		std::fstream dataFile;
		
		
		
		
};

template <typename T>
mpdFormat<T>::mpdFormat(S &systemBlob, std::string file, std::ios_base::openmode mode)
{

	
	S *sB=&systemBlob;
	
	this->open(file,mode);
}

template <typename T>
mpdFormat<T>::mpdFormat(S &systemBlob)
{
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
	commands[GAMMATYPE]="gammaType";
	
	//map the commands to their indexes to make them easy to find
	for(int i=0;i<nCommands;i++)
	{
		commandMap[commands[i]]=i;
		#ifdef NO_REQUIRED_COMMANDS
			required[i]=false;
		#else
			required[i]=true;
		#endif
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
	required[GAMMATYPE]=false;
	required[GAMMA]=false;
}

template <typename T>
void mpdFormat<T>::open(std::string file, std::ios_base::openmode mode)
{
	dataFile.open(file.c_str(),mode);
	if(!dataFile.is_open())
	{
		std::cout << file << " not found!\n";
		throw 0;
	}
}

template <typename T>
void mpdFormat<T>::close()
{
	if(dataFile.is_open())
	{
		dataFile.close();
	}
}

template <typename T>
mpdFormat<T>::~mpdFormat()
{
	this->close();
	
	delete[] commands;//the mysterious dynamic deallocator of arrays [], actually, this calls each destructor properly
	delete required;
	delete nAllocated;
}

template <typename T>
void mpdFormat<T>::store()
{
	if((*p)==NULL)
	{
		std::cout << "Null pointer in mpdFormat class, returning before segmentation fault!\n";
		throw 0;
	}
	if(dataFile.is_open())
	{
		dataFile << (*nP)-nSolvent << '\n' << "test\n";
		
		for(int i=0;i<(*nP);i++)
		#ifdef SOLVENT_FLAG
			if((*p)[i].type!=SOLVENT_FLAG || !noSolventFlag)
		#endif
		#ifdef SOLVENT
			if((*p)[i].type!=SOLVENT || !noSolventFlag)
		#endif
				dataFile << (*p)[i].type << '\t' << (*p)[i].x << '\t' << (*p)[i].y << '\t' << (*p)[i].z << '\n';
		if(dataFile.eof())
		{
			std::cout << "Something went horribly wrong with the file loaded in mpdFormat, end of file encountered while storing.\n";
			throw 0;
		}
		current++;
	}
	else
	{
		std::cout << "No file open.\n";
		throw 0;
	}
}

template <typename T>
bool mpdFormat<T>::load()
{
	if(dataFile.is_open())
	{
		char name[256];
		dataFile >> (*nP) >> name;
		if((*p)==NULL)
		{
			if(!allocated)
			{
				allocated=true;
				(*p)=new position<T>[(*nP)];
				std::cout << "Allocated position memory with mpdFormat...\n";
			}
			else
			{
				std::cout << "You broke mpdFormat memory allocation, set particle position pointer to NULL before creating class!\n";
				throw 0;
			}
		}
		for(int i=0;i<(*nP);i++)
			dataFile >> (*p)[i].type >> (*p)[i].x >> (*p)[i].y >> (*p)[i].z;
		if(dataFile.eof())
		{
			std::cout << "End of XYZ file!\n";
			return false;
		}
		current++;
	}
	else
	{
		std::cout << "No file open.\n";
		throw 0;
	}
	return true;
}

