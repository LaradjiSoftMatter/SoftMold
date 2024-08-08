//this is one set of file formats created by eric spangler
#include "../algorithms/functions.h"

#ifndef PDB_FORMAT
#define PDB_FORMAT

//this is actually allocating memory for you
template <typename T>
class PDBFormat {
public:
	PDBFormat(){allocated=false;};//this one does nothing
	//you kind of have to know what parameters you are using for the other algorithms
	PDBFormat(position<T> * &p, threeVector<T> * &v, threeVector<T> * &a, molecule<T> * &m,
		T &gamma, T &tempI, T &tempF, T * &twoBodyFconst, T * &twoBodyUconst, T * &twoBodyMFconst,
		T * &twoBodyMUconst, T * &threeBodyMFconst, T * &threeBodyMUconst, char *name, int &seed,
		int &types, int &molecules, int &particles, threeVector<bool> &wrap, T &cutoff,
		threeVector<T> &size, T &ti, T &duration, T &deltat, T &storeInterval, T &measureInterval, int &ntwoBodyFconst,
		int &ntwoBodyUconst, int &ntwoBodyMFconst, int &ntwoBodyMUconst, int &nthreeBodyMFconst, int &nthreeBodyMUconst);
	
	~PDBFormat();
	
	//store the system with known parameters
	bool store();
	
	//load a system from files, e.g. create a system with file parameters, allocates required memory for you
	bool load(char *toname);
	
	//store system to files
	bool store(char *toname);
	
	//load a system with these parameters, e.g. create a system with these, doesn't allocate anything
	bool load(position<T> * &p, threeVector<T> * &v, threeVector<T> * &a, molecule<T> * &m,
		T &gamma, T &tempI, T &tempF, T * &twoBodyFconst, T * &twoBodyUconst, T * &twoBodyMFconst,
		T * &twoBodyMUconst, T * &threeBodyMFconst, T * &threeBodyMUconst, char *name, int &seed,
		int &types, int &molecules, int &particles, threeVector<bool> &wrap, T &cutoff,
		threeVector<T> &size, T &ti, T &duration, T &deltat, T &storeInterval, T &measureInterval, int &ntwoBodyFconst,
		int &ntwoBodyUconst, int &ntwoBodyMFconst, int &ntwoBodyMUconst, int &nthreeBodyMFconst, int &nthreeBodyMUconst);
	
private:
	//this are the important parameters
	position<T> **p;
	threeVector<T> **v;
	threeVector<T> **a;
	molecule<T> **m;
	
	//constants loaded at run time
	T *gamma;
	T *tempI;
	T *tempF;
	T **twoBodyFconst;
	T **twoBodyUconst;
	T **twoBodyMFconst;
	T **twoBodyMUconst;
	T **threeBodyMFconst;
	T **threeBodyMUconst;
	
	//basic properties of the particle system
	char *name;    //name of particle system
	int *seed;      //random number seed
	int *types;     //number of particle types
	int *molecules; //number of molecules
	int *particles; //number of particles
	threeVector<bool> *wrap;//periodic boundaries
	T *cutoff;  //cutoff for algorithms
	threeVector<T> *size; //component sizes of system
	T *ti;      //initial time
	T *duration;//length of run
	T *deltat;  //time step size, depends on algorithm in process
	T *storeInterval;   //period for storing data for continuation, float?
	T *measureInterval; //measure and gather data at time steps, float?
	int *ntwoBodyFconst;      //number of two body constant types
	int *ntwoBodyUconst;      //number of two body molecule constant types
	int *ntwoBodyMFconst;     //number of three body molecule constant types
	int *ntwoBodyMUconst;     //number of two body constant types
	int *nthreeBodyMFconst;   //number of two body molecule constant types
	int *nthreeBodyMUconst;   //number of three body molecule constant types
	
	//a flag telling the object if it allocated it's own memory
	bool allocated;
};

template <typename T>
PDBFormat<T>::PDBFormat(position<T> * &p, threeVector<T> * &v, threeVector<T> * &a, molecule<T> * &m,
	T &gamma, T &tempI, T &tempF, T * &twoBodyFconst, T * &twoBodyUconst, T * &twoBodyMFconst,
	T * &twoBodyMUconst, T * &threeBodyMFconst, T * &threeBodyMUconst, char *name, int &seed,
	int &types, int &molecules, int &particles, threeVector<bool> &wrap, T &cutoff,
	threeVector<T> &size, T &ti, T &duration, T &deltat, T &storeInterval, T &measureInterval, int &ntwoBodyFconst,
	int &ntwoBodyUconst, int &ntwoBodyMFconst, int &ntwoBodyMUconst, int &nthreeBodyMFconst, int &nthreeBodyMUconst)
{
	this->p=&p;
	this->v=&v;
	this->a=&a;
	this->m=&m;
	this->gamma=&gamma;
	this->tempI=&tempI;
	this->tempF=&tempF;
	this->twoBodyFconst=&twoBodyFconst;
	this->twoBodyUconst=&twoBodyUconst;
	this->twoBodyMFconst=&twoBodyMFconst;
	this->twoBodyMUconst=&twoBodyMUconst;
	this->threeBodyMFconst=&threeBodyMFconst;
	this->threeBodyMUconst=&threeBodyMUconst;
	this->name=name;
	this->seed=&seed;
	this->types=&types;
	this->molecules=&molecules;
	this->particles=&particles;
	this->wrap=&wrap;
	this->cutoff=&cutoff;
	this->size=&size;
	this->ti=&ti;
	this->duration=&duration;
	this->deltat=&deltat;
	this->storeInterval=&storeInterval;
	this->measureInterval=&measureInterval;
	this->ntwoBodyFconst=&ntwoBodyFconst;
	this->ntwoBodyUconst=&ntwoBodyUconst;
	this->ntwoBodyMFconst=&ntwoBodyMFconst;
	this->ntwoBodyMUconst=&ntwoBodyMUconst;
	this->nthreeBodyMFconst=&nthreeBodyMFconst;
	this->nthreeBodyMUconst=&nthreeBodyMUconst;
	this->allocated=false;
}

template <typename T>
bool PDBFormat<T>::store()
{
	int i,j,periodic=0;
	char buf[256];
	std::fstream parametersfile;
	std::fstream positionsfile;
	std::fstream velocitiesfile;
	std::fstream forcesfile;
	std::fstream moleculesfile;

	//temporary
	if((*wrap).x)
		periodic=1;
	
	//store the parameter file
	sprintf(buf,"parameter_%s.p",this->name);
	parametersfile.open(buf, std::ios::out);
	
	if(!parametersfile.is_open())
	{
		std::cout << buf << " not found!\n";
		return false;
	}
	
	parametersfile << *seed << std::endl;
	parametersfile << *types << std::endl;
	parametersfile << *molecules << std::endl;
	parametersfile << *particles << std::endl;
	parametersfile << (*size).x << std::endl;
	parametersfile << (*size).y << std::endl;
	parametersfile << (*size).z << std::endl;
	parametersfile << periodic << std::endl;
	
	parametersfile << *cutoff << std::endl;
	parametersfile << *ti << std::endl;
	parametersfile << *duration << std::endl;
	parametersfile << *deltat << std::endl;
	parametersfile << *storeInterval << std::endl;
	parametersfile << *measureInterval << std::endl;
	parametersfile << *ntwoBodyFconst << std::endl;
	parametersfile << *ntwoBodyUconst << std::endl;
	parametersfile << *ntwoBodyMFconst << std::endl;
	parametersfile << *ntwoBodyMUconst << std::endl;
	parametersfile << *nthreeBodyMFconst << std::endl;
	parametersfile << *nthreeBodyMUconst << std::endl;
	
	if(parametersfile.eof())
	{
		std::cout << "Parameters file: End of file encountered!\n";
		parametersfile.close();
		return false;
	}
	
	parametersfile.close();
	
	//store the positions
	sprintf(buf,"position_%s.dat",this->name);
	positionsfile.open(buf, std::ios::out);
	if(!positionsfile.is_open())
	{
		std::cout << buf << " not found!\n";
		return false;
	}
	for(i=0;i<*particles;i++)
	{
		positionsfile << (*p)[i].type << '\t' << (*p)[i].x << '\t' << (*p)[i].y << '\t' << (*p)[i].z << std::endl;
	}
	if(positionsfile.eof())
	{
		std::cout << "Positions file: End of file encountered!\n";
		positionsfile.close();
		return false;
	}
	positionsfile.close();
	
	//store frames
	sprintf(buf,"frames_%s.xyz",this->name);
	positionsfile.open(buf, std::ios::out | std::ios::app);
	if(!positionsfile.is_open())
	{
		std::cout << buf << " not found!\n";
		return false;
	}
	positionsfile << *particles << std::endl << this->name << std::endl;
	for(i=0;i<*particles;i++)
		positionsfile << (*p)[i].type << '\t' << (*p)[i].x << '\t' << (*p)[i].y << '\t' << (*p)[i].z << std::endl;
	if(positionsfile.eof())
	{
		std::cout << "Frames file: End of file encountered!\n";
		positionsfile.close();
		return false;
	}
	positionsfile.close();
	
	//store the velocities
	sprintf(buf,"velocity_%s.dat",this->name);
	velocitiesfile.open(buf, std::ios::out);
	if(!velocitiesfile.is_open())
	{
		std::cout << buf << " not found!\n";
		return false;
	}
	for(i=0;i<*particles;i++)
	{
		velocitiesfile << (*v)[i].x << '\t' << (*v)[i].y << '\t' << (*v)[i].z << std::endl;
	}
	if(velocitiesfile.eof())
	{
		std::cout << "Velocity file: End of file encountered!\n";
		velocitiesfile.close();
		return false;
	}
	velocitiesfile.close();

	//store the molecules
	sprintf(buf,"molecules_%s.dat",this->name);
	moleculesfile.open(buf, std::ios::out);
	if(!moleculesfile.is_open())
	{
		std::cout << buf << " not found!\n";
		return false;
	}
	for(i=0;i<*molecules;i++)
	{
		moleculesfile << (*m)[i].type << '\t' << (*m)[i].bonded[0].group[0] << '\t' << (*m)[i].bonded[0].group[1] << '\t' << (*m)[i].bonded[0].group[2] << std::endl;
	}
	if(moleculesfile.eof())
	{
		std::cout << "Molecules file: End of file encountered!\n";
		moleculesfile.close();
		return false;
	}
	moleculesfile.close();

	//store the force constants
	sprintf(buf,"forces_%s.dat",this->name);
	forcesfile.open(buf, std::ios::out);
	if(!forcesfile.is_open())
	{
		std::cout << buf << " not found!\n";
		return false;
	}
	forcesfile << *gamma << std::endl << *tempI << std::endl << *tempF << std::endl;
	for(j=0;j<(*types)*(*types)*(*ntwoBodyFconst);j++)
	{
		forcesfile << (*twoBodyFconst)[j] <<'\n';
	}
	
	for(j=0;j<(*types)*(*types)*(*ntwoBodyUconst);j++)
	{
		forcesfile << (*twoBodyUconst)[j] <<'\n';
	}

	for(j=0;j<(*types)*(*types)*(*ntwoBodyMFconst);j++)
	{
		forcesfile << (*twoBodyMFconst)[j] <<'\n';
	}
	
	for(j=0;j<(*types)*(*types)*(*ntwoBodyMUconst);j++)
	{
		forcesfile << (*twoBodyMUconst)[j] <<'\n';
	}
	
	for(j=0;j<(*types)*(*types)*(*types)*(*nthreeBodyMFconst);j++)
	{
		forcesfile << (*threeBodyMFconst)[j] <<'\n';
	}
	
	for(j=0;j<(*types)*(*types)*(*types)*(*nthreeBodyMUconst);j++)
	{
		forcesfile << (*threeBodyMUconst)[j] <<'\n';
	}
	if(forcesfile.eof())
	{
		std::cout << "Forces file: End of file encountered!\n";
		forcesfile.close();
		return false;
	}
	forcesfile.close();
	
	return true;
}

template <typename T>
bool PDBFormat<T>::load(char *toname)
{
	char buf[256];
	int j,periodic;
	std::fstream parametersfile;
	std::fstream positionsfile;
	std::fstream velocitiesfile;
	std::fstream moleculesfile;
	std::fstream forcesfile;

	//read in parameter file
	this->name=toname;
	
	sprintf(buf,"parameter_%s.p",this->name);
	parametersfile.open(buf, std::fstream::in);
	if(parametersfile.is_open())
	{
		parametersfile >> *seed;
		parametersfile >> *types;
		parametersfile >> *molecules;
		parametersfile >> *particles;
		parametersfile >> (*size).x;
		parametersfile >> (*size).y;
		parametersfile >> (*size).z;
		parametersfile >> periodic;
		//I broke the periodic part in another release, defaulting it now
		(*wrap).x=true;
		(*wrap).y=true;
		(*wrap).z=true;
		parametersfile >> *cutoff;
		parametersfile >> *ti;
		parametersfile >> *duration;
		parametersfile >> *deltat;
		parametersfile >> *storeInterval;
		parametersfile >> *measureInterval;
		parametersfile >> *ntwoBodyFconst;
		parametersfile >> *ntwoBodyUconst;
		parametersfile >> *ntwoBodyMFconst;
		parametersfile >> *ntwoBodyMUconst;
		parametersfile >> *nthreeBodyMFconst;
		parametersfile >> *nthreeBodyMUconst;
		parametersfile.close();
	}
	else
	{
		std::cout << "No parameters file!\n";
		return false;
	}

	//allocate memory for particles and molecules
	*p=new position<T>[*particles];
	*v=new threeVector<T>[*particles];
	*a=new threeVector<T>[*particles];
	*m=new molecule<T>[*molecules];
	
	//load positions
	sprintf(buf,"position_%s.dat",this->name);
	positionsfile.open(buf, std::fstream::in);

	if(positionsfile.is_open())
	{
		for(j=0;j<*particles;j++)
		{
			if(positionsfile.eof())
			{
				std::cout << "Error reading position file: not enough places to go!\n";
				return false;
			}
			positionsfile >> (*p)[j].type >> (*p)[j].x >> (*p)[j].y >> (*p)[j].z;
		}
	}
	else
	{
		std::cout << "No positions file!\n";
		return false;
	}
	positionsfile.close();

	//load in velocities
	sprintf(buf,"velocity_%s.dat",this->name);
	velocitiesfile.open(buf, std::fstream::in);

	if(velocitiesfile.is_open())
	{
		for(j=0;j<*particles;j++)
		{
			if(velocitiesfile.eof())
			{
				std::cout << "Error reading velocity file: not enough velocities!\n";
				return false;
			}
			velocitiesfile >> (*v)[j].x >> (*v)[j].y >> (*v)[j].z;
		}		
	}
	else
	{
		std::cout << "No velocities file!\n";
		return false;
	}
	velocitiesfile.close();

	//initialize accelerations
	for(j=0;j<*particles;j++)
	{
		(*a)[j].x=0;
		(*a)[j].y=0;
		(*a)[j].z=0;
	}
	
	//load molecules
	sprintf(buf,"molecules_%s.dat",this->name);
	moleculesfile.open(buf, std::fstream::in);
	if(moleculesfile.is_open())
	{
		for(j=0;j<*molecules;j++)
		{
			if(moleculesfile.eof())
			{
				std::cout << "Error reading molecule file: not enough molecules!\n";
				return false;
			}
			(*m)[j].bonded=new bonds<T>[1];
			(*m)[j].nBonded=1;
			moleculesfile >> (*m)[j].type >> (*m)[j].bonded[0].group[0] >> (*m)[j].bonded[0].group[1] >> (*m)[j].bonded[0].group[2];
		}
		moleculesfile.close();
	}
	else
	{
		std::cout << "Error, no molecule file.\n";
		return false;
	}
	
	//allocate constant's memories
	*twoBodyFconst=new T[(*types)*(*types)*(*ntwoBodyFconst)];
	*twoBodyUconst=new T[(*types)*(*types)*(*ntwoBodyUconst)];
	*twoBodyMFconst=new T[(*types)*(*types)*(*ntwoBodyMFconst)];
	*twoBodyMUconst=new T[(*types)*(*types)*(*ntwoBodyMUconst)];
	*threeBodyMFconst=new T[(*types)*(*types)*(*types)*(*nthreeBodyMFconst)];
	*threeBodyMUconst=new T[(*types)*(*types)*(*types)*(*nthreeBodyMUconst)];
	
	//load constants
	sprintf(buf,"forces_%s.dat",this->name);
	forcesfile.open(buf, std::fstream::in);
	if(forcesfile.is_open())
	{
		forcesfile >> *gamma >> *tempI >> *tempF;
		
		for(j=0;j<(*types)*(*types)*(*ntwoBodyFconst);j++)
		{
			if(forcesfile.eof())
			{
				std::cout << "Error in force file: not enough forces!";
				return false;
			}
			forcesfile >> (*twoBodyFconst)[j];
		}
	
		for(j=0;j<(*types)*(*types)*(*ntwoBodyUconst);j++)
		{
			if(forcesfile.eof())
			{
				std::cout << "Error in force file: not enough forces!";
				return false;
			}
			forcesfile >> (*twoBodyUconst)[j];
		}

		for(j=0;j<(*types)*(*types)*(*ntwoBodyMFconst);j++)
		{
			if(forcesfile.eof())
			{
				std::cout << "Error in force file: not enough forces!";
				return false;
			}
			forcesfile >> (*twoBodyMFconst)[j];
		}
	
		for(j=0;j<(*types)*(*types)*(*ntwoBodyMUconst);j++)
		{
			if(forcesfile.eof())
			{
				std::cout << "Error in force file: not enough forces!";
				return false;
			}
			forcesfile >> (*twoBodyMUconst)[j];
		}
	
		for(j=0;j<(*types)*(*types)*(*types)*(*nthreeBodyMFconst);j++)
		{
			if(forcesfile.eof())
			{
				std::cout << "Error in force file: not enough forces!";
				return false;
			}
			forcesfile >> (*threeBodyMFconst)[j];
		}
	
		for(j=0;j<(*types)*(*types)*(*types)*(*nthreeBodyMUconst);j++)
		{
			if(forcesfile.eof())
			{
				std::cout << "Error in force file: not enough forces!";
				return false;
			}
			forcesfile >> (*threeBodyMUconst)[j];
		}
		forcesfile.close();
	}
	else
	{
		std::cout << "No forces file!";
		return false;
	}
	allocated=true;
	return true;
}

template <typename T>
bool PDBFormat<T>::load(position<T> * &p, threeVector<T> * &v, threeVector<T> * &a, molecule<T> * &m,
	T &gamma, T &tempI, T &tempF, T * &twoBodyFconst, T * &twoBodyUconst, T * &twoBodyMFconst,
	T * &twoBodyMUconst, T * &threeBodyMFconst, T * &threeBodyMUconst, char *name, int &seed,
	int &types, int &molecules, int &particles, threeVector<bool> &wrap, T &cutoff,
	threeVector<T> &size, T &ti, T &duration, T &deltat, T &storeInterval, T &measureInterval, int &ntwoBodyFconst,
	int &ntwoBodyUconst, int &ntwoBodyMFconst, int &ntwoBodyMUconst, int &nthreeBodyMFconst, int &nthreeBodyMUconst)
{
	p=&p;
	v=&v;
	a=&a;
	m=&m;
	gamma=&gamma;
	tempI=&tempI;
	tempF=&tempF;
	twoBodyFconst=&twoBodyFconst;
	twoBodyUconst=&twoBodyUconst;
	twoBodyMFconst=&twoBodyMFconst;
	twoBodyMUconst=&twoBodyMUconst;
	threeBodyMFconst=&threeBodyMFconst;
	threeBodyMUconst=&threeBodyMUconst;
	name=name;
	seed=&seed;
	types=&types;
	molecules=&molecules;
	particles=&particles;
	wrap=&wrap;
	cutoff=&cutoff;
	size=&size;
	ti=&ti;
	duration=&duration;
	deltat=&deltat;
	storeInterval=&storeInterval;
	measureInterval=&measureInterval;
	ntwoBodyFconst=&ntwoBodyFconst;
	ntwoBodyUconst=&ntwoBodyUconst;
	ntwoBodyMFconst=&ntwoBodyMFconst;
	ntwoBodyMUconst=&ntwoBodyMUconst;
	nthreeBodyMFconst=&nthreeBodyMFconst;
	nthreeBodyMUconst=&nthreeBodyMUconst;
}

template <typename T>
PDBFormat<T>::~PDBFormat()
{
	if(allocated)
	{
		delete *p;
		delete *v;
		delete *a;
		delete *twoBodyFconst;
		delete *twoBodyUconst;
		delete *twoBodyMFconst;
		delete *twoBodyMUconst;
		delete *threeBodyMFconst;
		delete *threeBodyMUconst;
		for(int i=0;i<*molecules;i++)
			delete (*m)[i].bonded;
		delete *m;
	}
	//add stuff here if needed
}

#endif

