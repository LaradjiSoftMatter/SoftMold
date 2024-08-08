#include "../algorithms/functions.h"
//vmd specific formating

#ifndef FILE_BUFFER_BYTES
#define FILE_BUFFER_BYTES 65536
#endif

//Need to improve frame seek speed
//allocates memory when you give it the pointers
template <typename T>
class xyzFormat {
	public:
		xyzFormat(position<T> * particles, int nParticles, std::string file, std::ios_base::openmode mode);
		xyzFormat(position<T> * particles, int nParticles);
		~xyzFormat();
	
		void open(std::string file, std::ios_base::openmode mode);
		void close();
		void store();
		bool load();
		bool load(int frame);
		bool loadHeader();
		bool loadLast();
		bool fillFrameBuffer();
		void setNoSolventFlag(bool flag);
		bool readNoSolventFlag();
		int readNParticles();
		position<T> *getPositions();
		
	private:
		int current;
		position<T> *p;
		int nP;
		std::vector< std::vector< position<T> > > frameBuffer;
		int nSolvent;
		bool noSolventFlag;
		bool allocated;
		std::fstream dataFile;
};

template <typename T>
void xyzFormat<T>::setNoSolventFlag(bool flag)
{
	noSolventFlag=flag;
	#ifdef SOLVENT_FLAG
		for(int i=0;i<nP;i++)
			if(p[i].type==SOLVENT_FLAG)
				nSolvent++;
	#endif
	#ifdef SOLVENT
		for(int i=0;i<nP;i++)
			if(p[i].type==SOLVENT)
				nSolvent++;
	#endif
}

template <typename T>
xyzFormat<T>::xyzFormat(position<T> * particles, int nParticles, std::string file, std::ios_base::openmode mode)
{
	p=particles;
	nP=nParticles;
	allocated=false;
	current=0;
	this->open(file,mode);
	nSolvent=0;
	noSolventFlag=false;
}

template <typename T>
xyzFormat<T>::xyzFormat(position<T> * particles, int nParticles)
{
	p=particles;
	nP=nParticles;
	allocated=false;
	current=0;
	nSolvent=0;
	noSolventFlag=false;
}

template <typename T>
void xyzFormat<T>::open(std::string file, std::ios_base::openmode mode)
{
	dataFile.open(file.c_str(),mode);
	if(!dataFile.is_open())
	{
		std::cerr << file << " not found!\n";
		throw 0;
	}
}

template <typename T>
void xyzFormat<T>::close()
{
	if(dataFile.is_open())
	{
		dataFile.close();
	}
}

template <typename T>
xyzFormat<T>::~xyzFormat()
{
	if(allocated)
	{
		std::cerr << "Deleting position memory!\n";
		delete p;
	}
	this->close();
}

template <typename T>
void xyzFormat<T>::store()
{
	if(p==NULL)
	{
		std::cerr << "Null pointer in xyzFormat class, returning before segmentation fault!\n";
		throw 0;
	}
	if(dataFile.is_open())
	{
		dataFile << nP-nSolvent << '\n' << "test\n";
		
		for(int i=0;i<nP;i++)
		#ifdef SOLVENT_FLAG
			if(p[i].type!=SOLVENT_FLAG || !noSolventFlag)
		#endif
		#ifdef SOLVENT
			if(p[i].type!=SOLVENT || !noSolventFlag)
		#endif
				dataFile << p[i].type << '\t' << p[i].x << '\t' << p[i].y << '\t' << p[i].z << '\n';
		if(dataFile.eof())
		{
			std::cerr << "Something went horribly wrong with the file loaded in xyzFormat, end of file encountered while storing.\n";
			throw 0;
		}
		current++;
	}
	else
	{
		std::cerr << "No file open.\n";
		throw 0;
	}
}

template <typename T>
bool xyzFormat<T>::load()
{
	if(dataFile.is_open())
	{
		char name[256];
		dataFile >> nP >> name;
		if(p==NULL)
		{
			if(!allocated)
			{
				allocated=true;
				p=new position<T>[nP];
				std::cerr << "Allocated position memory with xyzFormat...\n";
			}
			else
			{
				std::cerr << "You broke xyzFormat memory allocation, set particle position pointer to NULL before creating class!\n";
				throw 0;
			}
		}
		for(int i=0;i<nP;i++)
			dataFile >> p[i].type >> p[i].x >> p[i].y >> p[i].z;
		if(dataFile.eof())
		{
			std::cerr << "End of XYZ file!\n";
			return false;
		}
		current++;
	}
	else
	{
		std::cerr << "No file open.\n";
		throw 0;
	}
	return true;
}
/*
template <typename T>
bool xyzFormat<T>::fillFrameBuffer()
{
	if(dataFile.is_open())
	{
		if(frameBuffer.size()>0)
		{
			for(int i=0;i<frameBuffer.size();i++)
				frameBuffer[i].clear();
			frameBuffer.clear();
		}
		current=0;
		dataFile.seekg(std::ios::beg);
		
		std::vector< position<T> > singleFrame;
		char bufferA[FILE_BUFFER_BYTES], bufferB[FILE_BUFFER_BYTES];
		
		while(!dataFile.eof())
		{
			std::stringstream input;
			dataFile.read(bufferA, FILE_BUFFER_BYTES);
		}
		
		dataFile >> nP >> name;
		
		for(int i=0;i<nP;i++)
			dataFile >> p[i].type >> p[i].x >> p[i].y >> p[i].z;
		if(dataFile.eof())
		{
			std::cerr << "End of XYZ file!\n";
			return false;
		}
		current++;
	}
	else
	{
		std::cerr << "No file open.\n";
		throw 0;
	}
	return true;
}*/

template <typename T>
bool xyzFormat<T>::load(int frame)
{
	if(dataFile.is_open())
	{
		if(frame<current)
		{
			current=0;
			dataFile.seekg(std::ios::beg);
		}
		
		int tempT,tempNP;
		T tempX, tempY, tempZ;
		char name[256];
		
		for(;current<frame;current++)
		{
			dataFile >> tempNP >> name;
			if(!dataFile.eof())
			{
				for(int i=0;i<tempNP;i++)
					dataFile >> tempT >> tempX >> tempY >> tempZ;
			}
		}
		
		if(dataFile.eof())//this is repeated so that time isn't wasted loading nothing
		{
			std::cerr << "Frame does not exist!\n";
			return false;
		}
		
		dataFile >> nP >> name;
		
		if(p==NULL)
		{
			if(!allocated)
			{
				allocated=true;
				p=new position<T>[nP];
			}
			else
			{
				std::cerr << "You broke xyzFormat memory allocation, set particle position pointer to NULL!\n";
				throw 0;
			}
		}
		
		for(int i=0;i<nP;i++)
			dataFile >> p[i].type >> p[i].x >> p[i].y >> p[i].z;
		
		if(dataFile.eof())
		{
			std::cerr << "Frame does not exist!\n";
			return false;
		}
	}
	else
	{
		std::cerr << "No file open.\n";
		throw 0;
	}
	return true;
}

template <typename T>
bool xyzFormat<T>::loadLast()
{
	if(dataFile.is_open())
	{
		long filePtr=0;
		int tempT,tempNP;
		T tempX, tempY, tempZ;
		char name[256];
		
		//dataFile.seekg(std::ios::end);
		
		//std::cerr << filePtr << '\t' << dataFile.tellg() << '\n';
		
		if(dataFile.eof())
		{
			dataFile.seekg(std::ios::beg);
			dataFile.clear();
		}
		
		
		
		for(filePtr=dataFile.tellg();!dataFile.eof();current++)
		{
			//std::cerr << filePtr << '\n';
			dataFile >> tempNP >> name;
			if(!dataFile.eof())
			{
				for(int i=0;i<tempNP;i++)
					dataFile >> tempT >> tempX >> tempY >> tempZ;
			}
		}
		
		dataFile.clear();
		dataFile.seekg(filePtr);
		
		dataFile >> nP >> name;
		
		if(p==NULL)
		{
			if(!allocated)
			{
				allocated=true;
				p=new position<T>[nP];
			}
			else
			{
				std::cerr << "You broke xyzFormat memory allocation, set particle position pointer to NULL!\n";
				throw 0;
			}
		}
		
		if(!dataFile.eof())
		{
			for(int i=0;i<nP;i++)
				dataFile >> p[i].type >> p[i].x >> p[i].y >> p[i].z;
		}
		else
		{
			std::cerr << "Frame does not exist!\n";
			return false;
		}
	}
	else
	{
		std::cerr << "No file open.\n";
		throw 0;
	}
	return true;
}

template <typename T>
int xyzFormat<T>::readNParticles()
{
	return nP;
}

template <typename T>
position<T> *xyzFormat<T>::getPositions()
{
	return p;
}
