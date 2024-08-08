//This is a very short and simple program that shows you how to modify a variable already present in an mpd script.
//The reason this is used is because some variables (forces, potentials, geometry, etc...) are not easily modified
// in the mpd script. Using a file that opens the mpd script, modifies, and then saves it (as below) is far easier
// than modifying the actual variable. Most of it looks like the run.cpp file, it just replaces the "run loop"
// with the modification.

//For molecular dynamics forces and potentials
#include "../include/MD.h"

//For the molecular dynamics variables
#include "../include/system.h"

#include <unordered_map>
#include <forward_list>
#include <stack>

int main(int argc, char* argv[])
{
	if(argc!=5)
	{
		//so simple, this can't possibly mess it up
		std::cout << "usage: " << argv[0] << " name cutoff type tolerance\n";
		return 0;
	}
	
	///read in options
	std::string name;
	double cutoff,tolerance;
	int type;
	std::stringstream cmdArg;
	for(int i=1;i<argc;i++)
		cmdArg << argv[i] << ' ';
	cmdArg >> name >> cutoff >> type >> tolerance;
	
	///Configuration variables
	Blob<double> System;
	
	///load old variables, then initialize them
	Script<double, Blob <double> > fileIO(name.c_str(),std::ios::in,&System);
	fileIO.read();
	fileIO.close();
	
	///Initialize variables and density map
	std::vector<int> indices;
	for(int i=0;i<System.readNParticles();i++)
		if(System.getPositions()[i].type==type)
			indices.push_back(i);
	
	//our current size
	threeVector<double> size=System.readSize();
	
	//check if size varies
	std::string newName("size_");
	newName+=name;
	newName+=".dat";
	
	std::fstream sizeFile;
	sizeFile.open(newName.c_str(), std::ios::in);
	if(sizeFile.is_open())
	{
		double sTime=-1;
		threeVector<double> sCurrent;
		//for the first iteration, I expect this to exit the loop after one read,
		// otherwise, this is a template for subsequent reads in the loop that is
		// coming up.
		while(sTime<0 && !sizeFile.eof())
			sizeFile >> sTime >> sCurrent.x >> sCurrent.y >> sCurrent.z;
		size=sCurrent;
	}
	
	//Get our frames file
	std::string framesName("frames_");
	framesName+=name;
	framesName+=".xyz";
	
	position<double> *p=System.getPositions();
	int nParticles=System.readNParticles();
	
	xyzFormat<double> xyzFile(p,nParticles);
	xyzFile.open(framesName.c_str(), std::ios::in);
	
	double time=0;
	while(xyzFile.load())
	{
		//make sure the size hasn't changed
		if(sizeFile.is_open())
		{
			double sTime=0;
			threeVector<double> sCurrent=size;//just in case it is current
			//for the first iteration, I expect this to exit the loop after one read,
			// otherwise, this is a template for subsequent reads in the loop that is
			// coming up.
			while(sTime<time && !sizeFile.eof())
				sizeFile >> sTime >> sCurrent.x >> sCurrent.y >> sCurrent.z;
			size=sCurrent;
		}
		threeVector<int> nCells;
		nCells.x=size.x/cutoff;
		nCells.y=size.y/cutoff;
		nCells.z=size.z/cutoff;
		threeVector<double> cSize;
		cSize.x=size.x/nCells.x;
		cSize.y=size.y/nCells.y;
		cSize.z=size.z/nCells.z;
		std::vector<position<double>> q;
		for(auto &i:indices)
			q.push_back(p[i]);
		umapfwd<position<double>> cMap=createCMap(&(*q.begin()), &(*q.end()), nCells, cSize);
		
		std::vector<std::vector<position<double>>> surfaces=createSurfaces(cMap, nCells, cSize, size, cutoff);
		std::cerr << time << ' ' << surfaces.size() << ": ";
		for(auto &surface:surfaces)
			std::cerr << surface.size() << ' ';
		std::cerr << std::endl;
		
		time+=System.readStoreInterval();
	}
	return 0;
}

