//Loading an MD simulation "using systemMD.h".
//Everything needed to do a simulation is already included in systemMD.h.
//Even standard libraries like iostream and fstream are inlcuded.
//But, if any extra libraries are needed you can include them here or in systemMD.h.
//Permenent libraries (included with standard library or with this package) should
// go in systemMD.h, and anything you need that is external should go in this file.

//For molecular dynamics forces and potentials
#include "../include/MD.h"

//For the molecular dynamics variables
#include "../include/system.h"

//For radiusOfGyrationAxisIndex
#include "radiusOfGyrationInclude/radiusOfGyration.h"

int main(int argc, char **argv)
{
	if(argc!=4)
	{
		std::cerr << argv[0] << " name fromTime molIndex \n";
		return 0;
	}
	
	std::stringstream cmdArg;
	for(int i=2;i<argc;i++)
		cmdArg << argv[i] << '\t';
	
	double fromTime=0;
	threeVector<int> n;
	cmdArg >> fromTime;
	
	int molIndex;
	cmdArg >> molIndex;
	/*
	std::vector<int> types;
	
	for(int i=4;i<argc;i++)
	{
		int type;
		cmdArg >> type;
		types.push_back(type);
	}*/
	
	//the variables for the simulation, call with empty constructor
	Blob<double> System;
	
	//file input-output script engine
	//name part (argv[1] in this case) is used with no extension (*.mpd)
	Script<double, Blob <double> > fileIO(argv[1],std::ios::in,&System);
	
	//read the variables in
	fileIO.read();
	
	//close the file that it is reading from
	fileIO.close();
	
	//get the chains
	molecule<double,fourVector<int>> *m=System.getMolecule();
	int nMol=System.readNMolecules();
	if(molIndex<0 || molIndex>=nMol)
	{
		std::cerr << "molIndex " << molIndex << " is out of range [0," << nMol <<")!" << std::endl;
		return -1;
	}
	if(m[molIndex].readType()!=CHAIN)
	{
		std::cerr << "Molecule " << molIndex << " is not CHAIN type! (" 
			<< m[molIndex].readType() << "!=" << CHAIN << ")" << std::endl;
		return -1;
	}
	std::vector<std::vector<int>> chains;
	for(int bond=0;bond<m[molIndex].readNBond();bond++)
	{
		int start=m[molIndex].getBonds()[bond].s[START];
		int nChains=m[molIndex].getBonds()[bond].s[NCHAINS];
		int length=m[molIndex].getBonds()[bond].s[CHAINLENGTH];
		for(int i=start;i<nChains*length;i+=length)
		{
			std::vector<int> chain;
			for(int j=i;j<i+length;j++)
				chain.push_back(j);
			chains.push_back(chain);
		}
		std::cerr << "Found chain " << molIndex << " starting at " << start << " with " << nChains << " chains of length " << length << ".\n";
	}
	
	//set time
	double time=0;

	//grab sizes
	threeVector<double> size=System.readSize();
	
	std::string newName("size_");
	newName+=argv[1];
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
		time=sTime;
		size=sCurrent;
		//std::cout << size.x << '\t' << size.y << '\t' << size.z << '\n';
	}
		
	//Get our frames file
	std::string framesName("frames_");
	framesName+=argv[1];
	framesName+=".xyz";
	
	position<double> *p=System.getPositions();
	int nParticles=System.readNParticles();
	
	xyzFormat<double> xyzFile(p,nParticles);
	xyzFile.open(framesName.c_str(), std::ios::in);
	
	int frames=0;
	
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
			{
				sizeFile >> sTime >> sCurrent.x >> sCurrent.y >> sCurrent.z;
				//std::cout << time << '\t' << sTime << '\t' << size.x << '\t' << size.y << '\t' << size.z << std::endl;
				
			}
			size=sCurrent;
		}
		if(time>=fromTime)
		{
			std::vector<double> rogs;
			for(auto chain:chains)
				rogs.push_back(radiusOfGyrationIndex(chain.begin(),chain.end(),p,size));
			
			double rogAvg=0;
			for(auto rog:rogs)
				rogAvg+=rog;
			if(rogs.size()>0)
				rogAvg/=rogs.size();
			
			frames++;
			std::cerr << time << ' ' << rogAvg << std::endl;
			std::cout << time << ' ' << rogAvg << std::endl;
		}
		time+=System.readStoreInterval();
	}
	
	return 0;
}

