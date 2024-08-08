#define ANCHOR_DATA

//For molecular dynamics forces and potentials
#include "../include/MD.h"

//For the molecular dynamics variables
#include "../include/system.h"

//For extraction of molecular dynamics data (measurements)
#include "../dataExtractionLimited.h"

void help(std::string cmd)
{
	std::cerr << "Usage: " << cmd << " name options\n";
	std::cerr << "options:\n";
	std::cerr << "-t [integer]\n\tCenter by type.\n";
	std::cerr << "-i [integer]\n\tCenter by particle index.\n";
	std::cerr << "-T [floating]\n\tAssume start is at some initial time. For size file alignment.\n";
	std::cerr << "-m [integer]\n\tCenter by molecule index.\n";
	std::cerr << "-h \n\tList this help.\n";
	std::cerr << "Note: You can only have one name!\n";
	std::cerr << "Note: name must precede 3 options!\n";
}

int main(int argc, char* argv[])
{
	//encode arguments
	std::stringstream cmdArg;
	for(int i=1;i<argc;i++)
		cmdArg << argv[i] << ' ';
	
	///Configuration variables
	std::vector<std::vector<int> > iStack;
	double initialTime=0.0;
	Blob<double> System;
	std::string name;
	
	//parse arguments
	for(int i=1;i<argc;i++)
	{
		std::string toParse;
		cmdArg >> toParse;
		
		if(toParse.c_str()[0]=='-')
		{
			bool found=false;
			if(toParse=="-h")
			{
				help(argv[0]);
				found=true;
			}
			if(toParse=="-t")
			{
				if(name.size()==0)
				{
					std::cerr << "No system name defined before option!" << std::endl;
					return 2;
				}
				unsigned int tBuf;
				found=true;
				i++;
				cmdArg >> tBuf;
				if(tBuf>=System.readNTypes())
				{
					std::cerr << "Type out of bounds!" << std::endl;
					return 1;
				}
				position<double> *p=System.getPositions();
				int nParticles=System.readNParticles();
				std::vector<int> pIndex;
				for(int j=0;j<nParticles;j++)
					if(p[j].type==tBuf)
						pIndex.push_back(j);
				iStack.push_back(pIndex);
			}
			if(toParse=="-i")
			{
				if(name.size()==0)
				{
					std::cerr << "No system name defined before option!" << std::endl;
					return 2;
				}
				unsigned int pBuf;
				found=true;
				i++;
				cmdArg >> pBuf;
				if(pBuf>=System.readNParticles())
				{
					std::cerr << "Index out of bounds!" << std::endl;
					return 1;
				}
				std::vector<int> pIndex;
				pIndex.push_back(pBuf);
				iStack.push_back(pIndex);
			}
			if(toParse=="-m")
			{
				if(name.size()==0)
				{
					std::cerr << "No system name defined before option!" << std::endl;
					return 2;
				}
				unsigned int mol;
				found=true;
				i++;
				cmdArg >> mol;
				molecule<double,fourVector<int> > *m=System.getMolecule();
				if(mol>=System.readNMolecules())
				{
					std::cerr << "Molecule out of bounds!" << std::endl;
					return 1;
				}
				std::vector<int> pIndex;
				if(m[mol].readType()==CHAIN)
				{
					//go through each group of chains
					for(int bond=0;bond<m[mol].readNBond();bond++)
					{
						int start=m[mol].getBonds()[bond].s[START];
						int length=m[mol].getBonds()[bond].s[CHAINLENGTH];
						int nChains=m[mol].getBonds()[bond].s[NCHAINS];
						for(int j=start;j<start+length*nChains;j++)
							pIndex.push_back(j);
					}
				}
				if(m[mol].readType()==BEAD)
				{
					for(int bond=0;bond<m[mol].readNBond();bond++)
						pIndex.push_back(m[mol].getBonds()[bond].s[0]);
				}
				iStack.push_back(pIndex);
			}
			if(toParse=="-T")
			{
				if(name.size()==0)
				{
					std::cerr << "No system name defined before option!" << std::endl;
					return 2;
				}
				found=true;
				i++;
				cmdArg >> initialTime;
			}
			if(!found)
			{
				std::cerr << "Unrecognized option \"" << argv[i] << "\".\n\n";
				help(argv[0]);
				return -1;
			}
		}
		else
		{
			name=toParse;
			//load old variables, then initialize them
			Script<double, Blob <double> > fileIO(name.c_str(),std::ios::in,&System);
			fileIO.read();
			fileIO.close();
		}
	}
	
	if(iStack.size()!=3)
	{
		std::cerr << "Need 3 options!\n";
		return 0;
	}
	
	std::string framesName("frames_");
	framesName+=name;
	framesName+=".xyz";
	
	position<double> *p=System.getPositions();
	int nParticles=System.readNParticles();
	
	xyzFormat<double> framesFile(p,nParticles);
	framesFile.open(framesName.c_str(), std::ios::in);
	
	System.setInitialTime(initialTime);
	
	std::string sizeName("size_");
	sizeName+=name;
	sizeName+=".dat";
	
	std::fstream sizeFile;
	sizeFile.open(sizeName.c_str(), std::ios::in);
	if(sizeFile.is_open())
	{
		double sTime=-1;
		threeVector<double> sCurrent=System.readSize();
		//for the first iteration, I expect this to exit the loop after one read,
		// otherwise, this is a template for subsequent reads in the loop that is
		// coming up.
		while(sTime<0 && !sizeFile.eof())
			sizeFile >> sTime >> sCurrent.x >> sCurrent.y >> sCurrent.z;
		System.setSize(sCurrent);
		//std::cerr << size.x << '\t' << size.y << '\t' << size.z << '\n';
	}
	
	double simTime=0;
	for(int i=0;framesFile.load();i++)
	{
		std::cerr << static_cast<double>(i)*System.readStoreInterval() << std::endl;
		System.setInitialTime(static_cast<double>(i)*System.readStoreInterval()+initialTime);
		
		simTime=System.readInitialTime();
		
		//make sure the size hasn't changed
		if(sizeFile.is_open())
		{
			double sTime=0;
			threeVector<double> sCurrent=System.readSize();//just in case it is current
			while(sTime<simTime && !sizeFile.eof())
				sizeFile >> sTime >> sCurrent.x >> sCurrent.y >> sCurrent.z;
			System.setSize(sCurrent);
		}
		threeVector<double> size=System.readSize();
		
		std::vector<threeVector<double> > comStack;
		//find center of mass
		for(auto &comParticles:iStack)
		{
			threeVector<double> com=0;
			//position<double> offset=p[comParticles.front()];
			position<double> offset=p[0];
			for(int &j:comParticles)
			{
				threeVector<double> d,pAdj;
				d.x=p[j].x-offset.x;
				d.y=p[j].y-offset.y;
				d.z=p[j].z-offset.z;
				pAdj.x=p[j].x;
				pAdj.y=p[j].y;
				pAdj.z=p[j].z;
				pAdj.x-=d.x>=size.x/2.0?size.x:0;
				pAdj.y-=d.y>=size.y/2.0?size.y:0;
				pAdj.z-=d.z>=size.z/2.0?size.z:0;
				pAdj.x+=d.x<=-size.x/2.0?size.x:0;
				pAdj.y+=d.y<=-size.y/2.0?size.y:0;
				pAdj.z+=d.z<=-size.z/2.0?size.z:0;
				com.x+=pAdj.x;
				com.y+=pAdj.y;
				com.z+=pAdj.z;
			}
			if(comParticles.size()>0)
			{
				com.x/=static_cast<double>(comParticles.size());
				com.y/=static_cast<double>(comParticles.size());
				com.z/=static_cast<double>(comParticles.size());
			}
			comStack.push_back(com);
		}
		threeVector<double> da,db;
		da.x=comStack[0].x-comStack[1].x;
		da.y=comStack[0].y-comStack[1].y;
		da.z=comStack[0].z-comStack[1].z;
		db.x=comStack[1].x-comStack[2].x;
		db.y=comStack[1].y-comStack[2].y;
		db.z=comStack[1].z-comStack[2].z;
		double ra=sqrt(da.x*da.x+da.y*da.y+da.z*da.z);
		double rb=sqrt(db.x*db.x+db.y*db.y+db.z*db.z);
		std::cout << simTime << '\t' << -(da.x*db.x+da.y*db.y+da.z*db.z)/(ra*rb) << std::endl;
		//std::cout << "3\ntest\n";
		//std::cout << 1 << '\t' << comStack[0].x << '\t' << comStack[0].y << '\t' << comStack[0].z << std::endl;
		//std::cout << 2 << '\t' << comStack[1].x << '\t' << comStack[1].y << '\t' << comStack[1].z << std::endl;
		//std::cout << 3 << '\t' << comStack[2].x << '\t' << comStack[2].y << '\t' << comStack[2].z << std::endl;
	}
	return 0;
}
