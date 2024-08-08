#define ANCHOR_DATA

//For molecular dynamics forces and potentials
#include "../include/MD.h"

//For the molecular dynamics variables
#include "../include/system.h"

//For extraction of molecular dynamics data (measurements)
#include "../dataExtractionLimited.h"

int main(int argc, char* argv[])
{
	//encode arguments
	std::stringstream cmdArg;
	for(int i=1;i<argc;i++)
		cmdArg << argv[i] << ' ';
	
	std::vector<int> type, pIndex, mIndex;
	double initialTime=0.0;
	std::vector<std::string> names;
	
	//parse arguments
	for(int i=1;i<argc;i++)
	{
		std::string toParse;
		cmdArg >> toParse;
		
		if(toParse.c_str()[0]=='-')
		{
			bool found=false;
			if(toParse=="-h" || toParse=="-help")
			{
				found=true;
				std::cerr << "Usage: " << argv[0] << " name Options\n";
				std::cerr << "Options:\n";
				std::cerr << "-type [integer]\n\tCenter by type.\n";
				std::cerr << "-index [integer]\n\tCenter by particle index.\n";
				std::cerr << "-initialTime [floating]\n\tAssume start is at some initial time. For size file alignment.\n";
				//std::cerr << "-mol [integer]\n\tCenter by molecule index.\n";//Not implemented!
				std::cerr << "Note: You can have many names for a list!\n";
			}
			if(toParse=="-type")
			{
				int tBuf;
				found=true;
				i++;
				cmdArg >> tBuf;
				type.push_back(tBuf);
			}
			if(toParse=="-index")
			{
				int pBuf;
				found=true;
				i++;
				cmdArg >> pBuf;
				pIndex.push_back(pBuf);
			}
			/*
			if(toParse=="-mol")
			{
				found=true;
				i++;
				cmdArg >> mIndex;
			}
			*/
			if(toParse=="-initialTime")
			{
				found=true;
				i++;
				cmdArg >> initialTime;
			}
			
			if(!found)
			{
				std::cerr << "Unrecognized option \"" << argv[i] << "\".\n\n";
				
				std::cerr << "Usage: " << argv[0] << " name Options\n";
				std::cerr << "Options:\n";
				std::cerr << "-type [integer]\n\tCenter by type.\n";
				std::cerr << "-index [integer]\n\tCenter by particle index.\n";
				std::cerr << "-initialTime [floating]\n\tAssume start is at some initial time. For size file alignment.\n";
				//std::cerr << "-mol [integer]\n\tCenter by molecule index.\n";
				std::cerr << "Note: You can have many names for a list!\n";
				
				return -1;
			}
		}
		else
		{
			names.push_back(toParse);
		}
		
		
	}
	
	for(int nameIndex=0;nameIndex<names.size();nameIndex++)
	{
		///Configuration variables
		Blob<double> System;
		
		///load old variables, then initialize them
		Script<double, Blob <double> > fileIO(names[nameIndex].c_str(),std::ios::in,&System);
		
		fileIO.read();
		
		fileIO.close();
		
		position<double> *p=System.getPositions();
		int nParticles=System.readNParticles();
		
		std::string oldName("frames_");
		oldName+=names[nameIndex];
		oldName+=".xyz";
		
		xyzFormat<double> oldXYZFile(p,nParticles);
		oldXYZFile.open(oldName.c_str(), std::ios::in);
		
		std::string newName("framesCentered_");
		newName+=names[nameIndex];
		newName+=".xyz";
		
		System.setInitialTime(initialTime);
		
		std::string sizeName("size_");
		sizeName+=names[nameIndex];
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
		
		//particles to recenter system around
		std::vector<int> comParticles;
		if(type.size()>0)
		{
			for(int i=0;i<nParticles;i++)
				for(int j=0;j<type.size();j++)
					if(type[j]==p[i].type)
						comParticles.push_back(i);
		}
		for(int i=0;i<pIndex.size();i++)
			comParticles.push_back(pIndex[i]);
		
		
		double simTime=0;
		for(int i=0;oldXYZFile.load();i++)
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
			
			//find center of mass
			threeVector<double> com=0;
			for(int i=0;i<comParticles.size();i++)
			{
				com.x+=p[comParticles[i]].x;
				com.y+=p[comParticles[i]].y;
				com.z+=p[comParticles[i]].z;
			}
			if(comParticles.size()>0)
			{
				com.x/=static_cast<double>(comParticles.size());
				com.y/=static_cast<double>(comParticles.size());
				com.z/=static_cast<double>(comParticles.size());
			}
			
			//shift the coordinates to the center of the simulation box
			threeVector<double> s=System.readSize();
			com.x+=s.x/2.0;
			com.y+=s.y/2.0;
			com.z+=s.z/2.0;
			
			for(int i=0;i<nParticles;i++)
			{
				p[i].x-=com.x;
				p[i].y-=com.y;
				p[i].z-=com.z;
				while(p[i].x<0)p[i].x+=s.x;
				while(p[i].y<0)p[i].y+=s.y;
				while(p[i].z<0)p[i].z+=s.z;
			}
			xyzFormat<double> newXYZFile(p,nParticles);
			newXYZFile.open(newName.c_str(), std::ios::out | std::ios::app);
			newXYZFile.store();
			newXYZFile.close();
		}
		oldXYZFile.close();
	}
	return 0;
}
