//This is a very short and simple program that shows you how to modify a variable already present in an mpd script.
//The reason this is used is because some variables (forces, potentials, geometry, etc...) are not easily modified
// in the mpd script. Using a file that opens the mpd script, modifies, and then saves it (as below) is far easier
// than modifying the actual variable. Most of it looks like the run.cpp file, it just replaces the "run loop"
// with the modification.

//For molecular dynamics forces and potentials
#include "../include/MD.h"

//For the molecular dynamics variables
#include "../include/system.h"

int main(int argc, char* argv[])
{
	if(argc!=2 && argc!=5)
	{
		//so simple, this can't possibly mess it up
		std::cout << "usage: " << argv[0] << " name\n";
		std::cout << "\tOutputs (to terminal) columns:\n\tmolIndex molType molSize\n";
		//std::cout << "usage: " << argv[0] << " name molIndex\n";
		//std::cout << "\tOutputs (to name_2dDistribution.dat) body centered 2 dimensional histogram of spacial probability\n";
		std::cout << "usage: " << argv[0] << " name firstMolIndex secondMolIndex thirdMolIndex\n";
		std::cout << "\tOutputs center of mass distances between molecule indices.\n";
		//std::cout << "\tOutputs (to name_2dCOM.dat) body centered (of mass) 2 dimensional histogram of spacial probability\n";
		//std::cout << "\tOutputs (to name_1dCOM.dat) body centered (of mass) 1 dimensional histogram of spacial and time probability\n";
		return 0;
	}
	
	///read in options
	char *name=argv[1];
	//int index=atoi(argv[2]);
	
	///Configuration variables
	Blob<double> System;
	
	///load old variables, then initialize them
	Script<double, Blob <double> > fileIO(name,std::ios::in,&System);
	fileIO.read();
	fileIO.close();
	
	molecule<double, fourVector<int> > *m=System.getMolecule();
	
	if(argc==2)
	{
		for(int i=0;i<System.readNMolecules();i++)
		{
			if(m[i].readType()==CHAIN)
			{
				int length=m[i].getBonds()[0].s[CHAINLENGTH];
				int nChains=m[i].getBonds()[0].s[NCHAINS];
				std::cout << "Chain " << i << " with " << nChains << " chains of length " << length << ".\n";
			}
			if(m[i].readType()==BOND)
			{
				int length=m[i].getBonds()[0].s[CHAINLENGTH];
				int nChains=m[i].getBonds()[0].s[NCHAINS];
				std::cout << "Bond " << i << " with " << m[i].readNBond() << " bonds.\n";
			}
			if(m[i].readType()==BEND)
			{
				int length=m[i].getBonds()[0].s[CHAINLENGTH];
				int nChains=m[i].getBonds()[0].s[NCHAINS];
				std::cout << "Bend " << i << " with " << m[i].readNBond() << " bends.\n";
			}
		}
	}
	
	if(argc==5)
	{
		std::string newName("frames_");
		newName+=name;
		newName+=".xyz";
		
		position<double> *p=System.getPositions();
		int nParticles=System.readNParticles();
		
		xyzFormat<double> xyzFile(p,nParticles);
		xyzFile.open(newName.c_str(), std::ios::in);
		
		std::vector<int> molIndex;
		std::stringstream cmdArg;
		for(int i=2;i<argc;i++)
			cmdArg << argv[i] << '\t';
		for(int i=2;i<argc;i++)
		{
			int mol;
			cmdArg >> mol;
			if(mol>=System.readNMolecules() || mol<0)
			{
				std::cerr << "Molecule index " << mol << " is out of bounds [0, " << System.readNMolecules() << ").\n";
				return 0;
			}
			molIndex.push_back(mol);
		}
		
		std::vector<std::vector<int> > pIndex;
		
		for(int i=0;i<molIndex.size();i++)
		{
			std::vector<int> pIndexI;
			for(int bond=0;bond<m[molIndex[i]].readNBond();bond++)
			{
				switch(m[molIndex[i]].readType())
				{
					case CHAIN:
					{
						int start=m[molIndex[i]].getBonds()[bond].s[START];
						int length=m[molIndex[i]].getBonds()[bond].s[CHAINLENGTH];
						int nChains=m[molIndex[i]].getBonds()[bond].s[NCHAINS];
						
						for(int i=start;i<start+nChains*length;i+=length)
							for(int j=i;j<i+length-1;j++)
								pIndexI.push_back(j);
					}
					break;
					case BOND:
					{
						pIndexI.push_back(m[molIndex[i]].getBonds()[bond].s[0]);
						pIndexI.push_back(m[molIndex[i]].getBonds()[bond].s[1]);
					}
					break;
					case BEND:
					{
						pIndexI.push_back(m[molIndex[i]].getBonds()[bond].s[0]);
						pIndexI.push_back(m[molIndex[i]].getBonds()[bond].s[1]);
						pIndexI.push_back(m[molIndex[i]].getBonds()[bond].s[2]);
					}
					break;
				}
			}
			std::sort(pIndexI.begin(), pIndexI.end());
			std::vector<int>::iterator it;
			it=std::unique(pIndexI.begin(), pIndexI.end());
			pIndexI.resize(std::distance(pIndexI.begin(),it-1));
			pIndex.push_back(pIndexI);
		}
		
		double time=0;
		
		threeVector<double> size=System.readSize();
		
		//check if size varies
		std::string sizeName("size_");
		sizeName+=name;
		sizeName+=".dat";
		
		std::fstream sizeFile;
		sizeFile.open(sizeName.c_str(), std::ios::in);
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
			//std::cout << size.x << '\t' << size.y << '\t' << size.z << '\n';
		}
		
		for(int i=0;xyzFile.load();i++)
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
			
			std::vector< threeVector<double> > comMol;
			for(int j=0;j<pIndex.size();j++)
			{
				threeVector<double> comMolJ;
				for(int i=1;i<pIndex[j].size();i++)
				{
					threeVector<double> translated;
					translated.x=p[pIndex[j][i]].x;
					translated.y=p[pIndex[j][i]].y;
					translated.z=p[pIndex[j][i]].z;
					
					threeVector<double> d;
					d.x=translated.x-p[pIndex[j][0]].x;
					d.y=translated.y-p[pIndex[j][0]].y;
					d.z=translated.z-p[pIndex[j][0]].z;
					translated.x-=(d.x>size.x/2.0)?size.x:0;
					translated.y-=(d.y>size.y/2.0)?size.y:0;
					translated.z-=(d.z>size.z/2.0)?size.z:0;
					translated.x+=(d.x<-size.x/2.0)?size.x:0;
					translated.y+=(d.y<-size.y/2.0)?size.y:0;
					translated.z+=(d.z<-size.z/2.0)?size.z:0;
					
					comMolJ.x+=translated.x;
					comMolJ.y+=translated.y;
					comMolJ.z+=translated.z;
				}
				comMolJ.x/=static_cast<double>(pIndex[j].size());
				comMolJ.y/=static_cast<double>(pIndex[j].size());
				comMolJ.z/=static_cast<double>(pIndex[j].size());
				comMolJ.x-=(comMolJ.x>size.x)?size.x:0;
				comMolJ.y-=(comMolJ.y>size.y)?size.y:0;
				comMolJ.z-=(comMolJ.z>size.z)?size.z:0;
				comMolJ.x+=(comMolJ.x<0.0)?size.x:0;
				comMolJ.y+=(comMolJ.y<0.0)?size.y:0;
				comMolJ.z+=(comMolJ.z<0.0)?size.z:0;
				
				comMol.push_back(comMolJ);
			}
			
			std::vector<threeVector<double> > dist;
			for(int i=0;i<comMol.size()-1;i++)
			{
				dist[i].x=comMol[i].x-comMol[i+1].x;
				dist[i].y=comMol[i].y-comMol[i+1].y;
				dist[i].z=comMol[i].z-comMol[i+1].z;
				
				dist[i].x-=(dist[i].x>size.x/2.0)?size.x:0;
				dist[i].y-=(dist[i].y>size.y/2.0)?size.y:0;
				dist[i].z-=(dist[i].z>size.z/2.0)?size.z:0;
				dist[i].x+=(dist[i].x<-size.x/2.0)?size.x:0;
				dist[i].y+=(dist[i].y<-size.y/2.0)?size.y:0;
				dist[i].z+=(dist[i].z<-size.z/2.0)?size.z:0;
			}
			
			double dra=sqrt(dist[0].x*dist[0].x+dist[0].y*dist[0].y+dist[0].z*dist[0].z);
			double drb=sqrt(dist[1].x*dist[1].x+dist[1].y*dist[1].y+dist[1].z*dist[1].z);
			
			dist[0].x/=dra;
			dist[0].y/=dra;
			dist[0].z/=dra;
			
			dist[1].x/=drb;
			dist[1].y/=drb;
			dist[1].z/=drb;
			
			double cosTheta=(dist[0].x*dist[1].x)+(dist[0].y*dist[1].y)+(dist[0].z*dist[1].z);
			
			std::cout << time << '\t' << cosTheta << std::endl;
			std::cerr << time << std::endl;
			time+=System.readStoreInterval();
		}
	}
	return 0;
}

