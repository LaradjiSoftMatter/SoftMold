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
	if(argc!=2 && argc!=4)
	{
		//so simple, this can't possibly mess it up
		std::cout << "usage: " << argv[0] << " name\n";
		std::cout << "\tOutputs (to terminal) columns:\n\tmolIndex molType molSize\n";
		//std::cout << "usage: " << argv[0] << " name molIndex\n";
		//std::cout << "\tOutputs (to name_2dDistribution.dat) body centered 2 dimensional histogram of spacial probability\n";
		std::cout << "usage: " << argv[0] << " name firstMolIndex secondMolIndex\n";
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
	
	if(argc==4)
	{
		std::string newName("frames_");
		newName+=name;
		newName+=".xyz";
		
		position<double> *p=System.getPositions();
		int nParticles=System.readNParticles();
		
		xyzFormat<double> xyzFile(p,nParticles);
		xyzFile.open(newName.c_str(), std::ios::in);
		
		int molIndexOne,molIndexTwo;
		std::stringstream cmdArg;
		for(int i=2;i<argc;i++)
			cmdArg << argv[i] << '\t';
		cmdArg >> molIndexOne >> molIndexTwo;
		
		if(molIndexOne>=System.readNMolecules() || molIndexOne<0)
		{
			std::cerr << "Molecule index " << molIndexOne << " is out of bounds [0, " << System.readNMolecules() << ").\n";
			return 0;
		}
		
		if(molIndexTwo>=System.readNMolecules() || molIndexTwo<0)
		{
			std::cerr << "Molecule index " << molIndexTwo << " is out of bounds [0, " << System.readNMolecules() << ").\n";
			return 0;
		}
		
		std::vector<int> pIndexOne, pIndexTwo;
		
		for(int bond=0;bond<m[molIndexOne].readNBond();bond++)
		{
			switch(m[molIndexOne].readType())
			{
				case CHAIN:
				{
					int start=m[molIndexOne].getBonds()[bond].s[START];
					int length=m[molIndexOne].getBonds()[bond].s[CHAINLENGTH];
					int nChains=m[molIndexOne].getBonds()[bond].s[NCHAINS];
					
					for(int i=start;i<start+nChains*length;i+=length)
						for(int j=i;j<i+length-1;j++)
							pIndexOne.push_back(j);
				}
				break;
				case BOND:
				{
					pIndexOne.push_back(m[molIndexOne].getBonds()[bond].s[0]);
					pIndexOne.push_back(m[molIndexOne].getBonds()[bond].s[1]);
				}
				break;
				case BEND:
				{
					pIndexOne.push_back(m[molIndexOne].getBonds()[bond].s[0]);
					pIndexOne.push_back(m[molIndexOne].getBonds()[bond].s[1]);
					pIndexOne.push_back(m[molIndexOne].getBonds()[bond].s[2]);
				}
				break;
			}
		}
		std::sort(pIndexOne.begin(), pIndexOne.end());
		std::vector<int>::iterator it;
		it=std::unique(pIndexOne.begin(), pIndexOne.end());
		pIndexOne.resize(std::distance(pIndexOne.begin(),it-1));
		
		for(int bond=0;bond<m[molIndexTwo].readNBond();bond++)
		{
			switch(m[molIndexTwo].readType())
			{
				case CHAIN:
				{
					int start=m[molIndexTwo].getBonds()[bond].s[START];
					int length=m[molIndexTwo].getBonds()[bond].s[CHAINLENGTH];
					int nChains=m[molIndexTwo].getBonds()[bond].s[NCHAINS];
					
					for(int i=start;i<start+nChains*length;i+=length)
						for(int j=i;j<i+length-1;j++)
							pIndexTwo.push_back(j);
				}
				break;
				case BOND:
				{
					pIndexTwo.push_back(m[molIndexTwo].getBonds()[bond].s[0]);
					pIndexTwo.push_back(m[molIndexTwo].getBonds()[bond].s[1]);
				}
				break;
				case BEND:
				{
					pIndexTwo.push_back(m[molIndexTwo].getBonds()[bond].s[0]);
					pIndexTwo.push_back(m[molIndexTwo].getBonds()[bond].s[1]);
					pIndexTwo.push_back(m[molIndexTwo].getBonds()[bond].s[2]);
				}
				break;
			}
		}
		std::sort(pIndexTwo.begin(), pIndexTwo.end());
		it=std::unique(pIndexTwo.begin(), pIndexTwo.end());
		pIndexTwo.resize(std::distance(pIndexTwo.begin(),it-1));
		
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
			
			threeVector<double> comMolOne, comMolTwo;
			for(int i=0;i<pIndexOne.size();i++)
			{
				threeVector<double> translated;
				translated.x=p[pIndexOne[i]].x;
				translated.y=p[pIndexOne[i]].y;
				translated.z=p[pIndexOne[i]].z;
				
				if(i!=0)
				{
					threeVector<double> d;
					d.x=translated.x-p[pIndexOne[0]].x;
					d.y=translated.y-p[pIndexOne[0]].y;
					d.z=translated.z-p[pIndexOne[0]].z;
					translated.x-=(d.x>size.x/2.0)?size.x:0;
					translated.y-=(d.y>size.y/2.0)?size.y:0;
					translated.z-=(d.z>size.z/2.0)?size.z:0;
					translated.x+=(d.x<-size.x/2.0)?size.x:0;
					translated.y+=(d.y<-size.y/2.0)?size.y:0;
					translated.z+=(d.z<-size.z/2.0)?size.z:0;
				}
				comMolOne.x+=translated.x;
				comMolOne.y+=translated.y;
				comMolOne.z+=translated.z;
			}
			comMolOne.x/=static_cast<double>(pIndexOne.size());
			comMolOne.y/=static_cast<double>(pIndexOne.size());
			comMolOne.z/=static_cast<double>(pIndexOne.size());
			comMolOne.x-=(comMolOne.x>size.x)?size.x:0;
			comMolOne.y-=(comMolOne.y>size.y)?size.y:0;
			comMolOne.z-=(comMolOne.z>size.z)?size.z:0;
			comMolOne.x+=(comMolOne.x<0.0)?size.x:0;
			comMolOne.y+=(comMolOne.y<0.0)?size.y:0;
			comMolOne.z+=(comMolOne.z<0.0)?size.z:0;
			
			for(int i=0;i<pIndexTwo.size();i++)
			{
				threeVector<double> translated;
				translated.x=p[pIndexTwo[i]].x;
				translated.y=p[pIndexTwo[i]].y;
				translated.z=p[pIndexTwo[i]].z;
				
				if(i!=0)
				{
					threeVector<double> d;
					d.x=translated.x-p[pIndexTwo[0]].x;
					d.y=translated.y-p[pIndexTwo[0]].y;
					d.z=translated.z-p[pIndexTwo[0]].z;
					translated.x-=(d.x>size.x/2.0)?size.x:0;
					translated.y-=(d.y>size.y/2.0)?size.y:0;
					translated.z-=(d.z>size.z/2.0)?size.z:0;
					translated.x+=(d.x<-size.x/2.0)?size.x:0;
					translated.y+=(d.y<-size.y/2.0)?size.y:0;
					translated.z+=(d.z<-size.z/2.0)?size.z:0;
				}
				comMolTwo.x+=translated.x;
				comMolTwo.y+=translated.y;
				comMolTwo.z+=translated.z;
			}
			comMolTwo.x/=static_cast<double>(pIndexTwo.size());
			comMolTwo.y/=static_cast<double>(pIndexTwo.size());
			comMolTwo.z/=static_cast<double>(pIndexTwo.size());
			comMolTwo.x-=(comMolTwo.x>size.x)?size.x:0;
			comMolTwo.y-=(comMolTwo.y>size.y)?size.y:0;
			comMolTwo.z-=(comMolTwo.z>size.z)?size.z:0;
			comMolTwo.x+=(comMolTwo.x<0.0)?size.x:0;
			comMolTwo.y+=(comMolTwo.y<0.0)?size.y:0;
			comMolTwo.z+=(comMolTwo.z<0.0)?size.z:0;
			
			threeVector<double> dist;
			dist.x=comMolOne.x-comMolTwo.x;
			dist.y=comMolOne.y-comMolTwo.y;
			dist.z=comMolOne.z-comMolTwo.z;
			dist.x-=(dist.x>size.x/2.0)?size.x:0;
			dist.y-=(dist.y>size.y/2.0)?size.y:0;
			dist.z-=(dist.z>size.z/2.0)?size.z:0;
			dist.x+=(dist.x<-size.x/2.0)?size.x:0;
			dist.y+=(dist.y<-size.y/2.0)?size.y:0;
			dist.z+=(dist.z<-size.z/2.0)?size.z:0;
			
			std::cout << time << '\t' << sqrt(dist.x*dist.x+dist.y*dist.y+dist.z*dist.z) << std::endl;
			std::cerr << time << std::endl;
			time+=System.readStoreInterval();
		}
	}
	return 0;
}

