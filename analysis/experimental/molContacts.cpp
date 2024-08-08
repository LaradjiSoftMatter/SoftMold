//This is a very short and simple program that shows you how to modify a variable already present in an mpd script.
//The reason this is used is because some variables (forces, potentials, geometry, etc...) are not easily modified
// in the mpd script. Using a file that opens the mpd script, modifies, and then saves it (as below) is far easier
// than modifying the actual variable. Most of it looks like the run.cpp file, it just replaces the "run loop"
// with the modification.

//For molecular dynamics forces and potentials
#include "../include/MD.h"

//For the molecular dynamics variables
#include "../include/system.h"

//Tossed here for your inconvienence
class computePotential {
	public:
		//constructor
		computePotential(position<double> *particles, int nParticles, double *constants, double cutoff, int nTypes)
		{
			p=particles;
			nP=nParticles;
			c=constants;
			cutoffSquared=cutoff*cutoff;
			nT=nTypes;
			count=new int[nP];
			potentialI=new double[nP];
		};
		
		//destructor
		~computePotential()
		{
			if(nP>0)
			{
				delete count;
				delete potentialI;
			}
		};
		
		//set potential to 0
		void initialize()
		{
			potential=0;
			for(int i=0;i<nP;i++)
			{
				potentialI[i]=0;
				count[i]=0;
			}
		};
		
		//output the total potential
		double outputPotential()
		{
			return potential;
		};
		
		//get the calculated potential pointer
		double *getPotentials()
		{
			return potentialI;
		};
		
		//count
		int outputCount()
		{
			int total=0;
			for(int i=0;i<nP;i++)
				total+=count[i];
			return total;
		};
		
		//plain potential calculation
		void operator() (int &i, int &j)
		{
			double tPotential=Potential<double>(cutoffSquared, nT, p[i], p[j], c);
			if(tPotential!=0)
			{
				count[i]=1;
				count[j]=1;
			}
			potential+=tPotential;
		};
		
		//minimum image version
		void operator() (int &i, int &j, threeVector<double> &minImg)
		{
			position<double> pj=p[j];
			pj.x+=minImg.x;
			pj.y+=minImg.y;
			pj.z+=minImg.z;
			double tPotential=Potential<double>(cutoffSquared, nT, p[i], pj, c);
			if(tPotential!=0)
			{
				count[i]=1;
				count[j]=1;
			}
			potential+=tPotential;
		};
		
	private:
		double potential;
		position<double> *p;
		int nP;
		double *c, *potentialI;
		double cutoffSquared;
		int nT;
		int *count;
};

int main(int argc, char* argv[])
{
	if(argc<2)
	{
		//so simple, this can't possibly mess it up
		std::cout << "To see available chain type molecules:\n";
		std::cout << "\tusage: " << argv[0] << " name\n";
		
		std::cout << "To get the average length of chain type molecules:\n";
		std::cout << "\tusage: " << argv[0] << " name molType\n";
		
		std::cout << "To get the length of a particular chain:\n";
		std::cout << "\tusage: " << argv[0] << " name molIndex0 molIndex1 molIndex2 ...\n";
		return 0;
	}
	
	///read in options
	char *name=argv[1];
	
	std::stringstream cmdArg;
	for(int i=2;i<argc;i++)
		cmdArg << argv[i] << '\t';
	
	///Configuration variables
	Blob<double> System;
	
	///load old variables, then initialize them
	Script<double, Blob <double> > fileIO(name,std::ios::in,&System);
	fileIO.read();
	fileIO.close();
	
	position<double> *p=System.getPositions();
	int nParticles=System.readNParticles();
	molecules< fourVector<int> > *m=System.getMolecules();
	
	//grab sizes
	threeVector<double> size=System.readSize();
	
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
	}
	
	if(argc==2)//dump mol information
	{
		for(int i=0;i<System.readNMolecules();i++)
		{
			if(m[i].readType()==CHAIN)
			{
				int length=m[i].getBonds()[0].s[LENGTH];
				int nChains=m[i].getBonds()[0].s[NCHAINS];
				std::cout << "Chain (" << CHAIN << ") " << i << " with " << nChains << " chains of length " << length << ".\n";
			}
			if(m[i].readType()==BOND)
				std::cout << "Bond (" << BOND << ") " << i << " with " << m[i].readNBond() << " bonds.\n";
			if(m[i].readType()==BEND)
				std::cout << "Bend (" << BEND << ") " << i << " with " << m[i].readNBond() << " bends.\n";
		}
	}
	if(argc>=3)//check only molType
	{
		std::vector<std::vector<int> > molPI;
		if(argc==3)
		{
			int molType;
			cmdArg >> molType;
			if(molType!=BOND || molType!=BEND || molType!=CHAIN)
			{
				std::cerr << "molType " << molType << " is not one of " << BOND << ", " << BEND << ", or " << CHAIN << "!\n";
				return -1;
			}
			
			//unique set of molecule indices per molecule
			for(int i=0;i<System.readNMolecules();i++)
			{
				for(int bond=0;bond<m[i].readNBond() && m[i].readType()==molType;bond++)
				{
					std::vector<int> mol;
					switch(m[i].readType())
					{
						case CHAIN:
						{
							int start=m[i].getBonds()[bond].s[START];
							int length=m[i].getBonds()[bond].s[LENGTH];
							int nChains=m[i].getBonds()[bond].s[NCHAINS];
							
							for(int i=start;i<start+nChains*length;i+=length)
								for(int j=i;j<i+length-1;j++)
									mol.push_back(j);
						}
						break;
						case BOND:
						{
							mol.push_back(m[i].getBonds()[bond].s[0]);
							mol.push_back(m[i].getBonds()[bond].s[1]);
						}
						break;
						case BEND:
						{
							mol.push_back(m[i].getBonds()[bond].s[0]);
							mol.push_back(m[i].getBonds()[bond].s[1]);
							mol.push_back(m[i].getBonds()[bond].s[2]);
						}
						break;
					}
					if(mol.size()>0)
					{
						std::sort(mol.begin(), mol.end());
						std::vector<int>::iterator it;
						it=std::unique(mol.begin(), mol.end());
						mol.resize(std::distance(mol.begin(),it-1));
						molPI.push_back(mol);
					}
				}
			}
		}
		
		if(argc>3)//check contacts between these mols
		{
			std::vector<int> molIndices;
			do {
				int molIndex;
				cmdArg >> molIndex;
				molIndices.push_back(molIndex);
				if(molIndex>=System.readNMolecules() || molIndex<0)
				{
					std::cerr << "molIndex " << molIndex << " is out of range [0," << System.readNMolecules() << ")\n";
					return -1;
				}
			} while(!cmdArg.eof());
			
			//unique set of molecule indices per molecule
			for(int molI=0;molI<molIndices.size();molI++)
			{
				int i=molIndices[molI];
				for(int bond=0;bond<m[i].readNBond() && m[i].readType()==molType;bond++)
				{
					std::vector<int> mol;
					switch(m[i].readType())
					{
						case CHAIN:
						{
							int start=m[i].getBonds()[bond].s[START];
							int length=m[i].getBonds()[bond].s[LENGTH];
							int nChains=m[i].getBonds()[bond].s[NCHAINS];
							
							for(int i=start;i<start+nChains*length;i+=length)
								for(int j=i;j<i+length-1;j++)
									mol.push_back(j);
						}
						break;
						case BOND:
						{
							mol.push_back(m[i].getBonds()[bond].s[0]);
							mol.push_back(m[i].getBonds()[bond].s[1]);
						}
						break;
						case BEND:
						{
							mol.push_back(m[i].getBonds()[bond].s[0]);
							mol.push_back(m[i].getBonds()[bond].s[1]);
							mol.push_back(m[i].getBonds()[bond].s[2]);
						}
						break;
					}
					if(mol.size()>0)
					{
						std::sort(mol.begin(), mol.end());
						std::vector<int>::iterator it;
						it=std::unique(mol.begin(), mol.end());
						mol.resize(std::distance(mol.begin(),it-1));
						molPI.push_back(mol);
					}
				}
			}
		}
		//check for contacts between nanoparticles, using potential energy
		std::vector<int> combinedList;
		for(int i=0;i<molPI.size();i++)
			for(int j=0;j<molPI[i].size();j++)
				combinedList.push_back(molPI[i][j]);
		
		//Here???
		//Cell<double> pairInteractions(System.getPositions(), System.readNParticles(), System.readCutoff(), 
		//			      size, &combinedList[0], combinedList.size());
		
		std::string newName("frames_");
		newName+=name;
		newName+=".xyz";
		
		xyzFormat<double> xyzFile(p,nParticles);
		xyzFile.open(newName.c_str(), std::ios::in);
		
		double time=0;
		while(xyzFile.load())
		{
			//or is it stable enough to put here???
			Cell<double> pairInteractions(System.getPositions(), System.readNParticles(), System.readCutoff(), 
						      size, &combinedList[0], combinedList.size());
			
			computePotential cP(System.getPositions(), System.readNParticles(), 
					    System.getTwoBodyUconst(), System.readCutoff(), System.readNTypes());
			cP.initialize();//zero it out
			
			pairInteractions.build();
			
			//new and improved
			pairInteractions.twoWayComparePeriodicIndex(cP);
			
			//std::cout << System.readInitialTime() << '\t' << cP.readPotential() << '\t' << cP.outputCount(type1) << '\n';
			
			double avgPotential=cP.outputPotential()/static_cast<double>(cP.outputCount());
			int nContacts=0;
			double *potentialI=cP.getPotentials();
			
			//we are assuming contacts decrease the average energy
			for(int i=0;i<combinedList.size();i++)
				if(potentialI[i]<avgPotential)
					contacts++;
			
			std::cout << time << '\t' << contacts << std::endl;
				
			time=time+=System.readStoreInterval();
			
			if(sizeFile.is_open())
			{
				double sTime=-1;
				threeVector<double> sCurrent;
				//for the first iteration, I expect this to exit the loop after one read,
				// otherwise, this is a template for subsequent reads in the loop that is
				// coming up.
				while(sTime<time && !sizeFile.eof())
					sizeFile >> sTime >> sCurrent.x >> sCurrent.y >> sCurrent.z;
				time=sTime;//this might be self correcting
				size=sCurrent;
			}
		}
	}
	/*
	xyzFormat<double> xyzFile(p,nParticles);
	xyzFile.open(newName.c_str(), std::ios::in);
	
	while(xyzFile.load())
	{
		//Do data measurements in here
		int a=0;
		for(int i=0;i<nParticles;i++)
			if(p[i].x>System.readSize().x/2.0)
				a++;
		std::cout << a << std::endl;
	}
	*/
	return 0;
}

