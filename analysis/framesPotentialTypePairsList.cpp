/** \brief Extracts the potential between two types (interfacial, not internal) per frame.
 * Program that extracts the potential between two types. It does not extract internal potential
 * energy within the types (singlePotentialTypes.cpp does that).
 */

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
			potentialList=new double[nP];
		};
		
		//destructor
		~computePotential()
		{
			delete potentialList;
		};
		
		//set potential to 0
		void initialize()
		{
			potential=0;
			
			//Potential list per particle:
			for(int i=0;i<nP;i++)
				potentialList[i]=0;

		};
		
		//output the calculated potential
		double readPotential()
		{
			return potential;
		};

		double *getPotentialList()
		{
			return potentialList;
		};

		//plain potential calculation
		void operator() (int &i, int &j)
		{
			double tPotential=Potential<double>(cutoffSquared, nT, p[i], p[j], c);
			potential+=tPotential;
			potentialList[i]+=tPotential;
			potentialList[j]+=tPotential;
		};
		
		//minimum image version
		void operator() (int &i, int &j, threeVector<double> &minImg)
		{
			position<double> pj=p[j];
			pj.x+=minImg.x;
			pj.y+=minImg.y;
			pj.z+=minImg.z;
			double tPotential=Potential<double>(cutoffSquared, nT, p[i], pj, c);
			potentialList[i]+=tPotential;
			potentialList[j]+=tPotential;
		};
		
	private:
		double potential;
		position<double> *p;
		int nP;
		double *c;
		double cutoffSquared;
		int nT;
		double *potentialList;
};

int main(int argc, char* argv[])
{
	if(argc!=4)
	{
		//so simple, this can't possibly mess it up
		std::cout << "Extracts the potentials of a particular type from a name.mpd file, and lists them.\n";
		std::cout << "usage: " << argv[0] << " name type deltaP\n";
		return 0;
	}
	
	//read in options
	char *name=argv[1];
	int type;
	double deltaP;
	std::stringstream cmdInp;
	cmdInp << argv[2] << ' ' << argv[3];
	cmdInp >> type >> deltaP;
	
	//Configuration variables
	Blob<double> System;
	
	//load old variables, then initialize them
	Script<double, Blob <double> > fileIO(name,std::ios::in,&System);
	fileIO.read();
	fileIO.close();
	
	//set up our cell lists
	Cell<double> pairInteractions(System.getPositions(), System.readNParticles(), System.readCutoff(), System.readSize());
	
	//potential object
	computePotential cP(System.getPositions(), System.readNParticles(), System.getTwoBodyUconst(), System.readCutoff(), System.readNTypes());
	
	std::string newName("frames_");
	newName+=name;
	newName+=".xyz";
	
	position<double> *p=System.getPositions();
	int nParticles=System.readNParticles();
	
	xyzFormat<double> xyzFile(p,nParticles);
	xyzFile.open(newName.c_str(), std::ios::in);
	
	double time=0;
	while(xyzFile.load())
	{
		std::vector<double> hist;
		time+=System.readStoreInterval();
		//zero count and potential
		cP.initialize();
		
		//build the cell list, either one
		pairInteractions.build();
		
		//compute the potential and count
		pairInteractions.twoWayComparePeriodicIndex(cP);
		
		double min=0;
		double max=0;
		
		//output values
		for(int i=0;i<System.readNParticles();i++)
		{
			if(p[i].type==type)
			{
				if(cP.getPotentialList()[i]>max)
					max=cP.getPotentialList()[i];
				if(cP.getPotentialList()[i]<min)
					min=cP.getPotentialList()[i];
			}
		}
		
		for(int i=0;i<System.readNParticles();i++)
		{
			if(p[i].type==type)
			{
				int index=static_cast<int>((cP.getPotentialList()[i]-min)/deltaP);
				if(index<0)
				{
					std::cerr << index << " is out of bounds." << std::endl;
					std::cerr << min << '\t' << cP.getPotentialList()[i] << std::endl;
				}
				
				if(index>=hist.size())
				{
					hist.resize(index+1,0);
				}
				hist[index]++;
			}
		}
		for(int i=0;i<hist.size();i++)
			std::cout << static_cast<double>(i)*deltaP+min << '\t' << hist[i] << std::endl;
		std::cout << std::endl;
	}
	
	return 0;
}

