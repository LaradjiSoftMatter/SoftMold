/** \brief Extracts the potential between two types (interfacial, not internal).
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
		computePotential(position<double> *particles, int nParticles, double *constants, double cutoff, int nTypes, int type1, int type2)
		{
			p=particles;
			nP=nParticles;
			c=constants;
			cutoffSquared=cutoff*cutoff;
			nT=nTypes;
			t1=type1;
			t2=type2;
			count=new int[nP];
		};
		
		//destructor
		~computePotential()
		{
			delete count;
		};
		
		//set potential to 0
		void initialize()
		{
			potential=0;
			for(int i=0;i<nP;i++)
				count[i]=0;
		};
		
		//output the calculated potential
		double readPotential()
		{
			return potential;
		};
		
		//count by type
		int outputCount(int type)
		{
			int total=0;
			for(int i=0;i<nP;i++)
				if(p[i].type==type)
					total+=count[i];
			return total;
		};
		
		//plain potential calculation
		void operator() (int &i, int &j)
		{
			if((p[i].type==t1 && p[j].type==t2) || (p[i].type==t2 && p[j].type==t1))
			{
				double tPotential=Potential<double>(cutoffSquared, nT, p[i], p[j], c);
				if(tPotential!=0)
				{
					count[i]=1;
					count[j]=1;
				}
				potential+=tPotential;
			}
		};
		
		//minimum image version
		void operator() (int &i, int &j, threeVector<double> &minImg)
		{
			if((p[i].type==t1 && p[j].type==t2) || (p[i].type==t2 && p[j].type==t1))
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
			}
		};
		
	private:
		double potential;
		position<double> *p;
		int nP;
		double *c;
		double cutoffSquared;
		int nT, t1, t2;
		int *count;
};

int main(int argc, char* argv[])
{
	if(argc!=4)
	{
		//so simple, this can't possibly mess it up
		std::cout << "Extracts the potential between two types (interfacial, not internal) from a name.mpd file.\n";
		std::cout << "usage: " << argv[0] << " name type1 type2\n";
		return 0;
	}
	
	//read in options
	char *name=argv[1];
	int type1,type2;
	std::stringstream cmdInp;
	cmdInp << argv[2] << ' ' << argv[3] << '\n';
	cmdInp >> type1 >> type2;
	
	//Configuration variables
	Blob<double> System;
	
	//load old variables, then initialize them
	Script<double, Blob <double> > fileIO(name,std::ios::in,&System);
	fileIO.read();
	fileIO.close();
	
	std::cout << "#time potential count_of_Type_" << type1 << '\n';
	
	//set up our cell lists
	//new and improved
	Cell<double> pairInteractions(System.getPositions(), System.readNParticles(), System.readCutoff(), System.readSize());
	
	//old one
	//CellOpt<double, Potential<double>, Force <double> > pairInteractions(System.getPositions(), System.getAccelerations(), 
	//	System.getTwoBodyFconst(), System.getTwoBodyUconst(), System.readNParticles(), System.readNTypes(), System.readSize(),
	//	System.readPeriodic(), System.readCutoff());
	
	//a comparison object, old one
	//compareTypes compareObj(type1,type2);
	
	//new one
	computePotential cP(System.getPositions(), System.readNParticles(), System.getTwoBodyUconst(), System.readCutoff(), System.readNTypes(), type1, type2);
	cP.initialize();//zero it out
	
	//build the cell list, either one
	pairInteractions.build();
	
	//compute the potential, old one
	//std::cout << pairInteractions.computePotential(compareObj) << '\n';
	
	//new and improved
	pairInteractions.twoWayComparePeriodicIndex(cP);
	
	std::cout << System.readInitialTime() << '\t' << cP.readPotential() << '\t' << cP.outputCount(type1) << '\n';
	
	return 0;
}

