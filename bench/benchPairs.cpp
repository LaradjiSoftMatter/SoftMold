/** \brief Extracts the potential between two types (interfacial, not internal).
 * Program that extracts the potential between two types. It does not extract internal potential
 * energy within the types (singlePotentialTypes.cpp does that).
 */

#include "../include/systemMD.h"

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
			potential=0;
		};
		
		position<double>* getPositions(){return p;};
		int readNParticles(){return nP;};
		double *getConstants(){return c;};
		double readCutoff(){return sqrt(cutoffSquared);};
		int readNTypes(){return nT;};
		
		//destructor
		~computePotential()
		{
		};
		
		//set potential to 0
		void initialize()
		{
			potential=0;
		};
		
		//output the calculated potential
		double outputPotential()
		{
			return potential;
		};
		
		
		inline computePotential operator= (computePotential &oldObj)
		{
			computePotential newObj(oldObj.getPositions(), oldObj.readNParticles(), oldObj.getConstants(), oldObj.readCutoff(), oldObj.readNTypes());
			return newObj;
		}
		
		
		//plain potential calculation
		inline void operator() (int &i, int &j)
		{
			potential+=Potential<double>(cutoffSquared, nT, p[i], p[j], c);
		};
		
		//plain potential calculation
		inline void operator() (position<double> &pi, position<double> &pj)
		{
			potential+=Potential<double>(cutoffSquared, nT, pi, pj, c);
		};
		
		//minimum image version
		inline void operator() (int &i, int &j, threeVector<double> &minImg)
		{
			position<double> pj=p[j];
			pj.x+=minImg.x;
			pj.y+=minImg.y;
			pj.z+=minImg.z;
			potential+=Potential<double>(cutoffSquared, nT, p[i], pj, c);
		};
		
	private:
		double potential;
		position<double> *p;
		int nP;
		double *c;
		double cutoffSquared;
		int nT;
};

int main(int argc, char* argv[])
{
	if(argc!=3)
	{
		//so simple, this can't possibly mess it up
		std::cout << "Determines potential from a name.mpd file. Also outputs the time in doing so. This is for benchmark purposes.\n";
		std::cout << "usage: " << argv[0] << " name repeat\n";
		return 0;
	}
	
	//read in options
	char *name=argv[1];
	std::stringstream cmdInp;
	cmdInp << argv[2] << '\n';
	int repeat;
	cmdInp >> repeat;
	
	//Configuration variables
	Blob<double> System;
	
	//load old variables, then initialize them
	Script<double, Blob <double> > fileIO(name,std::ios::in,&System);
	fileIO.read();
	fileIO.close();
	
	//new one
	computePotential cP(System.getPositions(), System.readNParticles(), System.getTwoBodyUconst(), System.readCutoff(), System.readNTypes());
	
	//set up our cell lists
	//new and improved
	Cell<double> pairInteractions(System.getPositions(), System.readNParticles(), System.readCutoff(), System.readSize());
	
	//build the cell list, either one
	pairInteractions.build();
	
	double someCrap=0;
	
	time_t start=time(NULL);
	
	for(int i=0;i<repeat;i++)
	{
		cP.initialize();//zero it out
		
		//new and improved
		//computePotential newObj
		someCrap+=pairInteractions.twoWayComparePeriodic<computePotential>(cP).outputPotential();
		//someCrap+=newObj.outputPotential();
	}
	
	time_t finish=time(NULL);
	
	computePotential newObj=pairInteractions.twoWayComparePeriodic<computePotential>(cP);
	
	std::cout << "New version time to complete " << repeat << " executions:" << (finish-start) << '\n';
	std::cout << "Potential calculated: " << newObj.outputPotential() << '\t' << (someCrap/double(repeat)) << '\n';
	
	//old one
	CellOpt<double, Potential<double>, Force <double> > pairInteractionsOld(System.getPositions(), System.getAccelerations(), 
		System.getTwoBodyFconst(), System.getTwoBodyUconst(), System.readNParticles(), System.readNTypes(), System.readSize(),
		System.readPeriodic(), System.readCutoff());
	
	//build the cell list, either one
	pairInteractionsOld.build();
	
	start=time(NULL);
	
	for(int i=0;i<repeat;i++)
	{
		//new and improved
		pairInteractionsOld.computePotential();
	}
	
	finish=time(NULL);
	
	std::cout << "Old version time to complete " << repeat << " executions:" << (finish-start) << '\n';
	std::cout << "Potential calculated: " << pairInteractionsOld.computePotential() << '\n';
	
	return 0;
}

