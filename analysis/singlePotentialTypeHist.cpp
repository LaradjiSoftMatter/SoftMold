/** \brief Extracts the potential from all particles and generates a histagram based on type.
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
			count=new int[nP];
			potentials=new double[nT*nT];
		};
		
		//destructor
		~computePotential()
		{
			delete count;
			delete potentials;
		};
		
		//set potential to 0
		void initialize()
		{
			potential=0;
			for(int i=0;i<nP;i++)
				count[i]=0;
			for(int i=0;i<nT*nT;i++)
				potentials[i]=0;
		};
		
		//output the calculated potential
		double readPotential(int i)
		{
			#ifdef WARNINGS_ENABLED
			if(i>=nT*nT)
			{
				std::cout << i << " (i) exceeds " << nT*nT-1 << " (nTypes*nTypes-1) " << std::endl;
				#ifdef ERRORS_ENABLED
					throw 0;
				#endif
			}
			#endif
			return potentials[i];
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
			double tPotential=Potential<double>(cutoffSquared, nT, p[i], p[j], c);
			if(tPotential!=0)
			{
				count[i]=1;
				count[j]=1;
			}
			potentials[p[i].type*nT+p[j].type]+=tPotential;
			potentials[p[j].type*nT+p[i].type]+=tPotential;
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
			potentials[p[i].type*nT+p[j].type]+=tPotential;
			potentials[p[j].type*nT+p[i].type]+=tPotential;
		};
		
	private:
		double potential;
		position<double> *p;
		int nP;
		double *c;
		double cutoffSquared;
		int nT;
		int *count;
		double *potentials;
};

int main(int argc, char* argv[])
{
	if(argc!=2 && argc!=3)
	{
		//so simple, this can't possibly mess it up
		std::cout << "Extracts the potential between all pairs and outputs the resulting type interface matrix.\n";
		std::cout << "usage: " << argv[0] << " name\n";
		std::cout << "usage: " << argv[0] << " name frame\n";
		return 0;
	}
	
	//read in options
	char *name=argv[1];
	
	int frame=-1;//last frame
	if(argc==3)
	{
		std::stringstream cmdArg;
		cmdArg << argv[2];
		cmdArg >> frame;
	}
	
	//Configuration variables
	Blob<double> System;
	
	//load old variables, then initialize them
	Script<double, Blob <double> > fileIO(name,std::ios::in,&System);
	fileIO.read();
	fileIO.close();
	
	//set up our cell lists
	//new and improved
	Cell<double> pairInteractions(System.getPositions(), System.readNParticles(), System.readCutoff(), System.readSize());
	
	//new one
	computePotential cP(System.getPositions(), System.readNParticles(), System.getTwoBodyUconst(), System.readCutoff(), System.readNTypes());
	cP.initialize();//zero it out
	
	//build the cell list, either one
	pairInteractions.build();
	
	//compute the potential, old one
	//std::cout << pairInteractions.computePotential(compareObj) << '\n';
	
	if(frame>-1)
	{
		std::string newName("frames_");
		newName+=name;
		newName+=".xyz";
		
		position<double> *p=System.getPositions();
		int nParticles=System.readNParticles();
		
		xyzFormat<double> xyzFile(p,nParticles);
		xyzFile.open(newName.c_str(), std::ios::in);
		
		xyzFile.load(frame);
		
		xyzFile.close();
	}
	
	//new and improved
	pairInteractions.twoWayComparePeriodicIndex(cP);
	
	for(int i=0;i<System.readNTypes()*System.readNTypes();i++)
	{
		std::cout << cP.readPotential(i) << "\t\t";
		if((i+1)%System.readNTypes()==0)
			std::cout << std::endl;
	}
	std::cout << std::endl;
	return 0;
}

