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
				double dx=p[i].x-p[j].x;
				double dy=p[i].y-p[j].y;
				double dz=p[i].z-p[j].z;
				if(dx*dx+dy*dy+dz*dz<cutoffSquared)
				{
					count[i]=1;
					count[j]=1;
				}
				double tPotential=Potential<double>(cutoffSquared, nT, p[i], p[j], c);
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
				double dx=p[i].x-pj.x;
				double dy=p[i].y-pj.y;
				double dz=p[i].z-pj.z;
				if(dx*dx+dy*dy+dz*dz<cutoffSquared)
				{
					count[i]=1;
					count[j]=1;
				}
				double tPotential=Potential<double>(cutoffSquared, nT, p[i], pj, c);
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
	if(argc!=5)
	{
		//so simple, this can't possibly mess it up
		std::cout << "Extracts the potential between two types (interfacial, not internal) from a name.mpd file.\n";
		std::cout << "usage: " << argv[0] << " name type1 type2 cutoff\n";
		return 0;
	}
	
	//read in options
	char *name=argv[1];
	int type1,type2;
	std::stringstream cmdInp;
	double cutoff;
	cmdInp << argv[2] << ' ' << argv[3] << ' ' << argv[4] << '\n';
	cmdInp >> type1 >> type2 >> cutoff;
	
	//Configuration variables
	Blob<double> System;
	
	//load old variables, then initialize them
	Script<double, Blob <double> > fileIO(name,std::ios::in,&System);
	fileIO.read();
	fileIO.close();
	
	std::cout << "#time potential unique_count_of_Type_" << type1 << '\n';
	
	
	std::string newName("frames_");
	newName+=name;
	newName+=".xyz";
	position<double> *p=System.getPositions();
	int nParticles=System.readNParticles();
	xyzFormat<double> xyzFile(p,nParticles);
	xyzFile.open(newName.c_str(), std::ios::in);
	
	//grab sizes
	double time=0;
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
		time=sTime;
		size=sCurrent;
		//std::cout << size.x << '\t' << size.y << '\t' << size.z << '\n';
	}
	
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
				sizeFile >> sTime >> sCurrent.x >> sCurrent.y >> sCurrent.z;
			size=sCurrent;
		}
		//set up our cell lists
		Cell<double> pairInteractions(System.getPositions(), System.readNParticles(), cutoff, size);
		
		//potential object
		computePotential cP(System.getPositions(), System.readNParticles(), System.getTwoBodyUconst(), cutoff, System.readNTypes(), type1, type2);
		//zero count and potential
		cP.initialize();
		
		//build the cell list, either one
		pairInteractions.build();
		
		//compute the potential and count
		pairInteractions.twoWayComparePeriodicIndex(cP);
		
		//output values
		std::cout << time << '\t' << cP.readPotential() << '\t' << cP.outputCount(type1) << std::endl;
		time+=System.readStoreInterval();
	}
	
	return 0;
}

