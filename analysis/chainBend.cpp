//This is a very short and simple program that shows you how to modify a variable already present in an mpd script.
//The reason this is used is because some variables (forces, potentials, geometry, etc...) are not easily modified
// in the mpd script. Using a file that opens the mpd script, modifies, and then saves it (as below) is far easier
// than modifying the actual variable. Most of it looks like the run.cpp file, it just replaces the "run loop"
// with the modification.

//For molecular dynamics forces and potentials
#include "../include/MD.h"

//For the molecular dynamics variables
#include "../include/system.h"

bool compareAngles(fourVector<double> a, fourVector<double> b)
{
	return a.t>b.t;
}

class bendRange {
	public:
		bendRange(std::vector<fourVector<double> > bendMap, double cutoff)
		{
			bM=bendMap;
			cutoffSqr=cutoff*cutoff;
			for(int i=0;i<bM.size();i++)
			{
				avgBend.push_back(bM[i].t);
				nNeighbors.push_back(1);
			}
		};
		
		~bendRange()
		{};
		
		const std::vector<double> getBendAverages()
		{
			return avgBend;
		};
		
		void finalize()
		{
			for(int i=0;i<avgBend.size();i++)
				avgBend[i]/=static_cast<double>(nNeighbors[i]);
		};
		
		void operator () (int &i, int &j, threeVector<double> minImg)
		{
			threeVector<double> d;
			d.x=bM[i].x-bM[j].x+minImg.x;
			d.y=bM[i].y-bM[j].y+minImg.y;
			d.z=bM[i].z-bM[j].z+minImg.z;
			if(d.x*d.x+d.y*d.y+d.z*d.z<cutoffSqr)
			{
				nNeighbors[i]++;
				nNeighbors[j]++;
				avgBend[i]+=bM[j].t;
				avgBend[j]+=bM[i].t;
			}	
		};
		
		void operator () (int &i, int &j)
		{
			threeVector<double> d;
			d.x=bM[i].x-bM[j].x;
			d.y=bM[i].y-bM[j].y;
			d.z=bM[i].z-bM[j].z;
			if(d.x*d.x+d.y*d.y+d.z*d.z<cutoffSqr)
			{
				nNeighbors[i]++;
				nNeighbors[j]++;
				avgBend[i]+=bM[j].t;
				avgBend[j]+=bM[i].t;
			}	
		};
	private:
		std::vector<fourVector<double> > bM;
		std::vector<double> avgBend;
		std::vector<int> nNeighbors;
		double cutoffSqr;
};

int main(int argc, char* argv[])
{
	if(argc!=2 && argc!=3 && argc!=4 && argc!=5)
	{
		//so simple, this can't possibly mess it up
		std::cout << "To see available chain type molecules:\n";
		std::cout << "\tusage: " << argv[0] << " name\n";
		
		std::cout << "To get the bend average of chain type molecules:\n";
		std::cout << "\tusage: " << argv[0] << " name molIndex\n";
		
		std::cout << "To get the bend map of chain type molecules (and bendHist_name.dat):\n";
		std::cout << "\tusage: " << argv[0] << " name molIndex bendAngleSteps\n";
		
		std::cout << "To get the spatial average bend map of chain type molecules:\n";
		std::cout << "\tusage: " << argv[0] << " name molIndex bendAngleSteps cutoff\n";
		
		return 0;
	}
	
	///read in options
	char *name=argv[1];
	
	///Configuration variables
	Blob<double> System;
	
	///load old variables, then initialize them
	Script<double, Blob <double> > fileIO(name,std::ios::in,&System);
	fileIO.read();
	fileIO.close();
	
	std::string framesName("frames_");
	framesName+=name;
	framesName+=".xyz";
	
	position<double> *p=System.getPositions();
	int nParticles=System.readNParticles();
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
		}
	}
	
	//xyzFormat<double> xyzFile(p,nParticles);
	//xyzFile.open(framesName.c_str(), std::ios::in);
	
	if(argc>=3)
	{
		threeVector<double> size=System.readSize();
		int molIndex, bendAngleSteps;
		double cutoff;
		std::stringstream cmdArg;
		for(int i=2;i<argc;i++)
			cmdArg << argv[i] << '\t';
		cmdArg >> molIndex;
		if(argc>=4)
			cmdArg >> bendAngleSteps;
		if(argc>=5)
			cmdArg >> cutoff;
		
		if(molIndex>=System.readNMolecules() || molIndex<0)
		{
			std::cerr << "Molecule index " << molIndex << " is out of bounds [0, " << System.readNMolecules() << ").\n";
			return 0;
		}
		if(m[molIndex].readType()!=CHAIN)
		{
			std::cerr << "Molecule type is wrong!\n";
			return 0;
		}
		double time=0;
		if(argc==3)
			std::cout << "#time mean stdDev" << std::endl;
		//while(xyzFile.load())
		{
			time+=System.readStoreInterval();
			std::vector< fourVector<double> > bendMap;
			std::vector< position<double> > pMap;
			
			for(int bond=0;bond<m[molIndex].readNBond();bond++)
			{
				int start=m[molIndex].getBonds()[bond].s[START];
				int length=m[molIndex].getBonds()[bond].s[CHAINLENGTH];
				int nChains=m[molIndex].getBonds()[bond].s[NCHAINS];
				
				for(int i=start;i<start+nChains*length;i+=length)
				{
					fourVector<double> mapValue;
					mapValue.x=0;
					mapValue.y=0;
					mapValue.z=0;
					mapValue.t=0;
					//center of mass
					for(int j=i;j<i+length;j++)
					{
						mapValue.x+=p[j].x;
						mapValue.y+=p[j].y;
						mapValue.z+=p[j].z;
					}
					mapValue.x/=static_cast<double>(length);
					mapValue.y/=static_cast<double>(length);
					mapValue.z/=static_cast<double>(length);
					
					//average bend
					for(int j=i;j<i+length-2;j++)
					{
						threeVector<double> da,db;
						da.x=p[j].x-p[j+1].x;
						if(da.x>size.x/2.0) da.x-=size.x;
						if(da.x<-size.x/2.0) da.x+=size.x;
						da.y=p[j].y-p[j+1].y;
						if(da.y>size.y/2.0) da.y-=size.y;
						if(da.y<-size.y/2.0) da.y+=size.y;
						da.z=p[j].z-p[j+1].z;
						if(da.z>size.z/2.0) da.z-=size.z;
						if(da.z<-size.z/2.0) da.z+=size.z;
						
						db.x=p[j+1].x-p[j+2].x;
						if(db.x>size.x/2.0) db.x-=size.x;
						if(db.x<-size.x/2.0) db.x+=size.x;
						db.y=p[j+1].y-p[j+2].y;
						if(db.y>size.y/2.0) db.y-=size.y;
						if(db.y<-size.y/2.0) db.y+=size.y;
						db.z=p[j+1].z-p[j+2].z;
						if(db.z>size.z/2.0) db.z-=size.z;
						if(db.z<-size.z/2.0) db.z+=size.z;
						
						da=unitVector(da);
						db=unitVector(db);
						
						//should be a value from 0 to 1
						double costheta=(da.x*db.x)+(da.y*db.y)+(da.z*db.z);
						mapValue.t+=costheta;
					}
					mapValue.t/=static_cast<double>(length-2);
					bendMap.push_back(mapValue);
					position<double> comP;
					comP.type=1;
					comP.x=mapValue.x;
					comP.y=mapValue.y;
					comP.z=mapValue.z;
					pMap.push_back(comP);
				}
			}
			double mean=0;
			for(int i=0;i<bendMap.size();i++)
				mean+=bendMap[i].t;
			mean/=static_cast<double>(bendMap.size());
			double stdDev=0;
			for(int i=0;i<bendMap.size();i++)
				stdDev+=(bendMap[i].t-mean)*(bendMap[i].t-mean);
			stdDev/=static_cast<double>(bendMap.size()-1);
			stdDev=sqrt(stdDev);
			if(argc==3)
				std::cout << time << '\t' << mean << '\t' << stdDev << std::endl;
				//std::cout << "Mean cos(theta)=" << mean << " and standard deviation=" << stdDev << std::endl;
			
			double min=1.0, max=0.0, stepSize=0.0;
			if(argc>=4)
			{
				for(int i=0;i<bendMap.size();i++)
				{
					min=(bendMap[i].t<min)?bendMap[i].t:min;
					max=(bendMap[i].t>max)?bendMap[i].t:max;
				}
				stepSize=(max-min)/static_cast<double>(bendAngleSteps);
			}
			if(argc==4)
			{	
				//std::cout << bendMap.size() << "\ntest\n";
				//for(int i=0;i<bendMap.size();i++)
				//	std::cout << floor((bendMap[i].t-min)/stepSize) << '\t' << bendMap[i].x << '\t' << bendMap[i].y << '\t' << bendMap[i].z << '\n';
				
				std::string bendHistName("bendHist_");
				bendHistName+=name;
				bendHistName+=".dat";
				
				std::fstream bendHistFile;
				bendHistFile.open(bendHistName.c_str(), std::ios::out | std::ios::app);
				if(bendHistFile.is_open())
				{
					std::vector<double> bendHist(bendAngleSteps+1,0);
					for(int i=0;i<bendMap.size();i++)
						bendHist[floor((bendMap[i].t-min)/stepSize)]++;
					for(int i=0;i<bendHist.size();i++)
					{
						bendHist[i]/=static_cast<double>(bendMap.size());
						bendHistFile << (static_cast<double>(i)*stepSize)+min << '\t' << bendHist[i] << '\n';
					}
					bendHistFile << std::endl;
					bendHistFile.close();
				}
			}
			
			//new one
			bendRange bR(bendMap,cutoff);
			
			if(argc>=5)
			{
				Cell<double> pairInteractions(&pMap[0], pMap.size(), cutoff, System.readSize());
				
				//build the cell list, either one
				pairInteractions.build();
				pairInteractions.twoWayComparePeriodicIndex(bR);
				bR.finalize();
			}
			if(argc==5)
			{
				std::cout << bendMap.size() << "\ntest\n";
				std::vector<double> bendAvg=bR.getBendAverages();
				for(int i=0;i<bendMap.size();i++)
				{
					std::cout << floor((bendAvg[i]-min)/stepSize) << '\t';
					std::cout << bendMap[i].x << '\t' << bendMap[i].y << '\t' << bendMap[i].z << '\n';
				}
			}
			
			std::vector< fourVector<double> > averageBendMap;
			
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

