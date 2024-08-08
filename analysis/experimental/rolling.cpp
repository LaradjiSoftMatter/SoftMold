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
	if(argc!=3)
	{
		//so simple, this can't possibly mess it up
		std::cout << "usage: " << argv[0] << " name molIndex\n";
		return 0;
	}
	
	///read in options
	char *name=argv[1];
	int molIndex=atoi(argv[2]);
	
	///Configuration variables
	Blob<double> System;
	
	
	///load old variables, then initialize them
	Script<double, Blob <double> > fileIO(name,std::ios::in,&System);
	fileIO.read();
	fileIO.close();
	
	std::string newName("frames_");
	newName+=name;
	newName+=".xyz";
	
	position<double> *p=System.getPositions();
	int nParticles=System.readNParticles();
	threeVector<double> s=System.readSize();
	
	//check if size varies
	std::string newName("size_");
	newName+=name;
	newName+=".dat";
	
	std::fstream sizeFile;
	sizeFile.open(newName.c_str(), std::ios::in);
	if(sizeFile.is_open())
	{
		double sTime=-1;
		threeVector<double> sCurrent;
		//for the first iteration, I expect this to exit the loop after one read,
		// otherwise, this is a template for subsequent reads in the loop that is
		// coming up.
		while(sTime<0 && !sizeFile.eof())
			sizeFile >> sTime >> sCurrent.x >> sCurrent.y >> sCurrent.z;
		s=sCurrent;
	}
	
	xyzFormat<double> xyzFile(p,nParticles);
	xyzFile.open(newName.c_str(), std::ios::in);
	
	if(System.getMolecules()[molIndex].type==BOND)
	{
		molecule<double, fourVector<int> > *m=System.getMolecules();
		position<double> previousBase=p[m[molIndex].getBonds()[0].s[0]];
		for(int i=0;xyzFile.load();i++)
		{
			double time=static_cast<double>(i)*System.readStoreInterval();
			if(sizeFile.is_open())
			{
				double sTime=-1;
				threeVector<double> sCurrent=s;
				//for the first iteration, I expect this to exit the loop after one read,
				// otherwise, this is a template for subsequent reads in the loop that is
				// coming up.
				while(sTime<time && !sizeFile.eof())
					sizeFile >> sTime >> sCurrent.x >> sCurrent.y >> sCurrent.z;
				s=sCurrent;
			}
			threeVector<double> rAvg(0), pAvg(0);
			
			for(int bond=0;bond<m[molIndex].readNBonds();i++)
			{
				int firstParticle=m[molIndex].getBonds()[bond].s[0];
				int secondParticle=m[molIndex].getBonds()[bond].s[1];
				threeVector<T> d;
				d.x=p[firstParticle].x-p[secondParticle].x;
				d.y=p[firstParticle].y-p[secondParticle].y;
				d.z=p[firstParticle].z-p[secondParticle].z;
				d.x-=(d.x>s.x/2.0)?s.x:0.0;
				d.x+=(d.x<-s.x/2.0)?s.x:0.0;
				d.y-=(d.y>s.y/2.0)?s.y:0.0;
				d.y+=(d.y<-s.y/2.0)?s.y:0.0;
				d.z-=(d.z>s.z/2.0)?s.z:0.0;
				d.z+=(d.z<-s.z/2.0)?s.z:0.0;
				d=unitVector(d);
				
				int baseParticle=m[molIndex].getBonds()[0].s[0];
				
			}
			
			
			
			std::cout << time << '\t' << avg.x << '\t' << avg.y << '\t' << avg.z << '\t';
			std::cout << 
		}
	}
	if(sizeFile.is_open())
		sizeFile.close();
	xyzFile.close();
	return 0;
}

