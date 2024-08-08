/** \brief Extracts the potential between two types (interfacial, not internal) per frame.
 * Program that extracts the potential between two types. It does not extract internal potential
 * energy within the types (singlePotentialTypes.cpp does that).
 */

//For molecular dynamics forces and potentials
#include "../include/MD.h"

//For the molecular dynamics variables
#include "../include/system.h"

int main(int argc, char* argv[])
{
	if(argc!=4)
	{
		//so simple, this can't possibly mess it up
		std::cout << "Extracts the distance between two indices from a name.mpd file.\n";
		std::cout << "usage: " << argv[0] << " name index1 index2\n";
		return 0;
	}
	
	//read in options
	char *name=argv[1];
	int index1,index2;
	std::stringstream cmdInp;
	cmdInp << argv[2] << ' ' << argv[3] << '\n';
	cmdInp >> index1 >> index2;
	
	//Configuration variables
	Blob<double> System;
	
	//load old variables, then initialize them
	Script<double, Blob <double> > fileIO(name,std::ios::in,&System);
	fileIO.read();
	fileIO.close();
	
	std::cout << "#time distance("<< index1 << ", " << index2 << ")" << std::endl;
	
	std::string newName("frames_");
	newName+=name;
	newName+=".xyz";
	
	position<double> *p=System.getPositions();
	int nParticles=System.readNParticles();
	
	xyzFormat<double> xyzFile(p,nParticles);
	xyzFile.open(newName.c_str(), std::ios::in);
	
	//grab sizes
	threeVector<double> size=System.readSize();
	
	std::string sizeName("size_");
	sizeName+=name;
	sizeName+=".dat";
	
	std::fstream sizeFile;
	sizeFile.open(sizeName.c_str(), std::ios::in);
	double time=0;
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
	
	int frame=0;
	while(xyzFile.load())
	{
		threeVector<double> d;
		d.x=p[index1].x-p[index2].x;
		d.y=p[index1].y-p[index2].y;
		d.z=p[index1].z-p[index2].z;
		d.x-=(d.x>size.x/2.0)?size.x:0.0;
		d.y-=(d.y>size.y/2.0)?size.y:0.0;
		d.z-=(d.z>size.z/2.0)?size.z:0.0;
		d.x+=(d.x<-size.x/2.0)?size.x:0.0;
		d.y+=(d.y<-size.y/2.0)?size.y:0.0;
		d.z+=(d.z<-size.z/2.0)?size.z:0.0;
		
		//output values
		time=double(frame)*System.readStoreInterval();
		std::cout << time << '\t' << sqrt(d.x*d.x+d.y*d.y+d.z*d.z) << std::endl;
		frame++;
		time=double(frame)*System.readStoreInterval();
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
	}
	
	return 0;
}

