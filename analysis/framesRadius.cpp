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
	if(argc<3)
	{
		//so simple, this can't possibly mess it up
		std::cerr << "usage: " << argv[0] << " name centerType1 centerType2 centerType3 ... centerTypeN\n";
		std::cerr << "Output of radial distribution is only to the first encountered edge.\n";
		std::cerr << "\"centerType\"s are the types to center about.\n";
		return 0;
	}
	
	///read in options
	char *name=argv[1];
	std::vector<int> centerType;
	
	std::stringstream cmdArg;
	for(int i=2;i<argc;i++)
		cmdArg << argv[i] << ' ';
	for(int i;cmdArg >> i;)
		centerType.push_back(i);
	
	///Configuration variables
	Blob<double> System;
	
	///load old variables, then initialize them
	Script<double, Blob <double> > fileIO(name,std::ios::in,&System);

	fileIO.read();

	fileIO.close();
	
	//our current size
	threeVector<double> size=System.readSize();
	
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
		size=sCurrent;
	}
	
	//Get our frames file
	std::string framesName("frames_");
	framesName+=name;
	framesName+=".xyz";
	
	position<double> *p=System.getPositions();
	int nParticles=System.readNParticles();
	
	xyzFormat<double> xyzFile(p,nParticles);
	xyzFile.open(framesName.c_str(), std::ios::in);
	
	double time=0;
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
		std::vector<int> pI;
		for(int i=0;i<nParticles;i++)
			for(int j=0;j<centerType.size();j++)
				if(centerType[j]==p[i].type)
					pI.push_back(i);
		
		
		///Locate center of mass for centerTypes
		threeVector<double> com=0;
		int nPCom=0;
		for(int i:pI)
		{
			threeVector<double> d;
			d.x=p[i].x-p[0].x;
			d.y=p[i].y-p[0].y;
			d.z=p[i].z-p[0].z;
			
			if(d.x>size.x/2.0) p[i].x-=size.x;
			if(d.y>size.y/2.0) p[i].y-=size.y;
			if(d.z>size.z/2.0) p[i].z-=size.z;
			if(d.x<-size.x/2.0) p[i].x+=size.x;
			if(d.y<-size.y/2.0) p[i].y+=size.y;
			if(d.z<-size.z/2.0) p[i].z+=size.z;
			
			com.x+=p[i].x;
			com.y+=p[i].y;
			com.z+=p[i].z;
			nPCom++;
		}
		if(nPCom==0)
		{
			std::cerr << "Types selected are not present!\n";
			return -1;
		}
		
		com.x/=(double)nPCom;
		com.y/=(double)nPCom;
		com.z/=(double)nPCom;
		
		while(com.x>size.x) com.x-=size.x;
		while(com.y>size.y) com.y-=size.y;
		while(com.z>size.z) com.z-=size.z;
		while(com.x<0) com.x+=size.x;
		while(com.y<0) com.y+=size.y;
		while(com.z<0) com.z+=size.z;
		
		///Locate center of mass for centerTypes
		double radius=0;
		for(int i:pI)
		{
			threeVector<double> d;
			d.x=p[i].x-com.x;
			d.y=p[i].y-com.y;
			d.z=p[i].z-com.z;
			
			if(d.x>size.x/2.0) d.x-=size.x;
			if(d.y>size.y/2.0) d.y-=size.y;
			if(d.z>size.z/2.0) d.z-=size.z;
			if(d.x<-size.x/2.0) d.x+=size.x;
			if(d.y<-size.y/2.0) d.y+=size.y;
			if(d.z<-size.z/2.0) d.z+=size.z;
			
			radius+=sqrt(d.x*d.x+d.y*d.y+d.z*d.z);
		}
		radius/=(double)nPCom;
		std::cout << time << ' ' << radius << std::endl;
		std::cerr << time << std::endl;
		time+=System.readStoreInterval();
	}
	
	
	return 0;
}

