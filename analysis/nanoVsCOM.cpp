/** \brief Calculates eigenvectors for moment of inertia tensor.
 */

#include "../include/MD.h"
#include "../include/system.h"
#include <ctime>

int main(int argc, char* argv[])
{
	if(argc!=5)
	{
		//so simple, this can't possibly mess it up
		std::cerr << "Calculates distance between NPs and vesicle center of mass\n";
		std::cerr << "\tusage: " << argv[0] << " name type nanoType cutoff\n";
		return 0;
	}
	
	//read in options
	char *name=argv[1];
	int type,nanoType;
	std::stringstream cmdInp;
	double cutoff;
	cmdInp << argv[2] << ' ' << argv[3] << ' ' << argv[4];
	cmdInp >> type >> nanoType >> cutoff;
	
	//Configuration variables
	Blob<double> System;
	
	//load old variables, then initialize them
	Script<double, Blob <double> > fileIO(name,std::ios::in,&System);
	fileIO.read();
	fileIO.close();
	
	//Get our frames file
	position<double> *p=System.getPositions();
	int nParticles=System.readNParticles();
	
	//This one is now handling allocation of the 'p' pointer.
	std::string newName("frames_");
	newName+=name;
	newName+=".xyz";
	xyzFormat<double> xyzFile(p,nParticles);
	xyzFile.open(newName.c_str(), std::ios::in);
	
	std::vector<int> typeI, nanoI;
	for(int i=0;i<nParticles;i++)
	{
		if(p[i].type==type)
			typeI.push_back(i);
		if(p[i].type==nanoType)
			nanoI.push_back(i);
	}
	
	time_t start;
	std::time(&start);
	int frame=0;
	double time=0;
	
	//grab sizes
	threeVector<double> size=System.readSize();
	
	std::string sizeName("size_");
	sizeName+=argv[1];
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
	
	std::vector<std::vector<double>> moiELast;
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
		
		double cutoffSqr=cutoff*cutoff;
		
		//center of mass
		threeVector<double> com=0;
		int iCount=0;
		for(auto& i:typeI)
		{
			bool exclude=false;
			for(auto& j:nanoI)
			{
				threeVector<double> d=0;
				d.x=p[i].x-p[j].x;
				d.y=p[i].y-p[j].y;
				d.z=p[i].z-p[j].z;
				d.x+=(d.x>size.x/2.0)?-size.x:0;
				d.y+=(d.y>size.y/2.0)?-size.y:0;
				d.z+=(d.z>size.z/2.0)?-size.z:0;
				d.x+=(d.x<-size.x/2.0)?size.x:0;
				d.y+=(d.y<-size.y/2.0)?size.y:0;
				d.z+=(d.z<-size.z/2.0)?size.z:0;
				if(d.x*d.x+d.y*d.y+d.z*d.z<cutoffSqr)
					exclude=true;
			}
			if(!exclude)
			{
				threeVector<double> minImg=0;
				minImg.x+=(p[i].x-p[0].x>size.x/2.0)?size.x:0;
				minImg.y+=(p[i].y-p[0].y>size.y/2.0)?size.y:0;
				minImg.z+=(p[i].z-p[0].z>size.z/2.0)?size.z:0;
				minImg.x+=(p[i].x-p[0].x<-size.x/2.0)?-size.x:0;
				minImg.y+=(p[i].y-p[0].y<-size.y/2.0)?-size.y:0;
				minImg.z+=(p[i].z-p[0].z<-size.z/2.0)?-size.z:0;
				com.x+=p[i].x+minImg.x;
				com.y+=p[i].y+minImg.y;
				com.z+=p[i].z+minImg.z;
				iCount++;
			}
		}
		com.x/=iCount;
		com.y/=iCount;
		com.z/=iCount;
		while(com.x>size.x)com.x-=size.x;
		while(com.y>size.y)com.y-=size.y;
		while(com.z>size.z)com.z-=size.z;
		while(com.x<0)com.x+=size.x;
		while(com.y<0)com.y+=size.y;
		while(com.z<0)com.z+=size.z;
		
		std::cerr << time << std::endl;
		std::cout << time << ' ';
		for(auto& i:nanoI)
		{
			threeVector<double> d=0;
			d.x=p[i].x-com.x;
			d.y=p[i].y-com.y;
			d.z=p[i].z-com.z;
			d.x+=(d.x>size.x/2.0)?-size.x:0;
			d.y+=(d.y>size.y/2.0)?-size.y:0;
			d.z+=(d.z>size.z/2.0)?-size.z:0;
			d.x+=(d.x<-size.x/2.0)?size.x:0;
			d.y+=(d.y<-size.y/2.0)?size.y:0;
			d.z+=(d.z<-size.z/2.0)?size.z:0;
			std::cout << sqrt(d.x*d.x+d.y*d.y+d.z*d.z) << ' ';
		}
		std::cout << std::endl;
		time+=System.readStoreInterval();
		frame++;
	}
	time_t end;
	std::time(&end);
	std::cerr << "Loads at " << static_cast<double>(frame)/static_cast<double>(end-start) << " frames per second ";
	std::cerr << " in " << end-start << " seconds!\n";
	//close our old files
	xyzFile.close();
	
	return 0;
}

