/** \brief Loads each frame of a frames file and outputs frames in a new frames file.
 */

#include "../include/MD.h"
#include "../include/system.h"

int main(int argc, char* argv[])
{
	//check that options exist
	if(argc!=2)
	{
		//so simple, this can't possibly mess it up
		std::cerr << "Calculates radius of vesicle and frames per second loading\n";
		std::cerr << "\tusage: " << argv[0] << " name >> radius.dat\n";
		return 0;
	}
	
	//read in options
	char *name=argv[1];
	
	//Configuration variables
	Blob<double> System;
	
	//load old variables, then initialize them
	Script<double, Blob <double> > fileIO(name,std::ios::in,&System);
	fileIO.read();
	fileIO.close();
	
	//Get our frames file
	std::string frameName("frames_");
	frameName+=name;
	frameName+=".xyz";
	position<double> *p=System.getPositions();
	int nParticles=System.readNParticles();
	
	//This one is now handling allocation of the 'p' pointer.
	xyzFormat<double> xyzFile(p,nParticles);
	xyzFile.open(frameName.c_str(), std::ios::in);
	
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
		// otherwise, this is a template for subsequent reads in the loop that are
		// coming up.
		while(sTime<0 && !sizeFile.eof())
			sizeFile >> sTime >> sCurrent.x >> sCurrent.y >> sCurrent.z;
		time=sTime;
		size=sCurrent;
		//std::cerr << size.x << '\t' << size.y << '\t' << size.z << '\n';
	}
	
	//define the type we want to use for lipid head groups
	int lipidHead=2;
	
	//start reading frames
	time_t start=std::time(nullptr);
	int frame=0;
	while(xyzFile.load())
	{
		//Do stuff here, for this example, we will get the radius of a vesicle
		//initialize the center of mass to zero
		threeVector<double> com=0;
		
		//the number of head particles used to calculate center of mass
		double nLipidHeads=0;
		
		//go through all particles in system
		for(int i=0;i<nParticles;i++)
		{
			//choose only lipid head particles
			if(p[i].type==lipidHead)
			{
				//minimum image is an adjustment based on particle
				//distances.
				threeVector<double> minImg=0;
				
				//for minimum image:
				// We want to make sure the current particle is near
				// the first particle in the system. If it is farther
				// than half the size along each direction, we move 
				// it by one system distance closer to the first 
				// particle. All particles have x, y, z, and type
				// components.
				threeVector<double> d;
				d.x=p[i].x-p[0].x;
				d.y=p[i].y-p[0].y;
				d.z=p[i].z-p[0].z;
				
				if(d.x>size.x/2.0) minImg.x=size.x;
				if(d.y>size.y/2.0) minImg.y=size.y;
				if(d.z>size.z/2.0) minImg.z=size.z;
				
				if(d.x<-size.x/2.0) minImg.x=-size.x;
				if(d.y<-size.y/2.0) minImg.y=-size.y;
				if(d.z<-size.z/2.0) minImg.z=-size.z;
				
				//Now we sum the current position vector with 
				//the total calculated center of mass vector.
				com.x+=p[i].x-minImg.x;
				com.y+=p[i].y-minImg.y;
				com.z+=p[i].z-minImg.z;
				
				//increment n
				nLipidHeads++;
			}
		}
		
		//get the average center of mass vector:
		com.x/=nLipidHeads;
		com.y/=nLipidHeads;
		com.z/=nLipidHeads;
		
		//Now we want to the radius relative to the center of mass
		double radius=0;
		
		//go through all particles in system once more to obtain radius
		for(int i=0;i<nParticles;i++)
		{
			//choose only lipid head particles
			if(p[i].type==lipidHead)
			{
				//minimum image is an adjustment based on particle
				//distances.
				threeVector<double> minImg=0;
				
				//for minimum image:
				// We want to make sure the current particle is near
				// the center of mass in the system. If it is farther
				// than half the size along each direction, we move 
				// it by one system distance closer to the center of 
				// mass. All particles have x, y, z, and type
				// components.
				threeVector<double> d;
				d.x=p[i].x-com.x;
				d.y=p[i].y-com.y;
				d.z=p[i].z-com.z;
				
				if(d.x>size.x/2.0) minImg.x=size.x;
				if(d.y>size.y/2.0) minImg.y=size.y;
				if(d.z>size.z/2.0) minImg.z=size.z;
				
				if(d.x<-size.x/2.0) minImg.x=-size.x;
				if(d.y<-size.y/2.0) minImg.y=-size.y;
				if(d.z<-size.z/2.0) minImg.z=-size.z;
				
				//Now we recenter the current radial vector with 
				//the minimum image.
				d.x+=minImg.x;
				d.y+=minImg.y;
				d.z+=minImg.z;
				
				radius+=sqrt(d.x*d.x+d.y*d.y+d.z*d.z);
				
				//increment n, we already know what this is
				//nLipidHeads++;
			}
		}
		
		//average radius
		radius/=nLipidHeads;
		
		//output to "stdout"
		std::cout << time << ' ' << radius << std::endl;
		
		//output to "stderr"
		std::cerr << "Read frame " << frame << " with time " << time << " tau." << std::endl;
		
		//update time for next frame
		time+=System.readStoreInterval();
		
		//update size for next frame
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
		
		//update frame counter
		frame++;
	}
	time_t end=std::time(nullptr);
	std::cerr << "Loads at " << static_cast<double>(frame)/static_cast<double>(end-start) << " frames per second ";
	std::cerr << " in " << end-start << " seconds!\n";
	//close our old files
	xyzFile.close();
	
	return 0;
}

