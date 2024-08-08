/** \brief Loads each frame of a frames file and outputs frames in a new frames file.
 */

#include "../include/MD.h"
#include "../include/system.h"

int main(int argc, char* argv[])
{
	if(argc!=2)
	{
		//so simple, this can't possibly mess it up
		std::cerr << "\tUsage: " << argv[0] << " name\n";
		std::cerr << "Calculates frames per second loading an xyz file\n";
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
	int nParticles=0;
	
	//This one is now handling allocation of the 'p' pointer.
	xyzFormat<double> xyzFile(p,nParticles);
	xyzFile.open(frameName.c_str(), std::ios::in);
	
	double time=0;
	
	//start reading frames
	time_t start=std::time(nullptr);
	int frame=0;
	while(xyzFile.load())
	{
		//Do stuff here
		
		
		
		std::cerr << "Read frame " << frame << " with time " << time << " tau." << std::endl;
		
		//update time for next frame
		time+=System.readStoreInterval();
		
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

