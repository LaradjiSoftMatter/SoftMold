/** \brief Takes an mpd file 
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
		std::cout << "Takes frame and state data to generate a new system.\n";
		std::cout << "usage: " << argv[0] << " name frameNumber newName\n";
		return 0;
	}
	
	//read in options
	char *name=argv[1];
	int frameNumber;
	char *newName=argv[3];
	std::stringstream cmdInp;
	cmdInp << argv[2];
	cmdInp >> frameNumber;
	
	//Configuration variables
	Blob<double> System;
	
	//load old variables, then initialize them
	Script<double, Blob <double> > fileIO(name,std::ios::in,&System);
	fileIO.read();
	fileIO.close();
	
	//our current size
	threeVector<double> size=System.readSize();
	
	//check if size varies
	std::string sizeName("size_");
	sizeName+=name;
	sizeName+=".dat";
	
	std::fstream sizeFile;
	sizeFile.open(sizeName.c_str(), std::ios::in);
	
	//Get our frames file
	std::string framesName("frames_");
	framesName+=name;
	framesName+=".xyz";
	
	position<double> *p=System.getPositions();
	int nParticles=System.readNParticles();
	
	xyzFormat<double> xyzFile(p,nParticles);
	xyzFile.open(framesName.c_str(), std::ios::in);
	
	double time=0;
	int frame=0;
	while(xyzFile.load())
	{
		//make sure the size hasn't changed
		if(sizeFile.is_open())
		{
			double sTime=0;
			threeVector<double> sCurrent;
			//for the first iteration, I expect this to exit the loop after one read,
			// otherwise, this is a template for subsequent reads in the loop that is
			// coming up.
			while(sTime<time && !sizeFile.eof())
				sizeFile >> sTime >> sCurrent.x >> sCurrent.y >> sCurrent.z;
			size=sCurrent;
		}
		System.setSize(size);
		
		if(frame==frameNumber)
		{
			fileIO.open(newName,std::ios::out);
			fileIO.write();
			fileIO.close();
			
			break;
		}
		
		time+=System.readStoreInterval();
		frame++;
	}
	
	if(frame<frameNumber)
		std::cout << "Error(main): Frame " << frameNumber << " not found, there are only " << frame << " frames!\n";
	
	if(sizeFile.is_open())
		sizeFile.close();
	
	return 0;
}

