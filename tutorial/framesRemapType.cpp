/** \brief Loads each frame of a frames file and outputs frames in a new frames file.
 */

#include "../include/systemMD.h"

int main(int argc, char* argv[])
{
	if(argc!=6 && argc!=5)
	{
		//so simple, this can't possibly mess it up
		std::cout << "Loads each frame and outputs a new modified set.\n";
		std::cout << "\tusage: " << argv[0] << " oldFramesName.xyz newFramesName.xyz\n"
		return 0;
	}
	
	//read in options
	char *name=argv[1];
	char *newName=argv[2];
	
	//Get our frames file
	position<double> *p=NULL;
	int nParticles=0;
	
	//This one is now handling allocation of the 'p' pointer.
	xyzFormat<double> xyzFile(p,nParticles);
	xyzFile.open(name, std::ios::in);
	
	//Make sure there is something in the file
	if(!xyzFile.load())
	{
		std::cout << "Error(main): Frames file is empty!\n Exiting...\n";
		return 0;
	}
	if(p==NULL)
	{
		std::cout << "Error(main): p is NULL!\n Exiting...\n";
		return 0;
	}

	//Our new frames file, this one only reads p without changing it.
	xyzFormat<double> newXyzFile(p,nParticles);
	newXyzFile.open(newName, std::ios::out);
	
	while(xyzFile.load())
	{
		//Do stuff with frames here

		//save our new frame here
		newXyzFile.store();
	}

	//close our old files
	xyzFile.close();
	newXyzFile.close();
	
	return 0;
}

