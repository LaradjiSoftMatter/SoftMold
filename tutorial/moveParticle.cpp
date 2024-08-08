//Loading an MD simulation "using systemMD.h".
//Everything needed to do a simulation is already included in systemMD.h.
//Even standard libraries like iostream and fstream are inlcuded.
//But, if any extra libraries are needed you can include them here or in systemMD.h.
//Permenent libraries (included with standard library or with this package) should
// go in systemMD.h, and anything you need that is external should go in this file.

#include "../include/systemMD.h"

int main(int argc, char **argv)
{
	if(argc!=3)
	{
		std::cout << argv[0] << " name newName\n";
		return 0;
	}
	//the variables for the simulation, call with empty constructor
	Blob<double> System;
	
	//file input-output script engine
	//name part (argv[1] in this case) is used with no extension (*.mpd)
	Script<double, Blob <double> > fileIO(argv[1],std::ios::in,&System);
	
	//read the variables in
	fileIO.read();
	
	//close the file that it is reading from
	fileIO.close();
	
	//determine where the center of mass for particles with type 5
	int typeCount=0;
	threeVector<double> centerOfMass(0);
	for(int i=0;i<System.readNParticles();i++)
	{
		if(System.getPositions()[i].type==5)
		{
			typeCount++;
			centerOfMass.x+=System.getPositions()[i].x;
			centerOfMass.y+=System.getPositions()[i].y;
			centerOfMass.z+=System.getPositions()[i].z;
		}
	}
	centerOfMass/=(double)typeCount;
	
	//reposition it on the boundary
	for(int i=0;i<System.readNParticles();i++)
	{
		if(System.getPositions()[i].type==5)
		{
			position<double> temp=System.getPositions()[i];
			//temp.x+=(System.readSize().x-centerOfMass.x);
			temp.y+=(System.readSize().y-centerOfMass.y);
			//temp.z+=(System.readSize().z-centerOfMass.z);
			System.getPositions()[i].x=(temp.x>=System.readSize().x)?temp.x-System.readSize().x:temp.x;
			System.getPositions()[i].y=(temp.y>=System.readSize().y)?temp.y-System.readSize().y:temp.y;
			System.getPositions()[i].z=(temp.z>=System.readSize().z)?temp.z-System.readSize().z:temp.z;
		}
	}

	fileIO.open(argv[2],std::ios::out);
	
	fileIO.write();
	
	fileIO.close();
	
	return 0;
}