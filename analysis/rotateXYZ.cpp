#include "../include/fileFormats/xyzFormat.h"

int main(int argc, char **argv)
{
	if(argc!=5)
	{
		std::cout << "Usage: " << argv[0] << " file.xyz theta.x theta.y theta.z\n";
		return 0;
	}
	
	//Read in command line options:
	std::string xyzName(argv[1]);
	
	std::stringstream cmdInput;
	cmdInput << argv[2] << ' ' << argv[3] << ' ' << argv[4];
	threeVector<double> theta;
	cmdInput >> theta.x >> theta.y >> theta.z;
	
	//We need some containers
	position<double> *p=NULL;
	int nParticles=0;
	
	//Read in xyz file:
	xyzFormat<double> xyzFile(p, nParticles, xyzName, std::ios::in);
	xyzFile.load();
	xyzFile.close();
	
	//Locate center of mass:
	threeVector<double> com;
	com=0;
	
	for(int i=0;i<nParticles;i++)
	{
		com.x+=p[i].x;
		com.y+=p[i].y;
		com.z+=p[i].z;
	}
	
	if(nParticles!=0)
	{
		com.x/=static_cast<double>(nParticles);
		com.y/=static_cast<double>(nParticles);
		com.z/=static_cast<double>(nParticles);
	}
	
	//Perform rotation (one axis at a time, you could combine them into one operation):
	for(int i=0;i<nParticles;i++)
	{
		//Adjust system to place center of mass at the origin:
		position<double> newP;
		p[i].x-=com.x;
		p[i].y-=com.y;
		p[i].z-=com.z;
		
		//Each rotation must be represented seperately, unless you utilize a different
		// representation. Truthfully you only need 2 axis for the rotation, the third
		// one is just for the hell of it.
		
		//Rotate with respect to the x axis:
		newP.y=p[i].y*cos(theta.x)-p[i].z*sin(theta.x);
		newP.z=p[i].y*sin(theta.x)+p[i].z*cos(theta.x);
		p[i].y=newP.y;
		p[i].z=newP.z;
		
		//Rotate with respect to the y axis:
		newP.x=p[i].x*cos(theta.y)+p[i].z*sin(theta.y);
		newP.z=-p[i].x*sin(theta.y)+p[i].z*cos(theta.y);
		p[i].x=newP.x;
		p[i].z=newP.z;
		
		//Rotate with respect to the z axis:
		newP.x=p[i].x*cos(theta.z)-p[i].y*sin(theta.z);
		newP.y=p[i].x*sin(theta.z)+p[i].y*cos(theta.z);
		p[i].x=newP.x;
		p[i].y=newP.y;
		
		//A better version:
		//Rotate with respect to the x and y axis:
		//newP.x=p[i].x*cos(theta.y)+p[i].z*sin(theta.y);
		//p[i].y=p[i].y*cos(theta.x)-p[i].z*sin(theta.x);
		//newP.z=p[i].y*sin(theta.x)-p[i].x*sin(theta.y)+p[i].z*(cos(theta.x)+cos(theta.y));
		//p[i].x=newP.x;
		//p[i].z=newP.z;
		
		//Move center of mass back to original position:
		p[i].x+=com.x;
		p[i].y+=com.y;
		p[i].z+=com.z;
	}
	
	//Output the rotation:
	xyzFile.open(xyzName, std::ios::out | std::ios::app);
	xyzFile.store();
	xyzFile.close();
	
	return 0;
}