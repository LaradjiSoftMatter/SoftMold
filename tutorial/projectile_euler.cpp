//Projectile motion of a single body. Euler integration.
//Drag term means motion isn't conservative. Total energy decreases as projectile moves.


#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>

using namespace std;

struct twoVector {
	double x,y;
};

int main(int argc, char **argv)
{
	if(argc!=7)
	{
		cout << "Usage: " << argv[0] << " v.x v.y cF g dT fT\n";
		return 0;
	}

	//initial velocity
	twoVector v;
	v.x=atof(argv[1]);
	v.y=atof(argv[2]);
	
	//friction coefficient
	double cF=atof(argv[3]);
	
	//gravitational constant
	double g=atof(argv[4]);
	
	//time step
	double dT=atof(argv[5]);

	//final time
	double fT=atof(argv[6]);

	//initial positions and accelerations are 0
	twoVector p,a;
	p.x=0;
	p.y=0;
	a.x=0;
	a.y=0;

	//get the number of steps until completion
	int nSteps=fT/dT;

	//do an euler approximation of the equations of motion
	//stop when we hit the ground
	for(int i=0;i<nSteps && p.y>=0;i++)
	{
		//we assume mass is 1
		a.x=-v.x*v.x*cF;
		a.y=-g-v.y*v.y*cF;

		//velocity changes by a small increment of acceleration
		v.x+=a.x*dT;
		v.y+=a.y*dT;

		//position changes by a small increment of velocity
		p.x+=v.x*dT;
		p.y+=v.y*dT;

		//output positions and velocities
		cout << p.x << '\t' << p.y << '\t' << v.x << '\t' << v.y << '\n';
	}

	return 0;
}
