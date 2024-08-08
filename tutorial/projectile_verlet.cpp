//Projectile motion of single body. Velocity-Verlet integration.
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

	//initial positions are 0
	twoVector p;
	p.x=0;
	p.y=0;
	
	//initial accelerations are computed
	twoVector a;

	//accelerations are updated
	a.x=-v.x*v.x*cF;
	a.y=-g-v.y*v.y*cF;

	//get the number of steps until completion
	int nSteps=fT/dT;

	//do a verlet approximation of the equations of motion
	//stop when we hit the ground
	for(int i=0;i<nSteps && p.y>=0;i++)
	{
		//store half steps
		twoVector vHalf,aHalf;
		vHalf=v;
		aHalf=a;
		
		//velocity changes by a small increment of acceleration
		v.x+=a.x*dT*0.5;
		v.y+=a.y*dT*0.5;

		//position changes by a small increment of velocity
		p.x+=v.x*dT+a.x*dT*dT*0.5;
		p.y+=v.y*dT+a.y*dT*dT*0.5;

		
		//accelerations are updated
		a.x=-v.x*v.x*cF;
		a.y=-g-v.y*v.y*cF;
		
		//half step velocity update
		v.x=vHalf.x+(aHalf.x+a.x)*dT*0.5;
		v.y=vHalf.y+(aHalf.y+a.y)*dT*0.5;
		
		//output positions and velocities
		cout << p.x << '\t' << p.y << '\t';
		cout << v.x << '\t' << v.y << '\n';
	}

	return 0;
}
