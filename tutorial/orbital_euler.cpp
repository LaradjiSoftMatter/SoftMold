//Orbital motion of 2 bodies. Euler integration.

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
	if(argc!=13)
	{
		cout << "Usage: " << argv[0] << " p0.x p0.y p1.x p1.y v0.x v0.y v1.x v1.y cD G dT fT\n";
		return 0;
	}

	//initial positions
	twoVector p[2];
	p[0].x=atof(argv[1]);
	p[0].y=atof(argv[2]);
	p[1].x=atof(argv[3]);
	p[1].y=atof(argv[4]);

	//initial velocities
	twoVector v[2];
	v[0].x=atof(argv[5]);
	v[0].y=atof(argv[6]);
	v[1].x=atof(argv[7]);
	v[1].y=atof(argv[8]);

	//friction coefficient
	double cD=atof(argv[9]);
	
	//gravitational constant, note that this is for a central force,
	// units are different from the projectile motion problem
	double G=atof(argv[10]);
	
	//time step
	double dT=atof(argv[11]);

	//final time
	double fT=atof(argv[12]);

	twoVector d;
	d.x=p[0].x-p[1].x;
	d.y=p[0].y-p[1].y;
	double r=sqrt(d.x*d.x+d.y*d.y);
	//initial positions and accelerations are 0
	twoVector a[2];
	a[0].x=-d.x*G/(r*r*r);
	a[0].y=-d.y*G/(r*r*r);
	a[1].x=d.x*G/(r*r*r);
	a[1].y=d.y*G/(r*r*r);

	//get the number of steps until completion
	int nSteps=fT/dT;

	//make 2 files containing the coordinates
	fstream output0, output1;

	output0.open("p_euler0.dat", ios::out);
	output1.open("p_euler1.dat", ios::out);

	if(!output0.is_open() || !output1.is_open())
	{
		cout << "Cannot open p_euler0.dat or p_euler1.dat!\n";
		output0.close();
		output1.close();
		return 0;
	}

	//do an euler approximation of the equations of motion
	//stop when the time runs out
	for(int i=0;i<nSteps;i++)
	{
		//we need the displacement vector between the masses
		d.x=p[0].x-p[1].x;
		d.y=p[0].y-p[1].y;
		r=sqrt(d.x*d.x+d.y*d.y);
		
		//This is a little difficult to handle, we need a repulsive force to prevent this
		//But other methods could be used to handle it as well:
		// Globs (masses stick and dissipate energy)
		// Bounce (hard core repulsion)
		// Soft core repulsion
		if(r==0)
			break;

		//gravitational force
		//mass is assumed to be 1 for each mass
		a[0].x=-d.x*G/(r*r*r);
		a[0].y=-d.y*G/(r*r*r);
		a[1].x=d.x*G/(r*r*r);
		a[1].y=d.y*G/(r*r*r);

		//coefficient of drag (cD)
		//If cD is non zero, these should fall into each other, and will cause a division by 0
		a[0].x+=v[0].x*v[0].x*cD;
		a[0].y+=v[0].y*v[0].y*cD;
		a[1].x+=v[1].x*v[1].x*cD;
		a[1].y+=v[1].y*v[1].y*cD;

		//velocity changes by a small increment of acceleration
		v[0].x+=a[0].x*dT;
		v[0].y+=a[0].y*dT;
		v[1].x+=a[1].x*dT;
		v[1].y+=a[1].y*dT;

		//position changes by a small increment of velocity
		p[0].x+=v[0].x*dT;
		p[0].y+=v[0].y*dT;
		p[1].x+=v[1].x*dT;
		p[1].y+=v[1].y*dT;

		//output positions and velocities
		output0 << p[0].x << '\t' << p[0].y << '\n';
		output1 << p[1].x << '\t' << p[1].y << '\n';

		//Some information about the kinetic and potential energy
		double kE=(v[0].x*v[0].x+v[0].y*v[0].y)/2.0;
		kE+=(v[1].x*v[1].x+v[1].y*v[1].y)/2.0;
		//This is the gravitational potential, note the negative sign for attractive potentials
		double pE=-G/r;//both masses contribute
		cout << (double)i*dT << '\t' << kE << '\t' << pE << '\t' << kE+pE << endl;
	}

	output0.close();
	output1.close();

	return 0;
}
