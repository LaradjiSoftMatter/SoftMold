#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <vector>

struct point {
	double x,y;
};

struct parameters {
	double k, a;
	double temperature, volume, energy;
	double timeStep;
	int numberParticles;
};

//vector between points
point diff(point a, point b)
{
	point d;
	d.x=a.x-b.x;
	d.y=a.y-b.y;
	return d;
}

//magnitude of vector
double mag(point d)
{
	return sqrt(d.x*d.x+d.y*d.y);
}

//distance
double dist(point a, point b)
{
	return mag(diff(a,b));
}

//energy of a spring
double potential(point a, point b, parameters p)
{
	double d=dist(a,b);
	double E=0.5*p.k*(d-p.a)*(d-p.a);
	return E;
}

//force of spring
point force(point a, point b, parameters p)
{
	//vector between two points
	point D=diff(a,b);
	
	//distance between two points
	double d=mag(D);
	
	//magnitude of F=-grad(U)
	point F;
	F.x=-(d-p.a)*D.x/d;
	F.y=-(d-p.a)*D.y/d;
	return F;
}

//main program
int main()
{
	//initialize parameters
	parameters param;
	param.k=100;
	param.a=10;
	param.timeStep=0.00001;
	param.numberParticles=2;
	
	//build force, velocities, positions, and masses
	std::vector<point> forces, vel, pos;
	std::vector<double> masses;
	
	//this is pretty common
	point zero;
	zero.x=0;
	zero.y=0;
	
	for(int i=0;i<param.numberParticles;i++)
	{
		//initialize acceleration and velocity to zero
		//probably want to initialize velocities to something else
		forces.push_back(zero);
		vel.push_back(zero);
		//random position on a lattice
		point position;
		position.x=rand()%10;
		position.y=rand()%10;
		pos.push_back(position);
		
		//masses are just one
		masses.push_back(1.0);
	}
	
	//output for an xyz file
	std::fstream xyzFile("out.xyz", std::ios::out);
	
	//Do molecular dynamics here
	for(int step=0;step<1000000;step++)
	{
		//Calculate forces between particles i and j
		for(int i=0;i<param.numberParticles;i++)
		{
			for(int j=i+1;j<param.numberParticles;j++)
			{
				//if(calculate force?)
				//{
					//magnitude of the force
					point f=force(pos[i],pos[j],param);
					//apply to particle i
					forces[i].x+=f.x;
					forces[i].y+=f.y;
					//newton's third law for particle j
					forces[j].x-=f.x;
					forces[j].y-=f.y;
				//}
			}
		}
		
		//Update positions and velocities on particle i, "Euler integration"
		for(int i=0;i<param.numberParticles;i++)
		{
			vel[i].x+=(forces[i].x/masses[i])*param.timeStep;
			vel[i].y+=(forces[i].y/masses[i])*param.timeStep;
			
			pos[i].x+=vel[i].x*param.timeStep;
			pos[i].y+=vel[i].y*param.timeStep;
			
			forces[i]=zero;
		}
		
		//output particles' positions every 100th step
		if(step%100==0)
		{
			xyzFile << param.numberParticles << "\ntest\n";//vmd xyz header
			//vmd usually wants the "atomicnumber x y z", we set atomicnumber to 1 and z to 0
			for(int i=0;i<param.numberParticles;i++)
				xyzFile << "1\t" << pos[i].x << '\t' << pos[i].y << "\t0" << std::endl;
		}
	}
	
	return 0;
}
//hwen@memphis.edu	
