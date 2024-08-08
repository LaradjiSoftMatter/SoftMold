/** \brief Program that checks the diffusivity of a single particle.
 *  This program seems a little strange, but when the langevin thermostat is applied to a particle,
 *  there are little bumps in it's direction. This program outputs the diffusivity of a single 
 *  particle by repeating a diffusivity test on system of one particle.
 *
 *  The theory:
 *  A langevin thermostat utilizes the following relation F=-Fc-gamma*p
 *
 */

#include "../include/MD.h"
#include "../include/system.h"

using namespace std;

#define NHIST 1000


int main(int argc, char **argv)
{
	if(argc!=7)
	{
		std::cout << "Usage: " << argv[0] << " temperature deltaT totalTime gamma seedDisplacement nParticles\n";
		std::cout << "Output: temp gamma diffusionConstant\n";
		return 0;
	}
	
	//Read in command args
	std::stringstream cmdArg;
	cmdArg << argv[1] << ' ' << argv[2] << ' ' << argv[3] << ' ' << argv[4] << ' ' << argv[5] << ' ' << argv[6];
	double temp, deltaT, totalTime,gamma;
	int seedDisplacement,nParticles;
	cmdArg >> temp >> deltaT >> totalTime >> gamma >> seedDisplacement >> nParticles;
	
	//Set up the simulation parameters
	position<double> p;
	threeVector<double> a, v;
	int nSteps=totalTime/deltaT;
	threeVector<bool> wrap=false;
	threeVector<double> size=1000;
	Verlet<double> integrate(&p, &a, &v, 1, size, deltaT, wrap);
	Langevin<double> thermostat(&a, &v, 1, gamma, deltaT, time(NULL)+seedDisplacement);
	MTRand randNum(time(NULL)+seedDisplacement);
	
	//Initialize averages
	double maxD=0;
	vector<double> avgMS;
	for(int i=0;i<nSteps;i++)
		avgMS.push_back(0);
	

	//Perform simulation accross NPART particles
	for(int n=0;n<nParticles;n++)
	{
		//Initialize parameters
		p.x=0;
		p.y=0;
		p.z=0;
		p.type=1;
		a=0;
		double Vrms=sqrt(3.0*temp);
		double theta=M_PI*randNum.rand53();
		double phi=M_PI*2.0*randNum.rand53();
		v.x=Vrms*cos(phi)*sin(theta);
		v.y=Vrms*sin(phi)*sin(theta);
		v.z=Vrms*cos(theta);
		
		//Step through a simulation
		for(int i=0;i<nSteps;i++)
		{
			thermostat.compute(temp);
			integrate.first();
			a=0;
			//a.y=-0.01;//falling ball

			//Accumulation for ensemble average
			double ms=p.x*p.x+p.y*p.y+p.z*p.z;
			if(sqrt(ms)>maxD)
				maxD=sqrt(ms);
			avgMS[i]+=(ms/static_cast<double>(nParticles));
			
			//std::cout << i*deltaT << '\t' << ms << '\n';

			integrate.second();
		}
	}
	
	double xAvg=totalTime/2.0;
	double yAvg=0;
	
	//Perform ensemble average
	for(int i=0;i<nSteps;i++)
	{
		//std::cout << i*deltaT << '\t' << avgMS[i] << '\n';
		yAvg+=avgMS[i];
	}
	yAvg/=static_cast<double>(nSteps);
	
	
	//Regress the slope: y=yIntercept+dAvg*x
	double Sxy=0, Sxx=0;
	for(int i=0;i<nSteps;i++)
	{
		double dx=static_cast<double>(i)*deltaT-xAvg;
		double dy=avgMS[i]-yAvg;
		Sxy+=dx*dy;
		Sxx+=dx*dx;
	}
	double dAvg=Sxy/Sxx;
	double yIntercept=yAvg-dAvg*xAvg;
	
	//Find the standard error
	double Syy=0;
	for(int i=0;i<nSteps;i++)
	{
		double x=static_cast<double>(i)*deltaT;
		double y=dAvg*x+yIntercept;
		double dy=avgMS[i]-y;
		Syy+=(dy*dy);
	}
	double stdError=sqrt(Syy/static_cast<double>(nSteps))/totalTime;

	//Output the diffusion coefficient
	std::cout << temp << '\t' << gamma << '\t' << (dAvg/6.0) << '\t' << stdError << '\n';
	
	
	//std::cout << temp << '\t' << gamma << '\t' << (Sxy/(Sxx*6.0)) << '\n';
	//paths.push_back(path);

	/*
	int *nD=new int[NHIST];
	for(int i=0;i<NHIST;i++)
		nD[i]=0;
	
	double dS=maxD/NHIST;
	
	for(int i=0;i<paths.size();i++)
	{
		for(int j=0;j<paths[i].size();j++)
		{
			threeVector<double> d;
			d.x=paths[i][j].x;
			d.y=paths[i][j].y;
			d.z=paths[i][j].z;
			double dr=sqrt(d.x*d.x+d.y*d.y+d.z*d.z);
			int pos=int(dr/dS);
			pos=(pos==NHIST)?NHIST-1:pos;
			nD[pos]++;
		}
	}
	
	for(int i=0;i<NHIST;i++)
		cout << i*dS << '\t' << nD[i] << '\n';
	delete nD;
	*/
	return 0;
}
