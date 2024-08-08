#include "mpdCuda.h"

#include <fstream>
#include <iostream>
#include <random>
#include <omp.h>

int main()
{
	threeVector<float> size(800.0,800.0,50.0);
	int nLipids=size.x*size.y*3.11;
	int lGrid=sqrt(nLipids/2.0);
	float dlGrid=size.x/(lGrid+1);
	
	uint nParticles=3*nLipids,nTypes=4;
	float gamma=1.0, temperature=3.0, deltaT=0.02;
	
	mpd::randomAdaptor random(1234);
	std::vector<position<float>> p;
	std::vector<threeVector<float>> a,v;
	mpd::interaction<float,3,2> bendList(nParticles);
	mpd::interaction<float,2,2> bondList(nParticles);
	mpd::dataCollection<float> dataCollection(nParticles);
	
	//Bending and bonding geometry
	std::vector<int> chains={0,nLipids,3};//{START,NCHAINS,CHAINLENGTH}
	std::vector<float> constants={0.7,100.0,1,100.0};//{aBond,kBond,-cos(t0),kBend}
	
	
	//initialize constants (force and potential constants)
	//generate force constants
	std::vector<float> Umax(nTypes*nTypes);
	std::vector<float> Umin(nTypes*nTypes);
	
	//could be one loop, but it's good for an example
	//of how to set constants
	for(int i=0;i<nTypes;i++)
	{
		for(int j=0;j<nTypes;j++)
		{
			//i=first type, j=second type, indexed grid
			int k=i+j*nTypes;
			Umin[k]=0;
			Umax[k]=100;
		}
	}
	
	Umax[3*nTypes+3]=200;
	Umin[3*nTypes+3]=-6.0;
	
	std::vector<float> NBFconstants=mpd::laradjiRevaleeFC(Umax, Umin, 2.0f, 1.0f);
	std::vector<float> NBPconstants=mpd::laradjiRevaleePC(Umax, Umin, 2.0f, 1.0f);
	
	
	//Initialize bending geometry
	for(int i=chains[0];i<chains[1]*chains[2];i+=chains[2])
	{
		std::vector<int> bend;
		for(int j=i;j<i+chains[2];j++)
			bend.push_back(j);
		for(int j=0;j<chains[2]-2;j++)
		{
			bendList.addInteraction(bend[j],0,  &(bend[j]),&(constants[2]));
			bendList.addInteraction(bend[j+1],1,&(bend[j]),&(constants[2]));
			bendList.addInteraction(bend[j+2],2,&(bend[j]),&(constants[2]));
		}
	}
	
	//Initialize bonding
	for(int i=chains[0];i<chains[1]*chains[2];i+=chains[2])
	{
		std::vector<int> bond;
		for(int j=i;j<i+chains[2];j++)
			bond.push_back(j);
		for(int j=0;j<chains[2]-1;j++)
		{
			bondList.addInteraction(bond[j],0,  &(bond[j]),&(constants[0]));
			bondList.addInteraction(bond[j+1],1,&(bond[j]),&(constants[0]));
		}
	}
	
	//Initialize positions
	for(float x=0;x<size.x;x+=dlGrid)
	{
		for(float y=0;y<size.y;y+=dlGrid)
		{
			p.emplace_back(x,y,size.z/2.0-1.7,2);
			p.emplace_back(x,y,size.z/2.0-1.0,3);
			p.emplace_back(x,y,size.z/2.0-0.35,3);
			p.emplace_back(x,y,size.z/2.0+1.7,2);
			p.emplace_back(x,y,size.z/2.0+1.0,3);
			p.emplace_back(x,y,size.z/2.0+0.35,3);
		}
	}
	
	//Initialize velocities and accelerations
	for(int i=0;i<nParticles;i++)
	{
		//threeVector<float> vel(float(i%3)-1.0,float((i+5)%3)-1.0,float((i+7)%3)-1.0);
		threeVector<float> vel(0.0f,0.0f,0.0f);
		threeVector<float> acc(0.0f,0.0f,0.0f);
		v.push_back(vel);
		a.push_back(acc);
	}
	
	//Initialize state
	mpd::state<float> state(p.data(),v.data(),a.data(),NBFconstants.data(),NBPconstants.data(),
					nParticles,nTypes,deltaT,gamma,temperature,size);
	
	//Initialize grid
	mpd::cell<float> cData(nParticles, 2.0f, size);
	
	//Send to device
	bondList.toDevice();
	bendList.toDevice();
	state.toDevice();
	
	//Perform molecular dynamics algorithm
	mpd::bondForces_device(bondList.deviceInteraction(),state.deviceState());
	mpd::bendForces_device(bendList.deviceInteraction(),state.deviceState());
	
	//mpd::testCount.resize(nParticles,0);
	mpd::cellComputeForce_device(cData.deviceCell(), state.deviceState());
		
	{
		std::fstream xyzFile("out.xyz",std::ios::out);
		xyzFile << nParticles << "\ntest\n";
		for(int i=0;i<nParticles;i++)
			xyzFile << p[i].type << ' ' << p[i].x << ' ' << p[i].y << ' ' << p[i].z << std::endl;
		
		mpd::zeroKinetic_device(dataCollection.deviceState());
		mpd::kinetic_device(state.deviceState(),dataCollection.deviceState());
		float kinetic=mpd::reduceKinetic_device(dataCollection.deviceState());
		
		mpd::zeroPotential_device(dataCollection.deviceState());
		mpd::bondPotential_device(bondList.deviceInteraction(),state.deviceState(),
				dataCollection.deviceState());
		mpd::bendPotential_device(bendList.deviceInteraction(),state.deviceState(),
				dataCollection.deviceState());
		mpd::cellComputePotential_device(cData.deviceCell(), state.deviceState(),
				dataCollection.deviceState());
		float potential=mpd::reducePotential_device(dataCollection.deviceState());
		
		std::fstream energyFile("energy.dat", std::ios::out);
		energyFile << 0.0 << ' ' << potential << ' ' << kinetic << std::endl;
	}
	
	float begin=0.0,end=100.0;
	int nSteps=std::floor((end-begin)/deltaT);
	for(int step=0;step<nSteps+1;step++)
	{
		mpd::verletFirst_device(state.deviceState());
		mpd::zeroAccelerations_device(state.deviceState());
		mpd::langevin_device(state.deviceState(),random.deviceState());
		mpd::bondForces_device(bondList.deviceInteraction(),state.deviceState());
		mpd::bendForces_device(bendList.deviceInteraction(),state.deviceState());
		mpd::cellComputeForce_device(cData.deviceCell(), state.deviceState());
		mpd::verletSecond_device(state.deviceState());
		
		if(step%5000==0 && step!=0)
		{
			state.toHost();
			std::fstream xyzFile("out.xyz",std::ios::out | std::ios::app);
			xyzFile << nParticles << "\ntest\n";
			for(int i=0;i<nParticles;i++)
				xyzFile << p[i].type << ' ' << p[i].x << ' ' << p[i].y << ' ' << p[i].z << std::endl;
			
			mpd::zeroKinetic_device(dataCollection.deviceState());
			mpd::kinetic_device(state.deviceState(),dataCollection.deviceState());
			float kinetic=mpd::reduceKinetic_device(dataCollection.deviceState());
			
			mpd::zeroPotential_device(dataCollection.deviceState());
			mpd::bondPotential_device(bondList.deviceInteraction(),state.deviceState(),
					dataCollection.deviceState());
			mpd::bendPotential_device(bendList.deviceInteraction(),state.deviceState(),
					dataCollection.deviceState());
			mpd::cellComputePotential_device(cData.deviceCell(), state.deviceState(),
					dataCollection.deviceState());
			float potential=mpd::reducePotential_device(dataCollection.deviceState());
			
			std::fstream energyFile("energy.dat", std::ios::out | std::ios::app);
			energyFile << deltaT*step << ' ' << potential << ' ' << kinetic << std::endl;
		}
	}
	
	//dump data
	
	return 0;
}
	

	/*#pragma omp parallel for
	for(int i=0;i<nParticles;i++)
	{
		float cutoffSqr=2.0*2.0;
		for(int j=0;j<nParticles;j++)
		{
			if(i!=j)
			{
				threeVector<float> d;
				d.x=p[i].x-p[j].x;
				d.y=p[i].y-p[j].y;
				d.z=p[i].z-p[j].z;
				if(d.x>=size.x/2.0)d.x-=size.x;
				if(d.y>=size.y/2.0)d.y-=size.y;
				if(d.z>=size.z/2.0)d.z-=size.z;
				if(d.x<=-size.x/2.0)d.x+=size.x;
				if(d.y<=-size.y/2.0)d.y+=size.y;
				if(d.z<=-size.z/2.0)d.z+=size.z;
				threeVector<float> f= mpd::nonBondedF(d, cutoffSqr, NBconstants.data(),
						     p[i].type, p[j].type, nTypes);
				a[i].x+=f.x;
				a[i].y+=f.y;
				a[i].z+=f.z;
			}
		}
	}*/
