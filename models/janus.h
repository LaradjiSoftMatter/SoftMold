//Modified from Yu's Janus particle code.

#include "../include/algorithms/tesselaSphere.h"
#ifndef JANUS_MPD
#define JANUS_MPD

template <typename T>
void janus(Blob<T> &System, T vBond, T cBond, T kbend, T particlebond, T ratio, threeVector<T> pos, T particleradius, int tesselation)
{
	//some lists
	std::vector<position<T>> vertices;
	std::vector<std::vector<fourVector<int>>> edges;
	std::vector<T> lengths; //since the edge have different lengths,creat a list to save them
	std::vector<std::vector<T>> dMap; //adjanct map that save the distance between any two edges
	std::vector<T> cos_O;//different cos_O
	std::vector<std::vector<fourVector<int>>> classfiedcos_O;
	
	//std::vector<fourVector<int>> faces;
	T Radius;
	
	//Creat Sphere
	position<T> sphereCenter; //should always set the center at 0, because the relocation does not work right now
	sphereCenter.x = 0;
	sphereCenter.y = 0;
	sphereCenter.z = 0;
	
	//create the first Sphere
	tesselaSphere<T> FirstSphere(sphereCenter, particleradius, tesselation);
	FirstSphere.setcenterpointtype(5);
	FirstSphere.settype1(4); //adhesion part
	FirstSphere.settype2(5);
	FirstSphere.JanusRatioSet(ratio);
	
	vertices = FirstSphere.ReturnVertex();
	//Do not forget the center point of the sphere
	vertices.push_back(FirstSphere.ReturnCenterPoint());	
	lengths = FirstSphere.ReturnEdgeLengths();
	edges.resize(lengths.size());
	dMap = FirstSphere.ReturndistanceMap();
	edges = FirstSphere.ReturnClassfiedEdges();
	Radius = FirstSphere.ReturnRadius();
	cos_O = FirstSphere.Returncos_Os();
	classfiedcos_O = FirstSphere.ReturnClassfiedangleIndexList();

	std::cout<< "vertice number: "<<vertices.size()-1 << std::endl;
	
	//the offset for inserting particles, needed to remap molecules from model
	int newParticleOffset=System.readNParticles();
	
	//get the new number of particles
	int newNParticles=System.readNParticles()+vertices.size();
	System.allocParticle(newNParticles);
	
	//get the center of the system
	threeVector<T> center = System.readSize()/2;
	
	std::cout << "system center " << center.x << " " << center.y <<  " " << center.z << std::endl;
	
	//Locate center of mass, we are placing at pos
	position<T> comValue=com< position<T> >(&vertices[0], vertices.size());
	
	//Velocities
	MTRand *randNum=new MTRand(System.readSeed());
	T Vrms=sqrt(3.0*System.readInitialTemp());
	
	threeVector<T> zero=0;
	//put vertices in system, repositioning center of mass
	for(int i=0;i<vertices.size();i++)
	{
		vertices[i].x+=pos.x-comValue.x;
		vertices[i].y+=pos.y-comValue.y;
		vertices[i].z+=pos.z-comValue.z;
		threeVector<T> vel;
		T theta=M_PI*randNum->rand53();
		T phi=M_PI*2*randNum->rand53();
		vel.x=Vrms*cos(phi)*sin(theta);
		vel.y=Vrms*sin(phi)*sin(theta);
		vel.z=Vrms*cos(theta);
		
		System.addParticle(vertices[i],vel,zero);
	}
	
	int Nbefore = System.readNTypes();
	
	if(Nbefore<8)
		System.setNTypes(8);
	//set enough type of particles
	
	int Nafter = System.readNTypes();
	
	std::cout << Nbefore << " " << Nafter << std::endl;
	
	std::cout << "nMolecules:" << System.readNMolecules() << std::endl;
	std::cout << "lengths:" << lengths.size() << std::endl;
	//add the edge bonds into the system for the first system
	for(int i = 0; i < lengths.size();i++)
	{
		molecule<T,fourVector<int>> mBonds;
		mBonds.setType(BOND);
		mBonds.addConstant(lengths[i]);
		mBonds.addConstant(vBond);
		for(int j=0;j<edges[i].size();j++)
		{
			edges[i][j].s[0]+=newParticleOffset;
			edges[i][j].s[1]+=newParticleOffset;
			mBonds.addBond(edges[i][j]);
		}
		System.addMolecule(mBonds);
	}
	std::cout << "nMolecules:" << System.readNMolecules() << std::endl;
	
	//add the edge bonds between vertices and center point into the system for the first particle 
	molecule<T,fourVector<int>> cBonds;
	cBonds.setType(BOND);
	cBonds.addConstant(Radius);
	cBonds.addConstant(cBond);
	for(int i = 0; i < vertices.size()-1; i++)
	{
		fourVector<int> tempedge;
		tempedge.s[0] = i + newParticleOffset;
		tempedge.s[1] = vertices.size()-1 + newParticleOffset;
		cBonds.addBond(tempedge);
	}
	System.addMolecule(cBonds);
	
	std::cout << "nMolecules:" << System.readNMolecules() << std::endl;
	int firstParticleIndex = vertices.size() - 1 + newParticleOffset;

	//add the bend force for the first particle
	for(int i = 0; i < cos_O.size(); i++)
	{
		molecule<T,fourVector<int>> Bends;
		Bends.setType(BEND);
		Bends.addConstant(-cos_O[i]);
		Bends.addConstant(kbend);
		for(int j = 0; j< classfiedcos_O[i].size(); j++)
		{
			classfiedcos_O[i][j].s[0] += newParticleOffset;
			classfiedcos_O[i][j].s[1] += newParticleOffset;
			classfiedcos_O[i][j].s[2] += newParticleOffset;
			Bends.addBond(classfiedcos_O[i][j]);
		}
		System.addMolecule(Bends);
	}
	
	std::cout << "nMolecules:" << System.readNMolecules() << std::endl;
}

#endif
