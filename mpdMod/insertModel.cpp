//This is a very short and simple program that shows you how to modify a variable already present in an mpd script.
//The reason this is used is because some variables (forces, potentials, geometry, etc...) are not easily modified
// in the mpd script. Using a file that opens the mpd script, modifies, and then saves it (as below) is far easier
// than modifying the actual variable. Most of it looks like the run.cpp file, it just replaces the "run loop"
// with the modification.

#include "../include/MD.h"
#include "../include/system.h"

int main(int argc, char* argv[])
{
	if(argc!=10)
	{
		//so simple, this can't possibly mess it up
		std::cout << "Usage: " << argv[0] << " oldName newName vertices.dat edges.dat kBond lBond pos.x pos.y pos.z\n";
		return 0;
	}
	
	
	//read in options to parameters
	char *oldName=argv[1];
	char *newName=argv[2];
	char *vName=argv[3];
	char *eName=argv[4];
	//char *fName=argv[5];
	
	threeVector<double> pos;
	double kBond, lBond;
	//double kBend, cosThetaBend;
	std::stringstream cmdIn;
	for(int i=5;i<argc;i++)
		cmdIn << argv[i] << ' ';
	cmdIn >> kBond >> lBond >> pos.x >> pos.y >> pos.z;
	
	//some lists
	std::vector<position<double>> vertices;
	std::vector<fourVector<int>> edges;
	//std::vector<fourVector<int>> faces;
	
	std::fstream vIn(vName,std::ios::in);
	if(vIn.is_open())
	{
		position<double> vertex;
		while(vIn >> vertex.type >> vertex.x >> vertex.y >> vertex.z)
			vertices.push_back(vertex);
	}
	else
	{
		std::cerr << "Cannot open " << vName << "!" << std::endl;
		return -1;
	}
	
	std::fstream eIn(eName,std::ios::in);
	if(eIn.is_open())
	{
		fourVector<int> edge;
		while(eIn >> edge.x >> edge.y)
			edges.push_back(edge);
	}
	else
	{
		std::cerr << "Cannot open " << eName << "!" << std::endl;
		return -1;
	}
	
	//std::fstream fIn(fName,std::ios::in);
	//if(fIn.is_open())
	//{
	//	fourVector<int> face;
	//	while(eIn >> face.x >> face.y >> face.z)
	//		faces.push_back(face);
	//}
	//else
	//{
	//	std::cerr << "Cannot open " << fName << "!" << std::endl;
	//	return -1;
	//}
	
	//Configuration variables
	Blob<double> System;
	
	//load old variables, then initialize them
	Script<double, Blob <double> > fileIO(oldName,std::ios::in,&System);
	fileIO.read();
	fileIO.close();
	
	//the offset for inserting particles, needed to remap molecules from model
	int newParticleOffset=System.readNParticles();
	
	//get the new number of particles
	int newNParticles=System.readNParticles()+vertices.size();
	System.allocParticle(newNParticles);
	
	//Locate center of mass, we are placing at pos
	position<double> comValue=com< position<double> >(&vertices[0], vertices.size());
	
	//Velocities
	MTRand *randNum=new MTRand(System.readSeed());
	double Vrms=sqrt(3.0*System.readInitialTemp());
	
	threeVector<double> zero=0;
	//put vertices in system, repositioning center of mass
	for(int i=0;i<vertices.size();i++)
	{
		vertices[i].x+=pos.x-comValue.x;
		vertices[i].y+=pos.y-comValue.y;
		vertices[i].z+=pos.z-comValue.z;
		threeVector<double> vel;
		double theta=M_PI*randNum->rand53();
		double phi=M_PI*2*randNum->rand53();
		vel.x=Vrms*cos(phi)*sin(theta);
		vel.y=Vrms*sin(phi)*sin(theta);
		vel.z=Vrms*cos(theta);
		
		System.addParticle(vertices[i],vel,zero);
	}
	
	//add the edge bonds into the system
	molecule<double,fourVector<int>> mBonds;
	mBonds.setType(BOND);
	mBonds.addConstant(lBond);
	mBonds.addConstant(kBond);
	for(int i=0;i<edges.size();i++)
	{
		edges[i].s[0]+=newParticleOffset;
		edges[i].s[1]+=newParticleOffset;
		mBonds.addBond(edges[i]);
	}
	
	System.addMolecule(mBonds);
	
	//molecule<double,fourVector<int>> mBends;
	//mBends.setType(BEND);
	//mBends.addConstant(cosThetaBend);
	//mBends.addConstant(kBend);
	//for(int i=0;i<faces.size();i++)
	//{
	//	faces[i].s[0]+=newParticleOffset;
	//	faces[i].s[1]+=newParticleOffset;
	//	faces[i].s[2]+=newParticleOffset;
	//	mBends.addBond(faces[i]);
	//}
	
	//Save system
	fileIO.open(newName,std::ios::out);
	fileIO.write();
	fileIO.close();
	
	return 0;
}

