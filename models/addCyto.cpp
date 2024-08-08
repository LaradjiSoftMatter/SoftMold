//This is a very short and simple program that shows you how to modify a variable already present in an mpd script.
//The reason this is used is because some variables (forces, potentials, geometry, etc...) are not easily modified
// in the mpd script. Using a file that opens the mpd script, modifies, and then saves it (as below) is far easier
// than modifying the actual variable. Most of it looks like the run.cpp file, it just replaces the "run loop"
// with the modification.

#include "../include/systemMD.h"
#include <ctime>
#include <cstdlib>
#include <cstring>
#define THRESHOLD 1.0

//searches for midpoints in a walk around a polyhedra, derived from someone elses code for generating tesselated polyhedra
int search_midpoint(position<double> *vertices, threeVector<int> *edge, int &edge_walk, int &nVertices, int index_start, int index_end) 
{ 
	///Note these values for the edges
	//edge[i].x=START;
	//edge[i].y=MIDDLE;
	//edge[i].z=END;
	
	
	for (int i=0; i<edge_walk; i++) 
		if ((edge[i].x == index_start && edge[i].z == index_end) || 
			(edge[i].x == index_end && edge[i].z == index_start)) 
		{
			int res = edge[i].y;
			
			/* update the arrays */
			edge[i].x=edge[edge_walk-1].x;
			edge[i].z=edge[edge_walk-1].z;
			edge[i].y=edge[edge_walk-1].y;
			edge_walk--;
			return res; 
		}

	/* vertex not in the list, so we add it */
	edge[edge_walk].x = index_start;
	edge[edge_walk].z = index_end; 
	edge[edge_walk].y = nVertices; 
	
	/* create new vertex */ 
	vertices[nVertices].x = (vertices[index_start].x + vertices[index_end].x) / 2.0;
	vertices[nVertices].y = (vertices[index_start].y + vertices[index_end].y) / 2.0;
	vertices[nVertices].z = (vertices[index_start].z + vertices[index_end].z) / 2.0;
	
	/* normalize the new vertex */ 
	double length = sqrt (vertices[nVertices].x * vertices[nVertices].x +
				vertices[nVertices].y * vertices[nVertices].y +
				vertices[nVertices].z * vertices[nVertices].z);
	length = 1/length;
	vertices[nVertices].x *= length;
	vertices[nVertices].y *= length;
	vertices[nVertices].z *= length;
	
	nVertices++;
	edge_walk++;
	
	return edge[edge_walk-1].y;
}

int main(int argc, char* argv[])
{
	if(argc!=8)
	{
		//so simple, this can't possibly mess it up
		std::cout << "usage: command name nMonomers posX posY posZ radius nTess\n";
		return 0;
	}
	
	threeVector<double> pos,norm;
	char *name=argv[1];
	int nMonomers=atoi(argv[2]);
	pos.x=atoi(argv[3]);
	pos.y=atoi(argv[4]);
	pos.z=atoi(argv[5]);
	double radius=atof(argv[6]);
	int nTess=atoi(argv[7]);
	
	///Previous Configuration variables
	Blob<double> System;
	
	///load old variables, then initialize them
	Script<double, Blob <double> > fileIO(name,std::ios::in,&System);
	
	fileIO.read();
	
	fileIO.close();
	
	std::string newName("");
	newName+=name;
	newName+="_wC";
	
	///This initializes the new variables
	//adjust size of system if it is out of bounds
	threeVector<double> size=System.readSize();
	if(size.x==0 && size.y==0 && size.z==0)
	{
		size.x=pos.x;
		size.y=pos.y;
		size.z=pos.z;
		pos.x/=2.0;
		pos.y/=2.0;
		pos.z/=2.0;
	}
	else
	{
		size.x=size.x<radius+pos.x+5.0?radius+pos.x:size.x;
		size.y=size.y<radius+pos.y+5.0?radius+pos.y:size.y;
		size.z=size.z<radius+pos.z+5.0?radius+pos.z:size.z;
	}
	//set new size, if at all
	System.setSize(size);
	
	//tesselate an icosahedron
	int nVertices = 12, nEdges=30, nFaces = 20;
	
	//these are normalized constants
	double t = (1+(double)sqrt(5.0))/2;
	double tau = t/(double)sqrt(1.0+t*t);
	double one = 1/(double)sqrt(1.0+t*t);
	
	position <double> *anc = (position <double> *)malloc(nVertices*sizeof(position<double>));
	if(anc==NULL)
	{
		std::cout << "Not enough memory!\n";
		return 0;
	}
	
	//vertices of an icosahedron
	anc[0].x=tau;
	anc[0].y=one;
	anc[0].z=0.0;
	anc[1].x=-tau;
	anc[1].y=one;
	anc[1].z=0.0;
	anc[2].x=-tau;
	anc[2].y=-one;
	anc[2].z=0.0;
	anc[3].x=tau;
	anc[3].y=-one;
	anc[3].z=0.0;
	anc[4].x=one;
	anc[4].y=0.0;
	anc[4].z=tau;
	anc[5].x=one;
	anc[5].y=0.0;
	anc[5].z=-tau;
	anc[6].x=-one;
	anc[6].y=0.0;
	anc[6].z=-tau;
	anc[7].x=-one;
	anc[7].y=0.0;
	anc[7].z=tau;
	anc[8].x=0.0;
	anc[8].y=tau;
	anc[8].z=one;
	anc[9].x=0.0;
	anc[9].y=-tau;
	anc[9].z=one;
	anc[10].x=0.0;
	anc[10].y=-tau;
	anc[10].z=-one;
	anc[11].x=0.0;
	anc[11].y=tau;
	anc[11].z=-one;
	
	threeVector<int> *faces = (threeVector<int>*)malloc(nFaces*sizeof(threeVector<int>));
	if(faces==NULL)
	{
		std::cout << "Not enough memory!\n";
		return 0;
	}
	
	faces[0].s[0]=4;
	faces[0].s[1]=8;
	faces[0].s[2]=7;
	faces[1].s[0]=4;
	faces[1].s[1]=7;
	faces[1].s[2]=9;
	faces[2].s[0]=5;
	faces[2].s[1]=6;
	faces[2].s[2]=11;	
	faces[3].s[0]=5;
	faces[3].s[1]=10;
	faces[3].s[2]=6;
	faces[4].s[0]=0;
	faces[4].s[1]=4;
	faces[4].s[2]=3;
	faces[5].s[0]=0;
	faces[5].s[1]=3;
	faces[5].s[2]=5;
	faces[6].s[0]=2;
	faces[6].s[1]=7;
	faces[6].s[2]=1;
	faces[7].s[0]=2;
	faces[7].s[1]=1;
	faces[7].s[2]=6;
	faces[8].s[0]=8;
	faces[8].s[1]=0;
	faces[8].s[2]=11;
	faces[9].s[0]=8;
	faces[9].s[1]=11;
	faces[9].s[2]=1;
	faces[10].s[0]=9;
	faces[10].s[1]=10;
	faces[10].s[2]=3;
	faces[11].s[0]=9;
	faces[11].s[1]=2;
	faces[11].s[2]=10;
	faces[12].s[0]=8;
	faces[12].s[1]=4;
	faces[12].s[2]=0;
	faces[13].s[0]=11;
	faces[13].s[1]=0;
	faces[13].s[2]=5;
	faces[14].s[0]=4;
	faces[14].s[1]=9;
	faces[14].s[2]=3;
	faces[15].s[0]=5;
	faces[15].s[1]=3;
	faces[15].s[2]=10;
	faces[16].s[0]=7;
	faces[16].s[1]=8;
	faces[16].s[2]=1;
	faces[17].s[0]=6;
	faces[17].s[1]=1;
	faces[17].s[2]=11;
	faces[18].s[0]=7;
	faces[18].s[1]=2;
	faces[18].s[2]=9;
	faces[19].s[0]=6;
	faces[19].s[1]=10;
	faces[19].s[2]=2;
	
	for (int i=0; i<nTess; i++)
	{
		int n_a_new = nVertices+2*nEdges; 
		int n_faces_new = 4*nFaces; 
		
		int edge_walk = 0; 
		nEdges = 2*nVertices + 3*nFaces; 
		threeVector<int> *edge = new threeVector<int>[nEdges]; 
		
		for(int j=0;j<nEdges;j++)
		{
			edge[j].s[0]=-1;
			edge[j].s[1]=-1;
			edge[j].s[2]=-1;
		}
		
		threeVector<int> *faces_old=new threeVector<int>[nFaces]; 
		for(int j=0;j<nFaces;j++)
			faces_old[j]=faces[j];
		anc=(position<double>*)realloc ((void*)anc, n_a_new*sizeof(position<double>));
		if(anc==NULL)
		{
			std::cout << "Not enough memory!\n";
			return 0;
		}
		
		faces=(threeVector<int>*)realloc ((void*)faces, n_faces_new*sizeof(threeVector<int>));
		if(faces==NULL)
		{
			std::cout << "Not enough memory!\n";
			return 0;
		}
		
		n_faces_new=0; 
		
		for (int j=0; j<nFaces; j++) 
		{ 
			int xa = faces_old[j].s[0]; 
			int xb = faces_old[j].s[1]; 
			int xc = faces_old[j].s[2]; 
			
			int ab_midpoint = search_midpoint (anc, edge, edge_walk, nVertices, xb, xa); 
			int bc_midpoint = search_midpoint (anc, edge, edge_walk, nVertices, xc, xb); 
			int ca_midpoint = search_midpoint (anc, edge, edge_walk, nVertices, xa, xc); 
			
			faces[n_faces_new].s[0] = xa; 
			faces[n_faces_new].s[1] = ab_midpoint; 
			faces[n_faces_new].s[2] = ca_midpoint; 
			n_faces_new++;
			
			faces[n_faces_new].s[0] = ca_midpoint; 
			faces[n_faces_new].s[1] = ab_midpoint; 
			faces[n_faces_new].s[2] = bc_midpoint; 
			n_faces_new++;
			
			faces[n_faces_new].s[0] = ca_midpoint; 
			faces[n_faces_new].s[1] = bc_midpoint; 
			faces[n_faces_new].s[2] = xc;
			n_faces_new++;
			
			faces[n_faces_new].s[0] = ab_midpoint; 
			faces[n_faces_new].s[1] = xb; 
			faces[n_faces_new].s[2] = bc_midpoint; 
			n_faces_new++; 
		}
		
		nFaces = n_faces_new;
		
		delete faces_old;
		delete edge;
	}
	
	std::cout << "nFaces: " << nFaces << " nVertices: " << nVertices << '\n';
	
	//expand to radius-1.0 (just below radius of a liposome), assumes it is normalized
	
	for(int i=0;i<nVertices;i++)
	{
		anc[i].x*=(radius-1.0);
		anc[i].y*=(radius-1.0);
		anc[i].z*=(radius-1.0);
	} 
	
	///Adding new molecules to system
	fourVector<int> bond;
	molecule<double,fourVector<int> > m;
	
	m.setType(CHAIN);
	
	//constants for CHAIN
	m.addConstant(0.7);
	m.addConstant(100);
	m.addConstant(1.0);
	m.addConstant(100);
	
	//this is used for the anchors as well
	int startOffset=System.readNParticles();
	
	//bond for a chain type
	bond.s[START]=startOffset+nVertices;
	bond.s[NCHAINS]=(nFaces+nVertices)-2;
	bond.s[LENGTH]=nMonomers;
	m.addBond(bond);
	System.addMolecule(m);
	
	m.~molecule();
	
	///This is important to reduce excessive reallocations
	System.allocParticle(System.readNParticles()+nVertices+nMonomers*((nFaces+nVertices)-2));
	
	//Velocities
	MTRand *randNum=new MTRand(System.readSeed());
	double Vrms=sqrt(3*System.readInitialTemp());
	
	//put vertices in system
	for(int i=0;i<nVertices;i++)
	{
		position<double> p;
		threeVector<double> v,a;
		a.x=0;
		a.y=0;
		a.z=0;
		p.x=anc[i].x+pos.x;
		p.y=anc[i].y+pos.y;
		p.z=anc[i].z+pos.z;
		p.type=ANCHOR;
		//velocity
		double theta=M_PI*randNum->rand53();
		double phi=M_PI*2*randNum->rand53();
		v.x=Vrms*cos(phi)*sin(theta);
		v.y=Vrms*sin(phi)*sin(theta);
		v.z=Vrms*cos(theta);
		
		System.addParticle(p,v,a);
	}
	
	m.setType(BOND);
	
	//constants for BOND
	m.addConstant(0.7);
	m.addConstant(100);
	
	//Just allocate more than you expect
	m.allocBonds(nVertices*6);
	
	
	
	//place monomers between vertices
	for(int i=0;i<nFaces;i++)
	{
		//indices to points on the face
		int xa=faces[i].s[0];
		int xb=faces[i].s[1];
		int xc=faces[i].s[2];
		
		//checking for edge's existance on previous faces
		bool ab=true;
		bool bc=true;
		bool ca=true;
		
		//null vector
		threeVector<double> a;
		a.x=0;
		a.y=0;
		a.z=0;
		
		//check if edge pairs exist on other faces
		for(int j=0;j<i;j++)
		{
			int ya=faces[j].s[0];
			int yb=faces[j].s[1];
			int yc=faces[j].s[2];

			if((xa==ya || xa==yb || xa==yc) && (xb==ya || xb==yb || xb==yc))
				ab=false;
			if((xb==ya || xb==yb || xb==yc) && (xc==ya || xc==yb || xc==yc))
				bc=false;
			if((xc==ya || xc==yb || xc==yc) && (xa==ya || xa==yb || xa==yc))
				ca=false;
		}
		
		//place monomers between vertices if edge isn't a duplicate
		if(ab)
		{
			for(int j=1;j<nMonomers+1;j++)
			{
				position<double> p;
				threeVector<double> v;
				double phi,theta;
				
				//position
				p.x=(double)j*((anc[xa].x-anc[xb].x)/((double)nMonomers+1.0))+anc[xb].x+pos.x;
				p.y=(double)j*((anc[xa].y-anc[xb].y)/((double)nMonomers+1.0))+anc[xb].y+pos.y;
				p.z=(double)j*((anc[xa].z-anc[xb].z)/((double)nMonomers+1.0))+anc[xb].z+pos.z;
				p.type=j==1||j==nMonomers?CYTO:MONOMER;
				
				//velocity
				theta=M_PI*randNum->rand53();
				phi=M_PI*2*randNum->rand53();
				v.x=Vrms*cos(phi)*sin(theta);
				v.y=Vrms*sin(phi)*sin(theta);
				v.z=Vrms*cos(theta);
		
				System.addParticle(p,v,a);
			}
			
			if(nMonomers>2)
			{
				//bond for a bond type
				bond.s[0]=startOffset+xa;
				bond.s[1]=System.readNParticles()-1;
				m.addBond(bond);
				
				//bond for a bond type
				bond.s[0]=startOffset+xb;
				bond.s[1]=System.readNParticles()-nMonomers;
				m.addBond(bond);
			}
		}
		
		if(bc)
		{
			for(int j=1;j<nMonomers+1;j++)
			{
				position<double> p;
				threeVector<double> v;
				double phi,theta;
				
				//position
				p.x=(double)j*((anc[xb].x-anc[xc].x)/((double)nMonomers+1.0))+anc[xc].x+pos.x;
				p.y=(double)j*((anc[xb].y-anc[xc].y)/((double)nMonomers+1.0))+anc[xc].y+pos.y;
				p.z=(double)j*((anc[xb].z-anc[xc].z)/((double)nMonomers+1.0))+anc[xc].z+pos.z;
				p.type=j==1||j==nMonomers?CYTO:MONOMER;
				
				//velocity
				theta=M_PI*randNum->rand53();
				phi=M_PI*2*randNum->rand53();
				v.x=Vrms*cos(phi)*sin(theta);
				v.y=Vrms*sin(phi)*sin(theta);
				v.z=Vrms*cos(theta);
				
				System.addParticle(p,v,a);
			}
			if(nMonomers>2)
			{
				//bond for a bond type
				bond.s[0]=startOffset+xb;
				bond.s[1]=System.readNParticles()-1;
				m.addBond(bond);
				
				//bond for a bond type
				bond.s[0]=startOffset+xc;
				bond.s[1]=System.readNParticles()-nMonomers;
				m.addBond(bond);
			}
		}
		
		if(ca)
		{
			for(int j=1;j<nMonomers+1;j++)
			{
				position<double> p;
				threeVector<double> v;
				double phi,theta;
				
				//position
				p.x=(double)j*((anc[xc].x-anc[xa].x)/((double)nMonomers+1.0))+anc[xa].x+pos.x;
				p.y=(double)j*((anc[xc].y-anc[xa].y)/((double)nMonomers+1.0))+anc[xa].y+pos.y;
				p.z=(double)j*((anc[xc].z-anc[xa].z)/((double)nMonomers+1.0))+anc[xa].z+pos.z;
				p.type=j==1||j==nMonomers?CYTO:MONOMER;
				
				//velocity
				theta=M_PI*randNum->rand53();
				phi=M_PI*2*randNum->rand53();
				v.x=Vrms*cos(phi)*sin(theta);
				v.y=Vrms*sin(phi)*sin(theta);
				v.z=Vrms*cos(theta);
				
				System.addParticle(p,v,a);
			}
			if(nMonomers>2)
			{
				//bond for a bond type
				bond.s[0]=startOffset+xc;
				bond.s[1]=System.readNParticles()-1;
				m.addBond(bond);
				
				//bond for a bond type
				bond.s[0]=startOffset+xa;
				bond.s[1]=System.readNParticles()-nMonomers;
				m.addBond(bond);
			}
		}
	}
	
	//add the vertex bonds to the system
	System.addMolecule(m);
	m.~molecule();
	
	///This stores the new configuration
	Script<double, Blob <double> > fileIO_new(newName.c_str(),std::ios::out,&System);
	fileIO_new.write();
	fileIO_new.close();
	
	//This stores the xyz file to check new configuration
	position<double> *p=System.getPositions();
	int nParticles=System.readNParticles();
	
	xyzFormat<double> xyzFile(p, nParticles);
	newName+=".xyz";
	xyzFile.open(newName.c_str(), std::ios::out);
	xyzFile.store();
	xyzFile.close();
	
	free(faces);
	free(anc);
	return 0;
}

