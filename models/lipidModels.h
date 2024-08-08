/**
 * \brief Various model geometries associated with lipid systems.
 * We have flat and spherical cytoskeletons, solvent fills, and flat and spherical lipid bilayers. 
 * Many options for connectivity as well. Assumes your system blob has some functions to utilize 
 * with molecules, particles, and temperature.
 */

#include "../include/algorithms/functions.h"
#include "../include/algorithms/cellAlgorithm.h"

#ifndef R_THRESHOLD
#define R_THRESHOLD 1.0
#endif

/** \brief Searches for neighbors.
 * Generic search for neighbors.
 */
template <typename T>
class searchNeighbors {
	public:
		//constructor
		searchNeighbors(position<T> *particles, int nParticles, int offset, T threshold)
		{
			p=particles;
			nP=nParticles;
			this->threshold=threshold;
			this->offset=offset;
			myInterference=new bool[nP];
			for(int i=0;i<nP;i++)
				myInterference[i]=false;
		};
		
		//destructor
		~searchNeighbors()
		{
			if(myInterference!=NULL)
				delete myInterference;
		};
		
		//output the molecule
		bool interference(int i)
		{
			return myInterference[i];
		};
		
		//check if any neighbor fits the constraints
		void operator() (int &i, position<T> &nearby)
		{
			threeVector<T> d;
			
			d.x=p[i].x-nearby.x;
			d.y=p[i].y-nearby.y;
			d.z=p[i].z-nearby.z;
			
			if(sqrt(d.x*d.x+d.y*d.y+d.z*d.z)<threshold)
				interference[i]=true;
		};
		
	private:
		position<T> *p;
		int nP;
		int offset;
		T threshold;
		bool *myInterference;
};

//searches for midpoints in a walk around a polyhedra, derived from someone elses code for generating tesselated polyhedra
template <typename T>
int search_midpoint(position<T> *vertices, threeVector<int> *edge, int &edge_walk, int &nVertices, int index_start, int index_end) 
{ 
	//Note these values for the edges
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
	T length = sqrt (vertices[nVertices].x * vertices[nVertices].x +
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


template <typename T>
T sphericalCyto(Blob<T> &System, int nMonomers, threeVector<T> pos, T radius, int nTess, T *constants, int nConstants)
{
	//Error checking
	if(nConstants!=4)
	{
		std::cerr << "Error(sphericalCyto): nConstants is not 4!\n";
		return 0;
	}
	if(radius<=0)
	{
		std::cerr << "Error(sphericalCyto): No radius!\n";
		return 0;
	}
	if(nMonomers<3)
	{
		std::cerr << "Error(sphericalCyto): Not enough monomers!\n";
		return 0;
	}
	if(nTess<0)
	{
		std::cerr << "Error(sphericalCyto): No tesselation!\n";
		return 0;
	}
	
	//tesselate an icosahedron
	int nVertices = 12, nEdges=30, nFaces = 20;
	
	//these are normalized constants
	T t = (1+(T)sqrt(5.0))/2;
	T tau = t/(T)sqrt(1.0+t*t);
	T one = 1/(T)sqrt(1.0+t*t);
	
	position <T> *anc = (position <T> *)malloc(nVertices*sizeof(position<T>));
	if(anc==NULL)
	{
		std::cerr << "Not enough memory!\n";
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
		std::cerr << "Not enough memory!\n";
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
		anc=(position<T>*)realloc ((void*)anc, n_a_new*sizeof(position<T>));
		if(anc==NULL)
		{
			std::cerr << "Not enough memory!\n";
			return 0;
		}
		
		faces=(threeVector<int>*)realloc ((void*)faces, n_faces_new*sizeof(threeVector<int>));
		if(faces==NULL)
		{
			std::cerr << "Not enough memory!\n";
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
	
	std::cerr << "nFaces: " << nFaces << " nVertices: " << nVertices << '\n';
	
	//expand to radius-1.0 (just below radius of a liposome), assumes it is normalized
	
	for(int i=0;i<nVertices;i++)
	{
		anc[i].x*=(radius-constants[0]);
		anc[i].y*=(radius-constants[0]);
		anc[i].z*=(radius-constants[0]);
	} 
	
	//Adding new molecules to system
	fourVector<int> bond;
	molecule<T,fourVector<int> > chainMol;
	
	chainMol.setType(CHAIN);
	
	//constants for CHAIN
	chainMol.addConstant(constants[0]);
	chainMol.addConstant(constants[1]);
	chainMol.addConstant(constants[2]);
	chainMol.addConstant(constants[3]);
	
	//this is used for the anchors as well
	int startOffset=System.readNParticles();
	
	//bond for a chain type
	bond.s[START]=startOffset+nVertices;
	bond.s[NCHAINS]=(nFaces+nVertices)-2;
	bond.s[CHAINLENGTH]=nMonomers;
	chainMol.addBond(bond);
	System.addMolecule(chainMol);
	
	//This is important to reduce excessive reallocations
	System.allocParticle(System.readNParticles()+nVertices+nMonomers*((nFaces+nVertices)-2));
	
	//Velocities
	MTRand *randNum=new MTRand(System.readSeed());
	T Vrms=sqrt(3.0*System.readInitialTemp());
	
	molecule<T,fourVector<int> > anchorMol;
	
	//put vertices in system
	for(int i=0;i<nVertices;i++)
	{
		position<T> p;
		threeVector<T> v,a;
		a.x=0;
		a.y=0;
		a.z=0;
		p.x=anc[i].x+pos.x;
		p.y=anc[i].y+pos.y;
		p.z=anc[i].z+pos.z;
		p.type=ANCHOR;
		//velocity
		T theta=M_PI*randNum->rand53();
		T phi=M_PI*2*randNum->rand53();
		v.x=Vrms*cos(phi)*sin(theta);
		v.y=Vrms*sin(phi)*sin(theta);
		v.z=Vrms*cos(theta);
		
		System.addParticle(p,v,a);
	}
	
	anchorMol.setType(BOND);
	
	//constants for BOND
	anchorMol.addConstant(constants[0]);
	anchorMol.addConstant(constants[1]);
	
	//Just allocate more than you expect (12 of them have 5 each, but midpoints have 6)
	anchorMol.allocBonds(nVertices*6);
	
	
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
		threeVector<T> a;
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
				position<T> p;
				threeVector<T> v;
				T phi,theta;
				
				//position
				p.x=(T)j*((anc[xa].x-anc[xb].x)/((T)nMonomers+1.0))+anc[xb].x+pos.x;
				p.y=(T)j*((anc[xa].y-anc[xb].y)/((T)nMonomers+1.0))+anc[xb].y+pos.y;
				p.z=(T)j*((anc[xa].z-anc[xb].z)/((T)nMonomers+1.0))+anc[xb].z+pos.z;
				p.type=j==1||j==nMonomers?CYTO:MONOMER;
				
				//velocity
				theta=M_PI*randNum->rand53();
				phi=M_PI*2*randNum->rand53();
				v.x=Vrms*cos(phi)*sin(theta);
				v.y=Vrms*sin(phi)*sin(theta);
				v.z=Vrms*cos(theta);
		
				System.addParticle(p,v,a);
			}
			
			//bond for a bond type
			bond.s[0]=startOffset+xa;
			bond.s[1]=System.readNParticles()-1;
			anchorMol.addBond(bond);
			
			//bond for a bond type
			bond.s[0]=startOffset+xb;
			bond.s[1]=System.readNParticles()-nMonomers;
			anchorMol.addBond(bond);
		}
		
		if(bc)
		{
			for(int j=1;j<nMonomers+1;j++)
			{
				position<T> p;
				threeVector<T> v;
				T phi,theta;
				
				//position
				p.x=(T)j*((anc[xb].x-anc[xc].x)/((T)nMonomers+1.0))+anc[xc].x+pos.x;
				p.y=(T)j*((anc[xb].y-anc[xc].y)/((T)nMonomers+1.0))+anc[xc].y+pos.y;
				p.z=(T)j*((anc[xb].z-anc[xc].z)/((T)nMonomers+1.0))+anc[xc].z+pos.z;
				p.type=j==1||j==nMonomers?CYTO:MONOMER;
				
				//velocity
				theta=M_PI*randNum->rand53();
				phi=M_PI*2*randNum->rand53();
				v.x=Vrms*cos(phi)*sin(theta);
				v.y=Vrms*sin(phi)*sin(theta);
				v.z=Vrms*cos(theta);
				
				System.addParticle(p,v,a);
			}
			
			//bond for a bond type
			bond.s[0]=startOffset+xb;
			bond.s[1]=System.readNParticles()-1;
			anchorMol.addBond(bond);
			
			//bond for a bond type
			bond.s[0]=startOffset+xc;
			bond.s[1]=System.readNParticles()-nMonomers;
			anchorMol.addBond(bond);
		}
		
		if(ca)
		{
			for(int j=1;j<nMonomers+1;j++)
			{
				position<T> p;
				threeVector<T> v;
				T phi,theta;
				
				//position
				p.x=(T)j*((anc[xc].x-anc[xa].x)/((T)nMonomers+1.0))+anc[xa].x+pos.x;
				p.y=(T)j*((anc[xc].y-anc[xa].y)/((T)nMonomers+1.0))+anc[xa].y+pos.y;
				p.z=(T)j*((anc[xc].z-anc[xa].z)/((T)nMonomers+1.0))+anc[xa].z+pos.z;
				p.type=j==1||j==nMonomers?CYTO:MONOMER;
			
				//velocity
				theta=M_PI*randNum->rand53();
				phi=M_PI*2*randNum->rand53();
				v.x=Vrms*cos(phi)*sin(theta);
				v.y=Vrms*sin(phi)*sin(theta);
				v.z=Vrms*cos(theta);
				
				System.addParticle(p,v,a);
			}
			
			//bond for a bond type
			bond.s[0]=startOffset+xc;
			bond.s[1]=System.readNParticles()-1;
			anchorMol.addBond(bond);
			
			//bond for a bond type
			bond.s[0]=startOffset+xa;
			bond.s[1]=System.readNParticles()-nMonomers;
			anchorMol.addBond(bond);
		}
	}
	
	//add the vertex bonds to the system
	System.addMolecule(anchorMol);
	
	free(faces);
	free(anc);
	return radius;
}

//Same as above but with lipid anchors, must already have a bilayer
template <typename T>
T sphericalCyto(Blob<T> &System, int nMonomers, threeVector<T> pos, T radius, int nTess, T *constants, int nConstants, int lipidList)
{
	//Error checking
	if(nConstants!=4)
	{
		std::cerr << "Error(sphericalCyto): nConstants is not 4!\n";
		return 0;
	}
	if(radius<=0)
	{
		std::cerr << "Error(sphericalCyto): No radius!\n";
		return 0;
	}
	if(nMonomers<3)
	{
		std::cerr << "Error(sphericalCyto): Not enough monomers!\n";
		return 0;
	}
	if(nTess<0)
	{
		std::cerr << "Error(sphericalCyto): No tesselation!\n";
		return 0;
	}
	if(System.readNParticles()<=0)
	{
		std::cerr << "Error(sphericalCyto): No particles for search!\n";
		return 0;
	}
	if(System.readNMolecules()<=0 || System.readNMolecules()<=lipidList)
	{
		std::cerr << "Error(sphericalCyto): No lipidList present!\n";
		return 0;
	}
	if(System.getMolecule()[lipidList].readType()!=CHAIN)
	{
		std::cerr << "Error(sphericalCyto): Molecule at lipidList is not a CHAIN type!\n";
		return 0;
	}
	
	//tesselate an icosahedron
	int nVertices = 12, nEdges=30, nFaces = 20;
	
	//these are normalized constants
	T t = (1+(T)sqrt(5.0))/2;
	T tau = t/(T)sqrt(1.0+t*t);
	T one = 1/(T)sqrt(1.0+t*t);
	
	position <T> *anc = (position <T> *)malloc(nVertices*sizeof(position<T>));
	if(anc==NULL)
	{
		std::cerr << "Not enough memory!\n";
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
		std::cerr << "Not enough memory!\n";
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
		anc=(position<T>*)realloc ((void*)anc, n_a_new*sizeof(position<T>));
		if(anc==NULL)
		{
			std::cerr << "Not enough memory!\n";
			return 0;
		}
		
		faces=(threeVector<int>*)realloc ((void*)faces, n_faces_new*sizeof(threeVector<int>));
		if(faces==NULL)
		{
			std::cerr << "Not enough memory!\n";
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
	
	std::cerr << "nFaces: " << nFaces << " nVertices: " << nVertices << '\n';
	
	//expand to radius-1.0 (just below radius of a liposome), assumes it is normalized
	
	for(int i=0;i<nVertices;i++)
	{
		anc[i].x*=(radius-constants[0]);
		anc[i].y*=(radius-constants[0]);
		anc[i].z*=(radius-constants[0]);
	}
	
	//Adding new molecules to system
	fourVector<int> bond;
	molecule<T,fourVector<int> > chainMol;
	
	chainMol.setType(CHAIN);
	
	//constants for CHAIN
	chainMol.addConstant(constants[0]);
	chainMol.addConstant(constants[1]);
	chainMol.addConstant(constants[2]);
	chainMol.addConstant(constants[3]);
	
	//this is used for the anchors as well
	int startOffset=System.readNParticles();
	
	//bond for a chain type
	bond.s[START]=startOffset+nVertices;
	bond.s[NCHAINS]=(nFaces+nVertices)-2;
	bond.s[CHAINLENGTH]=nMonomers;
	chainMol.addBond(bond);
	System.addMolecule(chainMol);
	
	//This is important to reduce excessive reallocations
	System.allocParticle(System.readNParticles()+nVertices+nMonomers*((nFaces+nVertices)-2));
	
	//Velocities
	MTRand *randNum=new MTRand(System.readSeed());
	T Vrms=sqrt(3.0*System.readInitialTemp());
	
	//for transmembrane proteins
	int vertOffset=System.readNParticles();
	
	//put vertices in system
	for(int i=0;i<nVertices;i++)
	{
		position<T> p;
		threeVector<T> v,a;
		a.x=0;
		a.y=0;
		a.z=0;
		p.x=anc[i].x+pos.x;
		p.y=anc[i].y+pos.y;
		p.z=anc[i].z+pos.z;
		p.type=ANCHOR;
		//velocity
		T theta=M_PI*randNum->rand53();
		T phi=M_PI*2*randNum->rand53();
		v.x=Vrms*cos(phi)*sin(theta);
		v.y=Vrms*sin(phi)*sin(theta);
		v.z=Vrms*cos(theta);
		
		System.addParticle(p,v,a);
		
		//for(int j=0;j<System.getMolecule().readNBond();j++)
		//{
		//	bond=System.getMolecule().getBonds()[j];
		//}
	}
	
	//transmembrane bonding and bending
	molecule<T,fourVector<int> > transMemBond;
	molecule<T,fourVector<int> > transMemBend;
	
	transMemBond.setType(BOND);
	transMemBend.setType(BEND);
	
	//constants for transmembrane BOND
	transMemBond.addConstant(System.getMolecule()[lipidList].getConstants()[CHAINBOND]);
	transMemBond.addConstant(System.getMolecule()[lipidList].getConstants()[CHAINBOND+1]);
	
	//constants for transmembrane BEND
	transMemBend.addConstant(System.getMolecule()[lipidList].getConstants()[CHAINBEND]);
	transMemBend.addConstant(System.getMolecule()[lipidList].getConstants()[CHAINBEND+1]);
	
	for(int i=0;i<nVertices;i++)
	{
		threeVector<T> d;
		position<T> *p=System.getPositions();
		int start=System.getMolecule()[lipidList].getBonds()[0].s[START];
		int length=System.getMolecule()[lipidList].getBonds()[0].s[CHAINLENGTH];
		int nChains=System.getMolecule()[lipidList].getBonds()[0].s[NCHAINS];
		threeVector<T> size=System.readSize();
		int nearestHeadChain=-1;
		int nearestHeadMonomer=-1;
		double nearestHeadD=size.x*size.x+size.y*size.y+size.z*size.z;
		
		//Find the nearest head lipid
		for(int chain=0;chain<nChains;chain++)
		{
			for(int j=0;j<length;j++)
			{
				int index=chain*length+start+j;
				if(p[index].type==HEAD)
				{
					d.x=p[index].x-p[vertOffset+i].x;
					d.y=p[index].y-p[vertOffset+i].y;
					d.z=p[index].z-p[vertOffset+i].z;
					
					d.x-=(d.x>size.x/2.0)?size.x:0;
					d.x+=(d.x<-size.x/2.0)?size.x:0;
					
					d.y-=(d.y>size.y/2.0)?size.y:0;
					d.y+=(d.y<-size.y/2.0)?size.y:0;
					
					d.z-=(d.z>size.z/2.0)?size.z:0;
					d.z+=(d.z<-size.z/2.0)?size.z:0;
					
					double dr=d.x*d.x+d.y*d.y+d.z*d.z;
					
					if(dr<nearestHeadD)
					{
						nearestHeadD=dr;
						nearestHeadChain=chain;
						nearestHeadMonomer=j;
					}
				}
			}
		}
		
		//Set the new bonds and bends for the cytoskeleton membrane interface
		bond.s[0]=vertOffset+i;//our common anchor
		bond.s[1]=nearestHeadChain*length+start+nearestHeadMonomer;//our nearest head
		
		int tailEnd=0;
		
		//lexical orientation
		if(nearestHeadMonomer==0)//The monomer subsequently after our current monomer
		{
			bond.s[2]=bond.s[1]+1;
			tailEnd=bond.s[1]+length-1;
			//lets change our lipid now to an anchor type
			System.getPositions()[nearestHeadChain*length+start+nearestHeadMonomer]=HEAD_ANCHOR;
			for(int k=1;k<length;k++)
				System.getPositions()[nearestHeadChain*length+start+nearestHeadMonomer+k]=TAIL_ANCHOR;
			
		}
		else//or the monomer before our current monomer
		{
			bond.s[2]=bond.s[1]-1;
			tailEnd=bond.s[1]-length+1;
			//lets change our lipid now to an anchor type
			System.getPositions()[nearestHeadChain*length+start+nearestHeadMonomer]=HEAD_ANCHOR;
			for(int k=1;k<length;k++)
				System.getPositions()[nearestHeadChain*length+start+nearestHeadMonomer-k]=TAIL_ANCHOR;
		}
		transMemBond.addBond(bond);//We only worry about s[0] and s[1] in this case
		transMemBend.addBond(bond);//We already have information about s[2] from lexical orientation, completing our chain
		
		
		
		//Now locate an opposing lipid to add and convert convert
		for(int chain=0;chain<nChains;chain++)
		{
			threeVector<double> orient;
			for(int j=0;j<length;j++)
			{
				int index=chain*length+start+j;
				if(p[index].type==HEAD)
				{
					d.x=p[index].x-p[vertOffset+i].x;
					d.y=p[index].y-p[vertOffset+i].y;
					d.z=p[index].z-p[vertOffset+i].z;
					
					d.x-=(d.x>size.x/2.0)?size.x:0;
					d.x+=(d.x<-size.x/2.0)?size.x:0;
					
					d.y-=(d.y>size.y/2.0)?size.y:0;
					d.y+=(d.y<-size.y/2.0)?size.y:0;
					
					d.z-=(d.z>size.z/2.0)?size.z:0;
					d.z+=(d.z<-size.z/2.0)?size.z:0;
					
					double dr=d.x*d.x+d.y*d.y+d.z*d.z;
					
					if(dr<nearestHeadD)
					{
						nearestHeadD=dr;
						nearestHeadChain=chain;
						nearestHeadMonomer=j;
					}
				}
			}
		}
			
		
	}
	
	System.addMolecule(transMemBond);
	System.addMolecule(transMemBend);
	
	molecule<T,fourVector<int> > anchorMol;
	
	anchorMol.setType(BOND);
	
	//constants for BOND
	anchorMol.addConstant(constants[0]);
	anchorMol.addConstant(constants[1]);
	
	//Just allocate more than you expect (12 of them have 5 each, but midpoints have 6)
	anchorMol.allocBonds(nVertices*6);
	

	
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
		threeVector<T> a;
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
				position<T> p;
				threeVector<T> v;
				T phi,theta;
				
				//position
				p.x=(T)j*((anc[xa].x-anc[xb].x)/((T)nMonomers+1.0))+anc[xb].x+pos.x;
				p.y=(T)j*((anc[xa].y-anc[xb].y)/((T)nMonomers+1.0))+anc[xb].y+pos.y;
				p.z=(T)j*((anc[xa].z-anc[xb].z)/((T)nMonomers+1.0))+anc[xb].z+pos.z;
				p.type=j==1||j==nMonomers?CYTO:MONOMER;
				
				//velocity
				theta=M_PI*randNum->rand53();
				phi=M_PI*2*randNum->rand53();
				v.x=Vrms*cos(phi)*sin(theta);
				v.y=Vrms*sin(phi)*sin(theta);
				v.z=Vrms*cos(theta);
				
				System.addParticle(p,v,a);
			}
			
			//bond for a bond type
			bond.s[0]=startOffset+xa;
			bond.s[1]=System.readNParticles()-1;
			anchorMol.addBond(bond);
			
			//bond for a bond type
			bond.s[0]=startOffset+xb;
			bond.s[1]=System.readNParticles()-nMonomers;
			anchorMol.addBond(bond);
			
		}
		
		if(bc)
		{
			for(int j=1;j<nMonomers+1;j++)
			{
				position<T> p;
				threeVector<T> v;
				T phi,theta;
				
				//position
				p.x=(T)j*((anc[xb].x-anc[xc].x)/((T)nMonomers+1.0))+anc[xc].x+pos.x;
				p.y=(T)j*((anc[xb].y-anc[xc].y)/((T)nMonomers+1.0))+anc[xc].y+pos.y;
				p.z=(T)j*((anc[xb].z-anc[xc].z)/((T)nMonomers+1.0))+anc[xc].z+pos.z;
				p.type=j==1||j==nMonomers?CYTO:MONOMER;
				
				//velocity
				theta=M_PI*randNum->rand53();
				phi=M_PI*2*randNum->rand53();
				v.x=Vrms*cos(phi)*sin(theta);
				v.y=Vrms*sin(phi)*sin(theta);
				v.z=Vrms*cos(theta);
				
				System.addParticle(p,v,a);
			}
			
			//bond for a bond type
			bond.s[0]=startOffset+xb;
			bond.s[1]=System.readNParticles()-1;
			anchorMol.addBond(bond);
			
			//bond for a bond type
			bond.s[0]=startOffset+xc;
			bond.s[1]=System.readNParticles()-nMonomers;
			anchorMol.addBond(bond);
		}
		
		if(ca)
		{
			for(int j=1;j<nMonomers+1;j++)
			{
				position<T> p;
				threeVector<T> v;
				T phi,theta;
				
				//position
				p.x=(T)j*((anc[xc].x-anc[xa].x)/((T)nMonomers+1.0))+anc[xa].x+pos.x;
				p.y=(T)j*((anc[xc].y-anc[xa].y)/((T)nMonomers+1.0))+anc[xa].y+pos.y;
				p.z=(T)j*((anc[xc].z-anc[xa].z)/((T)nMonomers+1.0))+anc[xa].z+pos.z;
				p.type=j==1||j==nMonomers?CYTO:MONOMER;
			
				//velocity
				theta=M_PI*randNum->rand53();
				phi=M_PI*2*randNum->rand53();
				v.x=Vrms*cos(phi)*sin(theta);
				v.y=Vrms*sin(phi)*sin(theta);
				v.z=Vrms*cos(theta);
				
				System.addParticle(p,v,a);
			}
			
			//bond for a bond type
			bond.s[0]=startOffset+xc;
			bond.s[1]=System.readNParticles()-1;
			anchorMol.addBond(bond);
			
			//bond for a bond type
			bond.s[0]=startOffset+xa;
			bond.s[1]=System.readNParticles()-nMonomers;
			anchorMol.addBond(bond);
		}
	}
	
	//add the vertex bonds to the system
	System.addMolecule(anchorMol);
	
	free(faces);
	free(anc);
	return radius;
}


template <typename T>
T tube(Blob<T> &System, int nLipids, int lipidLength, threeVector<T> pos, T bondLength, T arealDensity, T *constants, int nConstants, T LRratio)
{
	//Error checking
	if(nConstants==4)
		std::cerr << "nConstants is 4, assuming CHAIN type...\n";
	if(nConstants==2)
		std::cerr << "nConstants is 2, assuming BOND type...\n";
	
	if(nConstants!=2 && nConstants!=4)
	{
		std::cerr << "Error: nConstants is not 2 or 4!(liposome)\n";
		return 0;
	}
	
	//Radius of tube
	T radius=sqrt((T)nLipids/(arealDensity*2.0*M_PI*LRratio));//density adjusted
	T length=sqrt((T)nLipids*LRratio/(arealDensity*2.0*M_PI));//density adjusted
	
	
	//adjust size of system if it is out of bounds
	threeVector<T> size=System.readSize();
	
	//This is oriented along the x axis, radiating along the y-z axii
	size.x=size.x<length?length:size.x;
	size.y=size.y<4.0*(radius+lipidLength*bondLength)?length:size.y;
	size.z=size.z<4.0*(radius+lipidLength*bondLength)?length:size.z;
	
	threeVector<T> offset;
	offset.x=0;
	offset.y=size.y/2.0;
	offset.z=size.z/2.0;
	
	System.setSize(size);
	
	//getting a fixed interval for the number of lipids along the radius and length.
	T latticeLength=sqrt(1.0/arealDensity);
	int nRLipids=floor(2*M_PI*radius/latticeLength);
	int nLLipids=floor(length/latticeLength);
	
	//This almost certainly changed the number of lipids.
	nLipids=nRLipids*nLLipids;
	
	int inner=floor(T(nRLipids)*(radius-(lipidLength*bondLength)/2.0)/radius/2.0);
	//int outer=floor(T(nRLipids)*(radius+(lipidLength*bondLength)/2.0)/radius/2.0);
	int outer=nRLipids-inner;
	
	//Initialize molecule
	fourVector<int> bond;
	molecule<T,fourVector<int> > m;
	
	if(nConstants==4)
	{
		m.setType(CHAIN);
		
		//bond for a chain type
		bond.s[START]=System.readNParticles();
		bond.s[NCHAINS]=nLipids;
		bond.s[CHAINLENGTH]=lipidLength;
		m.addBond(bond);
	}
	
	if(nConstants==2)
	{
		m.setType(BOND);
		m.allocBonds(nLipids*lipidLength);
	}
	
	for(int i=0;i<nConstants;i++)
		m.addConstant(constants[i]);
	
	std::cerr << "Number of inner lipids: " << inner << '\n';
	std::cerr << "Number of outer lipids: " << outer << '\n';
	
	
	MTRand *randNum=new MTRand(System.readSeed());
	//Velocities
	T Vrms=sqrt(3.0*System.readInitialTemp());
	
	//This is important to reduce excessive reallocations
	System.allocParticle(System.readNParticles()+nLipids*lipidLength);
	for(int i=0;i<nLLipids;i++)
	{
		T xPos=latticeLength*static_cast<T>(i);
		for(int j=0;j<inner;j++)
		{
			T thetaPos=static_cast<T>(j)*2.0*M_PI/static_cast<T>(inner);
			position <T> p;
			threeVector<T> v,a;
			//a is just a null vector
			a.x=0;
			a.y=0;
			a.z=0;
			
			for(int k=0;k<lipidLength;k++)
			{
				//inner monomers to outer monomers
				p.x=pos.x+offset.x+xPos;
				p.y=pos.y+offset.y+cos(thetaPos)*(radius-bondLength*(T)k-bondLength/2.0);
				p.z=pos.y+offset.z+sin(thetaPos)*(radius-bondLength*(T)k-bondLength/2.0);
				p.type=(k==lipidLength-1)?HEAD:TAIL;
				//velocity
				T theta=M_PI*randNum->rand53();
				T phi=M_PI*2.0*randNum->rand53();
				v.x=Vrms*cos(phi)*sin(theta);
				v.y=Vrms*sin(phi)*sin(theta);
				v.z=Vrms*cos(theta);
				System.addParticle(p,v,a);
				
				if(nConstants==2 && j!=lipidLength-1)
				{
					bond.s[0]=System.readNParticles();
					bond.s[1]=System.readNParticles()+1;
					m.addBond(bond);
				}
			}
			
		}
		for(int j=0;j<outer;j++)
		{
			T thetaPos=static_cast<T>(j)*2.0*M_PI/static_cast<T>(outer);
			position <T> p;
			threeVector<T> v,a;
			//a is just a null vector
			a.x=0;
			a.y=0;
			a.z=0;
			
			for(int k=0;k<lipidLength;k++)
			{
				//inner monomers to outer monomers
				p.x=pos.x+offset.x+xPos;
				p.y=pos.y+offset.y+cos(thetaPos)*(radius+bondLength*(T)k+bondLength/2.0);
				p.z=pos.y+offset.z+sin(thetaPos)*(radius+bondLength*(T)k+bondLength/2.0);
				p.type=(k==lipidLength-1)?HEAD:TAIL;
				//velocity
				T theta=M_PI*randNum->rand53();
				T phi=M_PI*2.0*randNum->rand53();
				v.x=Vrms*cos(phi)*sin(theta);
				v.y=Vrms*sin(phi)*sin(theta);
				v.z=Vrms*cos(theta);
				System.addParticle(p,v,a);
				
				if(nConstants==2 && j!=lipidLength-1)
				{
					bond.s[0]=System.readNParticles();
					bond.s[1]=System.readNParticles()+1;
					m.addBond(bond);
				}
			}
		}	
	}
	System.addMolecule(m);
	
	return radius;
}

template <typename T>
T liposome(Blob<T> &System, int nLipids, int lipidLength, threeVector<T> pos, T bondLength, T arealDensity, T *constants, int nConstants)
{
	//Error checking
	if(nConstants==4)
		std::cerr << "nConstants is 4, assuming CHAIN type...\n";
	if(nConstants==2)
		std::cerr << "nConstants is 2, assuming BOND type...\n";
	
	if(nConstants!=2 && nConstants!=4)
	{
		std::cerr << "Error: nConstants is not 2 or 4!(liposome)\n";
		return 0;
	}
	
	//Initialize molecule
	fourVector<int> bond;
	molecule<T,fourVector<int> > m;
	
	if(nConstants==4)
	{
		m.setType(CHAIN);
		
		//bond for a chain type
		bond.s[START]=System.readNParticles();
		bond.s[NCHAINS]=nLipids;
		bond.s[CHAINLENGTH]=lipidLength;
		m.addBond(bond);
	}
	
	if(nConstants==2)
	{
		m.setType(BOND);
		m.allocBonds(nLipids*lipidLength);
	}
	
	for(int i=0;i<nConstants;i++)
		m.addConstant(constants[i]);
	
	//Add lipids
	T radius=sqrt(((T)nLipids/arealDensity)/(4.0*M_PI));//density adjusted
	
	//adjust size of system if it is out of bounds
	threeVector<T> size=System.readSize();
	
	size.x=size.x<radius+pos.x+lipidLength*bondLength?pos.x*2.0:size.x;
	size.y=size.y<radius+pos.y+lipidLength*bondLength?pos.y*2.0:size.y;
	size.z=size.z<radius+pos.z+lipidLength*bondLength?pos.z*2.0:size.z;
	
	System.setSize(size);
	
	int inner=nLipids*pow(radius,2.0)/(pow(radius,2.0)+pow(radius+(T)(lipidLength*2-1)*bondLength,2.0));
	int outer=nLipids-inner;
	
	std::cerr << "Number of inner lipids: " << inner << '\n';
	std::cerr << "Number of outer lipids: " << outer << '\n';
	
	T s=3.6/sqrt(T(inner));
	T length=0;
	T dz=2.0/T(inner);
	T z=1.0-dz/2.0;
	
	MTRand *randNum=new MTRand(System.readSeed());
	//Velocities
	T Vrms=sqrt(3.0*System.readInitialTemp());
	
	//This is important to reduce excessive reallocations
	System.allocParticle(System.readNParticles()+nLipids*lipidLength);
	
	for(int i=0;i<inner;i++)
	{
		T phi,theta;
		position <T> p;
		threeVector<T> v,a;
		//a is just a null vector
		a.x=0;
		a.y=0;
		a.z=0;
		
		T r=sqrt(1.0-z*z);
		
		for(int j=0;j<lipidLength;j++)
		{
			//inner monomers to outer monomers
			p.x=pos.x+r*cos(length)*(radius+bondLength*(T)j);
			p.y=pos.y+r*sin(length)*(radius+bondLength*(T)j);
			p.z=pos.z+z*(radius+bondLength*(T)j);
			p.type=(j==0)?HEAD:TAIL;
			//velocity
			theta=M_PI*randNum->rand53();
			phi=M_PI*2.0*randNum->rand53();
			v.x=Vrms*cos(phi)*sin(theta);
			v.y=Vrms*sin(phi)*sin(theta);
			v.z=Vrms*cos(theta);
			System.addParticle(p,v,a);
			
			if(nConstants==2 && j!=lipidLength-1)
			{
				bond.s[0]=System.readNParticles();
				bond.s[1]=System.readNParticles()+1;
				m.addBond(bond);
			}
		}
		
		z=z-dz;
		length=length+s/r;
	}

	
	s=3.6/sqrt(T(outer));
	length=0;
	dz=2.0/T(outer);
	z=1.0-dz/2.0;
	
	for(int i=0;i<outer;i++)
	{
		T phi,theta;
		position <T> p;
		threeVector<T> v,a;
		//a is just a null vector
		a.x=0;
		a.y=0;
		a.z=0;
		
		T r=sqrt(1.0-z*z);
		
		for(int j=0;j<lipidLength;j++)
		{
			//inner monomers to outer monomers
			p.x=pos.x+r*cos(length)*(radius+(bondLength*(T)(lipidLength+j)));
			p.y=pos.y+r*sin(length)*(radius+(bondLength*(T)(lipidLength+j)));
			p.z=pos.z+z*(radius+(bondLength*(T)(lipidLength+j)));
			p.type=(j==lipidLength-1)?HEAD:TAIL;
			//velocity
			theta=M_PI*randNum->rand53();
			phi=M_PI*2.0*randNum->rand53();
			v.x=Vrms*cos(phi)*sin(theta);
			v.y=Vrms*sin(phi)*sin(theta);
			v.z=Vrms*cos(theta);
			System.addParticle(p,v,a);
			
			if(nConstants==2 && j!=lipidLength-1)
			{
				bond.s[0]=System.readNParticles();
				bond.s[1]=System.readNParticles()+1;
				m.addBond(bond);
			}
		}
		
		z=z-dz;
		length=length+s/r;
	}
	
	System.addMolecule(m);
	
	return radius;
}

template <typename T>
std::vector<position<T> > spiralShell(int nPoints, T radius, int type)
{
	T s=3.6/sqrt(T(nPoints));
	T length=0;
	T dz=2.0/T(nPoints);
	T z=1.0-dz/2.0;
	
	std::vector<position<T> > shell;
	
	for(int i=0;i<nPoints;i++)
	{
		T r=sqrt(1.0-z*z);
		
		position<T> p;
		p.x=r*cos(length)*radius;
		p.y=r*sin(length)*radius;
		p.z=z*radius;
		p.type=type;
		
		shell.push_back(p);
		
		z=z-dz;
		length=length+s/r;
	}
	
	return shell;
}

template <typename T>
T liposome(Blob<T> &System, int nLipids, int lipidLength, threeVector<T> pos, T bondLength, T arealDensity,
	   T iPercent, T oPercent, int *types, T *aConstants, T *bConstants, int nConstants)
{
	//Error checking
	if(nConstants==4)
		std::cerr << "nConstants is 4, assuming BOND and BEND types...\n";
	
	if(lipidLength<3)
	{
		std::cerr << "Error(liposome): lipidLength is less than 3\n";
		return 0;
	}
	
	if(nConstants!=4)
	{
		std::cerr << "Error(liposome): nConstants is not 4!\n";
		return 0;
	}
	
	//Initialize molecule
	fourVector<int> bond;
	molecule<T,fourVector<int> > mBondA;
	molecule<T,fourVector<int> > mBendA;
	molecule<T,fourVector<int> > mBondB;
	molecule<T,fourVector<int> > mBendB;
	
	//allocating way more than I need lol
	mBondA.setType(BOND);
	mBondA.allocBonds(nLipids*lipidLength);
	mBendA.setType(BEND);
	mBendA.allocBonds(nLipids*lipidLength);
	
	mBondB.setType(BOND);
	mBondB.allocBonds(nLipids*lipidLength);
	mBendB.setType(BEND);
	mBendB.allocBonds(nLipids*lipidLength);
	
	//We aren't sure of the order until it is checked
	mBondA.addConstant(aConstants[CHAINBOND]);
	mBondA.addConstant(aConstants[CHAINBOND+1]);
	mBendA.addConstant(aConstants[CHAINBEND]);
	mBendA.addConstant(aConstants[CHAINBEND+1]);
	
	mBondB.addConstant(bConstants[CHAINBOND]);
	mBondB.addConstant(bConstants[CHAINBOND+1]);
	mBendB.addConstant(bConstants[CHAINBEND]);
	mBendB.addConstant(bConstants[CHAINBEND+1]);
	
	//Add lipids
	T radius=sqrt(((T)nLipids/arealDensity)/(4.0*M_PI));//density adjusted
	
	//adjust size of system if it is out of bounds
	threeVector<T> size=System.readSize();
	
	size.x=size.x<radius+pos.x+lipidLength*bondLength?pos.x*2.0:size.x;
	size.y=size.y<radius+pos.y+lipidLength*bondLength?pos.y*2.0:size.y;
	size.z=size.z<radius+pos.z+lipidLength*bondLength?pos.z*2.0:size.z;
	
	System.setSize(size);
	
	//with inner and outer surface area compensation
	int inner=nLipids*pow(radius,2.0)/(pow(radius,2.0)+pow(radius+(T)(lipidLength*2-1)*bondLength,2.0));
	int outer=nLipids-inner;
	
	std::cerr << "Number of inner lipids: " << inner << '\n';
	std::cerr << "Number of outer lipids: " << outer << '\n';
	
	//We need a shuffle list to ensure that the ratio is correct
	T *iShuffle=new T[inner];
	T *oShuffle=new T[outer];
	
	for(int i=0;i<inner;i++)
		iShuffle[i]=T(i)/T(inner);
	for(int i=0;i<outer;i++)
		oShuffle[i]=T(i)/T(outer);
	
	MTRand randNum(System.readSeed());
	for(int i=0;i<inner;i++)
	{
		int otherLipid=randNum.randInt(inner-1);//this is range inclusive!
		T buf=iShuffle[i];
		iShuffle[i]=iShuffle[otherLipid];
		iShuffle[otherLipid]=buf;
	}
	for(int i=0;i<outer;i++)
	{
		int otherLipid=randNum.randInt(outer-1);//this is range inclusive
		T buf=oShuffle[i];
		oShuffle[i]=oShuffle[otherLipid];
		oShuffle[otherLipid]=buf;
	}
	
	//some constants
	T s=3.6/sqrt(T(inner));
	T length=0;
	T dz=2.0/T(inner);
	T z=1.0-dz/2.0;
	
	//Velocities
	T Vrms=sqrt(3.0*System.readInitialTemp());
	
	//This is important to reduce excessive reallocations
	System.allocParticle(System.readNParticles()+nLipids*lipidLength);
	
	for(int i=0;i<inner;i++)
	{
		T phi,theta;
		position <T> p;
		threeVector<T> v,a;
		//a is just a null vector
		a.x=0;
		a.y=0;
		a.z=0;
		
		T r=sqrt(1.0-z*z);
		
		for(int j=0;j<lipidLength;j++)
		{
			//inner monomers to outer monomers
			p.x=pos.x+r*cos(length)*(radius+bondLength*(T)j);
			p.y=pos.y+r*sin(length)*(radius+bondLength*(T)j);
			p.z=pos.z+z*(radius+bondLength*(T)j);
			if(iShuffle[i]>iPercent)
				p.type=(j!=0)?types[0]:types[1];
			else
				p.type=(j!=0)?types[2]:types[3];
			
			//velocity
			theta=M_PI*randNum.rand53();
			phi=M_PI*2.0*randNum.rand53();
			v.x=Vrms*cos(phi)*sin(theta);
			v.y=Vrms*sin(phi)*sin(theta);
			v.z=Vrms*cos(theta);
			System.addParticle(p,v,a);
			
			if(j<lipidLength-1)
			{
				bond.s[0]=System.readNParticles()-1;
				bond.s[1]=System.readNParticles();
				if(iShuffle[i]>iPercent)
					mBondA.addBond(bond);
				else
					mBondB.addBond(bond);
				
			}
			if(j<lipidLength-2)
			{
				bond.s[0]=System.readNParticles()-1;
				bond.s[1]=System.readNParticles();
				bond.s[2]=System.readNParticles()+1;
				if(iShuffle[i]>iPercent)
					mBendA.addBond(bond);
				else
					mBendB.addBond(bond);
				
			}
		}
		
		z=z-dz;
		length=length+s/r;
	}
	
	s=3.6/sqrt(T(outer));
	length=0;
	dz=2.0/T(outer);
	z=1.0-dz/2.0;
	
	for(int i=0;i<outer;i++)
	{
		T phi,theta;
		position <T> p;
		threeVector<T> v,a;
		//a is just a null vector
		a.x=0;
		a.y=0;
		a.z=0;
		
		T r=sqrt(1.0-z*z);
		
		for(int j=0;j<lipidLength;j++)
		{
			//inner monomers to outer monomers
			p.x=pos.x+r*cos(length)*(radius+(bondLength*(T)(lipidLength+j)));
			p.y=pos.y+r*sin(length)*(radius+(bondLength*(T)(lipidLength+j)));
			p.z=pos.z+z*(radius+(bondLength*(T)(lipidLength+j)));
			if(oShuffle[i]>oPercent)
				p.type=(j!=lipidLength-1)?types[0]:types[1];
			else
				p.type=(j!=lipidLength-1)?types[2]:types[3];
			//velocity
			theta=M_PI*randNum.rand53();
			phi=M_PI*2.0*randNum.rand53();
			v.x=Vrms*cos(phi)*sin(theta);
			v.y=Vrms*sin(phi)*sin(theta);
			v.z=Vrms*cos(theta);
			System.addParticle(p,v,a);
			
			if(j<lipidLength-1)
			{
				bond.s[0]=System.readNParticles()-1;
				bond.s[1]=System.readNParticles();
				if(oShuffle[i]>oPercent)
					mBondA.addBond(bond);
				else
					mBondB.addBond(bond);
				
			}
			if(j<lipidLength-2)
			{
				bond.s[0]=System.readNParticles()-1;
				bond.s[1]=System.readNParticles();
				bond.s[2]=System.readNParticles()+1;
				if(oShuffle[i]>oPercent)
					mBendA.addBond(bond);
				else
					mBendB.addBond(bond);
				
			}
		}
		
		z=z-dz;
		length=length+s/r;
	}
	if(mBondA.readNBond()>0)
		System.addMolecule(mBondA);
	if(mBendA.readNBond()>0)
		System.addMolecule(mBendA);
	if(mBondB.readNBond()>0)
		System.addMolecule(mBondB);
	if(mBendB.readNBond()>0)
		System.addMolecule(mBendB);
	
	delete iShuffle;
	delete oShuffle;
	
	return radius;
}

//aspectRatio is x/y
template <typename T>
bool bilayer(Blob<T> &System, int nLipids, int lipidLength, threeVector<T> pos, T bondLength, T arealDensity, T *constants, int nConstants, double aspectRatio)
{
return bilayer<T>(System, nLipids, lipidLength, pos, bondLength, arealDensity, constants, nConstants, aspectRatio, HEAD, TAIL, HEAD, TAIL);
}

//aspectRatio is x/y
template <typename T>
bool bilayer(Blob<T> &System, int nLipids, int lipidLength, threeVector<T> pos, T bondLength, T arealDensity, T *constants, int nConstants, double aspectRatio, int headTypeTop, int tailTypeTop, int headTypeBot, int tailTypeBot)
{
	//Error checking
	if(nConstants==4)
		std::cerr << "nConstants is 4, assuming CHAIN type (bilayer)...\n";
	if(nConstants==2)
		std::cerr << "nConstants is 2, assuming BOND type (bilayer)...\n";
	
	if(nConstants!=2 && nConstants!=4)
	{
		std::cerr << "Error: nConstants is not 2 or 4![bool bilayer()]\n";
		return true;
	}
	
	//Initialize molecule
	fourVector<int> bond;
	molecule<T,fourVector<int> > m;
	
	if(nConstants==4)
	{
		m.setType(CHAIN);
		
		//bond for a chain type
		bond.s[START]=System.readNParticles();
		bond.s[NCHAINS]=0;//nLipids;
		bond.s[CHAINLENGTH]=lipidLength;
		//m.addBond(bond);
	}
	
	if(nConstants==2)
	{
		m.setType(BOND);
		m.allocBonds(nLipids*lipidLength);
	}
	
	for(int i=0;i<nConstants;i++)
		m.addConstant(constants[i]);
	
	//Add lipids
	
	//adjust size of system if it is out of bounds
	threeVector<T> s;
	s.x=sqrt((T)nLipids/arealDensity)*sqrt(aspectRatio);
	s.y=sqrt((T)nLipids/arealDensity)/sqrt(aspectRatio);
	s.z=pos.z*2.0;
	std::cerr << "Minimum required system size is (bilayer): " << s.x << '\t' << s.y << '\t' << s.z << '\n';
	
	threeVector<T> size=System.readSize();
	
	std::cerr << "Current system size is (bilayer): " << size.x << '\t' << size.y << '\t' << size.z << '\n';
	
	size.x=(size.x<s.x)?s.x:size.x;
	size.y=(size.y<s.y)?s.y:size.y;
	size.z=(size.z<s.z)?s.z:size.z;
	
	System.setSize(size);
	
	//The worlds most pointless line. I don't think it helps anything, 
	// but I don't know why I did it...
	size=System.readSize();
	
	std::cerr << "Adjusted system size is (bilayer): " << size.x << '\t' << size.y << '\t' << size.z << '\n';	
	
	MTRand *randNum=new MTRand(System.readSeed());
	//Velocities
	T Vrms=sqrt(3.0*System.readInitialTemp());
	
	//This is important to reduce excessive reallocations
	System.allocParticle(System.readNParticles()+2*nLipids*lipidLength);
	
	threeVector<int> latticeSize;
	threeVector<T> latticeLength;
	
	//change the number of lipids if it doesn't produce the correct sized bilayer
	latticeLength.x=sqrt(2.0/arealDensity);
	latticeLength.y=sqrt(2.0/arealDensity);
	latticeLength.z=sqrt(2.0/arealDensity);
	
	latticeSize.x=static_cast<int>(size.x)/latticeLength.x+1;
	latticeSize.y=static_cast<int>(size.y)/latticeLength.y+1;
	latticeSize.z=static_cast<int>(size.z)/latticeLength.z;
	
	//initialize lipid positions
	for(int i=0;i<latticeSize.x;i++)
	{
		for(int j=0;j<latticeSize.y;j++)
		{
			int lipidIndex=2*(j*latticeSize.x+i);
			//if(lipidIndex%1000==0)
			//	std::cerr << lipidIndex << std::endl;
			T theta,phi;
			
			position<T> p;
			threeVector<T> v;
			threeVector<T> a;
			//set all acceleration to 0 initially
			a.x=0;
			a.y=0;
			a.z=0;
			
			if(lipidIndex<nLipids)
			{
				for(int k=0;k<lipidLength;k++)
				{
					//position
					p.x=i*latticeLength.x+0.001;
					p.y=j*latticeLength.y+0.001;
					p.z=pos.z+bondLength*(lipidLength-k);
					p.type=(k==lipidLength-1)?headTypeTop:tailTypeTop;
					
					//velocity
					theta=M_PI*randNum->rand53();
					phi=M_PI*2*randNum->rand53();
					v.x=Vrms*cos(phi)*sin(theta);
					v.y=Vrms*sin(phi)*sin(theta);
					v.z=Vrms*cos(theta);
					
					//save it
					System.addParticle(p,v,a);
					
					//if it is of BOND type and not the last monomer
					if(nConstants==2 && k!=lipidLength-1)
					{
						//bond
						bond.s[0]=System.readNParticles()-1;
						bond.s[1]=System.readNParticles();
						m.addBond(bond);
					}
				}
				if(nConstants==4)
				{
					bond.s[NCHAINS]++;
				}
			}
			if(lipidIndex+1<nLipids)
			{
				for(int k=0;k<lipidLength;k++)
				{
					//position
					p.x=i*latticeLength.x+0.001;
					p.y=j*latticeLength.y+0.001;
					p.z=pos.z+bondLength*(lipidLength+k+1);
					p.type=(k==lipidLength-1)?headTypeBot:tailTypeTop;
					
					//velocity
					theta=M_PI*randNum->rand53();
					phi=M_PI*2*randNum->rand53();
					v.x=Vrms*cos(phi)*sin(theta);
					v.y=Vrms*sin(phi)*sin(theta);
					v.z=Vrms*cos(theta);
					
					//save it
					System.addParticle(p,v,a);
					
					//if it is of BOND type and not the last monomer
					if(nConstants==2 && k!=lipidLength-1)
					{
						//bond
						bond.s[0]=System.readNParticles()-1;
						bond.s[1]=System.readNParticles();
						m.addBond(bond);
					}
				}
				if(nConstants==4)
				{
					bond.s[NCHAINS]++;
				}
			}
		}
	}
	if(nConstants==4)
	{
		m.addBond(bond);
	}
	System.addMolecule(m);
	return false;
}

/*
//aspectRatio is x/y
template <typename T>
bool bilayer(Blob<T> &System, int nLipids, int lipidLength, threeVector<T> pos, T bondLength, T arealDensity, 
	     T lPercent, T uPercent, int *types, T *aConstants, T *bConstants, int nConstants, int nConstants, double aspectRatio)
{
	//Error checking
	if(nConstants==4)
		std::cerr << "nConstants is 4, assuming CHAIN type (bilayer)...\n";
	if(nConstants==2)
		std::cerr << "nConstants is 2, assuming BOND type (bilayer)...\n";
	
	if(nConstants!=2 && nConstants!=4)
	{
		std::cerr << "Error: nConstants is not 2 or 4![bool bilayer()]\n";
		return true;
	}
	
	//Initialize molecule
	fourVector<int> bond;
	molecule<T,fourVector<int> > m;
	
	if(nConstants==4)
	{
		m.setType(CHAIN);
		
		//bond for a chain type
		bond.s[START]=System.readNParticles();
		bond.s[NCHAINS]=0;//nLipids;
		bond.s[CHAINLENGTH]=lipidLength;
		//m.addBond(bond);
	}
	
	if(nConstants==2)
	{
		m.setType(BOND);
		m.allocBonds(nLipids*lipidLength);
	}
	
	for(int i=0;i<nConstants;i++)
		m.addConstant(constants[i]);
	
	//Add lipids
	
	//adjust size of system if it is out of bounds
	threeVector<T> s;
	s.x=sqrt((T)nLipids/arealDensity)*sqrt(aspectRatio);
	s.y=sqrt((T)nLipids/arealDensity)/sqrt(aspectRatio);
	s.z=pos.z*2.0;
	std::cerr << "Minimum required system size is (bilayer): " << s.x << '\t' << s.y << '\t' << s.z << '\n';
	
	threeVector<T> size=System.readSize();
	
	std::cerr << "Current system size is (bilayer): " << size.x << '\t' << size.y << '\t' << size.z << '\n';
	
	size.x=(size.x<s.x)?s.x:size.x;
	size.y=(size.y<s.y)?s.y:size.y;
	size.z=(size.z<s.z)?s.z:size.z;
	
	System.setSize(size);
	
	//The worlds most pointless line. I don't think it helps anything, 
	// but I don't know why I did it...
	size=System.readSize();
	
	std::cerr << "Adjusted system size is (bilayer): " << size.x << '\t' << size.y << '\t' << size.z << '\n';	
	
	MTRand *randNum=new MTRand(System.readSeed());
	//Velocities
	T Vrms=sqrt(3.0*System.readInitialTemp());
	
	//This is important to reduce excessive reallocations
	System.allocParticle(System.readNParticles()+2*nLipids*lipidLength);
	
	threeVector<int> latticeSize;
	threeVector<T> latticeLength;
	
	//change the number of lipids if it doesn't produce the correct sized bilayer
	latticeLength.x=sqrt(2.0/arealDensity);
	latticeLength.y=sqrt(2.0/arealDensity);
	latticeLength.z=sqrt(2.0/arealDensity);
	
	latticeSize.x=static_cast<int>(size.x)/latticeLength.x+1;
	latticeSize.y=static_cast<int>(size.y)/latticeLength.y+1;
	latticeSize.z=static_cast<int>(size.z)/latticeLength.z;
	/*
	//readjust lattice size to minimize the lattice length difference
	for(int i=0;i<500;i++)
	{
		if(latticeLength.x-latticeLength.y>0.0001)
		{
			latticeSize.x++;
			latticeSize.y--;
			latticeLength.x=size.x/latticeSize.x;
			latticeLength.y=size.y/latticeSize.y;
		}
		if(latticeLength.x-latticeLength.y<-0.0001)
		{
			latticeSize.x--;
			latticeSize.y++;
			latticeLength.x=size.x/latticeSize.x;
			latticeLength.y=size.y/latticeSize.y;
		}
		//std::cerr << i << '\t' << latticeLength.x-latticeLength.y << '\t' << latticeSize.x << '\t' << latticeSize.y << std::endl;
	}
	*/
	/*
	std::vector<bool> lShuffle(latticeSize.x*latticeSize.y,false),uShuffle(latticeSize.x*latticeSize.y,false);
	
	for(int i=0;i<latticeLength.x*latticeLength.y*lPercent;i++)
		lShuffle[i]=true;
	for(int i=0;i<latticeLength.x*latticeLength.y*uPercent;i++)
		uShuffle[i]=true;
	
	MTRand randNum(System.readSeed());
	for(int i=0;i<latticeLength.x*latticeLength.y;i++)
	{
		int otherLipid=randNum.randInt(latticeLength.x*latticeLength.y-1);//this is range inclusive!
		bool buf=lShuffle[i];
		lShuffle[i]=lShuffle[otherLipid];
		lShuffle[otherLipid]=buf;
	}
	for(int i=0;i<latticeLength.x*latticeLength.y;i++)
	{
		int otherLipid=randNum.randInt(latticeLength.x*latticeLength.y-1);//this is range inclusive
		bool buf=uShuffle[i];
		uShuffle[i]=uShuffle[otherLipid];
		uShuffle[otherLipid]=buf;
	}
	
	//initialize lipid positions
	for(int i=0;i<latticeSize.x;i++)
	{
		for(int j=0;j<latticeSize.y;j++)
		{
			int lipidIndex=j*latticeSize.x+i;
			//if(lipidIndex%1000==0)
			//	std::cerr << lipidIndex << std::endl;
			T theta,phi;
			
			position<T> p;
			threeVector<T> v;
			threeVector<T> a;
			//set all acceleration to 0 initially
			a.x=0;
			a.y=0;
			a.z=0;
			
			//if(lipidIndex<nLipids)
			{
				for(int k=0;k<lipidLength;k++)
				{
					//position
					p.x=i*latticeLength.x+0.001;
					p.y=j*latticeLength.y+0.001;
					p.z=pos.z+bondLength*(lipidLength-k);
					if(lShuffle[lipidIndex])
						p.type=(k==lipidLength-1)?HEAD:TAIL;
					else
						p.type=(k==lipidLength-1)?HEAD2:TAIL2;
					
					//velocity
					theta=M_PI*randNum->rand53();
					phi=M_PI*2*randNum->rand53();
					v.x=Vrms*cos(phi)*sin(theta);
					v.y=Vrms*sin(phi)*sin(theta);
					v.z=Vrms*cos(theta);
					
					//save it
					System.addParticle(p,v,a);
					
					//if it is of BOND type and not the last monomer
					if(nConstants==2 && k!=lipidLength-1)
					{
						//bond
						bond.s[0]=System.readNParticles()-1;
						bond.s[1]=System.readNParticles();
						m.addBond(bond);
					}
				}
				if(nConstants==4)
				{
					bond.s[NCHAINS]++;
				}
			}
			//if(lipidIndex+1<nLipids)
			{
				for(int k=0;k<lipidLength;k++)
				{
					//position
					p.x=i*latticeLength.x+0.001;
					p.y=j*latticeLength.y+0.001;
					p.z=pos.z+bondLength*(lipidLength+k+1);
					if(uShuffle[lipidIndex])
						p.type=(k==lipidLength-1)?HEAD:TAIL;
					else
						p.type=(k==lipidLength-1)?HEAD2:TAIL2;
					
					//velocity
					theta=M_PI*randNum->rand53();
					phi=M_PI*2*randNum->rand53();
					v.x=Vrms*cos(phi)*sin(theta);
					v.y=Vrms*sin(phi)*sin(theta);
					v.z=Vrms*cos(theta);
					
					//save it
					System.addParticle(p,v,a);
					
					//if it is of BOND type and not the last monomer
					if(nConstants==2 && k!=lipidLength-1)
					{
						//bond
						bond.s[0]=System.readNParticles()-1;
						bond.s[1]=System.readNParticles();
						m.addBond(bond);
					}
				}
				if(nConstants==4)
				{
					bond.s[NCHAINS]++;
				}
			}
		}
	}
	if(nConstants==4)
	{
		m.addBond(bond);
	}
	System.addMolecule(m);
}
*/
#define HEXAGONAL_ASPECT_RATIO sqrt(3.0/2.0)

template <typename T>
bool flatHexagonalCyto(Blob<T> &System, twoVector<int> nAnchors, int nMonomers, threeVector<T> pos, T *constants, int nConstants)
{
	//Error Checking
	if(System.readSize().x<=0 || System.readSize().y<=0 || System.readSize().z<=0)
	{
		std::cerr << "System size not set for flatHexagonalCyto!\n";
		return true;
	}
	
	MTRand *randNum=new MTRand(System.readSeed());
	//Velocities
	T Vrms=sqrt(3.0*System.readInitialTemp());
	
	//This is important to reduce excessive reallocations, also to maintain the same pointer for the particle list
	System.allocParticle(System.readNParticles()+nAnchors.x*nAnchors.y*(1+3*nMonomers));
	
	twoVector<T> latticeLength;
	latticeLength.x=System.readSize().x/(T)nAnchors.x;
	latticeLength.y=System.readSize().y/(T)nAnchors.y;
	
	//Add Cytoskeleton
	//retrieve starting index for anchors
	int anchorStart=System.readNParticles();
	
	//set anchors
	for(int i=0;i<nAnchors.x;i++)
	{
		for(int j=0;j<nAnchors.y;j++)
		{
			T theta,phi;
			
			position<T> p;
			//add lattice length seperated anchors and half lengths every other row
			p.x=(T)i*latticeLength.x+pos.x+0.001;
			p.y=(T)j*latticeLength.y+pos.y+(i%2==0?latticeLength.y/2.0:0)+0.001;
			p.z=pos.z;
			p.type=ANCHOR;
			
			//make sure it is within box
			p.x-=(p.x>System.readSize().x)?System.readSize().x:0;
			p.y-=(p.y>System.readSize().y)?System.readSize().y:0;
			p.z-=(p.z>System.readSize().z)?System.readSize().z:0;
			
			//make some variable for velocity and acceleration
			threeVector<T> v;
			threeVector<T> a;
			//set all acceleration to 0 initially
			a.x=0;
			a.y=0;
			a.z=0;
			
			//velocity
			theta=M_PI*randNum->rand53();
			phi=M_PI*2*randNum->rand53();
			v.x=Vrms*cos(phi)*sin(theta);
			v.y=Vrms*sin(phi)*sin(theta);
			v.z=Vrms*cos(theta);
			
			//save it
			System.addParticle(p,v,a);
		}
	}
	
	//retrieve starting index for chains
	int chainStart=System.readNParticles();
	
	//Initialize molecule
	fourVector<int> bond;
	molecule<T,fourVector<int> > monomers,anchors;
	
	anchors.setType(BOND);
	
	if(nConstants==4)
	{
		monomers.setType(CHAIN);
		
		//bond for a chain type
		bond.s[START]=chainStart;
		bond.s[NCHAINS]=3*nAnchors.x*nAnchors.y;
		bond.s[CHAINLENGTH]=nMonomers;
		monomers.addBond(bond);
		
		//for anchor bonds
		anchors.addConstant(constants[CHAINBOND+ABOND]);
		anchors.addConstant(constants[CHAINBOND+KBOND]);
	}
	
	if(nConstants==2)
	{
		monomers.setType(BOND);
		monomers.allocBonds(nMonomers*nAnchors.x*nAnchors.y);
		
		//for anchor bonds
		anchors.addConstant(constants[ABOND]);
		anchors.addConstant(constants[KBOND]);
	}
	
	for(int i=0;i<nConstants;i++)
		monomers.addConstant(constants[i]);
	
	bool *flag=new bool[nAnchors.x*nAnchors.y];
	for(int i=0;i<nAnchors.x*nAnchors.y;i++)
		flag[i]=false;
	
	T rc=(latticeLength.x*latticeLength.x+latticeLength.y*latticeLength.y);
	
	//set monomers between anchors
	for(int i=0;i<nAnchors.x*nAnchors.y;i++)
	{
		//get an anchor position to create a starting position
		position<T> anchor=System.getPositions()[i+anchorStart];
		flag[i]=true;
		
		for(int j=0;j<nAnchors.x*nAnchors.y;j++)
		{
			threeVector<T> d;
			//vector between anchors
			d.x=System.getPositions()[j+anchorStart].x-System.getPositions()[i+anchorStart].x;
			d.y=System.getPositions()[j+anchorStart].y-System.getPositions()[i+anchorStart].y;
			d.z=System.getPositions()[j+anchorStart].z-System.getPositions()[i+anchorStart].z;
			
			//minimum image
			d.x-=(d.x>System.readSize().x/2.0)?System.readSize().x:0;
			d.y-=(d.y>System.readSize().y/2.0)?System.readSize().y:0;
			d.z-=(d.z>System.readSize().z/2.0)?System.readSize().z:0;
			d.x+=(d.x<-System.readSize().x/2.0)?System.readSize().x:0;
			d.y+=(d.y<-System.readSize().y/2.0)?System.readSize().y:0;
			d.z+=(d.z<-System.readSize().z/2.0)?System.readSize().z:0;
			
			if(d.x*d.x+d.y*d.y+d.z*d.z<rc*1.1 && !flag[j])
			{
				//now add monomers between anchors
				for(int k=0;k<nMonomers;k++)
				{
					position<T> monomer;
					
					//next monomer is an incremental distance away
					monomer.x=System.getPositions()[i+anchorStart].x+(1.0+k)*d.x/T(nMonomers+1);
					monomer.y=System.getPositions()[i+anchorStart].y+(1.0+k)*d.y/T(nMonomers+1);
					monomer.z=System.getPositions()[i+anchorStart].z+(1.0+k)*d.z/T(nMonomers+1);
					monomer.type=CYTO;
					
					//does it run over edge?
					monomer.x-=(monomer.x>System.readSize().x)?System.readSize().x:0;
					monomer.y-=(monomer.y>System.readSize().y)?System.readSize().y:0;
					monomer.z-=(monomer.z>System.readSize().z)?System.readSize().z:0;
					monomer.x+=(monomer.x<0)?System.readSize().x:0;
					monomer.y+=(monomer.y<0)?System.readSize().y:0;
					monomer.z+=(monomer.z<0)?System.readSize().z:0;
					
					threeVector<T> v,a;
					a.x=0;
					a.y=0;
					a.z=0;
					
					//velocity
					T theta=M_PI*randNum->rand53();
					T phi=M_PI*2*randNum->rand53();
					v.x=Vrms*cos(phi)*sin(theta);
					v.y=Vrms*sin(phi)*sin(theta);
					v.z=Vrms*cos(theta);
					
					//save it
					System.addParticle(monomer,v,a);
					
					//if it is of BOND type and not the last monomer
					if(nConstants==2 && k!=nMonomers-1)
					{
						//bond
						bond.s[0]=System.readNParticles()-1;
						bond.s[1]=System.readNParticles();
						monomers.addBond(bond);
					}
					//first monomer to current anchor
					if(k==0)
					{
						bond.s[0]=System.readNParticles()-1;
						bond.s[1]=i+anchorStart;
						anchors.addBond(bond);
					}
					//last monomer to endpoint anchor
					if(k==nMonomers-1)
					{
						bond.s[0]=System.readNParticles()-1;
						bond.s[1]=j+anchorStart;
						anchors.addBond(bond);
					}
				}
			}
		}
	}
	
	//put our new molecule structures into the system
	System.addMolecule(monomers);
	System.addMolecule(anchors);
	delete flag;
	
	return false;
}

template <typename T>
bool solventFill(Blob<T> &System, T density, int solventType)
{
	//Error Checking
	if(System.readSize().x<=0 || System.readSize().y<=0 || System.readSize().z<=0)
	{
		std::cerr << "System size not set for solventFill!\n";
		return true;
	}
	
	//Generate 3D cubic lattice to fill system
	T latticeLength=pow(1.0/density,1.0/3.0);
	T threshold=latticeLength-0.001;
	std::cerr << "Assuming no placement threshold: " << threshold << '\n';
	
	threeVector<int> latticeSize;
	latticeSize.x=System.readSize().x/latticeLength;
	latticeSize.y=System.readSize().y/latticeLength;
	latticeSize.z=System.readSize().z/latticeLength;
	int nSolvent=latticeSize.x*latticeSize.y*latticeSize.z;
	
	MTRand *randNum=new MTRand(System.readSeed());
	//Velocities
	T Vrms=sqrt(3.0*System.readInitialTemp());
	
	//This is important to reduce excessive reallocations
	System.allocParticle(System.readNParticles()+nSolvent);
	
	if(System.readNParticles()==0)
	{
		std::cerr << "Adding solvent...\n";
		
		//Add particles
		for(int i=0;i<latticeSize.x;i++)
		{
			for(int j=0;j<latticeSize.y;j++)
			{
				for(int k=0;k<latticeSize.z;k++)
				{
					int index=i+latticeSize.x*j+latticeSize.x*latticeSize.y*k;
					
					T theta,phi;
					
					position<T> p;
					p.x=(T)i*latticeLength+0.001;
					p.y=(T)j*latticeLength+0.001;
					p.z=(T)k*latticeLength+0.001;
					p.type=solventType;
					
					threeVector<T> v,a,minImg;
					//set all acceleration to 0 initially
					a.x=0;
					a.y=0;
					a.z=0;
					
					//velocity
					theta=M_PI*randNum->rand53();
					phi=M_PI*2*randNum->rand53();
					v.x=Vrms*cos(phi)*sin(theta);
					v.y=Vrms*sin(phi)*sin(theta);
					v.z=Vrms*cos(theta);
					//save it
					System.addParticle(p,v,a);
				}
			}
		}
			
		std::cerr << "Done adding solvent.\n";
	}
	else
	{
		bool *tooClose=new bool[nSolvent];
		
		//Place this after the step that allocates particles
		//Cell <T> nearby(System.getPositions(), System.readNParticles(), threshold, System.readSize());
		//nearby.build();
		
		std::cerr << "Checking distances...\n";
		
		for(int i=0;i<nSolvent;i++)
			tooClose[i]=false;
		
		int otherOffset=System.readNParticles();
		
		//Check to see if any other particles in the system are nearby
		for(int i=0;i<latticeSize.x;i++)
		{
			for(int j=0;j<latticeSize.y;j++)
			{
				for(int k=0;k<latticeSize.z;k++)
				{
					int index=i+latticeSize.x*j+latticeSize.x*latticeSize.y*k;
					
					position<T> p;
					p.x=(T)i*latticeLength+0.001;
					p.y=(T)j*latticeLength+0.001;
					p.z=(T)k*latticeLength+0.001;
					p.type=solventType;
					
					threeVector<T> minImg;
					minImg=0;
					
					//for(int m=nearby.query(p,minImg);m!=-1;m=nearby.query(p,minImg))
					for(int m=0;m<otherOffset;m++)
					{
						threeVector<T> d;
						d.x=p.x-System.getPositions()[m].x-minImg.x;
						d.y=p.y-System.getPositions()[m].y-minImg.y;
						d.z=p.z-System.getPositions()[m].z-minImg.z;
						//std::cerr << p.x << ' ' << p.y << ' ' << p.z << '\n';
						//std::cerr << d.x*d.x+d.y*d.y+d.z*d.z << '\n';
						//std::cin.get();
						if(d.x*d.x+d.y*d.y+d.z*d.z<threshold*threshold)
						{
							//std::cerr << "????????????????\n";
							//std::cin.get();
							tooClose[index]=true;
						}
					}
				}
			}
		}
		
		std::cerr << "Adding solvent...\n";
		
		//Add particles
		for(int i=0;i<latticeSize.x;i++)
		{
			for(int j=0;j<latticeSize.y;j++)
			{
				for(int k=0;k<latticeSize.z;k++)
				{
					int index=i+latticeSize.x*j+latticeSize.x*latticeSize.y*k;
					
					T theta,phi;
					
					position<T> p;
					p.x=(T)i*latticeLength+0.001;
					p.y=(T)j*latticeLength+0.001;
					p.z=(T)k*latticeLength+0.001;
					p.type=solventType;
					
					threeVector<T> v,a,minImg;
					//set all acceleration to 0 initially
					a.x=0;
					a.y=0;
					a.z=0;
					
					if(!tooClose[index])
					{
						//velocity
						theta=M_PI*randNum->rand53();
						phi=M_PI*2*randNum->rand53();
						v.x=Vrms*cos(phi)*sin(theta);
						v.y=Vrms*sin(phi)*sin(theta);
						v.z=Vrms*cos(theta);
						//save it
						//This behavior is undefined under certain circumstances
						System.addParticle(p,v,a);
					}
				}
			}
		}
			
		std::cerr << "Done adding solvent.\n";
			
		delete tooClose;
	}
	return true;
}
