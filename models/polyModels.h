//This is a very short and simple program that shows you how to modify a variable already present in an mpd script.
//The reason this is used is because some variables (forces, potentials, geometry, etc...) are not easily modified
// in the mpd script. Using a file that opens the mpd script, modifies, and then saves it (as below) is far easier
// than modifying the actual variable. Most of it looks like the run.cpp file, it just replaces the "run loop"
// with the modification.

#include "../include/algorithms/functions.h"
#include "../include/algorithms/cellAlgorithm.h"

#ifndef R_THRESHOLD
#define R_THRESHOLD 1.0
#endif


/**
 * \brief I wonder if this would be useful?
 * Kind of a build structure. Used for initialization of particles in the molecular or disspative
 * dynamics systems.
 */
template <typename T, typename S>
class buildOptions {
	///How many monomers, initially.
	int nMonomers;
	///Tesselate the monomers?
	int nTesselations;
	///A unit monomer.
	position<T> monomer;
	///Center of Mass.
	threeVector<T> avgMass;
	///Radius of curvature.
	T radius;
	///The constants for the molecules
	T **constants;
	int nMolecules;
	int *nConstants;
	
	
	
	///The system in which the options apply.
	S *System;
};


/**
 * \brief Builds a banana protein. Not an actual banana, just shaped like a banana.
 * Uses curvature to map a protein into a banana shape. nMonomers is how many monomers connect together.
 * monomer is a list of units that make up a monomer unit. nUnits is the number of units that make up a monomer.
 * It is assumed that the monomer is oriented at an infinite distance from the origin (flat, no curvature).
 * avgMass vector is the position of the protein's center of mass. orientation is a vector that points in the direction
 * of highest curvature (the center of a circle) and has a magnitude equal to the magnitude of the vector.
 * x-z profile (trapazoid, looking at a monomer):
 *      #------------#
 *      \           /
 *       \         /
 *        \       /
 *         #-----#
 * 
 * y-z profile (banana):
 *                ___
 *          ___---   ---___
 *    ___---               ---___
 *                   
 */
template <typename T>
T banana(Blob<T> &System, int nMonomers, position<T> *monomer, int nUnits, threeVector<T> avgMass, threeVector<T> orientation)
{
	// Error Checking
	if(orientation.x==0 && orientation.y==0 && orientation.z==0)
	{
		std::cout << "Error(banana): orientation vector is 0 length!\n";
		return 0;
	}
	
	if(nMonomers==0)
	{
		std::cout << "Error(banana): nMonomers is 0!\n";
		return 0;
	}
	
	if(nUnits==0)
	{
		std::cout << "Error(banana): nUnits is 0!\n";
		return 0;
	}
	
	if(monomer==NULL)
	{
		std::cout << "Error(banana): monomer is NULL (not allocated)!\n";
		return 0;
	}
	
	// Arrange particles and molecules nearly to equilibrium.
	fourVector<int> bond;
	molecule<T,fourVector<int> > mInner, mOuter, mIntermediate, mInnerHalf, mOuterHalf;
	
	m.setType(CHAIN);
	
	//constants for CHAIN
	m.addConstant(constants[0]);
	m.addConstant(constants[1]);
	m.addConstant(constants[2]);
	m.addConstant(constants[3]);
	
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
	T Vrms=sqrt(3.0*System.readInitialTemp());
	
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
	
	m.setType(BOND);
	
	//constants for BOND
	m.addConstant(constants[0]);
	m.addConstant(constants[1]);
	
	//Just allocate more than you expect (12 of them have 5 each, but midpoints have 6)
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
			m.addBond(bond);
			
			//bond for a bond type
			bond.s[0]=startOffset+xb;
			bond.s[1]=System.readNParticles()-nMonomers;
			m.addBond(bond);
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
			m.addBond(bond);
			
			//bond for a bond type
			bond.s[0]=startOffset+xc;
			bond.s[1]=System.readNParticles()-nMonomers;
			m.addBond(bond);
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
			m.addBond(bond);
			
			//bond for a bond type
			bond.s[0]=startOffset+xa;
			bond.s[1]=System.readNParticles()-nMonomers;
			m.addBond(bond);
		}
	}
	
	//add the vertex bonds to the system
	System.addMolecule(m);
	m.~molecule();
	
	free(faces);
	free(anc);
	return radius;
}

