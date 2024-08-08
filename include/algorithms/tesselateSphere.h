#include "functions.h"

#ifndef MD_TESSELATESPHERE
#define MD_TESSELATESPHERE


//this template could actually include other things, e.g. <typename T, typename S, typename BILL, typename PHYSICS, typename WHATEVER...>
template <typename T>
class TesselateSphere {
	public:
		//global functions that should usually be defined, their input varies from class to class
		TesselateSphere();
		~TesselateSphere();
		T build();
		T compute();
		//T output(T time, char *name);
	private:
		int search_midpoint(position *vertices, int3 *edge, int &edge_walk, int &nVertices, int index_start, int index_end);
};

//searches for midpoints in a walk around a polyhedra
template <typename T>
int TesselateSphere<T>::search_midpoint(position *vertices, int3 *edge, int &edge_walk, int &nVertices, int index_start, int index_end) 
{ 
	int i;
	for (i=0; i<edge_walk; i++) 
		if ((edge[i].dim[START] == index_start && edge[i].dim[END] == index_end) || 
			(edge[i].dim[START] == index_end && edge[i].dim[END] == index_start)) 
		{
			int res = edge[i].dim[MIDDLE];

			/* update the arrays */
			edge[i].dim[START]=edge[edge_walk-1].dim[START];
			edge[i].dim[END]=edge[edge_walk-1].dim[END];
			edge[i].dim[MIDDLE]=edge[edge_walk-1].dim[MIDDLE];
			edge_walk--;
	
			return res; 
		}

	/* vertex not in the list, so we add it */
	edge[edge_walk].dim[START] = index_start;
	edge[edge_walk].dim[END] = index_end; 
	edge[edge_walk].dim[MIDDLE] = nVertices; 
  
	/* create new vertex */ 
	vertices[nVertices].p[X] = (vertices[index_start].p[X] + vertices[index_end].p[X]) / 2.0;
	vertices[nVertices].p[Y] = (vertices[index_start].p[Y] + vertices[index_end].p[Y]) / 2.0;
	vertices[nVertices].p[Z] = (vertices[index_start].p[Z] + vertices[index_end].p[Z]) / 2.0;
  
	/* normalize the new vertex */ 
	double length = sqrt (vertices[nVertices].p[X] * vertices[nVertices].p[X] +
				vertices[nVertices].p[Y] * vertices[nVertices].p[Y] +
				vertices[nVertices].p[Z] * vertices[nVertices].p[Z]);
	length = 1/length;
	vertices[nVertices].p[X] *= length;
	vertices[nVertices].p[Y] *= length;
	vertices[nVertices].p[Z] *= length;
  
  nVertices++;
  edge_walk++;
  return edge[edge_walk-1].dim[MIDDLE];
}

Create::Create(int argc, char * argv[]):System()
{
	//some old variables
	int i, j, k, f, g, h, nearest, monomerChains, inner, outer, nearestTail;
	double *Umin,*Umax,*lbond,*kbond, *kbend;
	double radius, rms, density;
	double z,dz,length,s,r,theta,phi;
	//variables for the cytoskeleton chain
	int xa,xb,xc,ya,yb,yc,ab,bc,ca,n_a_new,n_faces_new,edge_walk,ab_midpoint,bc_midpoint,ca_midpoint;
	int nVertices = 12, nEdges=30, nFaces = 20;
	int3 *edge=NULL, *faces=NULL, *faces_old=NULL;
	int monomers=atoi(argv[9]);
	int depth=atoi(argv[10]);

	this->m=NULL;
	this->a=NULL;
	position *lip=NULL, *anc=NULL;//I really can't think of a better way to sort the particle list
	int lipids, anchors;//size relations for lip, anc, and mon.
	//lipids*CHAINLENGTH, anchors*(CHAINLENGTH*2+1), monomers*(nFaces+nVertices-2)

	//set up parameters
	this->name=argv[1];
	this->seed=atoi(argv[2]);
	this->ti=0.0;
	this->molecules=0;//the lipid chains for now
	lipids=atoi(argv[3]);
	this->duration=atof(argv[4]);
	this->deltat=atof(argv[5]);
	this->store=atof(argv[6]);
	this->measure=atof(argv[7]);
	density=sqrt((atoi(argv[3]))/atof(argv[8]));
	this->size[X]=sqrt((atoi(argv[3]))/atof(argv[8]))*5.0;
	this->size[Y]=this->size[X];
	this->size[Z]=this->size[X];
	this->tempI=atof(argv[11]);
	this->tempF=atof(argv[12]);
	this->periodic=PERIODIC;
	this->cutoff=CUTOFF;
	this->types=TYPES;
	this->ntwoBodyFconst=TWOBODYFCONST;
	this->ntwoBodyUconst=TWOBODYUCONST;
	this->ntwoBodyMFconst=TWOBODYMFCONST;
	this->ntwoBodyMUconst=TWOBODYMUCONST;
	this->nthreeBodyMFconst=THREEBODYMFCONST;
	this->nthreeBodyMUconst=THREEBODYMUCONST;
	
	this->randNum=new MTRand(this->seed);
	this->twoBodyFconst=(double *)malloc(this->types*this->types*this->ntwoBodyFconst*sizeof(double));
	this->twoBodyUconst=(double *)malloc(this->types*this->types*this->ntwoBodyUconst*sizeof(double));
	this->twoBodyMFconst=(double *)malloc(this->types*this->types*this->ntwoBodyMFconst*sizeof(double));
	this->twoBodyMUconst=(double *)malloc(this->types*this->types*this->ntwoBodyMUconst*sizeof(double));
	this->threeBodyMFconst=(double *)malloc(this->types*this->types*this->types*this->nthreeBodyMFconst*sizeof(double));
	this->threeBodyMUconst=(double *)malloc(this->types*this->types*this->types*this->nthreeBodyMUconst*sizeof(double));

	if(this->twoBodyFconst==NULL || this->twoBodyUconst==NULL || this->twoBodyMFconst==NULL || this->twoBodyMUconst==NULL || this->threeBodyMFconst==NULL || this->threeBodyMUconst==NULL)
	{
		cout << "No more memory.\n";
		throw "No more memory\n";
	}
	
	//make lipid list
	lip=(position *)malloc(lipids*CHAINLENGTH*sizeof(position));
	if(lip==NULL)
	{
		cout << "1:Not enough memory!" << endl;
		throw "No more memory\n";
	}

	radius=sqrt((density*density)/(4.0*M_PI));//density adjusted

	//with inner and outer surface area compensation
	inner=((4.0*M_PI*radius*radius)/((4.0*M_PI*radius*radius)+(4.0*M_PI*(radius+5.0*1.0)*(radius+5.0*1.0))))*lipids;
	outer=lipids-inner;
	s=3.6/sqrt(double(inner));
	length=0;
	dz=2.0/double(inner);
	z=1.0-dz/2.0;
	
	for(i=0;i<inner*CHAINLENGTH;i=i+CHAINLENGTH)
	{
		r=sqrt(1.0-z*z);
		//inner monomers to outer monomers
		lip[i].p[X]=r*cos(length)*radius;
		lip[i].p[Y]=r*sin(length)*radius;
		lip[i].p[Z]=z*radius;
		lip[i].type=HEAD;

		lip[i+1].p[X]=r*cos(length)*(radius+1.0);
		lip[i+1].p[Y]=r*sin(length)*(radius+1.0);
		lip[i+1].p[Z]=z*(radius+1.0);
		lip[i+1].type=TAIL;
			
		lip[i+2].p[X]=r*cos(length)*(radius+1.0*2.0);
		lip[i+2].p[Y]=r*sin(length)*(radius+1.0*2.0);
		lip[i+2].p[Z]=z*(radius+1.0*2.0);
		lip[i+2].type=TAIL;
		z=z-dz;
		length=length+s/r;
	}

	
	s=3.6/sqrt(double(outer));
	length=0;
	dz=2.0/double(outer);
	z=1.0-dz/2.0;
	
	for(i=inner*CHAINLENGTH;i<(outer+inner)*CHAINLENGTH;i=i+CHAINLENGTH)
	{
		r=sqrt(1.0-z*z);
		//outer monomers to inner monomers
		lip[i+2].p[X]=r*cos(length)*(radius+1.0*5.0);
		lip[i+2].p[Y]=r*sin(length)*(radius+1.0*5.0);
		lip[i+2].p[Z]=z*(radius+1.0*5.0);
		lip[i+2].type=HEAD;

		lip[i+1].p[X]=r*cos(length)*(radius+1.0*4.0);
		lip[i+1].p[Y]=r*sin(length)*(radius+1.0*4.0);
		lip[i+1].p[Z]=z*(radius+1.0*4.0);
		lip[i+1].type=TAIL;
		
		lip[i].p[X]=r*cos(length)*(radius+1.0*3.0);
		lip[i].p[Y]=r*sin(length)*(radius+1.0*3.0);
		lip[i].p[Z]=z*(radius+1.0*3.0);
		lip[i].type=TAIL;

		z=z-dz;
		length=length+s/r;
	}
	
	if(monomers!=0)
	{
		//cytoskeleton
		//these are normalized
		double t = (1+(double)sqrt(5.0))/2, tau = t/(double)sqrt(1.0+t*t), one = 1/(double)sqrt(1.0+t*t);
		anc = (position *)malloc(nVertices*sizeof(position));
		if(anc==NULL)
		{
			cout << "2:Not enough memory" << endl;
			throw "Not enough memory";
		}

		//vertices of an icosahedron
		anc[0].p[X]=tau;
		anc[0].p[Y]=one;
		anc[0].p[Z]=0.0;
		anc[1].p[X]=-tau;
		anc[1].p[Y]=one;
		anc[1].p[Z]=0.0;
		anc[2].p[X]=-tau;
		anc[2].p[Y]=-one;
		anc[2].p[Z]=0.0;
		anc[3].p[X]=tau;
		anc[3].p[Y]=-one;
		anc[3].p[Z]=0.0;
		anc[4].p[X]=one;
		anc[4].p[Y]=0.0;
		anc[4].p[Z]=tau;
		anc[5].p[X]=one;
		anc[5].p[Y]=0.0;
		anc[5].p[Z]=-tau;
		anc[6].p[X]=-one;
		anc[6].p[Y]=0.0;
		anc[6].p[Z]=-tau;
		anc[7].p[X]=-one;
		anc[7].p[Y]=0.0;
		anc[7].p[Z]=tau;
		anc[8].p[X]=0.0;
		anc[8].p[Y]=tau;
		anc[8].p[Z]=one;
		anc[9].p[X]=0.0;
		anc[9].p[Y]=-tau;
		anc[9].p[Z]=one;
		anc[10].p[X]=0.0;
		anc[10].p[Y]=-tau;
		anc[10].p[Z]=-one;
		anc[11].p[X]=0.0;
		anc[11].p[Y]=tau;
		anc[11].p[Z]=-one;
	
		faces = (int3*)malloc(nFaces*sizeof(int3));
	
		faces[0].dim[0]=4;
		faces[0].dim[1]=8;
		faces[0].dim[2]=7;
		faces[1].dim[0]=4;
		faces[1].dim[1]=7;
		faces[1].dim[2]=9;
		faces[2].dim[0]=5;
		faces[2].dim[1]=6;
		faces[2].dim[2]=11;	
		faces[3].dim[0]=5;
		faces[3].dim[1]=10;
		faces[3].dim[2]=6;
		faces[4].dim[0]=0;
		faces[4].dim[1]=4;
		faces[4].dim[2]=3;
		faces[5].dim[0]=0;
		faces[5].dim[1]=3;
		faces[5].dim[2]=5;
		faces[6].dim[0]=2;
		faces[6].dim[1]=7;
		faces[6].dim[2]=1;
		faces[7].dim[0]=2;
		faces[7].dim[1]=1;
		faces[7].dim[2]=6;
		faces[8].dim[0]=8;
		faces[8].dim[1]=0;
		faces[8].dim[2]=11;
		faces[9].dim[0]=8;
		faces[9].dim[1]=11;
		faces[9].dim[2]=1;
		faces[10].dim[0]=9;
		faces[10].dim[1]=10;
		faces[10].dim[2]=3;
		faces[11].dim[0]=9;
		faces[11].dim[1]=2;
		faces[11].dim[2]=10;
		faces[12].dim[0]=8;
		faces[12].dim[1]=4;
		faces[12].dim[2]=0;
		faces[13].dim[0]=11;
		faces[13].dim[1]=0;
		faces[13].dim[2]=5;
		faces[14].dim[0]=4;
		faces[14].dim[1]=9;
		faces[14].dim[2]=3;
		faces[15].dim[0]=5;
		faces[15].dim[1]=3;
		faces[15].dim[2]=10;
		faces[16].dim[0]=7;
		faces[16].dim[1]=8;
		faces[16].dim[2]=1;
		faces[17].dim[0]=6;
		faces[17].dim[1]=1;
		faces[17].dim[2]=11;
		faces[18].dim[0]=7;
		faces[18].dim[1]=2;
		faces[18].dim[2]=9;
		faces[19].dim[0]=6;
		faces[19].dim[1]=10;
		faces[19].dim[2]=2;
	
		for (i=0; i<depth; i++)
		{
			n_a_new = nVertices+2*nEdges; 
			n_faces_new = 4*nFaces; 
	
			edge_walk = 0; 
			nEdges = 2*nVertices + 3*nFaces; 
			edge = (int3 *)malloc(nEdges*sizeof(int3)); 
			for(j=0;j<nEdges;j++)
			{
			edge[j].dim[0]=-1;
			edge[j].dim[1]=-1;
			edge[j].dim[2]=-1;
			}
		faces_old = (int3*)malloc (nFaces*sizeof(int3)); 
		faces_old = (int3*)memcpy((void*)faces_old, (void*)faces, nFaces*sizeof(int3)); 
		anc = (position *)realloc ((void*)anc, n_a_new*sizeof(position)); 
		faces = (int3*)realloc ((void*)faces, n_faces_new*sizeof(int3)); 
		n_faces_new = 0; 

		for (j=0; j<nFaces; j++) 
		{ 
			xa = faces_old[j].dim[0]; 
			xb = faces_old[j].dim[1]; 
			xc = faces_old[j].dim[2]; 

			ab_midpoint = search_midpoint (anc, edge, edge_walk, nVertices, xb, xa); 
			bc_midpoint = search_midpoint (anc, edge, edge_walk, nVertices, xc, xb); 
			ca_midpoint = search_midpoint (anc, edge, edge_walk, nVertices, xa, xc); 

			faces[n_faces_new].dim[0] = xa; 
			faces[n_faces_new].dim[1] = ab_midpoint; 
			faces[n_faces_new].dim[2] = ca_midpoint; 
			n_faces_new++; 
			faces[n_faces_new].dim[0] = ca_midpoint; 
			faces[n_faces_new].dim[1] = ab_midpoint; 
			faces[n_faces_new].dim[2] = bc_midpoint; 
			n_faces_new++; 
			faces[n_faces_new].dim[0] = ca_midpoint; 
			faces[n_faces_new].dim[1] = bc_midpoint; 
			faces[n_faces_new].dim[2] = xc; 
			n_faces_new++; 
			faces[n_faces_new].dim[0] = ab_midpoint; 
			faces[n_faces_new].dim[1] = xb; 
			faces[n_faces_new].dim[2] = bc_midpoint; 
			n_faces_new++; 
			} 
			nFaces = n_faces_new;
			free(edge);
			free(faces_old);
		}

	
		//expand the radius, assumes it is normalized
		for(i=0;i<nVertices;i++)
		{
			anc[i].p[X]*=(radius-1);
			anc[i].p[Y]*=(radius-1);
			anc[i].p[Z]*=(radius-1);
		}


		//find nearest lipids to anchors and shuffle them from start 
		anc = (position *)realloc ((void*)anc, nVertices*(CHAINLENGTH*2+1)*sizeof(position)); 

		for(i=0;i<nVertices;i++)
		{	
			//locate nearest point to i, this assumes the head is always closer than the tail
			nearest=0;
			for(j=1,nearest=0;j<inner*CHAINLENGTH;j++)
			{
				if(((lip[nearest].p[X]-anc[i].p[X])*(lip[nearest].p[X]-anc[i].p[X])+
					(lip[nearest].p[Y]-anc[i].p[Y])*(lip[nearest].p[Y]-anc[i].p[Y])+
					(lip[nearest].p[Z]-anc[i].p[Z])*(lip[nearest].p[Z]-anc[i].p[Z]))>
					((lip[j].p[X]-anc[i].p[X])*(lip[j].p[X]-anc[i].p[X])+
					(lip[j].p[Y]-anc[i].p[Y])*(lip[j].p[Y]-anc[i].p[Y])+
					(lip[j].p[Z]-anc[i].p[Z])*(lip[j].p[Z]-anc[i].p[Z])))
				{
					nearest=j;
				}
			}
			for(j=0;j<CHAINLENGTH;j++)
			{
				anc[nVertices+(CHAINLENGTH*2*i)+j].p[X]=lip[nearest+j].p[X];
				anc[nVertices+(CHAINLENGTH*2*i)+j].p[Y]=lip[nearest+j].p[Y];
				anc[nVertices+(CHAINLENGTH*2*i)+j].p[Z]=lip[nearest+j].p[Z];
				if(lip[nearest+j].type==TAIL)//convert types
					anc[nVertices+(CHAINLENGTH*2*i)+j].type=TAILANCH;
				if(lip[nearest+j].type==HEAD)
					anc[nVertices+(CHAINLENGTH*2*i)+j].type=HEADANCH;
				lip[nearest+j].type=-1;//mark the particle for deletion
			}

			//locate nearest point to nearest, this assumes the tail is always closer than the head
			nearest=nVertices+(CHAINLENGTH*2*i)+CHAINLENGTH-1;
			nearestTail=inner*CHAINLENGTH;
			for(j=nearestTail+1;j<(inner+outer)*CHAINLENGTH;j++)
			{
				if(((lip[nearestTail].p[X]-anc[nearest].p[X])*(lip[nearestTail].p[X]-anc[nearest].p[X])+
					(lip[nearestTail].p[Y]-anc[nearest].p[Y])*(lip[nearestTail].p[Y]-anc[nearest].p[Y])+
					(lip[nearestTail].p[Z]-anc[nearest].p[Z])*(lip[nearestTail].p[Z]-anc[nearest].p[Z]))>
					((lip[j].p[X]-anc[nearest].p[X])*(lip[j].p[X]-anc[nearest].p[X])+
					(lip[j].p[Y]-anc[nearest].p[Y])*(lip[j].p[Y]-anc[nearest].p[Y])+
					(lip[j].p[Z]-anc[nearest].p[Z])*(lip[j].p[Z]-anc[nearest].p[Z])))
				{
					nearestTail=j;
				}
			}

			for(j=0;j<CHAINLENGTH;j++)
			{
				anc[nVertices+(CHAINLENGTH*2*i)+j+CHAINLENGTH].p[X]=lip[nearestTail+j].p[X];
				anc[nVertices+(CHAINLENGTH*2*i)+j+CHAINLENGTH].p[Y]=lip[nearestTail+j].p[Y];
				anc[nVertices+(CHAINLENGTH*2*i)+j+CHAINLENGTH].p[Z]=lip[nearestTail+j].p[Z];
				if(lip[nearestTail+j].type==TAIL)//convert types
					anc[nVertices+(CHAINLENGTH*2*i)+j+CHAINLENGTH].type=TAILANCH;
				if(lip[nearestTail+j].type==HEAD)
					anc[nVertices+(CHAINLENGTH*2*i)+j+CHAINLENGTH].type=HEADANCH;
				lip[nearestTail+j].type=-1;//mark the particle for deletion
			}
		}

	anchors=nVertices;
	nVertices+=(nVertices*CHAINLENGTH*2);
	monomerChains=0;
	//place monomers between
	anc = (position *)realloc ((void*)anc, (nVertices+monomers*((nFaces+nVertices)-2))*sizeof(position));
	for(i=0;i<nFaces;i++)
	{
		xa=faces[i].dim[0];
		xb=faces[i].dim[1];
		xc=faces[i].dim[2];
		ab=1;
		bc=1;
		ca=1;
		//check if ab pair exists on other faces
		for(j=0;j<i;j++)
		{
			ya=faces[j].dim[0];
			yb=faces[j].dim[1];
			yc=faces[j].dim[2];

			if((xa==ya || xa==yb || xa==yc) && (xb==ya || xb==yb || xb==yc))
				ab=0;
			if((xb==ya || xb==yb || xb==yc) && (xc==ya || xc==yb || xc==yc))
				bc=0;
			if((xc==ya || xc==yb || xc==yc) && (xa==ya || xa==yb || xa==yc))
				ca=0;
		}
		//place monomers between anchors
		if(ab)
		{
			monomerChains++;
			for(j=1;j<monomers+1;j++)
			{
				anc[nVertices].p[X]=j*((anc[xa].p[X]-anc[xb].p[X])/(monomers+1))+anc[xb].p[X];
				anc[nVertices].p[Y]=j*((anc[xa].p[Y]-anc[xb].p[Y])/(monomers+1))+anc[xb].p[Y];
				anc[nVertices].p[Z]=j*((anc[xa].p[Z]-anc[xb].p[Z])/(monomers+1))+anc[xb].p[Z];
				nVertices++;
			}
			
			if(monomers>2)
			{

					this->m=(molecule *)realloc ((void*)this->m, (this->molecules+1)*sizeof(molecule));
					this->m[this->molecules].type=BOND;
					this->m[this->molecules].info[0]=xa*(CHAINLENGTH*2+1);//normally every 7
					this->m[this->molecules].info[1]=nVertices-1;
					this->m[this->molecules].info[2]=0;//nothing on this part of BOND
					this->molecules++;

					this->m=(molecule *)realloc ((void*)this->m, (this->molecules+1)*sizeof(molecule));
					this->m[this->molecules].type=BOND;
					this->m[this->molecules].info[0]=xb*(CHAINLENGTH*2+1);//normally every 7
					this->m[this->molecules].info[1]=nVertices-monomers;
					this->m[this->molecules].info[2]=0;//nothing on this part of BOND
					this->molecules++;
			}
			else
			{
				cout << "3:Not enough monomers!" << endl;
				throw "Not enough monomers!";
			}
			
			
			}
			if(bc)
			{
				monomerChains++;
				for(j=1;j<monomers+1;j++)
				{
					anc[nVertices].p[X]=j*((anc[xb].p[X]-anc[xc].p[X])/(monomers+1))+anc[xc].p[X];
					anc[nVertices].p[Y]=j*((anc[xb].p[Y]-anc[xc].p[Y])/(monomers+1))+anc[xc].p[Y];
					anc[nVertices].p[Z]=j*((anc[xb].p[Z]-anc[xc].p[Z])/(monomers+1))+anc[xc].p[Z];
					nVertices++;
				}
			
				if(monomers>2)
				{

					this->m=(molecule *)realloc ((void*)this->m, (this->molecules+1)*sizeof(molecule));
					this->m[this->molecules].type=BOND;
					this->m[this->molecules].info[0]=xb*(CHAINLENGTH*2+1);//normally every 7
					this->m[this->molecules].info[1]=nVertices-1;
					this->m[this->molecules].info[2]=0;//nothing on this part of BOND
					this->molecules++;

					this->m=(molecule *)realloc ((void*)this->m, (this->molecules+1)*sizeof(molecule));
					this->m[this->molecules].type=BOND;
					this->m[this->molecules].info[0]=xc*(CHAINLENGTH*2+1);//normally every 7
					this->m[this->molecules].info[1]=nVertices-monomers;
					this->m[this->molecules].info[2]=0;//nothing on this part of BOND
					this->molecules++;
				}
				else
				{
					cout << "4:Not enough monomers!" << endl;
					throw "Not enough monomers!";
				}
			}
			if(ca)
			{
				monomerChains++;
				for(j=1;j<monomers+1;j++)
				{
					anc[nVertices].p[X]=j*((anc[xc].p[X]-anc[xa].p[X])/(monomers+1))+anc[xa].p[X];
					anc[nVertices].p[Y]=j*((anc[xc].p[Y]-anc[xa].p[Y])/(monomers+1))+anc[xa].p[Y];
					anc[nVertices].p[Z]=j*((anc[xc].p[Z]-anc[xa].p[Z])/(monomers+1))+anc[xa].p[Z];
					nVertices++;
				}
				if(monomers>2)
				{
	
						this->m=(molecule *)realloc ((void*)this->m, (this->molecules+1)*sizeof(molecule));
						this->m[this->molecules].type=BOND;
						this->m[this->molecules].info[0]=xc*(CHAINLENGTH*2+1);//normally every 7
						this->m[this->molecules].info[1]=nVertices-1;
						this->m[this->molecules].info[2]=0;//nothing on this part of BOND
						this->molecules++;
	
						this->m=(molecule *)realloc ((void*)this->m, (this->molecules+1)*sizeof(molecule));
						this->m[this->molecules].type=BOND;
						this->m[this->molecules].info[0]=xa*(CHAINLENGTH*2+1);//normally every 7
						this->m[this->molecules].info[1]=nVertices-monomers;
						this->m[this->molecules].info[2]=0;//nothing on this part of BOND
						this->molecules++;
				}
				else
				{
					cout << "5:Not enough monomers!" << endl;
					throw "Not enough monomers!";
				}
			}
		}

		for(i=0;i<anchors;i++)
			anc[i].type=CYTOANCH;
		for(i=anchors*(1+CHAINLENGTH*2);i<nVertices;i++)
			anc[i].type=CYTO;
	}
	else
	{
		nVertices=0;
		anchors=0;
		monomers=0;
		monomerChains=0;
	}
	this->particles=0;
	p=(position *)malloc((nVertices+lipids*CHAINLENGTH)*sizeof(position));
	if(p==NULL)
	{
		cout << "6:Not enough memory!" << endl;
		throw "Not enough memory!";
	}
	if(monomers!=0)
	{
		//put each vertice/anchor at the beginning of each molecule
		this->m=(molecule *)realloc((void *)this->m, (this->molecules+1)*sizeof(molecule));
		this->m[this->molecules].info[START]=0;
		this->m[this->molecules].info[LENGTH]=1;//+(CHAINLENGTH*2);
		this->m[this->molecules].info[CHAINS]=anchors;
		this->m[this->molecules].type=LIST;
		this->molecules++;
	}

	for(i=0;i<anchors;i++)
	{
		p[this->particles].p[X]=anc[i].p[X];
		p[this->particles].p[Y]=anc[i].p[Y];
		p[this->particles].p[Z]=anc[i].p[Z];
		p[this->particles].type=anc[i].type;
		this->particles++;
		
		for(j=0;j<CHAINLENGTH*2;j++)
		{
			p[this->particles].p[X]=anc[anchors+(CHAINLENGTH*2*i)+j].p[X];
			p[this->particles].p[Y]=anc[anchors+(CHAINLENGTH*2*i)+j].p[Y];
			p[this->particles].p[Z]=anc[anchors+(CHAINLENGTH*2*i)+j].p[Z];
			p[this->particles].type=anc[anchors+(CHAINLENGTH*2*i)+j].type;
			this->particles++;
		}
		
	}

	if(monomers!=0)
	{
		//start placing monomers in lists
		this->m=(molecule *)realloc((void*)this->m, (this->molecules+1)*sizeof(molecule));
		this->m[this->molecules].info[START]=this->particles;
		this->m[this->molecules].info[LENGTH]=monomers;
		this->m[this->molecules].info[CHAINS]=monomerChains;//I don't know how many edges there are...
		this->m[this->molecules].type=LIST;
		this->molecules++;
	}

	for(i=anchors*(1+CHAINLENGTH*2);i<nVertices;i++)//put the monomers at the end of the anchor list
	{
		p[this->particles].p[X]=anc[i].p[X];
		p[this->particles].p[Y]=anc[i].p[Y];
		p[this->particles].p[Z]=anc[i].p[Z];
		p[this->particles].type=anc[i].type;
		this->particles++;
	}
		
	//start placing lipids in lists
	this->m=(molecule *)realloc((void*)this->m, (this->molecules+1)*sizeof(molecule));
	this->m[this->molecules].info[START]=this->particles;
	this->m[this->molecules].info[LENGTH]=CHAINLENGTH;
	this->m[this->molecules].info[CHAINS]=lipids-(anchors*2);//adjust for lipids that are added to anchors
	this->m[this->molecules].type=LIST;
	this->molecules++;
	
	for(i=0;i<lipids*CHAINLENGTH;i++)//put the lipids at the end
	{
		if(lip[i].type!=-1)
		{
			p[this->particles].p[X]=lip[i].p[X];
			p[this->particles].p[Y]=lip[i].p[Y];
			p[this->particles].p[Z]=lip[i].p[Z];
			p[this->particles].type=lip[i].type;
			this->particles++;
		}
	}
	
	if(lip)
		free(lip);
	if(anc)
		free(anc);
	if (faces)
		free(faces); 

	
	//center this at center of system rather than origin
	for(i=0;i<this->particles;i++)
	{
		p[i].p[X]=p[i].p[X]+this->size[X]/2;
		p[i].p[Y]=p[i].p[Y]+this->size[Y]/2;
		p[i].p[Z]=p[i].p[Z]+this->size[Z]/2;
	}
	
	//two body interactions
	this->m=(molecule *)realloc((void*)this->m, (this->molecules+1)*sizeof(molecule));
	this->m[this->molecules].info[PARTICLES]=this->particles;//number of particles in this interaction
	this->m[this->molecules].info[INTERACTING]=((TAIL+1)<<1) | ((HEAD+1)<<1) | ((CYTO+1)<<1);//| any other types
	this->m[this->molecules].type=TWOBODY;
	this->molecules++;

	//generate force constants
	Umin=new double[this->types*this->types];
	Umax=new double[this->types*this->types];
	lbond=new double[this->types*this->types];
	kbond=new double[this->types*this->types*this->types];
	kbend=new double[this->types*this->types*this->types];
	
	if(Umin==0 || Umax==0 || lbond==0 || kbond==0 || kbend==0)
	{
		cout << "Out of memory\n";
		throw "Out of memory\n";
	}

	//could be one loop, but it's good for an example
	//of how to set constants
	for(i=0;i<this->types;i++)
	{
		for(j=0;j<this->types;j++)
		{
			//i=first type, j=second type, indexed grid
			g=h=k=i+j*this->types;
			Umin[g]=0;
			Umax[h]=100;
			lbond[k]=0.7;
			kbond[g]=100;
		}
	}
	
	//umin and umax exceptions
	Umin[TAIL+TAIL*this->types]=-6;
	Umax[TAIL+TAIL*this->types]=200;
	
	//TAIL anchors
	Umin[TAILANCH+TAIL*this->types]=-6.5;
	Umax[TAILANCH+TAIL*this->types]=200;
	Umin[TAIL+TAILANCH*this->types]=-6.5;
	Umax[TAIL+TAILANCH*this->types]=200;
	Umin[TAILANCH+TAILANCH*this->types]=-6.5;
	Umax[TAILANCH+TAILANCH*this->types]=200;
	
	//how to do it with one loop
	for(i=0;i<this->types*this->types*this->types;i++)
	{
		kbend[i]=100;
	}
	
	this->gamma=GAMMA;
//	this->tempI=ITEMP;
//	this->tempF=FTEMP;

	//two body constants
	for(i=0;i<this->types;i++)
	{
		for(j=0;j<this->types;j++)
		{
			g=i+j*this->types;
			h=this->ntwoBodyFconst*g;
			k=this->ntwoBodyUconst*g;
			//U[r<=rmin]=((Umax-Umin)*(rmin-r)^2/rmin^2)+Umin
			//U[rmin<r<=rc]=(-2*Umin*(rc-r)^3/(rc-rmin)^3)+(3*Umin*(rc-r)^2/(rc-rm)^2)
			//F[r<=rmin]=(-2*(Umax-Umin)/rmin^2)*(rmin-r)/r
			//F[rmin<r<=rc]=(6*Umin/(rc-rmin)^3)*(rc-r)^2/r-(6*Umin/(rc-rmin)^2)*(rc-r)/r
			
			//constants[U[rmin<r<=rc]]: 0:-2*Umin/(rc-rmin)^3 1:3*Umin/(rc-rmin)^2
			//constants[F[rmin<r<=rc]]: 2:-6*Umin/(rc-rmin)^3  3:6*Umin/(rc-rmin)^2
			
			//constants[U[r<=rmin]]: 4:(Umax-Umin)/rmin^2  5:Umin
			//constants[F[r<=rmin]]: 6:2*(Umax-Umin)/rmin^2
			//constants[general]: 7:rc 8:rmin 9:rc^2 10:rmin^2
			
			//F constants, force constants
			this->twoBodyFconst[h+0]=RMIN;//C8
			this->twoBodyFconst[h+1]=(2.0*(Umax[g]-Umin[g]))/(RMIN*RMIN);//C6
			this->twoBodyFconst[h+2]=0;//part of index trick
			this->twoBodyFconst[h+3]=CUTOFF;//C7
			this->twoBodyFconst[h+4]=(6.0*Umin[g])/((CUTOFF-RMIN)*(CUTOFF-RMIN));//C3
			this->twoBodyFconst[h+5]=(6.0*Umin[g])/((CUTOFF-RMIN)*(CUTOFF-RMIN)*(CUTOFF-RMIN));//C2
			
			//U constants, potential constants
			this->twoBodyUconst[k+0]=RMIN;//C8
			this->twoBodyUconst[k+1]=(Umax[g]-Umin[g])/(RMIN*RMIN);//C4
			this->twoBodyUconst[k+2]=Umin[g];//C5,no index trick
			this->twoBodyUconst[k+3]=CUTOFF;//C7
			this->twoBodyUconst[k+4]=(3.0*Umin[g])/((CUTOFF-RMIN)*(CUTOFF-RMIN));//C1
			this->twoBodyUconst[k+5]=(2.0*Umin[g])/((CUTOFF-RMIN)*(CUTOFF-RMIN)*(CUTOFF-RMIN));//C0
		}
	}

	//bond interactions for three body constants
	for(i=0;i<this->types;i++)
	{
		for(j=0;j<this->types;j++)
		{
			g=i+j*this->types;
			h=this->ntwoBodyMFconst*g;
			k=this->ntwoBodyMUconst*g;
			//U=0.5*kbond*(distance-lbond)^2
			//constants: 0:0.5*kbond 1:-kbond 2:-lbond
			//F=-kbond*(distance-lbond)
			//constants: 1:-kbond 2:-lbond
			this->twoBodyMFconst[h+0]=-lbond[g];//C2
			this->twoBodyMFconst[h+1]=-kbond[g];//C1
			
			this->twoBodyMUconst[k+0]=-lbond[g];//C2
			this->twoBodyMUconst[k+1]=kbond[g]*0.5;//C0
		}
	}

	//three body interactions
	for(i=0;i<this->types;i++)
	{
		for(j=0;j<this->types;j++)
		{
			for(k=0;k<this->types;k++)
			{
				g=i+j*this->types+k*this->types*this->types;
				h=this->nthreeBodyMFconst*g;
				f=this->nthreeBodyMUconst*g;
				//U=0.5*kbend*(1+cos(theta))^2
				//constants: 0:kbend
				//F=?!?!?!?!?!?! Total blasphomy!  I stole this part!  Bad Eric! Look at this in more depth.
				//constants: 0:kbend
				this->threeBodyMFconst[h]=kbend[g];
				this->threeBodyMUconst[f]=kbend[g];
			}
		}
	}

	this->v=new velocity[this->particles];
	
	//Velocities
	rms=sqrt(3*this->tempI);
	for(i=0;i<this->particles;i++)
	{
		theta=M_PI*this->randNum->rand53();
		phi=M_PI*2*this->randNum->rand53();
		this->v[i].v[X]=rms*cos(phi)*sin(theta);
		this->v[i].v[Y]=rms*sin(phi)*sin(theta);
		this->v[i].v[Z]=rms*cos(theta);
	}
	
	delete Umin;
	delete Umax;
	delete lbond;
	delete kbond;
	delete kbend;
}

Create::~Create()
{
	delete v;
	free(p);
	free(m);
	free(twoBodyFconst);
	free(twoBodyUconst);
	free(twoBodyMFconst);
	free(twoBodyMUconst);
	free(threeBodyMFconst);
	free(threeBodyMUconst);
}
