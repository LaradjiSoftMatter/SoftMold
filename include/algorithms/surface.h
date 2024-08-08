//This algorithm is used to find and correlate volumes from different configurations

#ifndef SURFACE_MD
#define SURFACE_MD

#include "cellAlgorithm.h"

#ifndef BSPLINE_ORDER

template <typename T>
class SurfaceExtraction {
	public:
		SurfaceExtraction();
		//cutoff may need to be quite large to accomodate thickness
		SurfaceExtraction(position<T> *particles, int nParticles, threeVector<double> size, 
				 T cutoff, int surfaceType);
		~SurfaceExtraction();
		
		void build();//build a surfaceMap from scratch
		void correlate();//correlate previous surfaceMap to a new surfaceMap
		
		//exports surface map to VMD's XYZ
		void exportMap(const char *name, std::ios_base::openmode mode);
		int readNSurfaces(){return nSurfaces;};//returns 0 if build hasn't been called
		
		int readSurfaceIndex(int index)//figures out which volume the position index belongs to, -1 being exclusions
		{
			#ifdef ERRORS_ENABLED
				if(index>nP)
				{
					std::cout << "Error(SurfaceExtraction): Index is out of bounds!\n";
					throw 0;
				}
			#endif
			return surfaceIndex[index];
		}
		
		T readSurfaceArea(int surface)
		{
			#ifdef ERRORS_ENABLED
				if(surface>=nSurfaces || surface<0)
				{
					std::cout << "Error(SurfaceExtraction): Surface index is out of bounds!\n";
					throw 0;
				}
			#endif
			return surfaceSize[surface];
		}
		
	private:
		T rc;//cutoff
		fourVector<int> nCells;//number of volume elements
		threeVector<T> cellSize;//size of each volume element
		<T> *surfaceMap;//grid and spacial map hash
		<T> *oldSurfaceMap;//previous surfaceMap
		
		position<T> *p;//particles
		int nP;//number of particles
		threeVector<T> s;//size
		int *eT;//excluded types
		int nET;//number of excluded types
		
		threeVector<int> *surfaceIndex;//whatever volume a particular particle belongs to
		int nVolumes;//number of volumes
		T *volumeSize;
		T *surfaceArea;
		int outer;
		int nInner;
		int nOuterMap;
		
		int *stack;//stack for volume exploration
		MTRand *randNum;
		int *indices;
};

template <typename T>
SurfaceExtraction<T>::SurfaceExtraction()
{
	p=NULL;
	nP=0;
	s.x=0;
	s.y=0;
	s.z=0;
	rc=0;
	nET=0;
	eT=NULL;
	nCells.x=0;
	nCells.y=0;
	nCells.z=0;
	nCells.t=0;
	cellSize.x=0;
	cellSize.y=0;
	cellSize.z=0;
	volumeIndex=NULL;
	volumeMap=NULL;
	oldVolumeMap=NULL;
	stack=NULL;
	volumeSize=NULL;
	surfaceArea=NULL;
	outer=0;
	nVolumes=0;
	randNum=NULL;
	indices=NULL;
}

template <typename T>
SurfaceExtraction<T>::SurfaceExtraction(position<T> *particles, int nParticles, threeVector<double> size,
				      T cutoff, int *excludeType, int nExcludeType, int seed)
{
	this->p=particles;
	this->nP=nParticles;
	this->s=size;
	this->rc=cutoff;
	this->nET=nExcludeType;
	
	this->eT=new int[nET];
	for(int i=0;i<nET;i++)
		eT[i]=excludeType[i];
	
	if(rc<=0)
	{
		std::cout << "No cutoff in Volume!\n";
		throw 0;
	}
	this->nCells.x=int(s.x/rc);
	this->nCells.y=int(s.y/rc);
	this->nCells.z=int(s.z/rc);
	this->nCells.t=nCells.x*nCells.y*nCells.z;
	
	this->cellSize.x=s.x/nCells.x;
	this->cellSize.y=s.y/nCells.y;
	this->cellSize.z=s.z/nCells.z;
	
	this->volumeIndex=new threeVector<int>[nP];
	for(int i=0;i<nP;i++)
	{
		volumeIndex[i].s[0]=-1;
		volumeIndex[i].s[1]=-1;
	}
	this->volumeMap=new fourVector<T>[nCells.t];
	for(int i=0;i<nCells.x;i++)
	{
		for(int j=0;j<nCells.y;j++)
		{
			for(int k=0;k<nCells.z;k++)
			{
				int index=i+j*nCells.x+k*nCells.x*nCells.y;
				volumeMap[index].x=cellSize.x*((double)i+0.5);
				volumeMap[index].y=cellSize.y*((double)j+0.5);
				volumeMap[index].z=cellSize.z*((double)k+0.5);
			}
		}
	}
	
	this->oldVolumeMap=new fourVector<T>[nCells.t];
	this->stack=new int[nCells.t];
	this->volumeSize=NULL;
	this->surfaceArea=NULL;
	this->outer=0;
	this->nVolumes=0;
	this->randNum=new MTRand(seed);
	this->indices=new int[nP];
}

template <typename T>
SurfaceExtraction<T>::~SurfaceExtraction()
{
	if(volumeMap!=NULL)
		delete volumeMap;
	if(oldVolumeMap!=NULL)
		delete oldVolumeMap;
	if(eT!=NULL)
		delete eT;
	if(volumeIndex!=NULL)
		delete volumeIndex;
	if(stack!=NULL)
		delete stack;
	if(randNum!=NULL)
		delete randNum;
	if(indices!=NULL)
		delete indices;
}

template <typename T>
void SurfaceExtraction<T>::build()
{
	//initialize volume map
	#pragma omp parallel for
	for(int i=0;i<nCells.t;i++)
		volumeMap[i].t=-1;
	
	T rcSqr=rc*rc;
	
	
	//exlude volumes
	for(int i=0;i<nP;i++)
	{
		//Is this an exclusive type?
		bool exclude=false;
		for(int k=0;k<nET && !exclude;k++)
			exclude=(p[i].type==eT[k]);
			
		//Get the mapped value
		int x=p[i].x/cellSize.x;
		int y=p[i].y/cellSize.y;
		int z=p[i].z/cellSize.z;
		
		int index=x+y*nCells.x+z*nCells.x*nCells.y;
		//no?
		if(exclude)
		{
			//Get the mapped value
			int x=p[i].x/cellSize.x;
			int y=p[i].y/cellSize.y;
			int z=p[i].z/cellSize.z;
			
			int index=x+y*nCells.x+z*nCells.x*nCells.y;
			
			//it is close enough
			//if(volumeMap[index].t==-1)
				volumeMap[index].t=0;
			//else
			//	volumeMap[index].t-=(volumeMap[index].t*volumeMap[index].t)
		}
		//else if(volumeMap[index
		//{
		//	
		//}
	}
	
	nVolumes=0;
	
	///This is really fast
	//locate contiguous volumes
	for(int i=0;i<nCells.t;i++)
	{
		if(volumeMap[i].t==-1)
		{
			nVolumes++;
			stack[0]=i;
			for(int stackP=0;stackP>-1;)
			{
				int x=stack[stackP]%nCells.x;
				int y=int(stack[stackP]/nCells.x)%nCells.y;
				int z=int(stack[stackP]/(nCells.x*nCells.y));
				for(int a=-1;a<2;a++)
				{
					for(int b=-1;b<2;b++)
					{
						for(int c=-1;c<2;c++)
						{
							int f=a+x+nCells.x;
							int g=b+y+nCells.y;
							int h=c+z+nCells.z;
							
							int index=f%nCells.x+g%nCells.y*nCells.x+h%nCells.z*nCells.x*nCells.y;
							
							if(volumeMap[index].t==-1 && index!=stack[stackP])
							{
								volumeMap[index].t=nVolumes;
								#ifdef ERRORS_ENABLED
									if(stackP+1>nCells.t)
									{
										std::cout << "Error(SurfaceExtraction): stackP is out of bounds!\n";
										throw 0;
									}
								#endif
								stack[stackP++]=index;
							}
						}
					}
				}
				stackP--;
			}
			
			//only found one volume, excluding it anyway
			if(volumeMap[i].t==-1)
			{
				//Found a very small amount of mixing
				//#ifdef WARNINGS_ENABLED
				//	std::cout << "Warning(SurfaceExtraction): Possible permiation.\n";
				//#endif
				volumeMap[i].t=0;
				nVolumes--;
			}
		}
	}
	
	//also include the excluded volume volume
	nVolumes++;
	
	if(volumeSize!=NULL)
		delete volumeSize;
	volumeSize=new T[nVolumes];
	
	if(surfaceArea!=NULL)
		delete surfaceArea;
	surfaceArea=new T[nVolumes];
	
	for(int i=0;i<nVolumes;i++)
	{
		surfaceArea[i]=0;
		volumeSize[i]=0;
	}
	
	
	//This sum needs to be fixed, accumulation of very small values causes inaccuracies.
	//Should be good for ~10,000,000 points though.
	for(int i=0;i<nCells.t;i++)
		if(volumeMap[i].t>=0)
			volumeSize[volumeMap[i].t]++;
	
	for(int i=0;i<nVolumes;i++)
		volumeSize[i]*=cellSize.x*cellSize.y*cellSize.z;
		
	#ifdef STAT_OUT
		std::cout << "nVolumes found: " << nVolumes << '\n';
		std::cout << "Volume sizes are:\n";
		for(int i=0;i<nVolumes;i++)
			std::cout << '\t' << volumeSize[i];
		std::cout << '\n';
	#endif
	
	//determine which volume a particle belongs to using the grid hash
	#pragma omp parallel for
	for(int i=0;i<nP;i++)
	{
		//Is this an exclusive type?
		bool exclude=false;
		for(int k=0;k<nET && !exclude;k++)
			exclude=(p[i].type==eT[k]);
			
		//yes?
		if(exclude)
		{
			volumeIndex[i].s[0]=0;
			volumeIndex[i].s[1]=0;
		}
		else
		{
			//get map index
			int x=p[i].x/cellSize.x;
			int y=p[i].y/cellSize.y;
			int z=p[i].z/cellSize.z;
			
			int index=x+y*nCells.x+z*nCells.x*nCells.y;
			
			/*
			//check nearby cells
			if(vType==0)
			{
				for(int a=-1;a<2;a++)
				{
					for(int b=-1;b<2;b++)
					{
						for(int c=-1;c<2;c++)
						{
							int f=a+x+nCells.x;
							int g=b+y+nCells.y;
							int h=c+z+nCells.z;
							
							int k=f%nCells.x+g%nCells.y*nCells.x+h%nCells.z*nCells.x*nCells.y;
							if(volumeMap[k].t!=0)
								volumeMap[index].t=volumeMap[k].t;
						}
					}
				}
			}
			*/
			
			//shift if not the same, this is related to permiation
			if(volumeIndex[i].s[0]==-1 || volumeIndex[i].s[1]==-1)
			{
				//initialize
				volumeIndex[i].s[0]=volumeMap[index].t;
				volumeIndex[i].s[1]=volumeMap[index].t;
			}				
			else if(volumeMap[index].t!=0)//if(volumeIndex[i].s[0]!=volumeMap[index].t && volumeMap[index].t!=0)
			{
				//shift
				volumeIndex[i].s[1]=volumeIndex[i].s[0];
				volumeIndex[i].s[0]=volumeMap[index].t;
			}
		}
	}
	
	int *outerCount=new int[nVolumes];
	for(int i=0;i<nVolumes;i++)
		outerCount[i]=0;
	
	
	//determine inner or outer volumes and surface area
	for(int i=0;i<nCells.t;i++)
	{
		int x=i%nCells.x;
		int y=int(i/nCells.x)%nCells.y;
		int z=int(i/(nCells.x*nCells.y));
		
		if(x==0 || x==nCells.x-1)
			outerCount[volumeMap[i].t]++;
		if((y==0 || y==nCells.y-1) && !(x==0 || x==nCells.x-1))
			outerCount[volumeMap[i].t]++;
		if((z==0 || z==nCells.z-1) && !(x==0 || x==nCells.x-1) && !(y==0 || y==nCells.y-1))
			outerCount[volumeMap[i].t]++;
		
		//y-z planes
		for(int j=-1;j<2;j++)
		{
			int f=j+x+nCells.x;
			int g=y+nCells.y;
			int h=z+nCells.z;
			
			int index=f%nCells.x+g%nCells.y*nCells.x+h%nCells.z*nCells.x*nCells.y;
			if(volumeMap[i].t!=volumeMap[index].t)
			{
				surfaceArea[volumeMap[i].t]+=cellSize.z*cellSize.y;
			}
		}
		
		//x-z planes
		for(int j=-1;j<2;j++)
		{
			int f=x+nCells.x;
			int g=j+y+nCells.y;
			int h=z+nCells.z;
			
			int index=f%nCells.x+g%nCells.y*nCells.x+h%nCells.z*nCells.x*nCells.y;
			if(volumeMap[i].t!=volumeMap[index].t)
			{
				surfaceArea[volumeMap[i].t]+=cellSize.x*cellSize.z;
			}
		}
		
		//x-y planes
		for(int j=-1;j<2;j++)
		{
			int f=x+nCells.x;
			int g=y+nCells.y;
			int h=j+z+nCells.z;
			
			int index=f%nCells.x+g%nCells.y*nCells.x+h%nCells.z*nCells.x*nCells.y;
			if(volumeMap[i].t!=volumeMap[index].t)
			{
				surfaceArea[volumeMap[i].t]+=cellSize.x*cellSize.y;
			}
		}
	}
	
	//The one with the most volume on the outer edge is most likely to be outside.
	//Otherwise the system is composed of many open volumes (non-physical for a closed object).
	outer=0;
	for(int i=0;i<nVolumes;i++)
		if(outerCount[i]>outerCount[outer])
			outer=i;
	
	#ifdef STAT_OUT
		std::cout << "Outer volume found: " << outer << '\n';
		std::cout << "Surface areas are:\n";
		for(int i=0;i<nVolumes;i++)
			std::cout << '\t' << surfaceArea[i];
		std::cout << '\n';
	#endif
	
	nOuterMap=0;
	//collect
	for(int i=0;i<nCells.t;i++)
		if(volumeMap[i].t==outer)
			stack[nOuterMap++]=i;
			
	//everyday I'm shufflin...
	nInner=0;
	for(int i=0;i<nP;i++)
		if(volumeIndex[i].s[0]!=outer && volumeIndex[i].s[0]!=0)
			indices[nInner++]=i;
		
	delete outerCount;
}

template <typename T>
void SurfaceExtraction<T>::correlate()
{
	#ifdef ERRORS_ENABLED
	if(nVolumes==0)
	{
		std::cout << "Error(SurfaceExtraction): No previous volumes to correlate!\n";
		throw 0;
	}
	#endif
	/*
	//swap volume map pointers (better than copying whole set)
	position<T> *buf=oldVolumeMap;
	oldVolumeMap=volumeMap;
	volumeMap=buf;
	
	//copy other previous data
	int nOldVolumes=nVolumes;
	T *oldVolumeSize=new T[nOldVolumes];
	for(int i=0;i<nOldVolumes;i++)
		oldVolumeSize[i]=volumeSize[i];
	
	//initialize correlation matrix
	T *corrMatrix=new T[nVolumes*nOldVolumes];
	for(int i=0;i<nVolumes*nOldVolumes;i++)
	{
		corrMatrix[i]=0;
	}
	
	//build new map
	this->build();
	
	//build correlation matrix
	for(int i=0;i<nCells.t;i++)
	{
		//all of the shared volumes add twice
		//any loss or gain volume will add once
		//any with no overlap will be 0
		corrMatrix[volumeMap[i].t+nVolumes*oldVolumeMap[i].t]+=1.0/(volumeSize[volumeMap[i].t]+oldVolumeSize[oldVolumeMap[i].t]);
		corrMatrix[oldVolumeMap[i].t+nOldVolumes*volumeMap[i].t]+=1.0/(oldVolumeSize[oldVolumeMap[i].t]+volumeSize[volumeMap[i].t]);
	}
	
	//Analyze first (highest) correlations only, but
	//all correlations past that can determine if volumes split.
	//Furthur splits become more unlikely until it is the entire volume, 
	//which is a split from first configuration.
	//This all assumes that a volume moves slightly.
	//Spontaneous volumes will consume the last known volumes.
	if(nVolumes>nOldVolumes)
	{
		int *volume=new int[nVolumes];
		for(int i=0;i<nVolumes;i++)
			volume[i]=-1;
		for(int i=0;i<nOldVolumes;i++)
		{
			int oldVolume=0;
			for(int j=0;j<nVolumes;j++)
				if(corrMatrix[j+nVolumes*oldVolume]<corrMatrix[j+nVolumes*i])
					oldVolume=j;
			volume[i]=oldVolume;
		}
		
		//last unused volume value
		int newVal=nOldVolumes;
		
		//make sure it is full
		for(int i=0;i<nVolumes;i++)
			if(volume[i]==-1)
				volume[i]=newVal++;
			
		//swap volume values
		for(int i=0;i<nVolumes;i++)
		{
			oldVolumeSize[i]=volumeSize[volume[i]];
		}
	}
	else
	{
		int *volume=new int[nOldVolumes];
		for(int i=0;i<nOldVolumes;i++)
			volume[i]=-1;
		for(int i=0;i<
	}
	
	//set new volumes
	
	delete corrMatrix,oldVolume,oldVolumeSize;
	*/
}

template<typename T>
int SurfaceExtraction<T>::grabInner()
{
	#ifdef ERRORS_ENABLED
	if(outer==0)
	{
		std::cout << "Error(SurfaceExtraction): No outer volume defined!\n";
		throw 0;
	}
	if(nVolumes<=2)
	{
		std::cout << "Error(SurfaceExtraction): No inner volumes defined!\n";
		throw 0;
	}
	#endif
	
	int choose=nInner*randNum->rand53();
	
	int buf=indices[choose];
	indices[choose]=indices[nInner-1];
	indices[nInner-1]=buf;
	
	nInner--;
	
	return indices[nInner];
}

template<typename T>
void SurfaceExtraction<T>::moveToOuter(int index)
{
	#ifdef ERRORS_ENABLED
	if(index>=nP || index<0)
	{
		std::cout << "Error(SurfaceExtraction): Index out of bounds!\n";
		throw 0;
	}
	if(outer==0)
	{
		std::cout << "Error(SurfaceExtraction): No outer volume defined!\n";
		throw 0;
	}
	#endif
	
	//is it already an outer particle?
	if(volumeIndex[index].s[0]!=outer)
	{
		//pick one position
		int choose=randNum->rand53()*nOuterMap;
		
		//set particle at that position
		p[index].x=volumeMap[stack[choose]].x;
		p[index].y=volumeMap[stack[choose]].y;
		p[index].z=volumeMap[stack[choose]].z;
		//volumeIndex[index].s[1]=volumeIndex[index].s[0];
		//volumeIndex[index].s[0]=outer;
	}
	#ifdef WARNINGS_ENABLED
	else
	{
		std::cout << "Warning(SurfaceExtraction): Already an outer particle, doing nothing.\n";
	}
	#endif
}

template <typename T>
int SurfaceExtraction<T>::nExchanged()
{
	int n=0;
	for(int i=0;i<nP;i++)
		if(volumeIndex[i].s[0]!=volumeIndex[i].s[1] && volumeIndex[i].s[0]!=0)
			n++;
	return n;
}

template <typename T>
void SurfaceExtraction<T>::exportMap(const char *name, std::ios_base::openmode mode)
{
	std::fstream outFile;
	outFile.open(name,mode);
	if(outFile.is_open())
	{
		outFile << nCells.t << '\n' << "volume" << '\n';
		for(int i=0;i<nCells.t;i++)
			outFile << volumeMap[i].t << '\t' << volumeMap[i].x << '\t' << volumeMap[i].y << '\t' << volumeMap[i].z << '\n';
		outFile.close();
	}
	else
	{
		#ifdef WARNINGS_ENABLED
			std::cout << "Warning(SurfaceExtraction): Cannot open " << name << ".\n";
		#endif
	}
}

#endif