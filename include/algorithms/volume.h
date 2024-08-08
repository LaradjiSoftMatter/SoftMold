//This algorithm is used to find and correlate volumes from different configurations

#ifndef VOLUME_MD
#define VOLUME_MD

#include "cellAlgorithm.h"

#define ALL_EXCLUDE 0
#define ALL_SOLVENT 1

template <typename T>
class VolumeExtraction {
	public:
		VolumeExtraction();
		VolumeExtraction(position<T> *particles, int nParticles, threeVector<double> size, 
				 T cutoff, int *excludeType, int nExcludeType, int seed);
		
		void initialize(position<T> *particles, int nParticles, threeVector<double> size, 
			T cutoff, int *excludeType, int nExcludeType, int seed);
		
		~VolumeExtraction();
		
		void build();//build a volumeMap from scratch
		void correlate();//correlate previous volumeMap to a new volumeMap
		
		void exportMap(const char *name, std::ios_base::openmode mode);
		int readNVolumes(){return nVolumes;};//returns 0 if build hasn't been called
		int readOuter(){return outer;};//returns 0 if build hasn't been called
		int grabInner();
		void moveToOuter(int index);
		int nExchanged();
		
		int readVolumeIndex(int index)//figures out which volume the position index belongs to, -1 being exclusions
		{
			#ifdef ERRORS_ENABLED
				if(index>nP)
				{
					std::cout << "Error(VolumeExtraction): Index is out of bounds!\n";
					throw 0;
				}
			#endif
			return volumeIndex[index].s[1];
		}
		
		T readVolumeSize(int volume)
		{
			#ifdef ERRORS_ENABLED
				if(volume>=nVolumes || volume<0)
				{
					std::cout << "Error(VolumeExtraction): Volume index is out of bounds!\n";
					throw 0;
				}
			#endif
			return volumeSize[volume];
		}
		
		T readSurfaceArea(int volume)
		{
			#ifdef ERRORS_ENABLED
				if(volume>=nVolumes || volume<0)
				{
					std::cout << "Error(VolumeExtraction): Volume index is out of bounds!\n";
					throw 0;
				}
			#endif
			return surfaceArea[volume];
		}
		
	private:
		T rc;//cutoff
		fourVector<int> nCells;//number of volume elements
		threeVector<T> cellSize;//size of each volume element
		position<T> *volumeMap;//grid and spacial map hash
		position<T> *oldVolumeMap;//previous volumeMap
		
		position<T> *p;//particles
		int nP;//number of particles
		threeVector<T> s;//size
		int *eT;//excluded types
		int nET;//number of excluded types
		
		threeVector<int> *volumeIndex;//whatever volume a particular particle belongs to
		threeVector<bool *> path;
		int nVolumes;//number of volumes
		int nVolumesAlloc;//number of volumes allocated
		T *volumeSize;
		T *surfaceArea;
		int outer;
		int nInner;
		std::vector<int> outerMap;
		
		MTRand *randNum;
		int *indices;
		
		//range list, hash defined by 'vRanges[i][j]=xCell+yCell*nX+zCell*nX*nY',
		//where i=xCell+yCell*nX; a bin defined by a depth queue
		//and j is just an element in the bin
		//it is a depth queue
		std::vector< std::vector< twoVector<int> > > vRanges;
		
		//hash index same as above, but it stores the direction of the local curl
		std::vector< std::vector< threeVector< T > > > vCurl;
		
		//all the excluded indices in one convienient list.
		std::vector<int> excludedIndices;
		
		//sets vPicked, we don't even care what the value is, we just increment it
		unsigned int vCurrent;
		
		//Initialized at 0, every time a volume element is picked, it sets vPicked[element] to the new value
		std::vector<int> vPicked;
};

template <typename T>
VolumeExtraction<T>::VolumeExtraction()
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
	volumeSize=NULL;
	surfaceArea=NULL;
	outer=0;
	nVolumes=0;
	nVolumesAlloc=0;
	randNum=NULL;
	indices=NULL;
	path.x=NULL;
	path.y=NULL;
	path.z=NULL;
}

template <typename T>
VolumeExtraction<T>::VolumeExtraction(position<T> *particles, int nParticles, threeVector<double> size,
				      T cutoff, int *excludeType, int nExcludeType, int seed)
{
	this->initialize(particles, nParticles, size, cutoff, excludeType, nExcludeType, seed);
}

template <typename T>
void VolumeExtraction<T>::initialize(position<T> *particles, int nParticles, threeVector<double> size,
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
	this->volumeMap=new position<T>[nCells.t];
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
	
	this->oldVolumeMap=new position<T>[nCells.t];
	this->volumeSize=NULL;
	this->surfaceArea=NULL;
	this->outer=0;
	this->nVolumes=0;
	this->nVolumesAlloc=0;
	this->randNum=new MTRand(seed);
	this->indices=new int[nP];
	this->path.x=new bool[nCells.x];
	this->path.y=new bool[nCells.y];
	this->path.z=new bool[nCells.z];
	
	for(int i=0;i<nP;i++)
	{
		//Is this an exclusive type?
		bool exclude=false;
		for(int k=0;k<nET && !exclude;k++)
			exclude=(p[i].type==eT[k]);
			
		//yes?
		if(exclude)
			excludedIndices.push_back(i);
	}
	
	vCurrent=0;//could be anything, but we need to update this every 2^32 or 2^64 build calls!
	for(int i=0;i<nCells.t;i++)
		vPicked.push_back(vCurrent);
	
	//Reserve vRanges!
	for(int i=0;i<nCells.x*nCells.y;i++)
	{
		vRanges[i].reserve(nCells.z);//can be an arbitrary value if you are running out of memory!
		vCurl[i].reserve(nCells.z);//also arbitrary, but must be same as above line
	}
}

template <typename T>
VolumeExtraction<T>::~VolumeExtraction()
{
	if(volumeMap!=NULL)
		delete volumeMap;
	if(oldVolumeMap!=NULL)
		delete oldVolumeMap;
	if(eT!=NULL)
		delete eT;
	if(volumeIndex!=NULL)
		delete volumeIndex;
	if(randNum!=NULL)
		delete randNum;
	if(indices!=NULL)
		delete indices;
	if(path.x!=NULL)
		delete path.x,path.y,path.z;
}



template <typename T>
void VolumeExtraction<T>::build()
{
	//vector<int> stack;
	std::vector< twoVector<int> > stack;
	
	//Divergence version
	
	//We've reached the end of our integer, we have to reset the vPicked list
	if(vCurrent+1<vCurrent)
	{
		vCurrent=0;
		for(int i=0;i<vPicked.size();i++)
			vPicked[i]=vCurrent;
		std::cout << "resetting vCurrent!\n";
	}
	vCurrent++;
	
	//reset vRanges
	for(int i=0;i<nCells.x*nCells.y;i++)
	{
		vRanges[i].clear();
		vCurl[i].clear();
	}
	
	//generate vRanges, d^2 
	for(int i=0;i<excludedIndices.size();i++)
	{
		//Get the mapped value
		int x=p[excludedIndices[i]].x/cellSize.x;
		int y=p[excludedIndices[i]].y/cellSize.y;
		int z=p[excludedIndices[i]].z/cellSize.z;
		
		int index=x+y*nCells.x+z*nCells.x*nCells.y;
		if(vPicked[index]!=vCurrent)
		{
			twoVector<int> surf;
			surf.x=index;//the cell index
			surf.y=-1;//which surface it belongs to
			
			//set up our surface search
			vRanges[x+y*nCells.x].push_back(surf);
			
			//flag it as current element
			vPicked[index]=vCurrent;
		}
	}
	
	
	
	//n*log(n) sort, bilayer is d^2*log(d^2), can search in log(n) or log(d^2) time
	//std::sort(vRanges.begin(),vRanges.end(),compareX< threeVector<int> >);
	

	
	//Getting the curl for all of our surfaces, we want to avoid ones with curl=0 when we surf
	//Also, this can be done in parallel
	//#pragma omp parallel for
	for(int i=0;i<vRanges.size();i++)
	{
		for(int j=0;j<vRanges[i].size();j++)
		{
			//initialize our current position's rotation vector to 0
			threeVector<T> currentP(0);
			
			//get our cell index
			int cellHash=vRanges[i][j].x;
			
			//decode our cell x, y, and z indices.
			threeVector<int> cellI;
			cellI.x=cellHash%nCells.x;
			cellI.y=int(cellHash/nCells.x)%nCells.y;
			cellI.z=int(cellHash/(nCells.x*nCells.y));
			
			//I've unrolled the loop to be much quicker, 
			//also the compiler will probably make it faster 
			//because man of these values are known
			//at compile time.
			
			int cellN;//our next cell index
			
			//Get our cellI.x-1 index:
			cellN=((cellI.x-1+nCells.x)%nCells.x)+\
				cellI.y*nCells.x+\
				cellI.z*nCells.x*nCells.y;
			if(vPicked[cellN]==vCurrent)
			{
				currentP.y-=cellSize.y;
				currentP.z-=cellSize.z;
			}
			
			//Get our cellI.x+1 index:
			cellN=((cellI.x+1+nCells.x)%nCells.x)+\
				cellI.y*nCells.x+\
				cellI.z*nCells.x*nCells.y;
			if(vPicked[cellN]==vCurrent)
			{
				currentP.y+=cellSize.y;
				currentP.z+=cellSize.z;
			}
			
			//Get our cellI.y-1 index:
			cellN=cellI.x+\
				(((cellI.y-1+nCells.y)%nCells.y)*nCells.x)+\
				cellI.z*nCells.x*nCells.y;
			if(vPicked[cellN]==vCurrent)
			{
				currentP.x-=cellSize.x;
				currentP.z+=cellSize.z;
			}
			
			//Get our cellI.y+1 index:
			cellN=cellI.x+\
				(((cellI.y+1+nCells.y)%nCells.y)*nCells.x)+\
				cellI.z*nCells.x*nCells.y;
			if(vPicked[cellN]==vCurrent)
			{
				currentP.x+=cellSize.x;
				currentP.z-=cellSize.z;
			}
			
			//Get our cellI.z-1 index:
			cellN=cellI.x+\
				cellI.y*nCells.x+\
				(((cellI.z-1+nCells.z)%nCells.z)*nCells.x*nCells.y)*nCells.x*nCells.y;
			if(vPicked[cellN]==vCurrent)
			{
				currentP.x-=cellSize.x;
				currentP.y-=cellSize.y;
			}
			
			//Get our cellI.z+1 index:
			cellN=cellI.x+\
				cellI.y*nCells.x+\
				(((cellI.z+1+nCells.z)%nCells.z)*nCells.x*nCells.y)*nCells.x*nCells.y;
			if(vPicked[cellN]==vCurrent)
			{
				currentP.x+=cellSize.x;
				currentP.y+=cellSize.y;
			}
			
			vCurl[i].push_back(currentP);
		}
	}
	
	twoVector<int> stackI;//stackIndex
	int nSurfaces=0;//number of surfaces traversed
	
	//surf the surface, to determine inner and outer pivots
	for(stackI.x=0;stackI.x<vRanges.size();stackI.x++)
	{
		for(stackI.y=0;stackI.y<vRanges[stackI.x].size();stackI.y++)
		{
			//Make sure the current surface element hasn't been flagged and that it has some curl
			if(vRanges[stackI.x][stackI.y].y==-1 && magnitude(vCurl[stackI.x][stackI.y])>0.000001)
			{
				//check current cell to see if any nearby indices only have solvent
				stack.clear();
				stack.push_back(stackI);
				threeVector<T> crossPSum(0);
				vRanges[stackI.x][stackI.y].y=nSurfaces;
				while(stack.size()!=0)
				{
					twoVector<int> cIndex=*(stack.end()-1);//our current surface
					stack.pop_back();
					
					//get our surface depth index
					twoVector<int> surfI;
					surfI.x=cIndex.x%nCells.x;
					surfI.y=static_cast<int>(cIndex.x/nCells.x)%nCells.y;
					
					//get our current cell indices
					threeVector<T> cCIndex;
					cCIndex.x=vRanges[cIndex.x][cIndex.y].x%nCells.x;
					cCIndex.y=static_cast<int>(vRanges[cIndex.x][cIndex.y].x/nCells.x)%nCells.y;
					cCIndex.z=static_cast<int>(vRanges[cIndex.x][cIndex.y].x/(nCells.x*nCells.y))%nCells.z;
					
					//check nearby surfaces for connection
					for(int a=-1;a<2;a++)
					{
						for(int b=-1;b<2;b++)
						{
							//our next surface depth index
							twoVector<int> surfN;
							surfN.x=(a+surfI.x+nCells.x)%nCells.x;
							surfN.y=(b+surfI.y+nCells.y)%nCells.y;
							
							twoVector<int> nIndex;
							nIndex.x=surfN.x+surfN.y*nCells.x;
							
							for(nIndex.y=0;nIndex.y<vRanges[nIndex.x].size();nIndex.y++)
							{
								//get the distance to the next cell index
								threeVector<T> d;
								d.x=vRanges[nIndex.x][nIndex.y].x%nCells.x-cCIndex.x;
								d.y=static_cast<int>(vRanges[nIndex.x][nIndex.y].x/nCells.x)%nCells.y-cCIndex.y;
								d.z=static_cast<int>(vRanges[nIndex.x][nIndex.y].x/(nCells.x*nCells.y))%nCells.z-cCIndex.z;
								
								//periodic boundaries
								d.x-=(d.x>nCells.x/2)?nCells.x:0;
								d.x+=(d.x<-nCells.x/2)?nCells.x:0;
								
								d.y-=(d.y>nCells.y/2)?nCells.y:0;
								d.y+=(d.y<-nCells.y/2)?nCells.y:0;
								
								d.z-=(d.z>nCells.z/2)?nCells.z:0;
								d.z+=(d.z<-nCells.z/2)?nCells.z:0;
								
								//Does the surface connect? 1.8 chosen because the only acceptable distances are:
								//sqrt(1^2+0^2+0^2)=1
								//sqrt(1^2+1^2+0^2)=1.414
								//sqrt(1^2+1^2+1^2)=1.732
								if(vRanges[nIndex.x][nIndex.y].y==-1 && d.x*d.x+d.y*d.y+d.z*d.z<1.8)
								{
									//d is left in previous form because we are only interested in the sign
									//Also, this product is simplified by the fact that only two vectors contribute
									//Although we use all three to keep computations from diverging
									//Remember one of these is definately 0.
									crossPSum.x+=d.x*vCurl[cIndex.x][cIndex.y].y*vCurl[nIndex.x][nIndex.y].z;
									crossPSum.y-=d.y*vCurl[cIndex.x][cIndex.y].x*vCurl[nIndex.x][nIndex.y].z;
									crossPSum.z+=d.z*vCurl[cIndex.x][cIndex.y].x*vCurl[nIndex.x][nIndex.y].y;
									
									
									vRanges[nIndex.x][nIndex.y].y=nSurfaces;
									stack.push_back(nIndex);//push back our next index
								}
							}
							
						}
					}
				}
				nSurfaces++;
			}
		}
	}
	
	for(stackI.x=0;stackI.x<vRanges.size();stackI.x++)
	{
		for(stackI.y=0;stackI.y<vRanges[stackI.x].size();stackI.y++)
		{
			
		}
	}
	
	//std::vector<int> stack;
	
	//initialize volume map
	#pragma omp parallel for
	for(int i=0;i<nCells.t;i++)
		volumeMap[i].type=-1;
	
	//exlude volumes
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
			//Get the mapped value
			int x=p[i].x/cellSize.x;
			int y=p[i].y/cellSize.y;
			int z=p[i].z/cellSize.z;
			
			int index=x+y*nCells.x+z*nCells.x*nCells.y;
			volumeMap[index].type=0;
		}
	}
	
	
	nVolumes=0;
	outer=-1;
	
	//This is really fast
	//locate contiguous volumes
	for(int i=0;i<nCells.t;i++)
	{
		//it hasn't been catagorized
		if(volumeMap[i].type==-1)
		{
			//reset path test
			for(int j=0;j<nCells.x;j++)
				path.x[j]=false;
			for(int j=0;j<nCells.y;j++)
				path.y[j]=false;
			for(int j=0;j<nCells.z;j++)
				path.z[j]=false;
			
			//incrmement number of volumes
			nVolumes++;
			
			//reset stack
			stack.clear();
			stack.push_back(i);
			
			//flag and search expanse
			while(stack.size()!=0)
			{
				int currentIndex=(*(stack.end()-1)).x;
				stack.pop_back();
				//this volume's position hash
				int x=currentIndex%nCells.x;
				int y=int(currentIndex/nCells.x)%nCells.y;
				int z=int(currentIndex/(nCells.x*nCells.y));
				
				//mark the current position as traversed
				path.x[x]=true;
				path.y[y]=true;
				path.z[z]=true;
				
				//check the nearby cells for traversal
				for(int j=0;j<27;j++)
				{
					int a=j%3-1;
					int b=int(j/3)%3-1;
					int c=int(j/9)-1;
					int index=((x+a+nCells.x)%nCells.x)+\
					(((y+b+nCells.y)%nCells.y)*nCells.x)+\
					(((z+c+nCells.z)%nCells.z)*nCells.x*nCells.y);
					
					//has this volume been traversed?
					if(volumeMap[index].type==-1 && index!=currentIndex)
					{
						volumeMap[index].type=nVolumes;
						//place it on the stack
						stack.push_back(index);
					}
				}
			}
			
			//check if all lines have been traversed
			threeVector<bool> isOuter(true);//true until otherwise
			for(int j=0;j<nCells.x;j++)
				isOuter.x=(path.x[j] && isOuter.x);
			for(int j=0;j<nCells.y;j++)
				isOuter.y=(path.y[j] && isOuter.y);
			for(int j=0;j<nCells.z;j++)
				isOuter.z=(path.z[j] && isOuter.z);
			
			if(isOuter.x && isOuter.y && isOuter.z)
				outer=nVolumes;
			
			//infinite periodic surfaces
			/*
			if(!isOuter.x && isOuter.y && isOuter.z)
				YZinfinite=nVolumes;
			
			if(!isOuter.y && isOuter.x && isOuter.z)
				XZinfinite=nVolumes;
			
			if(!isOuter.z && isOuter.x && isOuter.y)
				XYinfinite=nVolumes;
			*/
			
			//only found one volume, excluding it anyway
			if(volumeMap[i].type==-1)
			{
				//Found a very small amount of mixing
				//#ifdef WARNINGS_ENABLED
				//	std::cout << "Warning(VolumeExtraction): Possible permiation.\n";
				//#endif
				volumeMap[i].type=0;
				nVolumes--;
			}
		}
	}
	
	//also include the excluded volume volume
	nVolumes++;
	
	if(nVolumes>nVolumesAlloc)
	{
		nVolumesAlloc=nVolumes;
		if(volumeSize!=NULL)
			delete volumeSize;
		volumeSize=new T[nVolumes];
		
		if(surfaceArea!=NULL)
			delete surfaceArea;
		surfaceArea=new T[nVolumes];
	}
	
	for(int i=0;i<nVolumes;i++)
	{
		surfaceArea[i]=0;
		volumeSize[i]=0;
	}
	
	
	//This sum needs to be fixed, accumulation of very small values causes inaccuracies.
	//Should be good for ~10,000,000 points though.
	for(int i=0;i<nCells.t;i++)
		if(volumeMap[i].type>=0)
			volumeSize[volumeMap[i].type]++;
	
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
			//check nearby volumes if there is an exclude type mismatch
			if(volumeMap[index].type==0)
			{
				for(int j=0;j<27;j++)
				{
					int a=j%3-1;
					int b=int(j/3)%3-1;
					int c=int(j/9)-1;
					int nearbyIndex=((x+a+nCells.x)%nCells.x)+\
					(((y+b+nCells.y)%nCells.y)*nCells.x)+\
					(((z+c+nCells.z)%nCells.z)*nCells.x*nCells.y);
				}
			}
			*/
			//shift if not the same, this is related to permiation
			if(volumeIndex[i].s[0]==-1 || volumeIndex[i].s[1]==-1)
			{
				//initialize
				volumeIndex[i].s[0]=volumeMap[index].type;
				volumeIndex[i].s[1]=volumeMap[index].type;
			}				
			else if(volumeMap[index].type!=0)//if(volumeIndex[i].s[0]!=volumeMap[index].type && volumeMap[index].type!=0)
			{
				//shift
				volumeIndex[i].s[1]=volumeIndex[i].s[0];
				volumeIndex[i].s[0]=volumeMap[index].type;
			}
		}
	}
	
	//determine surface area
	for(int i=0;i<nCells.t;i++)
	{
		int x=i%nCells.x;
		int y=int(i/nCells.x)%nCells.y;
		int z=int(i/(nCells.x*nCells.y));
		
		//y-z planes
		for(int j=-1;j<2;j++)
		{
			int f=j+x+nCells.x;
			int g=y+nCells.y;
			int h=z+nCells.z;
			
			int index=f%nCells.x+g%nCells.y*nCells.x+h%nCells.z*nCells.x*nCells.y;
			if(volumeMap[i].type!=volumeMap[index].type)
			{
				surfaceArea[volumeMap[i].type]+=cellSize.z*cellSize.y;
			}
		}
		
		//x-z planes
		for(int j=-1;j<2;j++)
		{
			int f=x+nCells.x;
			int g=j+y+nCells.y;
			int h=z+nCells.z;
			
			int index=f%nCells.x+g%nCells.y*nCells.x+h%nCells.z*nCells.x*nCells.y;
			if(volumeMap[i].type!=volumeMap[index].type)
			{
				surfaceArea[volumeMap[i].type]+=cellSize.x*cellSize.z;
			}
		}
		
		//x-y planes
		for(int j=-1;j<2;j++)
		{
			int f=x+nCells.x;
			int g=y+nCells.y;
			int h=j+z+nCells.z;
			
			int index=f%nCells.x+g%nCells.y*nCells.x+h%nCells.z*nCells.x*nCells.y;
			if(volumeMap[i].type!=volumeMap[index].type)
			{
				surfaceArea[volumeMap[i].type]+=cellSize.x*cellSize.y;
			}
		}
	}
	
	#ifdef STAT_OUT
		std::cout << "Outer volume found: " << outer << '\n';
		std::cout << "Surface areas are:\n";
		for(int i=0;i<nVolumes;i++)
			std::cout << '\t' << surfaceArea[i];
		std::cout << '\n';
	#endif
	
	outerMap.clear();
	//collect
	for(int i=0;i<nCells.t;i++)
		if(volumeMap[i].type==outer)
			outerMap.push_back(i);
	
	
	//everyday I'm shufflin...
	nInner=0;
	for(int i=0;i<nP;i++)
		if(volumeIndex[i].s[0]!=outer && volumeIndex[i].s[0]!=0)
			indices[nInner++]=i;
}

template <typename T>
void VolumeExtraction<T>::correlate()
{
	#ifdef ERRORS_ENABLED
	if(nVolumes==0)
	{
		std::cout << "Error(VolumeExtraction): No previous volumes to correlate!\n";
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
		corrMatrix[volumeMap[i].type+nVolumes*oldVolumeMap[i].type]+=1.0/(volumeSize[volumeMap[i].type]+oldVolumeSize[oldVolumeMap[i].type]);
		corrMatrix[oldVolumeMap[i].type+nOldVolumes*volumeMap[i].type]+=1.0/(oldVolumeSize[oldVolumeMap[i].type]+volumeSize[volumeMap[i].type]);
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
int VolumeExtraction<T>::grabInner()
{
	#ifdef ERRORS_ENABLED
	if(outer==0)
	{
		std::cout << "Error(VolumeExtraction): No outer volume defined!\n";
		throw 0;
	}
	if(nVolumes<=2)
	{
		std::cout << "Error(VolumeExtraction): No inner volumes defined!\n";
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
void VolumeExtraction<T>::moveToOuter(int index)
{
	#ifdef ERRORS_ENABLED
	if(index>=nP || index<0)
	{
		std::cout << "Error(VolumeExtraction): Index out of bounds!\n";
		throw 0;
	}
	if(outer==0)
	{
		std::cout << "Error(VolumeExtraction): No outer volume defined!\n";
		throw 0;
	}
	#endif
	
	//is it already an outer particle?
	if(volumeIndex[index].s[0]!=outer)
	{
		//pick one position
		int choose=randNum->rand53()*outerMap.size();
		
		//set particle at that position
		p[index].x=volumeMap[outerMap[choose]].x;
		p[index].y=volumeMap[outerMap[choose]].y;
		p[index].z=volumeMap[outerMap[choose]].z;
		//volumeIndex[index].s[1]=volumeIndex[index].s[0];
		//volumeIndex[index].s[0]=outer;
	}
	#ifdef WARNINGS_ENABLED
	else
	{
		std::cout << "Warning(VolumeExtraction): Already an outer particle, doing nothing.\n";
	}
	#endif
}

template <typename T>
int VolumeExtraction<T>::nExchanged()
{
	int n=0;
	for(int i=0;i<nP;i++)
		if(volumeIndex[i].s[0]!=volumeIndex[i].s[1] && volumeIndex[i].s[0]!=0)
			n++;
	return n;
}

template <typename T>
void VolumeExtraction<T>::exportMap(const char *name, std::ios_base::openmode mode)
{
	std::fstream outFile;
	outFile.open(name,mode);
	if(outFile.is_open())
	{
		outFile << nCells.t << '\n' << "volume" << '\n';
		for(int i=0;i<nCells.t;i++)
			outFile << volumeMap[i].type << '\t' << volumeMap[i].x << '\t' << volumeMap[i].y << '\t' << volumeMap[i].z << '\n';
		outFile.close();
	}
	else
	{
		#ifdef WARNINGS_ENABLED
			std::cout << "Warning(VolumeExtraction): Cannot open " << name << ".\n";
		#endif
	}
}

#endif
