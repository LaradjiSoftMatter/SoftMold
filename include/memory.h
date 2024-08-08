#include "algorithms/dataTypes.h"

class subdivide {
	public:
		subdivide(threeVector<double> size, int nDivisions, double overlap)
		{
			this->s=size;
			this->ol=overlap;
			if(nDivisions<1)
			{
				#ifdef WARNINGS_ENABLED
					std::cout << "Warning(subdivide): nDivisions is negative!\n";
				#endif
				nDivisions=1;
			}
			
			//is it cubic?
			vector<int> factors;//LCM
			while(nDivisions!=1)
			{
				for(int i=2;i<=nDivisions;i++)
				{
					if(nDivisions%i==0)
					{
						factors.push_back(i);
						nDivisions/=i;
						//cout << i << '\t';
						break;
					}
				}
			}
			
			switch(factors.size())
			{
				case 1:
					
					break;
				case 2:
					
					break;
				case 3:
					
					break;
				default:
					
					break;
			}
			
			//set the divider sizes
			divider.x=s.x/static_cast<double>(nDividers.x);
			divider.y=s.y/static_cast<double>(nDividers.y);
			divider.z=s.z/static_cast<double>(nDividers.z);
		}
		~subdivide();
		int getDivision(position<double> p)
		{
			return (p.x/divider.x)+(p.y/divider.y)+(p.z/divider.z)
		}
		int getOverlap(position<double> p)
		{
			return 
		}
	private:
		threeVector<double> s;//size
		double ol;//overlap
		threeVector<double> divider;
		threeVector<int> nDividers;
};

/**
 * \brief Memory handler.
 * Handles memory per thread. Divides memory to eliminate thread collisions.
 * Entire object is accessed via master thread (where the object is declared and resides) 
 * and information attained is sent to slave threads.
 */
template <typename T, typename subdivide>
class memory {
	public:
		memory(position<double> *pos, threeVector<double> *vel, threeVector<double> *acc,
			   int nParticles, molecule<double, fourVector<int> > *mol, int nMolecules,
			   double *TBFC, double *TBUC, int nTypes, threeVector<double> size, double cutoff, int nThreads)
		{
			//set host stuff
			this->p=pos;
			this->v=vel;
			this->a=acc;
			this->m=mol;
			this->nM=nMolecules;
			this->nP=nParticles;
			this->twoBodyFconst=TBFC;
			this->twoBodyUconst=TBUC;
			this->nT=nTypes;
			this->threeVector<double> s=size;
			this->cutoff=cutoff;
			deviceCount=nThreads;
			if(deviceCount<=0)
			{
				std::cout << "Error(memory): No devices found!\n";
				throw 0;
			}
			
			//initialize device variables
			this->pDevice=NULL;
			this->vDevice=NULL;
			this->aDevice=NULL;
			this->mDevice=NULL;
			this->twoBodyFconstDevice=NULL;
			this->index=NULL;
			
			//allocate device memory
			if(allocate())
			{
				std::cout << "Error(memory): Unable to allocate memory!\n";
				throw 0;
			}
			
			//allocate host memory
			this->iMDevice=new int[nM];
			this->iMCDevice=new int[nM];
			this->leftIndex=new threeVector<double>[deviceCount];
			threeVector<int> *ghostIndex;//ghost index is stored in device adjacent to leftIndex of next device
			position<double> *hostPBuf;
			threeVector<double> *hostVBuf;
			int *hostIBuf;
		};
		
		~memory()
		{
			//delete device memory
			#pragma omp parallel for
			for(int i=0;i<deviceCount;i++)
			{
				
			}
			delete pDevice;
			delete aDevice;
			delete vDevice;
			delete mDevice;
			delete twoBodyFconstDevice;
			delete twoBodyUconstDevice;
			delete index;
			delete iMDevice;
			delete iMCDevice;
			delete ghostIndex;
			delete hostPBuf;
			delete hostVBuf;
			delete hostIBuf;
			delete leftIndex;
		};
		
		//allocate device memory
		bool allocate()
		{
			//For better performance determine density gradient
			//ignore last comment for now
			
			//set partitioning index and allocate device memory
			leftIndex=new threeVector<double>[deviceCount];
			
			for(int i=0;i<deviceCount;i++)
			{
				leftIndex[i].x=(s.x/(double)deviceCount)*(double)i;
				if(pDevice[i]!=NULL)
					cudaMalloc(&pDevice[i],sizeof(position<double>)*nP);
				if(vDevice[i]!=NULL)
					cudaMalloc(&vDevice[i],sizeof(threeVector<double>)*nP);
				if(aDevice[i]!=NULL)
					cudaMalloc(&aDevice[i],sizeof(threeVector<double>)*nP);
				if(mDevice!=NULL)//Because I'm lazy!
				{
					nMDevice=0;
					for(int j=0;j<nM;j++)
					{
						nMDevice+=m[j].readNBond();
						iMDevice[j]=nMDevice;
					}
					cudaMalloc(&mDevice[i],nMDevice*sizeof(fourVector <int>));
				}
				if(mConstantsDevice!=NULL)
				{
					nMCDevice=0;
					for(int j=0;j<nM;j++)
					{
						nMDCevice+=m[j].readNConstant();
						iMDCevice[j]=nMDevice;
					}
					cudaMalloc(&mDevice[i],nMCDevice*sizeof(int));
				}
				if(twoBodyFconstDevice[i]!=NULL)
					cudaMalloc(&twoBodyFconstDevice[i], sizeof(double)*nC);
				if(index[i]!=NULL)
					cudaMalloc(&index[i],sizeof(twoVector<int>)*nP);
			}
			return false;
		};
		
		//reshape the partitioning
		bool reshape()
		{
			
		};
		
		bool hostToDevice()
		{
			//scatter
			
			return false;
		};
		
		bool deviceToHost()
		{
			//gather
			
			return false;
		};
		
		position<double> * getPositions(int device)
		{
			if(device<deviceCount && device>=0)
			{
				return pDevice[device];
			}
			std::cout << "Error(memory.getPositions): Device number " << device << " not available!\n";
			throw 0;
		};
		
		threeVector<double> * getVelocities(int device)
		{
			if(device<deviceCount && device>=0)
			{
				return vDevice[device];
			}
			std::cout << "Error(memory.getVelocities): Device number " << device << " not available!\n";
			throw 0;
		};
		
		threeVector<double> * getAccelerations(int device)
		{
			if(device<deviceCount && device>=0)
			{
				return aDevice[device];
			}
			std::cout << "Error(memory.getAccelerations): Device number " << device << " not available!\n";
			throw 0;
		};
		
		double * getTwoBodyFconst(int device)
		{
			if(device<deviceCount && device>=0)
			{
				return twoBodyFconstDevice[device];
			}
			std::cout << "Error(memory.getTwoBodyFconst): Device number " << device << " not available!\n";
			throw 0;
		};
		
		double * getTwoBodyUconst(int device)
		{
			if(device<deviceCount && device>=0)
			{
				return twoBodyUconstDevice[device];
			}
			std::cout << "Error(memory.getTwoBodyUconst): Device number " << device << " not available!\n";
			throw 0;
		};
		
		//moleculeIndex is how partition information is handled
		double * getMConstants(int device, int moleculeIndex)
		{
			if(device<deviceCount)
			{
				if(moleculeIndex>=nMCDevice || moleculeIndex<0)
				{
					std::cout << "Error(memory.getMConstants): Molecule number " << moleculeIndex << " not available!\n";
					throw 0;
				}
				return mConstants[device];
			}
			std::cout << "Error(memory.getMConstants): Device number " << device << " not available!\n";
			throw 0;
		};
		
		int *getIndex(int device)
		{
			
		}
		
		int *getAddress(int device)
		{
			
		}
		
		threeVector<double> readPartition(int device)
		{
			if(device<deviceCount && device>=0)
				return leftIndex;
			std::cout << "Error(memory.readPartition): Device number " << device << " not available!\n";
			throw 0;
		}
		
		void setPartition(int device, threeVector<double> newPartition)
		{
			if(device<deviceCount && device>=0)
			{
				leftIndex=newPartition;
			}
			std::cout << "Error(memory.setPartition): Device number " << device << " not available!\n";
			throw 0; 
		}
		
		void setGhost(int device, threeVector<int> newGhost)
		{
			if(device<deviceCount && device>=0)
			{
				ghostIndex=newGhost;
			}
			std::cout << "Error(memory.setGhost): Device number " << device << " not available!\n";
			throw 0; 
		}
		
	private:
		///Host memory
		position<double> *p;
		threeVector<double> *v;
		threeVector<double> *a;
		molecule<double, fourVector<int> > *m;
		int nM;//number of molecules
		int nP;//number of particles
		double *twoBodyFconst;
		double *twoBodyUconst;
		int nT;//number of types
		double cutoff;
		int deviceCount;
		
		///Device memory
		position<double> **pDevice;//16 bytes per pos
		threeVector<double> **vDevice;//12 bytes per vel
		threeVector<double> **aDevice;//12 bytes per acc
		fourVector<int> **mDevice;//just list all bonds, uses a separate list, 16 bytes per mol
		double **mConstantsDevice;//constants per molecule, 4 bytes per constant
		double **twoBodyFconstDevice;//4 bytes per const
		double **twoBodyUconstDevice;//4 bytes per const
		//not sure if potentials will be computed on device, could do them in parallel on host
		int **index;//4 bytes per pos, Part 1 of two way associative indices, points to current index
		int **address;//4 bytes per pos, Part 2 of two way associative indices, points to correctly ordered index
		//Total mem=(8+16+12+12)*nP+16*nM+4*nC
		//nC=nTBC*nT^2
		
		///Device partition related information
		int nMDevice;//number of device molecules
		int nMCDevice;//number of device molecule constants
		int *iMDevice;//partition indices of deviceMolecules
		int *iMCDevice;//partition indices of mConstants
		threeVector<double> *leftIndex;//spacial partitioning index per device
		threeVector<int> *ghostIndex;//ghost index is stored in device adjacent to leftIndex of next device
		position<double> *hostPBuf;
		threeVector<double> *hostVBuf;
		int *hostIBuf;
};