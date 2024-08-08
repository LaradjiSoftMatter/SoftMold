/**
 * \brief Basic molecular dynamics optimized routines. Used for single processor routines, like those in MPI.
 */

#define SIZE_BOUND

template <typename T, void POTENTIAL(T cutoffSquared, int nT, position<T> &p1, position<T> &p2, T *pC)>
T computeDPotential(position<T> *p, 
			int nParticles, 
			T cutoff, 
			T *TBPC, 
			int nTypes
			#ifdef SIZE_BOUND
				,threeVector<T> size
			#endif
			);

//this is sort of an all in one variant of a cell/verlet list
template <typename T, void FORCE(T cutoffSquared, int nT, position<T> &p1, position<T> &p2, threeVector<T> &a1, threeVector<T> &a2, T *fC)>
void computeForce(position<T> *p, 
			threeVector<T> *v, 
			threeVector<T> *a, 
			int nParticles, 
			T cutoff, 
			T *TBFC, 
			int nTypes
			#ifdef SIZE_BOUND
				,threeVector<T> size
			#endif
			)
{
	int nextBoxStatic[13][3]={{1,0,0},{0,1,0},{1,1,0},{0,0,1},{1,0,1},{0,1,1},{1,1,1},{-1,0,1},{0,-1,1},{-1,1,1},{1,-1,1},{-1,-1,1},{1,-1,0}};
	T cutoffSquared=cutoff*cutoff;
	fourVector<unsigned int> nCells;//housekeeping
	std::vector< keyVal<unsigned int,unsigned int> > hashIndices;//for setting up verlet lists
	#ifndef SIZE_BOUND
		//we can change the shape if need be, 
		// but this is the largest generic variant we can use with 32 bits
		nCells.x=2048;//2^11 }
		nCells.y=2048;//2^11 }  > 2^32
		nCells.z=1024;//2^10 }
		nCells.t=nCells.x*nCells.y*nCells.z;
	#else
		nCells.x=size.x/cutoff;
		nCells.y=size.y/cutoff;
		nCells.z=size.z/cutoff;
		nCells.t=nCells.x*nCells.y*nCells.z;
		threeVector<T> cellSize;
		cellSize.x=size.x/nCells.x;
		cellSize.y=size.y/nCells.y;
		cellSize.z=size.z/nCells.z;
	#endif
	
	std::map<unsigned int, unsigned int> cells;//accessed by (cell,head element)
	
	//find all particles' cell index
	for(unsigned int i=0;i<nParticles;i++)
	{
		#ifndef SIZE_BOUND
			unsigned int x=int(p[i].x/cutoff)%nCells.x;
			unsigned int y=int(p[i].y/cutoff)%nCells.y;
			unsigned int z=int(p[i].z/cutoff)%nCells.z;
		#else
			unsigned int x=int(p[i].x/cellSize.x)%nCells.x;
			unsigned int y=int(p[i].y/cellSize.y)%nCells.y;
			unsigned int z=int(p[i].z/cellSize.z)%nCells.z;
		#endif
		
		keyVal<unsigned int, unsigned int> hI;
		
		hI.value=i;
		hI.key=x+y*nCells.x+z*nCells.x*nCells.y;
		hashIndices.push_back(hI);
		
		std::map<unsigned int, unsigned int>::iterator it;
		//it=cells.find(hashIndices[i].key);
		it=cells.lower_bound(hashIndices[i].key);//supposedly faster (Effective STL)
		//cell linked list
		if(it!=cells.end() && !(cells.key_comp()(it->first, hashIndices[i].key)))
		{
			//push operation
			//unsigned int buf=cells[hashIndices[i].key];//use our old cell head as the next list tail
			//cells[hashIndices[i].key]=hashIndices[i].value;//new cell list head
			unsigned int buf=it->second;//use our old cell head as the next list tail
			it->second=hashIndices[i].value;//new cell list head
			hashIndices[i].value=buf;
		}
		else
		{
			//initialize operation
			//cells[hashIndices[i].key]=hashIndices[i].value;//without hint
			cells.insert(it, std::map<unsigned int, unsigned int>::value_type(hashIndices[i].key, hashIndices[i].value));
			hashIndices[i].value=-1;//it is now the last iterator
		}
	}
//std::cerr << cells.size() << '\t' << hashIndices.size() << '\n';
//std::cin.get();
	//calculate forces between nearby particles, in O(ln(cells.size())+nParticles*density) time
	for(std::map<unsigned int, unsigned int>::iterator current=cells.begin();
	    current!=cells.end();
	    ++current)
	{
		int m,ni;
		//this is for local sorting
		threeVector<T> tempA[MAX_CELL_SIZE];
		position<T> tempP[MAX_CELL_SIZE];
		threeVector<T> currentA;
		position<T> currentP;
//std::cout << current->first << '\n';
//std::cin.get();
		//load each element of current list while computing accelerations
		for(ni=0, m=current->second; m!=-1; ni++, m=hashIndices[m].value)
		{
			#ifdef CELL_SIZE_FAILURE
				if(ni>MAX_CELL_SIZE)
				{
					std::cout << "Error(twoWayCompare): Too many particles in cell!\n";
					throw 0;//I don't know if this even works multithreaded
				}
			#endif
			
			//load an element of current box and initialize
			tempP[ni]=p[m];
			tempA[ni].x=0;
			tempA[ni].y=0;
			tempA[ni].z=0;
			
//std::cout << tempP[ni].type << '\t' << tempP[ni].x << '\t' << tempP[ni].y << '\t' << tempP[ni].z << '\n';
			
			//comapre it to previously loaded elements
			for(int j=0;j<ni;j++)
				FORCE(cutoffSquared,nTypes, tempP[ni], tempP[j], tempA[ni], tempA[j], TBFC);
		}
		
		//unhash our keys
		int x=current->first%nCells.x;
		int y=int(current->first/nCells.x)%nCells.y;
		int z=int(current->first/(nCells.x*nCells.y));
//std::cout << "1\t" << x << '\t' << y << '\t' << z << '\n';
		//look at nearby keys, conservative search
		
		for(int j=0;j<13;j++)
		{
			//compute a nearby key, k
			int nextX=x+nextBoxStatic[j][0];
			int nextY=y+nextBoxStatic[j][1];
			int nextZ=z+nextBoxStatic[j][2];
			
			unsigned int nearbyKey=nextX%nCells.x+((nextY%nCells.y)*nCells.x)+((nextZ%nCells.z)*nCells.x*nCells.y);
			std::map<unsigned int, unsigned int>::iterator it;
			it=cells.find(nearbyKey);
			
			//check for existence of nearby key in O(ln(cells.size())) time
			//if(cells.find(nearbyKey)!=cells.end())
			if(it!=cells.end())
			{
				#ifdef SIZE_BOUND
					//minimum image, for size bound
					threeVector<T> minImg;
					minImg.x=0;
					if(nextX<0) minImg.x-=size.x;
					if(nextX>=nCells.x) minImg.x+=size.x;
					
					minImg.y=0;
					if(nextY<0) minImg.y-=size.y;
					if(nextY>=nCells.y) minImg.y+=size.y;
					
					minImg.z=0;
					if(nextZ<0) minImg.z-=size.z;
					if(nextZ>=nCells.z) minImg.z+=size.z;
				#endif
				//std::cout << "2\t" << x << '\t' << y << '\t' << z << '\n';
				//compare all elements in current cell to every element of next cell
				for(int n=it->second; n!=-1; n=hashIndices[n].value)
				{
					//this is confusing, currentP is actually an element of the next box
					currentP=p[n];
					#ifdef SIZE_BOUND
						currentP.x+=minImg.x;
						currentP.y+=minImg.y;
						currentP.z+=minImg.z;
					#endif
					currentA.x=0;
					currentA.y=0;
					currentA.z=0;
					//compare element to everything in current box
					for(m=0;m<ni;m++)
						FORCE(cutoffSquared, nTypes, tempP[m], currentP, tempA[m], currentA, TBFC);
						//input.force(tempP[m], currentP, tempA[m], currentA);
					
//std::cout << currentA.x << '\t' << currentA.y << '\t' << currentA.z << '\n';
//std::cout << currentP.type << '\t' << currentP.x << '\t' << currentP.y << '\t' << currentP.z << '\n';
					//reduce next acceleration element back to array
					a[n].x+=currentA.x;
					a[n].y+=currentA.y;
					a[n].z+=currentA.z;
//std::cout << a << '\n';
//std::cin.get();
				}
			}
		}
//std::cin.get();
		//reduce the current accelerations of cell back to array
		for(ni=0, m=current->second; m!=-1; ni++, m=hashIndices[m].value)
		{
			a[m].x+=tempA[ni].x;
			a[m].y+=tempA[ni].y;
			a[m].z+=tempA[ni].z;
		}
	}
	//return input;
}

//this one does clusters
template <typename T, void FORCE(T cutoffSquared, int nT, position<T> &p1, position<T> &p2, threeVector<T> &a1, threeVector<T> &a2, T *fC)>
void computeForce(position<T> *p, 
		  threeVector<T> *v, 
		  threeVector<T> *a, 
		  int nParticles, 
		  T cutoff, 
		  T *TBFC, 
		  int nTypes
		  #ifdef SIZE_BOUND
		  	,threeVector<T> size,
		  #endif
		  int *cluster,
		  int nCluster,
		  T clusterRadius)
{
	std::vector< std::vector< int > > clusterNeighbors(nCluster, std::vector<int>());
	
	#pragma omp parallel for
	for(int i=0;i<clusterNeighbors.size();i++)
	{
		int index=cluster[i];
		position<T> myPos=p[cluster[i]];
		for(int j=0;j<nParticles;j++)
		{
			if(j!=index)
			{
				threeVector<T> d;
				d.x=myPos.x-p[j].x;
				d.y=myPos.y-p[j].y;
				d.z=myPos.z-p[j].z;
				int key=
				if(d.x<
					clusterNeighbors[i].push_back(j);
				
			}
		}
	}
	
	int nextBoxStatic[13][3]={{1,0,0},{0,1,0},{1,1,0},{0,0,1},{1,0,1},{0,1,1},{1,1,1},{-1,0,1},{0,-1,1},{-1,1,1},{1,-1,1},{-1,-1,1},{1,-1,0}};
	T cutoffSquared=cutoff*cutoff;
	fourVector<unsigned int> nCells;//housekeeping
	std::vector< keyVal<unsigned int,unsigned int> > hashIndices;//for setting up verlet lists
	#ifndef SIZE_BOUND
		//we can change the shape if need be, 
		// but this is the largest generic variant we can use with 32 bits
		nCells.x=2048;//2^11 }
		nCells.y=2048;//2^11 }  > 2^32
		nCells.z=1024;//2^10 }
		nCells.t=nCells.x*nCells.y*nCells.z;
	#else
		nCells.x=size.x/cutoff;
		nCells.y=size.y/cutoff;
		nCells.z=size.z/cutoff;
		nCells.t=nCells.x*nCells.y*nCells.z;
		threeVector<T> cellSize;
		cellSize.x=size.x/nCells.x;
		cellSize.y=size.y/nCells.y;
		cellSize.z=size.z/nCells.z;
	#endif
	
	std::map<unsigned int, unsigned int> cells;//accessed by (cell,head element)
	
	//find all particles' cell index
	for(unsigned int i=0;i<nParticles;i++)
	{
		#ifndef SIZE_BOUND
			unsigned int x=int(p[i].x/cutoff)%nCells.x;
			unsigned int y=int(p[i].y/cutoff)%nCells.y;
			unsigned int z=int(p[i].z/cutoff)%nCells.z;
		#else
			unsigned int x=int(p[i].x/cellSize.x)%nCells.x;
			unsigned int y=int(p[i].y/cellSize.y)%nCells.y;
			unsigned int z=int(p[i].z/cellSize.z)%nCells.z;
		#endif
		
		keyVal<unsigned int, unsigned int> hI;
		
		hI.value=i;
		hI.key=x+y*nCells.x+z*nCells.x*nCells.y;
		hashIndices.push_back(hI);
		
		std::map<unsigned int, unsigned int>::iterator it;
		//it=cells.find(hashIndices[i].key);
		it=cells.lower_bound(hashIndices[i].key);//supposedly faster (Effective STL)
		//cell linked list
		if(it!=cells.end() && !(cells.key_comp()(it->first, hashIndices[i].key)))
		{
			//push operation
			//unsigned int buf=cells[hashIndices[i].key];//use our old cell head as the next list tail
			//cells[hashIndices[i].key]=hashIndices[i].value;//new cell list head
			unsigned int buf=it->second;//use our old cell head as the next list tail
			it->second=hashIndices[i].value;//new cell list head
			hashIndices[i].value=buf;
		}
		else
		{
			//initialize operation
			//cells[hashIndices[i].key]=hashIndices[i].value;//without hint
			cells.insert(it, std::map<unsigned int, unsigned int>::value_type(hashIndices[i].key, hashIndices[i].value));
			hashIndices[i].value=-1;//it is now the last iterator
		}
	}
//std::cerr << cells.size() << '\t' << hashIndices.size() << '\n';
//std::cin.get();
	//calculate forces between nearby particles, in O(ln(cells.size())+nParticles*density) time
	for(std::map<unsigned int, unsigned int>::iterator current=cells.begin();
	    current!=cells.end();
	    ++current)
	{
		int m,ni;
		//this is for local sorting
		threeVector<T> tempA[MAX_CELL_SIZE];
		position<T> tempP[MAX_CELL_SIZE];
		threeVector<T> currentA;
		position<T> currentP;
//std::cout << current->first << '\n';
//std::cin.get();
		//load each element of current list while computing accelerations
		for(ni=0, m=current->second; m!=-1; ni++, m=hashIndices[m].value)
		{
			#ifdef CELL_SIZE_FAILURE
				if(ni>MAX_CELL_SIZE)
				{
					std::cout << "Error(twoWayCompare): Too many particles in cell!\n";
					throw 0;//I don't know if this even works multithreaded
				}
			#endif
			
			//load an element of current box and initialize
			tempP[ni]=p[m];
			tempA[ni].x=0;
			tempA[ni].y=0;
			tempA[ni].z=0;
			
//std::cout << tempP[ni].type << '\t' << tempP[ni].x << '\t' << tempP[ni].y << '\t' << tempP[ni].z << '\n';
			
			//comapre it to previously loaded elements
			for(int j=0;j<ni;j++)
				FORCE(cutoffSquared,nTypes, tempP[ni], tempP[j], tempA[ni], tempA[j], TBFC);
		}
		
		//unhash our keys
		int x=current->first%nCells.x;
		int y=int(current->first/nCells.x)%nCells.y;
		int z=int(current->first/(nCells.x*nCells.y));
//std::cout << "1\t" << x << '\t' << y << '\t' << z << '\n';
		//look at nearby keys, conservative search
		
		for(int j=0;j<13;j++)
		{
			//compute a nearby key, k
			int nextX=x+nextBoxStatic[j][0];
			int nextY=y+nextBoxStatic[j][1];
			int nextZ=z+nextBoxStatic[j][2];
			
			unsigned int nearbyKey=nextX%nCells.x+((nextY%nCells.y)*nCells.x)+((nextZ%nCells.z)*nCells.x*nCells.y);
			std::map<unsigned int, unsigned int>::iterator it;
			it=cells.find(nearbyKey);
			
			//check for existence of nearby key in O(ln(cells.size())) time
			//if(cells.find(nearbyKey)!=cells.end())
			if(it!=cells.end())
			{
				#ifdef SIZE_BOUND
					//minimum image, for size bound
					threeVector<T> minImg;
					minImg.x=0;
					if(nextX<0) minImg.x-=size.x;
					if(nextX>=nCells.x) minImg.x+=size.x;
					
					minImg.y=0;
					if(nextY<0) minImg.y-=size.y;
					if(nextY>=nCells.y) minImg.y+=size.y;
					
					minImg.z=0;
					if(nextZ<0) minImg.z-=size.z;
					if(nextZ>=nCells.z) minImg.z+=size.z;
				#endif
				//std::cout << "2\t" << x << '\t' << y << '\t' << z << '\n';
				//compare all elements in current cell to every element of next cell
				for(int n=it->second; n!=-1; n=hashIndices[n].value)
				{
					//this is confusing, currentP is actually an element of the next box
					currentP=p[n];
					#ifdef SIZE_BOUND
						currentP.x+=minImg.x;
						currentP.y+=minImg.y;
						currentP.z+=minImg.z;
					#endif
					currentA.x=0;
					currentA.y=0;
					currentA.z=0;
					//compare element to everything in current box
					for(m=0;m<ni;m++)
						FORCE(cutoffSquared, nTypes, tempP[m], currentP, tempA[m], currentA, TBFC);
						//input.force(tempP[m], currentP, tempA[m], currentA);
					
//std::cout << currentA.x << '\t' << currentA.y << '\t' << currentA.z << '\n';
//std::cout << currentP.type << '\t' << currentP.x << '\t' << currentP.y << '\t' << currentP.z << '\n';
					//reduce next acceleration element back to array
					a[n].x+=currentA.x;
					a[n].y+=currentA.y;
					a[n].z+=currentA.z;
//std::cout << a << '\n';
//std::cin.get();
				}
			}
		}
//std::cin.get();
		//reduce the current accelerations of cell back to array
		for(ni=0, m=current->second; m!=-1; ni++, m=hashIndices[m].value)
		{
			a[m].x+=tempA[ni].x;
			a[m].y+=tempA[ni].y;
			a[m].z+=tempA[ni].z;
		}
	}
	//return input;
}

template <typename T>
void langevinThermostat(threeVector<T> *v, 
			threeVector<T> *a, 
			int nP, 
			T temperature,
			T dT,
			T g,
			MTRand *randNum //pre initialized
			)
{
	T sigma=sqrt((6.0*temperature*1.0)/dT);
	for(int i=0;i<nP;i++)
	{
		T psi;
		psi=2.0*randNum[0].rand53()-1.0;
		a[i].x+=(-g*v[i].x+sigma*psi);
		psi=2.0*randNum[0].rand53()-1.0;
		a[i].y+=(-g*v[i].y+sigma*psi);
		psi=2.0*randNum[0].rand53()-1.0;
		a[i].z+=(-g*v[i].z+sigma*psi);
	}
}


template <typename T>
void firstIntegration(position<T> *p,
		      threeVector<T> *v,
		      threeVector<T> *a, 
		      int nP, 
		      T dT
)
{
	T halfDt=0.5*dT;
	for(int i=0;i<nP;i++)
	{
		//increment velocity
		v[i].x+=(a[i].x*halfDt);
		v[i].y+=(a[i].y*halfDt);
		v[i].z+=(a[i].z*halfDt);
		//increment position
		p[i].x+=v[i].x*dT;
		p[i].y+=v[i].y*dT;
		p[i].z+=v[i].z*dT;
		
	}
}

template <typename T>
void secondIntegration(position<T> *p,
		      threeVector<T> *v,
		      threeVector<T> *a, 
		      int nP, 
		      T dT
)
{
	T halfDt=dT*0.5;
	for(int i=0;i<nP;i++)
	{
		v[i].x+=(a[i].x*halfDt);
		v[i].y+=(a[i].y*halfDt);
		v[i].z+=(a[i].z*halfDt);
	}
}
