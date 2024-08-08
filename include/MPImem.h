/**
 * \brief Thread memory handling object.
 * Handles memory per thread. Divides memory to eliminate thread collisions.
 * Entire object is accessed via master thread (where the object is declared and resides) 
 * and information attained is sent to slave threads.
 */

#include "algorithms/dataTypes.h"
#include <mpi.h>

template<typename T>
class positionKeyCompare {
	public:
		positionKeyCompare(position<T> *particles, int dimension)
			:p(particles), dim(dimension)
		{
		};
		bool operator () (int a, int b)
		{
			//the '<' operation will have a significant effect on
			// bound searching in the section
			// "consolidating pivots".
			return p[a].s[dim]<p[b].s[dim];
		};
	private:
		position<T> *p;
		int dim;
};

template<typename T>
class mpiMemory {
	public:
		mpiMemory()
			:masterP(NULL), masterV(NULL), masterA(NULL), masterNP(0),
			masterBondList(NULL), masterBendList(NULL), masterNBonds(0), masterNBends(0),
			masterBondConst(NULL), masterBendConst(NULL), masterBondConstIndex(NULL), masterBendConstIndex(NULL),
			masterNBondConst(0), masterNBendConst(0)
		{ 
			//#ifdef MPI_VERSION
			MPI_Type_contiguous(sizeof(position<T>), MPI_BYTE, &MPI_P_Type);
			MPI_Type_commit(&MPI_P_Type);
			MPI_Type_contiguous(sizeof(threeVector<T>), MPI_BYTE, &MPI_ThreeVectorFloat_Type);
			MPI_Type_commit(&MPI_ThreeVectorFloat_Type);
			MPI_Type_contiguous(sizeof(twoVector<T>), MPI_BYTE, &MPI_TwoVectorFloat_Type);
			MPI_Type_commit(&MPI_TwoVectorFloat_Type);
			MPI_Type_contiguous(sizeof(threeVector<int>), MPI_BYTE, &MPI_ThreeVectorInt_Type);
			MPI_Type_commit(&MPI_ThreeVectorInt_Type);
			MPI_Type_contiguous(sizeof(twoVector<int>), MPI_BYTE, &MPI_TwoVectorInt_Type);
			MPI_Type_commit(&MPI_TwoVectorInt_Type);
			//#endif
		};
		
		~mpiMemory()
		{ };
		
		//initialize communicators from slaves
		void initComm(position<T> *particles, threeVector<T> *velocities, threeVector<T> *accelerations,
			      int nParticles, twoVector<int> *bondList, int nBonds, threeVector<int> *bendList, int nBends,
			      twoVector<T> *bondConst, int *bondConstIndex, int nBondConst,
			      twoVector<T> *bendConst, int *bendConstIndex, int nBendConst,
			      threeVector<T> boundSize)
		{
			//#ifndef MPI_VERSION
			//#else
			int rank, nThreads;
			MPI_Group group;
			MPI_Status *status;
			MPI_Comm_rank (MPI_COMM_WORLD, &rank);
			MPI_Comm_size (MPI_COMM_WORLD, &nThreads);
			MPI_Comm_group(MPI_COMM_WORLD, &group);
			//#endif
			
			//buffers for the master thread
			std::vector<int> sendNP(nThreads,0);
			std::vector<int> sendPartitions(nThreads,0);
			std::vector<position<T> > sendP;
			std::vector<threeVector<T> > sendV;
			std::vector<threeVector<T> > sendA;
			std::vector<int> sendIP;
			std::vector<int> divisors;
			std::vector<threeVector<T> > lBounds;
			std::vector<threeVector<T> > rBounds;
			std::vector< std::vector<int> > partitions;
		
			if(rank==0)
			{
				//rank 0 is the only rank with the master lists
				// The other ranks should be NULL or 0
				masterP=particles;
				masterV=velocities;
				masterA=accelerations;
				masterNP=nParticles;
				masterSize=boundSize;
				
				masterBondList=bondList;
				masterNBonds=nBonds;
				masterBondConst=bondConst;
				masterBondConstIndex=bondConstIndex;
				masterNBondConst=nBondConst;
				
				masterBendList=bendList;
				masterNBends=nBends;
				masterBendConst=bendConst;
				masterBendConstIndex=bendConstIndex;
				masterNBendConst=nBendConst;
				
				//these are resized to the max size so they can act as an organization buffer
				p.resize(masterNP);
				v.resize(masterNP);
				a.resize(masterNP);
				iP.resize(masterNP);
				boL.resize(masterNBonds);
				beL.resize(masterNBends);
				boC.resize(masterNBondConst);
				beC.resize(masterNBendConst);
				boCI.resize(masterNBondConst);
				beCI.resize(masterNBendConst);
				
				//some extra initializations
				nBo=masterNBonds;
				nBe=masterNBends;
				nBoC=masterNBondConst;
				nBeC=masterNBendConst;
				size=masterSize;
				
				//initialize master's buffers
				for(int i=0;i<masterNP;i++)
					iP[i]=i;
				
				for(int i=0;i<masterNBonds;i++)
					boL[i]=masterBondList[i];
				
				for(int i=0;i<masterNBends;i++)
					beL[i]=masterBendList[i];
				
				for(int i=0;i<masterNBondConst;i++)
				{
					boC[i]=masterBondConst[i];
					boCI[i]=masterBondConstIndex[i];
				}
				
				for(int i=0;i<masterNBendConst;i++)
				{
					beC[i]=masterBendConst[i];
					beCI[i]=masterBendConstIndex[i];
				}
				
				//factoring primes!
				int remaining=nThreads;
				for(int i=2;i<=remaining;i++)
				{
					bool isPrime=true;
					for(int j=2;j<i && isPrime;j++)
						if(i%j==0)
							isPrime=false;
					while(remaining%i==0 & isPrime)
					{
						divisors.push_back(i);
						remaining/=i;
					}
				}
				
				//swap the order to get largest common denominators first
				for(int i=0;i<divisors.size()/2;i++)
				{
					int buf=divisors[i];
					divisors[i]=(*(divisors.end()-i-1));
					(*(divisors.end()-i-1))=buf;
				}
				
				//set up partitions using divisors, partitions are ordered x,y,z,x,y,z...
				for(int i=0;i<divisors.size();i++)
				{
					std::vector<int> partitionGroup;
					if(i==0)//initial pivots
					{
						for(int j=0;j<divisors[i];j++)
							partitionGroup.push_back(masterNP*j/divisors[i]);
						partitionGroup.push_back(masterNP);
					}
					else
					{
						for(int j=0;j<partitions[i-1].size()-1;j++)
						{
							int start=partitions[i-1][j];
							int end=partitions[i-1][j+1];
							int length=end-start;
							for(int k=0;k<divisors[i];k++)
								partitionGroup.push_back(start+length*k/divisors[i]);
						}
						partitionGroup.push_back(partitions[i-1][partitions[i-1].size()-1]);
					}
					partitions.push_back(partitionGroup);
				}
				
				//now that we have the partitions, sort the data using partitions, kdTree style
				if(partitions.size()>0)
				{
					positionKeyCompare<T> sortKeys(masterP, 0);
					for(int j=0;j<partitions[0].size()-1;j++)
					{
						int start;
						if(j==0)
							start=0;
						else
							start=partitions[0][j-1];
						int middle=partitions[0][j];
						int end=partitions[0][j+1];
						std::partial_sort(iP.begin()+start, iP.begin()+middle, iP.begin()+end, sortKeys);
					}
				}
				for(int i=1;i<partitions.size();i++)
				{
					int dim=i%3;
					positionKeyCompare<T> sortKeys(masterP, dim);
					for(int j=0;j<partitions[i].size()-1;j++)
					{
						int start;
						if(j==0)
							start=0;
						else
							start=partitions[i][j-1];
						int middle=partitions[i][j];
						int end=partitions[i][j+1];
						int partIndex=find(partitions[i-1].begin(),partitions[i-1].end(),middle)
								-partitions[i-1].begin();
						if(partIndex==partitions[i-1].size())
							std::partial_sort(iP.begin()+start, iP.begin()+middle, iP.begin()+end, sortKeys);
					}
				}
				
				//order the master's lists so they can be used as a buffer
				for(int i=0;i<(partitions.end()-1)->size()-1;i++)
				{
					for(int j=(*(partitions.end()-1))[i];j<(*(partitions.end()-1))[i+1];j++)
					{
						sendP.push_back(masterP[iP[j]]);
						sendV.push_back(masterV[iP[j]]);
						sendA.push_back(masterA[iP[j]]);
						sendIP.push_back(iP[j]);
					}
				}
				
				//locating our bounds, first assume everyone is inside the system box
				lBounds.resize(nThreads);
				rBounds.resize(nThreads);
				for(int i=0;i<lBounds.size();i++)
					lBounds[i]=0;
				for(int i=0;i<rBounds.size();i++)
					rBounds[i]=masterSize;
				
				//assuming {i++, f(i)} ordering in "i++,groupSize/=divisors[i]"
				for(int i=0, groupSize=nThreads/divisors[0];i<divisors.size();i++,groupSize/=divisors[i])
				{
					std::vector<T> boundLowPivots;
					std::vector<T> boundHighPivots;
					
					//find some nearby pivots
					for(int j=0;j<nThreads;j+=groupSize)
					{
						T boundLow=rBounds[j].s[i%3];//masterSize.s[i%3];
						T boundHigh=lBounds[j].s[i%3];
						for(int k=(*(partitions.end()-1))[j];k<(*(partitions.end()-1))[j+groupSize];k++)
						{
							boundLow=(masterP[iP[k]].s[i%3]<boundLow)?masterP[iP[k]].s[i%3]:boundLow;
							boundHigh=(masterP[iP[k]].s[i%3]>boundHigh)?masterP[iP[k]].s[i%3]:boundHigh;
						}
						boundLowPivots.push_back(boundLow);
						boundHighPivots.push_back(boundHigh);
					}
					
					//consolidate pivots, order determined from positionKeyCompare
					for(int j=0;j<boundLowPivots.size()-1;j++)
					{
						//it is a middle pivot in need of consolidation
						if((j+1)%divisors[i]!=0)
						{
							boundLowPivots[j+1]=(boundLowPivots[j+1]+boundHighPivots[j])/2.0;
							boundHighPivots[j]=boundLowPivots[j+1];
						}
					}
					
					//fix the ends to match the bounds
					for(int j=0;j<boundLowPivots.size();j++)
					{
						//left ends
						if(j%divisors[i]==0)
							boundLowPivots[j]=lBounds[j].s[i%3];
						//right ends
						if(j%divisors[i]==divisors[i]-1)
							boundHighPivots[j]=rBounds[j].s[i%3];
					}
					
					//redistribute pivots to sizes
					for(int j=0;j<nThreads;j++)
					{
						lBounds[j].s[i%3]=boundLowPivots[floor(j/groupSize)];
						rBounds[j].s[i%3]=boundHighPivots[floor(j/groupSize)];
					}
				}
				/*
				for(int i=0;i<lBounds.size();i++)
				{
					std::cout << i << '\t' << lBounds[i].s[0] << '\t' << lBounds[i].s[1] << '\t' << lBounds[i].s[2] << std::endl;
					std::cout << i << '\t' << rBounds[i].s[0] << '\t' << lBounds[i].s[1] << '\t' << lBounds[i].s[2] << std::endl;
					std::cout << i << '\t' << lBounds[i].s[0] << '\t' << rBounds[i].s[1] << '\t' << lBounds[i].s[2] << std::endl;
					std::cout << i << '\t' << lBounds[i].s[0] << '\t' << lBounds[i].s[1] << '\t' << rBounds[i].s[2] << std::endl;
					std::cout << i << '\t' << rBounds[i].s[0] << '\t' << rBounds[i].s[1] << '\t' << lBounds[i].s[2] << std::endl;
					std::cout << i << '\t' << lBounds[i].s[0] << '\t' << rBounds[i].s[1] << '\t' << rBounds[i].s[2] << std::endl;
					std::cout << i << '\t' << rBounds[i].s[0] << '\t' << lBounds[i].s[1] << '\t' << rBounds[i].s[2] << std::endl;
					std::cout << i << '\t' << rBounds[i].s[0] << '\t' << rBounds[i].s[1] << '\t' << rBounds[i].s[2] << std::endl;
				}*/
				
				//set everyones' nParticles
				for(int i=0;i<nThreads-1;i++)
					sendNP[i]=(*(partitions.end()-1))[i+1]-(*(partitions.end()-1))[i];
				sendNP[nThreads-1]=masterNP-(*(partitions.end()-1))[nThreads-1];
				
				//to maintain partitions through threads
				for(int i=0;i<(partitions.end()-1)->size();i++)
					sendPartitions[i]=(*(partitions.end()-1))[i];
				//sendPartitions[nThreads-1]=masterNP;
			}
			
			//set up groups
			int nDivisions=0;
			if(rank==0)
				nDivisions=divisors.size();
			
			MPI_Barrier(MPI_COMM_WORLD);
			
			//tell everyone how many groups there are
			MPI_Bcast(&nDivisions, 1, MPI_INT, 0, MPI_COMM_WORLD);
			if(rank!=0)
				divisors.resize(nDivisions);
			
			MPI_Barrier(MPI_COMM_WORLD);
			
			//send out divisors so everyone knows where the splits are
			MPI_Bcast(&divisors[0], nDivisions, MPI_INT, 0, MPI_COMM_WORLD);
			
			MPI_Barrier(MPI_COMM_WORLD);
			
			//very top level
			groupHandles.push_back(group);
			commHandles.push_back(MPI_COMM_WORLD);
			
			//initial divisions
			boundLevel.s[0]=0;
			boundLevel.s[1]=1;
			boundLevel.s[2]=2;
			
			for(int i=0;i<divisors.size();i++)
			{
				int commSize, commRank;
				MPI_Group newGroup;
				MPI_Comm newComm;
				
				//for current group
				MPI_Comm_rank(commHandles[i], &commRank); 
				MPI_Comm_size(commHandles[i], &commSize);
				
				//determine own group
				std::vector<int> groupMembers;
				int selfDiv=commSize/divisors[i];
				int self=commRank/selfDiv;
				
				//figure out who is nearby
				for(int j=0;j<nThreads;j++)
					if(j/selfDiv==self)
						groupMembers.push_back(j);
					
				//include everyone into group
				MPI_Group_incl(groupHandles[i], groupMembers.size(), &groupMembers[0], &newGroup);
				MPI_Comm_create(commHandles[i], newGroup, &newComm);
				
				//store the group for later
				groupHandles.push_back(newGroup);
				commHandles.push_back(newComm);
				
				//keeps track of where boundary came from
				if(commRank==0 && i>2)
					boundLevel.s[i%3]+=3;
				MPI_Bcast(&boundLevel.s[i%3], 1, MPI_INT, 0, newComm);
				
				MPI_Barrier(MPI_COMM_WORLD);
			}
			
			//T startTime=MPI_Wtime();
			
			//Send bounds and number of particles to everyone (top level)
			MPI_Scatter(&sendNP[0], 1, MPI_INT, &nP, 1, MPI_INT, 0, MPI_COMM_WORLD);
			MPI_Scatter(&lBounds[0], 1, MPI_ThreeVectorFloat_Type, &lBound, 1, MPI_ThreeVectorFloat_Type, 0, MPI_COMM_WORLD);
			MPI_Scatter(&rBounds[0], 1, MPI_ThreeVectorFloat_Type, &rBound, 1, MPI_ThreeVectorFloat_Type, 0, MPI_COMM_WORLD);
			nPSelf=nP;//unique never changes through edge scattering
			
			//make sure we have enough room for particle state information
			if(rank!=0)
			{
				p.resize(nP);
				v.resize(nP);
				a.resize(nP);
				iP.resize(nP);
			}
			MPI_Barrier(MPI_COMM_WORLD);
			
			//send particle state information
			MPI_Scatterv(&sendP[0], &sendNP[0], &sendPartitions[0], MPI_P_Type, 
				     &p[0], nP, MPI_P_Type, 0, MPI_COMM_WORLD);
			MPI_Scatterv(&sendV[0], &sendNP[0], &sendPartitions[0], MPI_ThreeVectorFloat_Type, 
				     &v[0], nP, MPI_ThreeVectorFloat_Type, 0, MPI_COMM_WORLD);
			MPI_Scatterv(&sendA[0], &sendNP[0], &sendPartitions[0], MPI_ThreeVectorFloat_Type, 
				     &a[0], nP, MPI_ThreeVectorFloat_Type, 0, MPI_COMM_WORLD);
			MPI_Scatterv(&sendIP[0], &sendNP[0], &sendPartitions[0], MPI_INT, 
				     &iP[0], nP, MPI_INT, 0, MPI_COMM_WORLD);
			MPI_Barrier(MPI_COMM_WORLD);
			
			//Send particle bonding and bending quantity (top level)
			MPI_Bcast(&nBo, 1, MPI_INT, 0, MPI_COMM_WORLD);
			MPI_Bcast(&nBe, 1, MPI_INT, 0, MPI_COMM_WORLD);
			MPI_Bcast(&nBoC, 1, MPI_INT, 0, MPI_COMM_WORLD);
			MPI_Bcast(&nBeC, 1, MPI_INT, 0, MPI_COMM_WORLD);
			
			//make sure we have room to hold the bonding and bending information
			if(rank!=0)
			{
				boL.resize(nBo);
				beL.resize(nBe);
				boC.resize(nBoC);
				beC.resize(nBeC);
				boCI.resize(nBoC);
				beCI.resize(nBeC);
			}
			MPI_Barrier(MPI_COMM_WORLD);
			
			//send bonding and bending states
			MPI_Bcast(&boL[0], nBo, MPI_TwoVectorInt_Type, 0, MPI_COMM_WORLD);
			MPI_Bcast(&beL[0], nBe, MPI_ThreeVectorInt_Type, 0, MPI_COMM_WORLD);
			MPI_Bcast(&boC[0], nBoC, MPI_TwoVectorFloat_Type, 0, MPI_COMM_WORLD);
			MPI_Bcast(&beC[0], nBeC, MPI_TwoVectorFloat_Type, 0, MPI_COMM_WORLD);
			MPI_Bcast(&boCI[0], nBoC, MPI_INT, 0, MPI_COMM_WORLD);
			MPI_Bcast(&beCI[0], nBeC, MPI_INT, 0, MPI_COMM_WORLD);
			
			//size information
			MPI_Bcast(&size, 1, MPI_ThreeVectorFloat_Type, 0, MPI_COMM_WORLD);
			MPI_Barrier(MPI_COMM_WORLD);
		};
		
		//takes ~nP*log((nP)^(2/3)) or nP*log(nP) time
		void passEdgePositionsRight(T cutoff)
		{
			//particles to distribute among communicators,
			// this only needs to exist long enough to perform the transfer!
			std::vector< position<T> > LXEdge;
			std::vector< position<T> > LYEdge;
			std::vector< position<T> > LZEdge;
			std::vector< position<T> > RXEdge;
			std::vector< position<T> > RYEdge;
			std::vector< position<T> > RZEdge;
			
			//nP is the overall number of particles
			//nPSelf is the groups own particles
			nP=nPSelf;
			
			//my particles' "local indices" to distribute are cleared since we are starting over
			mLXEdge.clear();
			mLYEdge.clear();
			mLZEdge.clear();
			mRXEdge.clear();
			mRYEdge.clear();
			mRZEdge.clear();
			//communicator source offset is cleared since we are starting over
			/*
			sLXEdge.clear();
			sLYEdge.clear();
			sLZEdge.clear();
			sRXEdge.clear();
			sRYEdge.clear();
			sRZEdge.clear();
			*/
			
			//gather the edges from own thread
			for(int i=0;i<nPSelf;i++)
			{
				if(p[i].x-lBound.x<cutoff)
				{
					mLXEdge.push_back(i);
					LXEdge.push_back(p[i]);
				}
				if(p[i].y-lBound.y<cutoff)
				{
					mLYEdge.push_back(i);
					LYEdge.push_back(p[i]);
				}
				if(p[i].z-lBound.z<cutoff)
				{
					mLZEdge.push_back(i);
					LZEdge.push_back(p[i]);
				}
				
				if(rBound.x-p[i].x<cutoff)
				{
					mRXEdge.push_back(i);
					RXEdge.push_back(p[i]);
				}
				if(rBound.y-p[i].y<cutoff)
				{
					mRYEdge.push_back(i);
					RYEdge.push_back(p[i]);
				}
				if(rBound.z-p[i].z<cutoff)
				{
					mRZEdge.push_back(i);
					RZEdge.push_back(p[i]);
				}
			}
			
			//communicator source offset is cleared since we are starting over
			sGroup.clear();
			
			//set up initial offset for own thread
			//sLXEdge.push_back(mLXEdge.size());
			//sLYEdge.push_back(mLYEdge.size());
			//sLZEdge.push_back(mLZEdge.size());
			//sRXEdge.push_back(mRXEdge.size());
			//sRYEdge.push_back(mRYEdge.size());
			//sRZEdge.push_back(mRZEdge.size());
			
			//for each level, pass up all gathered
			for(int i=commHandles.size()-1;i>=0;i--)
			{
				int rank, nThreads;
				MPI_Comm_rank(commHandles[i], &rank); 
				MPI_Comm_size(commHandles[i], &nThreads);
				
				//small buffer for left and right groups
				std::vector<int> left, right;
				
				for(int j=0;j<nThreads;j++)
				{
					int selfDiv=nThreads/divisors[i];
					int self=rank/selfDiv;
					int lOffset=self;
					for(int k=i+1;k<divisors.size();k+=3)
					{
						lOffset=lOffset%divisors[k];
					}
					if(lOffset==0)
						left.push_back(
				}
				
				//int smallestDivisor=divisors.size()-(divisors.size()-1)%3-(i%3);
				int lcm=1;
				for(int j=i;j<divisors.size();j++)
					lcm*=divisors[j];
				
				//gather together all the ends
				//for lcm=1, every end is used, and rGroup=lGroup
				//for lcm=2, every other end is used, and rGroup!=lGroup
				//for lcm=3, every other end is used, and rGroup!=lGroup
				//etc...
				for(int j=0;j<nThreads;j++)
				{
					if(j%lcm==lcm-1 || (j%lcm==0 && i<3))
						right.push_back(j);
					if(j%lcm==0)
						left.push_back(j);
				}
				
				int rightEdgeSize=0;
				int *rightEdge=NULL;
				position<T> *rightEdgeParticles;
				switch(i%3)
				{
					case 0:
						rightEdgeSize=mRXEdge.size();
						rightEdge=&mRXEdge[0];
						rightEdgeParticles=&RXEdge[0];
						break;
					case 1:
						rightEdgeSize=mRYEdge.size();
						rightEdge=&mRYEdge[0];
						rightEdgeParticles=&RYEdge[0];
						break;
					case 2:
						rightEdgeSize=mRZEdge.size();
						rightEdge=&mRZEdge[0];
						rightEdgeParticles=&RZEdge[0];
						break;
					default:
						break;
				}
				
				int largestGroup=(right.size()>left.size())?right.size():left.size();
				//graph overlap checking can be used to shorten communication time for huge numbers of processors
				//but we will just march through groups for now
				for(int j=0;j<lcm;j++)
				{
					for(int k=0;k<lcm
				}
				
				//periodic boundaries
				if(i<3)
				{
					//odd transfer
					for(int j=0;j<lGroup.size();j+=2)
					{
						for(int k=1;k<rGroup.size();k+=2)
						{
							if(rank==lGroup[j])
						{
							int rightEdgeSize=0;
							switch(i%3)
							{
								case 0:
									rightEdgeSize=mRXEdge.size();
									break;
								case 1:
									rightEdgeSize=mRYEdge.size();
									break;
								case 2:
									rightEdgeSize=mRZEdge.size();
									break;
								default:
									break;
							}
							
							int *rightEdge=NULL;
							switch(i%3)
							{
								case 0:
									rightEdge=&mRXEdge[0];
									break;
								case 1:
									rightEdge=&mRYEdge[0];
									break;
								case 2:
									rightEdge=&mRZEdge[0];
									break;
								default:
									break;
							}
							
							
						}
						
					}
				}
				
				//even transfer
				//odd transfer
				for(int j=0;j<lGroup.size();j+=2)
				{
					for(int k=1;k<rGroup.size();k+=2)
					{
						
					}
				}
				
				//periodic boundaries
				//odd transfer
				for(int j=0;j<lGroup.size();j+=2)
				{
					for(int k=1;k<rGroup.size();k+=2)
					{
						
					}
				}
				
				//small buffer for number of incoming indices
				std::vector<int> lGather(nThreads, 0);
				std::vector<int> rGather(nThreads, 0);
				
				int leftEdgeSize=0, rightEdgeSize=0;
				switch(i%3)
				{
					case 0:
						leftEdgeSize=mLXEdge.size();
						rightEdgeSize=mRXEdge.size();
						break;
					case 1:
						leftEdgeSize=mLYEdge.size();
						rightEdgeSize=mRYEdge.size();
						break;
					case 2:
						leftEdgeSize=mLZEdge.size();
						rightEdgeSize=mRZEdge.size();
						break;
					default:
						break;
				}
				MPI_Gather(&leftEdgeSize, 1, MPI_INT, &lGather[0], 1, MPI_INT, 0, commHandles[i]);
				MPI_Gather(&rightEdgeSize, 1, MPI_INT, &rGather[0], 1, MPI_INT, 0, commHandles[i]);
				
				//get some totals for this communicator (tree depth)
				int lGatherTotal=0, rGatherTotal=0;
				for(int j=0;j<lGather.size();j++)
				{
					//if(j!=0 || i<3)
						lGatherTotal+=lGather[i];
					//if(j!=lGather.size()-1 || i<3)
						rGatherTotal+=rGather[i];
				}
				
				//here is the new number of particles
				int nPOffset=nP;
				nP+=lGatherTotal+rGatherTotal;
				
				//make some room
				sGroup.push_back(lGatherTotal+rGatherTotal);
				if(iP.size()<nP)
					iP.resize(nP);
				if(p.size()<nP)
					p.resize(nP);
				if(v.size()<nP)
					v.resize(nP);
				if(a.size()<nP)
					a.resize(nP);
				int *leftEdge=NULL, *rightEdge=NULL;
				switch(i%3)
				{
					case 0:
						leftEdge=&mLXEdge[0];
						rightEdge=&mRXEdge[0];
						break;
					case 1:
						leftEdge=&mLYEdge[0];
						rightEdge=&mRYEdge[0];
						break;
					case 2:
						leftEdge=&mLZEdge[0];
						rightEdge=&mRZEdge[0];
						break;
					default:
						break;
				}
				MPI_Gatherv(leftEdge, mLXEdge.size(), MPI_INT, &iP[nPOffset], &lGather[0], MPI_INT, 0, commHandles[i]);
				MPI_Gatherv(rightEdge, mRXEdge.size(), MPI_INT, &iP[nPOffset+lGatherTotal], &rGather[0], MPI_INT, 0, commHandles[i]);
			}
		};
		
		//takes ~log((nP)^(2/3)) time
		void accumulateAccelerationsLeft()
		{
			//accelerations to distribute among communicators,
			// this only needs to exist long enough to perform the transfer!
			std::vector< threeVector<T> > LXEdge;
			std::vector< threeVector<T> > LYEdge;
			std::vector< threeVector<T> > LZEdge;
			std::vector< threeVector<T> > RXEdge;
			std::vector< threeVector<T> > RYEdge;
			std::vector< threeVector<T> > RZEdge;
			
			
		};
		
		void updatePivots()
		{
			//for each level
			for(int i=commHandles.size()-1;i>2;i++)
			{
				int commRank, commSize;
				MPI_Comm_rank(commHandles[i], &commRank); 
				MPI_Comm_size(commHandles[i], &commSize);
				
				int dim=i%3;
				
				if(commRank%2==0)//even ranks decide bounds 
				{
					if(commRank+1<commSize)//not the last communicator
					{
						
					}
					//last communicator doesn't change right bounds
				}
				else//odd ranks send bounds
				{
					//MPI_Send(&,
				}
			}
			
			//finish up blocks that have edges at boundaries
		};
		
		void scatter()
		{
			
		};
		
		void gather()
		{
			
		};
		
		position<T> *getPositions() {return &p[0];};
		threeVector<T> *getVelocities() {return &v[0];};
		threeVector<T> *getAccelerations() {return &a[0];};
		threeVector<T> readLBound() {return lBound;};//user might want this
		threeVector<T> readRBound() {return rBound;};//user might want this
		threeVector<T> readSize() {return size;};//overall size
		int *getIndices() {return &iP[0];};
		int readNParticles() {return nP;};
		int readNSelfParticles() {return nPSelf;};
		
		twoVector<int> *getBondList() {return &boL[0];};
		int readNBonds() {return nBo;};
		twoVector<T> *getBondConstants() {return &boC[0];};
		int *getBondConstantIndices() {return &boCI[0];};
		int readNBondConstants() {return nBoC;};
		
		threeVector<int> *getBendList() {return &beL[0];};
		int readNBends() {return nBe;};
		twoVector<T> *getBendConstants() {return &beC[0];};
		int *getBendConstantIndices() {return &beCI[0];};
		int readNBendConstants() {return nBeC;};
		
	private:
		//master thread list, these really should be const, but with initComm,
		// this won't happen.
		position<T> *masterP;
		threeVector<T> *masterV;
		threeVector<T> *masterA;
		threeVector<T> masterSize;
		int masterNP;
		
		twoVector<int> *masterBondList;
		int masterNBonds;
		twoVector<T> *masterBondConst;
		int *masterBondConstIndex;
		int masterNBondConst;
		
		threeVector<int> *masterBendList;
		int masterNBends;
		twoVector<T> *masterBendConst;
		int *masterBendConstIndex;
		int masterNBendConst;
		
		//current slave thread list
		std::vector<position<T> > p;
		std::vector<threeVector<T> > v;
		std::vector<threeVector<T> > a;
		threeVector<T> lBound;
		threeVector<T> rBound;
		threeVector<int> boundLevel;//each bound gets a comm bound level
		threeVector<T> size;
		std::vector<int> iP;
		int nP;
		int nPSelf;
		
		std::vector<twoVector<int> > boL;
		int nBo;
		std::vector<twoVector<T> > boC;
		std::vector<int> boCI;
		int nBoC;
		
		std::vector<threeVector<int> > beL;
		int nBe;
		std::vector<twoVector<T> > beC;
		std::vector<int> beCI;
		int nBeC;
		
		//we only need the source offset, the order never changes
		std::vector<int> sGroup;
		
		//communicator source offset, how much is up for edge call
		//it only looks like this because c++03 doesn't allow anonymous unions
		//std::vector<int> sLXEdge;
		//std::vector<int> sLYEdge;
		//std::vector<int> sLZEdge;
		//std::vector<int> sRXEdge;
		//std::vector<int> sRYEdge;
		//std::vector<int> sRZEdge;
		
		//which of my own particles is up for edge call
		//it only looks like this because c++03 doesn't allow anonymous unions
		std::vector<int> mLXEdge;
		std::vector<int> mLYEdge;
		std::vector<int> mLZEdge;
		std::vector<int> mRXEdge;
		std::vector<int> mRYEdge;
		std::vector<int> mRZEdge;
		
		//#ifdef MPI_VERSION
		//MPI_Group groupHandles[6];
		//#endif
		
		//#ifdef MPI_VERSION
		MPI_Datatype MPI_P_Type;
		MPI_Datatype MPI_ThreeVectorFloat_Type;
		MPI_Datatype MPI_TwoVectorFloat_Type;
		MPI_Datatype MPI_ThreeVectorInt_Type;
		MPI_Datatype MPI_TwoVectorInt_Type;
		
		std::vector<MPI_Group> groupHandles;
		std::vector<MPI_Comm> commHandles;
		std::vector<std::vector<int> > worldGroupCon;
		//#endif
};
