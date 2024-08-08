//Computes clusters int T defined by COMPARE given some sorted container S.

#include <iostream>

#ifndef MD_CLUSTER
#define MD_CLUSTER

template <typename T, typename S, bool COMPARE(T,T)>
class cluster {
	public:
		cluster()
		{
			p=NULL;
			c=NULL;
			clusterSize=NULL;
			clusterIndex=NULL;
			stack=NULL;
		};
		cluster(T *pos, S *container, int nPos)
		~cluster();
		void initialize(T *pos, S *container, int nPos);

		int compute();//returns nCluster
		int readNCluster();
		int readClusterSize(int index);
		int readClusterIndex(int index);
	private:
		T *p;
		int nP;
		S *c;
		int nCluster;
		int *clusterSize;
		int *clusterIndex;
		int *stack;
};

template <typename T, typename S, bool COMPARE(T,T)>
cluster<T,S,COMPARE>::cluster(T *pos, S *container, int nPos)
{
	this->initialize(pos,container,nPos);
}

template <typename T, typename S, bool COMPARE(T,T)>
cluster<T,S,COMPARE>::~cluster()
{
	if(clusterSize!=NULL)
		delete clusterSize;
	if(clusterIndex!=NULL)
		delete clusterIndex;
	if(stack!=NULL)
		delete stack;
}

template <typename T, typename S, bool COMPARE(T,T)>
void cluster<T,S,COMPARE>::initialize(T *pos, S *container, int nPos)
{
	p=pos;
	c=container;
	nCluster=0;
	nP=nPos;
	clusterIndex=new int[nP];
	clusterSize=new int[nP];
	stack=new int[nP];
}

template <typename T, typename S, bool COMPARE(T,T)>
int cluster<T,S,COMPARE>::compute()
{
	for(int i=0;i<nP;i++)
		clusterIndex[i]=-1;//null index
	
	for(int i=0;i<nP;i++)
	{
		for(int nStack=0;nStack>0;nStack--)
		{
			
		}
	}
}

template <typename T, typename S, bool COMPARE(T,T)>
int cluster<T,S,COMPARE>::readNCluster()
{
	return nCluster;
}

template <typename T, typename S, bool COMPARE(T,T)>
int cluster<T,S,COMPARE>::readClusterSize(int index)
{
	if(index>=nCluster || index<0)
	{
		std::cout << "Error (cluster.readClusterSize): index " << index << " not a cluster!\n";
		throw 0;
	}
	return clusterSize[index];
}

template <typename T, typename S, bool COMPARE(T,T)>
int cluster<T,S,COMPARE>::readClusterIndex(int index)
{
	if(index>=nP || index<0)
	{
		std::cout << "Error (cluster.readClusterIndex): index " << index << " not a position!\n";
		throw 0;
	}
	return clusterIndex[index];
}

#endif