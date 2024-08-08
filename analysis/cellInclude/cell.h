//fast version of cell list

//For the threeVector data types
#include "../../include/algorithms/dataTypes.h"

#ifndef ACELL_H
#define ACELL_H

#include <unordered_map>
#include <unordered_set>
#include <forward_list>
#include <stack>

//Cell map data structure, access is cMap[hash]->forward_list
template <typename T>
using umapfwd= typename std::unordered_map<int,std::forward_list<T>>;

//hash the coordinate indices to a cell index
int hash(const threeVector<int> &in, const threeVector<int> &s)
{
	return in.x+in.y*s.x+in.z*s.x*s.y;
}

//unhash the cell index to a coordinate indices
threeVector<int> unhash(const int &in, const threeVector<int> &s)
{
	threeVector<int> out;
	out.x=in%s.x;
	out.y=(in/s.x)%s.y;
	out.z=(in/(s.x*s.y))%s.z;
	return out;
}

//retrieve the neighbor indices
template <typename NCELLS>
std::vector<int> neighIndices(const int &index, const NCELLS &nCells)
{
	threeVector<int> cell=unhash(index,nCells);
	std::vector<int> out;
	//i, j, and k are all the nearby coordinates
	for(int i=-1;i<2;i++)
	{
		for(int j=-1;j<2;j++)
		{
			for(int k=-1;k<2;k++)
			{
				//don't count the central index, off for this sample
				//if(i!=0 && j!=0 && k!=0)
				{
					threeVector<int> key;
					key.x=(cell.x+i+nCells.x)%nCells.x;
					key.y=(cell.y+j+nCells.y)%nCells.y;
					key.z=(cell.z+k+nCells.z)%nCells.z;
					out.push_back(hash(key,nCells));
				}
			}
		}
	}
	return out;
}

//retrieve the forward neighbor indices
template <typename NCELLS>
std::vector<int> neighForwardIndices(const int &index, const NCELLS &nCells)
{
	int nextBox[13][3]={{1,0,0},{0,1,0},{1,1,0},{0,0,1},{1,0,1},{0,1,1},{1,1,1},{-1,0,1},{0,-1,1},{-1,1,1},{1,-1,1},{-1,-1,1},{1,-1,0}};
	threeVector<int> cell=unhash(index,nCells);
	std::vector<int> out;
	//i, j, and k are all the nearby coordinates
	for(int i=0;i<13;i++)
	{
		threeVector<int> key;
		key.x=(cell.x+nextBox[i][0]+nCells.x)%nCells.x;
		key.y=(cell.y+nextBox[i][1]+nCells.y)%nCells.y;
		key.z=(cell.z+nextBox[i][2]+nCells.z)%nCells.z;
		out.push_back(hash(key,nCells));
	}
	return out;
}

//get the coordinate index from a particle position
template <typename T, typename U>
threeVector<int> getCell(const T &p, const U cellSize)
{
	threeVector<int> cell;
	cell.x=std::floor(p.x/cellSize.x);
	cell.y=std::floor(p.y/cellSize.y);
	cell.z=std::floor(p.z/cellSize.z);
	return cell;
}

//create a cell map from some particles' positions
template <typename T, typename NCELLS, typename CELLSIZE>
umapfwd<T> createCMap(T *pBegin, T *pEnd, const NCELLS &nCells, const CELLSIZE &cellSize)
{
	umapfwd<T> cMap;
	//for all particles
	for(T *p=pBegin;p!=pEnd;++p)
	{
		NCELLS k=getCell(*p,cellSize);
		int h=hash(k,nCells);
		//locate an existing cell
		auto iter=cMap.find(h);
		if(iter!=cMap.end()) //if found, push particle to forward_list
			iter->second.push_front(*p);
		else //if not found, construct new forward list and push particle to forward list
			cMap[h].push_front(*p);
	}
	return cMap;
}

//create a cell map from some particles' positions' indices
template <typename T, typename NCELLS, typename CELLSIZE>
umapfwd<int> createIMap(T *p, int nP, const NCELLS &nCells, const CELLSIZE &cellSize)
{
	umapfwd<int> iMap;
	//for all particles
	for(int i=0;i<nP;i++)
	{
		NCELLS k=getCell(p[i],cellSize);
		int h=hash(k,nCells);
		//locate an existing cell
		auto iter=iMap.find(h);
		if(iter!=iMap.end()) //if found, push particle to forward_list
			iter->second.push_front(i);
		else //if not found, construct new forward list and push particle to forward list
			iMap[h].push_front(i);
	}
	return iMap;
}

//create a surface data structure, each surface is a collection of particles
template <typename T, typename NCELLS, typename CELLSIZE, typename U>
std::vector<std::vector<T>> 
createSurfaces(umapfwd<T> cMap, const NCELLS &nCells, const CELLSIZE &cellSize, const CELLSIZE &size, const U &cutoff)
{
	std::vector<std::vector<T>> surfaces;
	U cutoffSqr=cutoff*cutoff;
	//for all cells
	for(auto &c:cMap)
	{
		//make sure cell isn't empty
		while(!c.second.empty())
		{
			std::vector<T> surface;
			//create and initialize a stack for the particles' positions
			std::stack<T> stk;
			stk.push(c.second.front());
			//remove particle from cell
			c.second.pop_front();
			while(stk.size()!=0)
			{
				//too much data!
				if(stk.size()>10000000)
					throw 0;
				//grab the particle off the top
				T buf=stk.top();
				surface.push_back(buf);
				stk.pop();
				//get the cell's coordinates
				NCELLS k=getCell(buf,cellSize);
				int cHash=hash(k,nCells);
				//get the neighbors
				auto neigh=neighIndices(cHash, nCells);
				//for all neighbors
				for(auto n:neigh)
				{
					//get the forward_list
					auto &c=cMap[n];
					//for all particles in forward_list
					// last keeps track of the previous position for deletion
					for(auto buf2=c.begin(),last=c.before_begin();
						buf2!=c.end();last=buf2,buf2++)
					{
						//distance with minimum image
						threeVector<U> d;
						d.x=buf.x-buf2->x;
						d.y=buf.y-buf2->y;
						d.z=buf.z-buf2->z;
						d.x-=d.x>size.x/2.0?size.x:0.0;
						d.y-=d.y>size.y/2.0?size.y:0.0;
						d.z-=d.z>size.z/2.0?size.z:0.0;
						d.x+=d.x<-size.x/2.0?size.x:0.0;
						d.y+=d.y<-size.y/2.0?size.y:0.0;
						d.z+=d.z<-size.z/2.0?size.z:0.0;
						U dist=d.x*d.x+d.y*d.y+d.z*d.z;
						//if it is within range
						if(dist<cutoffSqr)
						{
							//push onto stack
							stk.push(*buf2);
							//delete current particle
							c.erase_after(last);
							//reset current particle to previous
							buf2=last;
						}
					}
				}
			}
			//accumulate the surfaces
			surfaces.push_back(surface);
		}
	}
	return surfaces;
}

//sort surface lists by size
template <typename T>
bool sortBySize(const std::vector<position<T>> &a, const std::vector<position<T>> &b)
{
	return a.size()>b.size();
}

template <typename T>
double getAngle(const T &in, const T &com, const int &axis)
{
	int x=(axis+1)%3;
	int y=(axis+2)%3;
	threeVector<double> dac;
	dac.s[x]=in.s[x]-com.s[x];
	dac.s[y]=in.s[y]-com.s[y];
	return atan2(dac.s[y],dac.s[x]);
}

template <typename T>
T centerOfMass(T *begin, T *end)
{
	T com;
	com.x=0;
	com.y=0;
	com.z=0;
	for(T *p=begin;p!=end;++p)
	{
		com.x+=p->x;
		com.y+=p->y;
		com.z+=p->z;
	}
	int n=end-begin;
	if(n>0)
	{
		com.x/=n;
		com.y/=n;
		com.z/=n;
	}
	return com;
}

template <typename T, typename U>
std::vector<std::vector<T>> createAxisSlices(T *begin, T *end, const int &axis, const U &sliceThickness,
	const int nAngles)
{
	T com=centerOfMass(begin,end);
	T zero;
	zero.x=0;
	zero.y=0;
	zero.z=0;
	std::vector<std::vector<T>> slices;
	for(T *p=begin;p!=end;++p)
	{
		int slice=p->s[axis]/sliceThickness;
		while(slices.size()<=slice)
			slices.push_back(std::vector<T>());
		slices[slice].push_back(*p);
	}
	U dAngle=2.0*M_PI/nAngles;
	for(auto &slice:slices)
	{
		std::vector<T> mergedSlice(nAngles,zero);
		std::vector<int> countSlice(nAngles,0);
		U cutoffSqr=sliceThickness*sliceThickness;
		for(T *p=&(*slice.begin());p!=&(*slice.end());++p)
		{
			U angle=getAngle(*p,com,axis);
			int iSlice=(angle+M_PI)/dAngle;
			if(iSlice==nAngles)
				iSlice=0;
			mergedSlice[iSlice].x+=p->x;
			mergedSlice[iSlice].y+=p->y;
			mergedSlice[iSlice].z+=p->z;
			countSlice[iSlice]++;
		}
		for(int i=0;i<mergedSlice.size();i++)
		{
			if(countSlice[i]!=0)
			{
				mergedSlice[i].x/=countSlice[i];
				mergedSlice[i].y/=countSlice[i];
				mergedSlice[i].z/=countSlice[i];
			}
		}
		slice=mergedSlice;
	}
	return slices;
}
#endif
