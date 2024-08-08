//this library is dependent on the following libraries (for their functions)
//multiprocessing
#ifdef _OPENMP
	#include <omp.h>
#endif
//stl, might as well include the whole damn thing
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <map>
#include <cstdlib>
#include <algorithm>
#include <vector>
#include <iomanip>
//#include <GL/glut.h>
//#include <GL/gl.h>
//#include <SDL/SDL.h>
//included in directory
#include "MersenneTwister.h"
//miscelanious molecular dynamics includes
#include "dataTypes.h"

#ifdef CUDA_TYPES
#include <thrust/sort.h>
#endif

#ifndef MD_FUNCTIONS
#define MD_FUNCTIONS


///Macros
#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

#ifndef minimumImage
#define minimumImage(X1,X2,DX,MAX) DX=X1-X2; if(DX>=MAX/2.0) DX=DX-MAX; if(DX<-MAX/2.0) DX=DX+MAX;
#endif

///inline functions
template <typename T>
inline double distance(T a, T b)
{
	T d;
	for(int i=0;i<d.nCoordinates();i++)
		d.s[i]=a.s[i]-b.s[i];
	double sum=0;
	for(int i=0;i<d.nCoordinates();i++)
		sum+=d.s[i]*d.s[i];
	return sqrt(sum);
}

//minimum image version
template <typename T>
inline double distance(T a, T b, T size)
{
	T d;
	for(int i=0;i<d.nCoordinates();i++)
	{
		d.s[i]=a.s[i]-b.s[i];
		d.s[i]-=(d.s[i]>size.s[i]/2.0)?size.s[i]:0.0;
		d.s[i]+=(d.s[i]<-size.s[i]/2.0)?size.s[i]:0.0;
	}
	double sum=0;
	for(int i=0;i<d.nCoordinates();i++)
		sum+=d.s[i]*d.s[i];
	return sqrt(sum);
}

inline int unique3DHash(threeVector<int> index, threeVector<int> limit)
{
	return index.x+index.y*limit.x+index.z*limit.y*limit.x;
}

inline int unique2DHash(twoVector<int> index, twoVector<int> limit)
{
	return index.x+index.y*limit.x;
}

template <typename T, typename S>
inline bool sortByKey(keyVal<T,S> a, keyVal<T,S> b)
{
	return a.key<b.key;
}

template <typename T, typename S>
inline bool sortByValue(keyVal<T,S> a, keyVal<T,S> b)
{
	return a.value<b.value;
}

template <typename T, int dim>
inline bool sortByDim(T a, T b)
{
	return a.s[dim]<b.s[dim];
}

///math functions
#ifndef PARTIAL_SUM_LIMIT
#define PARTIAL_SUM_LIMIT 100
#endif
template <typename T>
T pairwiseSum(T *list, int start, int end)
{
        if(end-start<=PARTIAL_SUM_LIMIT)
        {
                T sum=0;
                for(int i=start;i<end;i++)
                        sum+=list[i];
                return sum;
        }
        else
        {
                return pairwiseSum(list, start, (start+end)/2)+pairwiseSum<T>(list, (start+end)/2+1, end);
        }
}

template <typename T>
T partialSum(T *list, int start, int end)
{
	if(end==start)
		return list[end];
	else
		return partialSum(list,start,(start+end)/2)+partialSum(list,(start+end)/2+1,end);
}

template <typename T>
T average(T *list, int nElements)
{
	T avg=partialSum<T>(list,0,nElements-1);
	avg/=static_cast<T>(nElements);
	return avg;
}

template <typename T>
T quickAverage(T *list, int nElements)
{
	T avg=pairwiseSum<T>(list,0,nElements-1);
	avg/=static_cast<T>(nElements);
	return avg;
}

//vector functions
template <typename T>
T dotProduct(fourVector<T> a, fourVector<T> b)
{
	return (a.x*b.x)+(a.y*b.y)+(a.z*b.z)+(a.t*b.t);
}

template <typename T>
T dotProduct(threeVector<T> a, threeVector<T> b)
{
	return (a.x*b.x)+(a.y*b.y)+(a.z*b.z);
}

template <typename T>
T dotProduct(position<T> a, position<T> b)
{
	return (a.x*b.x)+(a.y*b.y)+(a.z*b.z);
}

template <typename T>
T dotProduct(position<T> a, threeVector<T> b)
{
	return (a.x*b.x)+(a.y*b.y)+(a.z*b.z);
}

template <typename T>
T dotProduct(threeVector<T> a, position<T> b)
{
	return (a.x*b.x)+(a.y*b.y)+(a.z*b.z);
}

template <typename T>
T dotProduct(twoVector<T> a, twoVector<T> b)
{
	return (a.x*b.x)+(a.y*b.y);
}

template <typename T>
double dotProduct(T a, T b)
{
	double sum=0;
	for(int i=0;i<a.nCoordinates();i++)
		sum+=(a.s[i]*b.s[i]);
	return sum;
}

template <typename T, typename S>
double dotProduct(T a, S b)
{
	double sum=0;
	for(int i=0;i<a.nCoordinates() && i<b.nCoordinates();i++)
		sum+=(a.s[i]*b.s[i]);
	return sum;
}
/*
template <typename T>
T exteriorProduct(T a, T b)
{
	T result;
	for(int i=0;i<a.nCoordinates();i++)
		result.s[i]=
*/
template <typename T>
threeVector<T> crossProduct(threeVector<T> a, threeVector<T> b)
{
	threeVector<T> result;
	result.x=(a.y*b.z)-(a.z*b.y);
	result.y=-((a.x*b.z)-(a.z*b.x));
	result.z=(a.y*b.x)-(a.x*b.y);
	return result;
}

template <typename T>
position<T> crossProduct(position<T> a, position<T> b)
{
	position<T> result;
	result.x=(a.y*b.z)-(a.z*b.y);
	result.y=-((a.x*b.z)-(a.z*b.x));
	result.z=(a.y*b.x)-(a.x*b.y);
	return result;
}
/*
template <typename T>
T crossProduct(T a, T b)
{
	T result;
	result.x=(a.y*b.z)-(a.z*b.y);
	result.y=-((a.x*b.z)-(a.z*b.x));
	result.z=(a.y*b.x)-(a.x*b.y);
	return result;
}*/

template <typename T, typename V>
V crossProduct(T a, T b)
{
	T result;
	result.x=(a.y*b.z)-(a.z*b.y);
	result.y=-((a.x*b.z)-(a.z*b.x));
	result.z=(a.y*b.x)-(a.x*b.y);
	return result;
}

/** \brief Returns the magnitude of a vector.
 * Takes a vector (T a) and returns a double
 */
template <typename T>
double magnitude(T a)
{
	double sum=0;
	for(int i=0;i<a.nCoordinates();i++)
		sum+=a.s[i]*a.s[i];
	return sqrt(sum);
}

//Generic unitVector
template <typename T>
T unitVector(T a)
{
	double m=magnitude<T>(a);
	for(int i=0;i<a.nCoordinates();i++)
		a.s[i]/=m;
	return a;
}




//these are some transformation of coordinate functions
/** \brief Doesn't do anything to transform the vector.
 * Null transformation on input (T a) to output upon return.
 */
template <typename T>
T noTransform(T a)
{
	//this kind of does nothing
	return a;
}

/** \brief Template function that converts cartesian to cylindrical coordinates.
 * This template function converts cartesian coordinates (T a) to cylindrical coordinates upon return.
 * Requires typename T to have x, y, and z members.
 */
template <typename T>
T cartesianToCylindrical(T a)
{
	//transform into cylindrical coordinates
	T c;
	c.x=sqrt(a.x*a.x+a.y+a.y);
	c.y=atan(abs(a.y/a.x));
	c.z=a.z;
	return c;
}

/** \brief Template function that converts cylindrical to cartesian coordinates.
 * This template function converts cylindrical coordinates (T a) to cartesian coordinates upon return.
 * Requires typename T to have x, y, and z members.
 */
template <typename T>
T cylindricalToCartesian(T a)
{
	//transform into cartesian coordinates
	T c;
	c.x=a.x*cos(a.y);
	c.y=a.x*sin(a.y);
	c.z=a.z;
	return c;
}

/** \brief Template function that converts cartesian to spherical coordinates.
 * This template function converts cartesian coordinates (T a) to spherical coordinates upon return.
 * Requires typename T to have x, y, and z members.
 */
template <typename T>
T cartesianToSpherical(T a)
{
	//transform into spherical coordinates
	T c;
	c.x=sqrt(a.x*a.x+a.y*a.y+a.z*a.z);
	c.y=atan(sqrt(a.x*a.x+a.y*a.y)/a.z);
	c.z=atan(a.y/a.x);
	return c;
}

/** \brief Template function that converts spherical to cartesian coordinates.
 * This template function converts spherical coordinates (T a) to cartesian coordinates upon return.
 * Requires typename T to have x, y, and z members.
 */
template <typename T>
T sphericalToCartesian(T a)
{
	T c;
	c.x=a.x*sin(a.x)*cos(a.z);
	c.y=a.x*sin(a.x)*sin(a.z);
	c.z=a.x*cos(a.x);
	return c;
}

/** \brief Template function that returns the average center of mass.
 * This template function evaluates the average center of mass from a list of
 * particles (T *p). Takes the number of particles (int nParticles) as a parameter as well.
 * Returns the vector for center of mass.
 */
template <typename T>
T com(T *p, int nParticles)
{
	T comValue;
	for(int j=0;j<comValue.nCoordinates();j++)
		comValue.s[j]=0;
	
	for(int i=0;i<nParticles;i++)
		for(int j=0;j<p[i].nCoordinates();j++)
			comValue.s[j]+=p[i].s[j];
	
	if(nParticles!=0)
		for(int j=0;j<p[0].nCoordinates();j++)
			comValue.s[j]/=static_cast<double>(nParticles);//should I static cast here???????!?????????!
	return comValue;
}

/** \brief Template function that returns the average distance from the center of mass.
 * This template function determines the average radius of a group of points. It takes
 * a list of points (T *p) and the number of particles (int nParticles) in the list. It
 * returns the average radius (as a double floating point value).
 */
template <typename T>
double radius(T *p, int nParticles)
{
	T comValue=com(p,nParticles);
	double radiusValue=0;
	for(int i=0;i<nParticles;i++)
		radiusValue+=distance(p[i],comValue);
	if(nParticles!=0)
		radiusValue/=static_cast<double>(nParticles);
	return radiusValue;
}

/** \brief This template function performs a gimbal like rotation along three axis.
 * This template function takes a template type, typically a position or threeVector, and
 * performs a rotation about the origin. It takes a position list (pointer offset),
 * number of particles in position list, and a vector (same template type as
 * the function) variable. The theta vector value is the theta rotation along each
 * axis. It performs them in an x, y, and z order. If you would like to do them
 * out of order, just perform a rotation along only the y or z axis, then call the
 * function again to perform a rotation along the x and/or y axis. Requires an object
 * with x, y, and z members (readable and writable).
 */
template <typename T>
void gimbalRotate(T *p, int nParticles, T theta)
{
	//Perform rotation (one axis at a time, you could combine them into one operation):
	for(int i=0;i<nParticles;i++)
	{
		//Temporary position value:
		T newP;
		
		//Each rotation must be represented seperately, unless you utilize a different
		// representation. Truthfully you only need 2 axis for the rotation, the third
		// one is just for the hell of it.
		
		//Rotate with respect to the x axis:
		newP.y=p[i].y*cos(theta.x)-p[i].z*sin(theta.x);
		newP.z=p[i].y*sin(theta.x)+p[i].z*cos(theta.x);
		p[i].y=newP.y;
		p[i].z=newP.z;
		
		//Rotate with respect to the y axis:
		newP.x=p[i].x*cos(theta.y)+p[i].z*sin(theta.y);
		newP.z=-p[i].x*sin(theta.y)+p[i].z*cos(theta.y);
		p[i].x=newP.x;
		p[i].z=newP.z;
		
		//Rotate with respect to the z axis:
		newP.x=p[i].x*cos(theta.z)-p[i].y*sin(theta.z);
		newP.y=p[i].x*sin(theta.z)+p[i].y*cos(theta.z);
		p[i].x=newP.x;
		p[i].y=newP.y;
		
		//A better version:
		//Rotate with respect to the x and y axis:
		//newP.x=p[i].x*cos(theta.y)+p[i].z*sin(theta.y);
		//p[i].y=p[i].y*cos(theta.x)-p[i].z*sin(theta.x);
		//newP.z=p[i].y*sin(theta.x)-p[i].x*sin(theta.y)+p[i].z*(cos(theta.x)+cos(theta.y));
		//p[i].x=newP.x;
		//p[i].z=newP.z;
	}
}

/** \brief This template function performs a gimbal like rotation along three axis on the center of mass.
 * This template function takes a template type, typically a position or threeVector, and
 * performs a rotation about the center of mass. It takes a particle list (T *p),
 * number of particles (int nParticles) in the particle list, and a vector (same template type as
 * the function) variable. The theta value is the theta rotation along each
 * axis. It performs them in an x, y, and z order. If you would like to do them
 * out of order, just perform a rotation along only the y or z axis, then call the
 * function again to perform a rotation along the x and/or y axis. This function 
 * returns the center of mass.  Requires an object with x, y, and z members (readable and writable).
 */
template <typename T>
T gimbalRotateCOM(T *p, int nParticles, T theta)
{
	T comValue=com<T>(p,nParticles);
	
	//Perform rotation (one axis at a time, you could combine them into one operation):
	for(int i=0;i<nParticles;i++)
	{
		//Adjust system to place center of mass at the origin:
		T newP;
		p[i].x-=comValue.x;
		p[i].y-=comValue.y;
		p[i].z-=comValue.z;
		
		//Each rotation must be represented seperately, unless you utilize a different
		// representation. Truthfully you only need 2 axis for the rotation, the third
		// one is just for the hell of it.
		
		//Rotate with respect to the x axis:
		newP.y=p[i].y*cos(theta.x)-p[i].z*sin(theta.x);
		newP.z=p[i].y*sin(theta.x)+p[i].z*cos(theta.x);
		p[i].y=newP.y;
		p[i].z=newP.z;
		
		//Rotate with respect to the y axis:
		newP.x=p[i].x*cos(theta.y)+p[i].z*sin(theta.y);
		newP.z=-p[i].x*sin(theta.y)+p[i].z*cos(theta.y);
		p[i].x=newP.x;
		p[i].z=newP.z;
		
		//Rotate with respect to the z axis:
		newP.x=p[i].x*cos(theta.z)-p[i].y*sin(theta.z);
		newP.y=p[i].x*sin(theta.z)+p[i].y*cos(theta.z);
		p[i].x=newP.x;
		p[i].y=newP.y;
		
		//A better version:
		//Rotate with respect to the x and y axis:
		//newP.x=p[i].x*cos(theta.y)+p[i].z*sin(theta.y);
		//p[i].y=p[i].y*cos(theta.x)-p[i].z*sin(theta.x);
		//newP.z=p[i].y*sin(theta.x)-p[i].x*sin(theta.y)+p[i].z*(cos(theta.x)+cos(theta.y));
		//p[i].x=newP.x;
		//p[i].z=newP.z;
		
		//Adjust the system back to the original center of mass position
		p[i].x+=comValue.x;
		p[i].y+=comValue.y;
		p[i].z+=comValue.z;
	}
	
	return comValue;
}

/** \brief Moves a group of particles. 
 * Moves a group (int nParticles) of particles (T *p) by a displacement (T dir).
 * 
 */
template <typename T>
void moveParticles(T *p, int nParticles, T dir)
{
	for(int i=0;i<nParticles;i++)
		for(int j=0;j<p[i].nCoordinates();j++)
			p[i].s[j]+=dir.s[j];
}

/** \brief Moves a group of particles. 
 * Moves a group (int nParticles) of particles (T *p) to a new position (T pos).
 * Returns center of mass vector.
 */
template <typename T>
T moveParticlesCOM(T *p, int nParticles, T pos)
{
	T comValue=com<T>(p, nParticles);
	
	for(int i=0;i<nParticles;i++)
		for(int j=0;j<p[i].nCoordinates();j++)
			p[i].s[j]+=(pos.s[j]-comValue.s[j]);
	return comValue;
}

/** \brief Check two values using X value as pivot.
 * Meant for sorting, something like 'threeVector *list; sort(list,list+nElements,compareX);'
 */
template <typename T>
bool compareX(T a, T b)
{
	return a.x<b.x;
}

/** \brief Check two values using Y value as pivot.
 * Meant for sorting, something like 'threeVector *list; sort(list,list+nElements,compareY);'
 */
template <typename T>
bool compareY(T a, T b)
{
	return a.y<b.y;
}

/** \brief Check two values using Z value as pivot.
 * Meant for sorting, something like 'threeVector *list; sort(list,list+nElements,compareZ);'
 */
template <typename T>
bool compareZ(T a, T b)
{
	return a.z<b.z;
}



/** \brief Finds minimum distance from a point to a plane.
 * Takes a point (with x,y, and z public members) and a polygon of sorts
 */
template <typename T>
double pointToPlaneDistance(T point, T *polygon, int nVertices)
{
	//find the closest vertex
	//int nearestVert=0;
	double nearestDist=distance(polygon[0],point);
	
	
	
	for(int i=0;i<nVertices;i++)
	{
		
	}
	return 0;
}

/** \brief Cluster stuff together.
 * stuff is a class with member functions nElements() and compare(int, int).
 * All integers are relative indices.
 */
template <typename S>
std::vector< std::vector<int> > cluster(S &stuff)
{
	//create a stack to draw from
	std::vector<int> pStack;
	
	//fill it with indices
	for(int i=0;i<stuff.nElements();++i)
		pStack.push_back(i);
	
	//initialize a group of clusters
	std::vector< std::vector<int> > cStack;
	
	//while elements of pStack are still present (we assume axiom of choice)
	while(pStack.size()>0)
	{
		//create a new stack for clustering
		std::vector<int> stack;
		
		//put the initial one onto the stack
		stack.push_back(pStack.back());
		
		//remove the initial one
		pStack.pop_back();
		
		//go through all elements of pStack, starting from the last one
		for(int i=0;i<stack.size();++i)
		{
			for(int j=pStack.size()-1;j>=0;--j)
			{
				//is it part of the cluster?
				if(i!=j && stuff.compare(stack[i],pStack[j]))
				{
					//swap it to the end
					int buf=pStack.back();
					pStack.back()=pStack[j];
					pStack[j]=buf;
					
					//push the element onto stack
					stack.push_back(pStack.back());
					
					//remove it
					pStack.pop_back();
				}
			}
		}
		
		//push the cluster onto the stack
		cStack.push_back(stack);
	}
	
	//return our clusters
	return cStack;
}

/** \brief Cluster stuff together with control of nearby elements.
 * stuff is a class with member functions nElements(), compare(int, int), and next(int).
 * All integers are relative indices. Requires next() to be the next element on circularly linked list.
 */
template <typename S>
std::vector< std::vector<int> > clusterNext(S &stuff)
{
	//flag's elements when they are clustered
	std::vector<bool> flag;
	
	//initialize all flags to false (unclustered classification)
	for(int i=0;i<stuff.nElements();++i)
		flag.push_back(false);
	
	//initialize a group of clusters
	std::vector< std::vector<int> > cStack;
	
	//go through all elements
	for(int i=0;i<stuff.nElements();++i)
	{
		//check if unclustered
		if(!flag[i])
		{
			//create a new stack for the cluster
			std::vector<int> stack;
			
			//push the initial index onto the stack
			stack.push_back(i);
			
			//go through nearby elements
			for(int j=0;j<stack.size();++j)
			{
				for(int k=stuff.next(stack[j]);k>=0;k=stuff.next(stack[j]))
				{
					if(k==1)
					{
						std::cout << k << ' ' << j << ' ' << flag[k] << ' ' << stack[j] << ' ' << stack.size() << std::endl;
						std::cin.get();
					}
					//is the next element unclustered?
					if(!flag[k])
					{
						//is it part of the cluster?
						if(stuff.compare(stack[j],k))
						{
							//push the element onto stack
							stack.push_back(k);
							
							//flag it as  clustered
							flag[k]=true;
						}
					}
				}
			}
			
			//push the cluster onto the stack
			cStack.push_back(stack);
		}
	}
	
	//return our clusters
	return cStack;
}

/*
//this thing can fail without warning
template <typename T>
void addIcosohedran(polyhedra<T> &a, threeVector<T> center)
{
	a.nVertices=12;
	a.nEdges=30;
	a.nFaces=20;
	if(a.edges==NULL)
		a.edges=new position<T> *[a.nEdges];
	if(a.vertices==NULL)
		a.vertices=new position<T> *[a.nVertices];
	if(a.nFaces==NULL)
		a.faces=new <T> *[a.nFaces];
	//these are normalized
	double t = (1+(double)sqrt(5.0))/2, tau = t/(double)sqrt(1.0+t*t), one = 1/(double)sqrt(1.0+t*t);

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
}
*/



//end of MD_FUNCTIONS
#endif
