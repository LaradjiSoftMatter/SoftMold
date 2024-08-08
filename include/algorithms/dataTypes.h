//data types used in common molecular dynamics algorithms


#ifndef MD_DATATYPES
#define MD_DATATYPES

//make some readDat.h for common column types e.g. 1 column, 2 column, 3 column etc...
//overload the functions for headers etc... Use templates etc...

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

//common definitions for all files
//definitions, DPD



///some useful vectors for various properties, most of the library is built on these

//std lib has a pair data type, but this one is super simple.
//Although, I also have twoVector, but it can't have two independent types.
//Although, twoVector should be move assign and move construct safe, as this is.
template <typename T, typename S>
struct keyVal {
	int key;
	int value;
};

//special position vector with a type included
template <typename T>
struct position {
	union {
		struct {T x,y,z;};
		T s[3];
	};

#ifdef __CUDACC__
	__device__ __host__
#endif
	position(T x, T y, T z, int t):x(x),y(y),z(z),type(t){};
	
#ifdef __CUDACC__
	__device__ __host__
#endif
	position(){};
	
#ifdef __CUDACC__
	__device__ __host__
#endif
	int nCoordinates(){return 3;};
	int type;
};

template <typename T>
struct twoVector {
	union {
		struct {T x,y;};
		T s[2];
	};
	
#ifdef __CUDACC__
	__device__ __host__
#endif
	twoVector(){};
#ifdef __CUDACC__
	__device__ __host__
#endif
	twoVector(T init){this->x=init;this->y=init;};
	
#ifdef __CUDACC__
	__device__ __host__
#endif
	T distanceSqr(twoVector vec)
	{
		twoVector d;
		d.x=this->x-vec.x;
		d.y=this->y-vec.y;
		return d.x*d.x+d.y*d.y;
	}

#ifdef __CUDACC__
	__device__ __host__
#endif
	twoVector operator += (twoVector vec)
	{
		this->x+=vec.x;
		this->y+=vec.y;
		return (*this);
	};

#ifdef __CUDACC__
	__device__ __host__
#endif
	twoVector operator + (twoVector vec)
	{
		twoVector buf;
		buf.x=vec.x+this->x;
		buf.y=vec.y+this->y;
		return buf;
	};
	
#ifdef __CUDACC__
	__device__ __host__
#endif
	twoVector operator -= (twoVector vec)
	{
		this->x-=vec.x;
		this->y-=vec.y;
		return (*this);
	};
	
#ifdef __CUDACC__
	__device__ __host__
#endif
	twoVector operator - (twoVector vec)
	{
		twoVector buf;
		buf.x=this->x-vec.x;
		buf.y=this->y-vec.y;
		return buf;
	};
	
#ifdef __CUDACC__
	__device__ __host__
#endif
	twoVector operator /= (twoVector vec)
	{
		this->x/=vec.x;
		this->y/=vec.y;
		return (*this);
	};
	
#ifdef __CUDACC__
	__device__ __host__
#endif
	twoVector operator / (twoVector vec)
	{
		twoVector buf;
		buf.x=this->x/vec.x;
		buf.y=this->y/vec.y;
		return buf;
	};
	
#ifdef __CUDACC__
	__device__ __host__
#endif
	twoVector operator *= (twoVector vec)
	{
		this->x*=vec.x;
		this->y*=vec.y;
		return (*this);
	};
	
#ifdef __CUDACC__
	__device__ __host__
#endif
	twoVector operator * (twoVector vec)
	{
		twoVector buf;
		buf.x=this->x*vec.x;
		buf.y=this->y*vec.y;
		return buf;
	};
	
#ifdef __CUDACC__
	__device__ __host__
#endif
	twoVector operator = (T val)
	{
		this->x=val;
		this->y=val;
		return (*this);
	};
	
	
#ifdef __CUDACC__
	__device__ __host__
#endif
	int nCoordinates(){return 2;};
};

template <typename T>
struct threeVector {
	union {
		struct {T x,y,z;};
		T s[3];
	};
	


#ifdef __CUDACC__
	__device__ __host__
#endif
	threeVector(){};
	
#ifdef __CUDACC__
	__device__ __host__
#endif
	threeVector(T init){this->x=init;this->y=init;this->z=init;};
	
#ifdef __CUDACC__
	__device__ __host__
#endif
	threeVector(T x, T y, T z):x(x),y(y),z(z){};
	
#ifdef __CUDACC__
	__device__ __host__
#endif
	T distanceSqr(threeVector vec)
	{
		threeVector d;
		d.x=this->x-vec.x;
		d.y=this->y-vec.y;
		d.z=this->z-vec.z;
		return d.x*d.x+d.y*d.y+d.z*d.z;
	}
	
#ifdef __CUDACC__
	__device__ __host__
#endif
	T distanceSqr(position<T> vec)
	{
		threeVector d;
		d.x=this->x-vec.x;
		d.y=this->y-vec.y;
		d.z=this->z-vec.z;
		return d.x*d.x+d.y*d.y+d.z*d.z;
	}
	
#ifdef __CUDACC__
	__device__ __host__
#endif
	threeVector operator += (threeVector vec)
	{
		this->x+=vec.x;
		this->y+=vec.y;
		this->z+=vec.z;
		return (*this);
	};
	
#ifdef __CUDACC__
	__device__ __host__
#endif
	threeVector operator + (threeVector vec)
	{
		threeVector buf;
		buf.x=vec.x+this->x;
		buf.y=vec.y+this->y;
		buf.z=vec.z+this->z;
		return buf;
	};

#ifdef __CUDACC__
	__device__ __host__
#endif
	threeVector operator -= (threeVector vec)
	{
		this->x-=vec.x;
		this->y-=vec.y;
		this->z-=vec.z;
		return (*this);
	};
	
#ifdef __CUDACC__
	__device__ __host__
#endif
	threeVector operator - (threeVector vec)
	{
		threeVector buf;
		buf.x=this->x-vec.x;
		buf.y=this->y-vec.y;
		buf.z=this->z-vec.z;
		return buf;
	};
	
#ifdef __CUDACC__
	__device__ __host__
#endif
	threeVector operator /= (threeVector vec)
	{
		this->x/=vec.x;
		this->y/=vec.y;
		this->z/=vec.z;
		return (*this);
	};
	
#ifdef __CUDACC__
	__device__ __host__
#endif
	threeVector operator / (threeVector vec)
	{
		threeVector buf;
		buf.x=this->x/vec.x;
		buf.y=this->y/vec.y;
		buf.z=this->z/vec.z;
		return buf;
	};
	
#ifdef __CUDACC__
	__device__ __host__
#endif
	threeVector operator *= (threeVector vec)
	{
		this->x*=vec.x;
		this->y*=vec.y;
		this->z*=vec.z;
		return (*this);
	};
	
#ifdef __CUDACC__
	__device__ __host__
#endif
	threeVector operator * (threeVector vec)
	{
		threeVector buf;
		buf.x=this->x*vec.x;
		buf.y=this->y*vec.y;
		buf.z=this->z*vec.z;
		return buf;
	};
	
#ifdef __CUDACC__
	__device__ __host__
#endif
	threeVector operator = (T val)
	{
		this->x=val;
		this->y=val;
		this->z=val;
		return (*this);
	};
	
#ifdef __CUDACC__
	__device__ __host__
#endif
	int nCoordinates(){return 3;};
};

template <typename T>
struct fourVector {
	union {
		struct {T x,y,z,t;};
		T s[4];
	};
	
#ifdef __CUDACC__
	__device__ __host__
#endif
	int nCoordinates(){return 4;};
};



///these are generic vectors, etc...
//this is a rigid model of a position vector
template <typename T, int N>
struct genericPosition {
	union {
		struct {T x,y,z;};
		T s[3+N];//N is the extra number of dimensions
	};
	int type;
	
#ifdef __CUDACC__
	__device__ __host__
#endif
	int nCoordinates(){return N;};
};

//this is a generic vector of N dimensions
template <typename T, int N>
struct genericVector {
	T s[N];
	
#ifdef __CUDACC__
	__device__ __host__
#endif
	int nCoordinates(){return N;};
};

//this is like the position vector, but can have more dimensions added to it
template <typename T, int N>
struct typeVector {
	T s[N];
	int type;
	
#ifdef __CUDACC__
	__device__ __host__
#endif
	int nCoordinates(){return N;};
};

/*
//for making polyhedra
template <typename T>
struct polyhedra {
	int nVertices, nEdges, nFaces;
	position<T> **vertices, **edges
	vector< position <T> > **faces;//the mesh is supposed to be a fixed polygon type
};
*/

//this is useful for measuring parameters, it is the ultimate generalized vector
template <typename T, int N>
struct parameter {
	T s[N];
	int i[N];
	
#ifdef __CUDACC__
	__device__ __host__
#endif
	int nCoordinates(){return N;};
};


///Some ideas for an anonymous data type
///Note the lack of type safety, it is passed to the handler
/*
enum datumType {DATUM_INT,DATUM_INT_PTR,DATUM_FLOAT,DATUM_FLOAT_PTR};

//data list type
struct datum
{
	void *var;//the variable
	datumType type;//the type, pass to the handler
};

struct datum
{
	union {
		int integer;
		int *integer_ptr;
		float floating;
		float *floating_ptr;
	};
	datumType type;
};
*/



/*
template <typename T>
class linkList {
	public:
		linkList(T *input, int nElements);
		~linkList();
		bool allocate() {nextElement=new T[n]; lastElement=new T[n];};
		bool next();
		bool next(T *input);
		bool last();
		bool last(T *input);
		T operator[] (int i) {return (*current);};
	private:
		T *current;//current iterator
		T *val;//the information that is present 
		T *nextElement;//pointer to the next iterator
		T *lastElement;//pointer to the last iterator
		int n;
};
*/



//End of MD_DATATYPES
#endif
