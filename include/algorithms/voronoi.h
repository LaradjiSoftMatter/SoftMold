//this includes dataTypes.h, the functions needed, and any other libraries with functions needed, like the STL
//Who knows, maybe openGL, boost (meh), CUDA, or whatever will be needed in the future
//conflicts may occur, and I cannot guarantee a conflict free program, but I can take suggestions and hopefully will fix them
#include "functions.h"

#ifndef MD_VORONOI
#define MD_VORONOI


//this template could actually include other things, e.g. <typename T, typename S, typename BILL, typename PHYSICS, typename WHATEVER...>
template <typename T>
class Voronoi {
	public:
		//global functions that should usually be defined, their input varies from class to class
		Voronoi();
		~Voronoi();
		T build();
		T compute();
		T output(T time, char *name);
	private:
		//variables required by this class for internal use
};

//the constructor should include setup, etc...
//so that when the user uses the functions, they don't need to do all the setup, just the input
template <typename T>
Voronoi<T>::Voronoi()
{

}

//the destructor should be defined, but usually used to deallocate memory
template <typename T>
Voronoi<T>::~Voronoi()
{

}

//for rebuilding run time data structures related to execution
//this isn't always needed, but is a good idea if you don't need to build all the time
//e.g. when one data structure lasts 20 time steps before needing to update
template <typename T>
T Voronoi<T>::build()
{

}

//This should just compute and output a value, but you could pass a pointer and find other values
template <typename T>
T Voronoi<T>::compute()
{

}

//This should output whatever compute() does, and takes time, name, and whatever compute needs for arguments
template <typename T>
T Voronoi<T>::output(T time, char *name)
{
	//Ideally, this should be the output, as a single value, and T may be a vector, but anything is possible
	T value=compute();
	
	//storage of the value (or values) should follow this convention
	std::fstream dataFile;
	std::string buf("prototype_");
	buf+=name;
	buf+=".dat";
	dataFile.open(buf.c_str(), std::ios::app | std::ios::out);
	dataFile << time << '\t' << value << std::endl;//potentially more values here, maybe even more files
	dataFile.close();
}

#endif