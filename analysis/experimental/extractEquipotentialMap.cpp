//This is a very short and simple program that shows you how to modify a variable already present in an mpd script.
//The reason this is used is because some variables (forces, potentials, geometry, etc...) are not easily modified
// in the mpd script. Using a file that opens the mpd script, modifies, and then saves it (as below) is far easier
// than modifying the actual variable. Most of it looks like the run.cpp file, it just replaces the "run loop"
// with the modification.

#define NO_REQUIRED_COMMANDS

//For molecular dynamics forces and potentials
#include "../include/MD.h"

//For the molecular dynamics variables
#include "../include/system.h"

#include <sstream>

int main(int argc, char* argv[])
{
	if(argc!=4 || argc!=7)
	{
		//so simple, this can't possibly mess it up
		std::cout << "usage: " << argv[0] << " name type density\n";
		std::cout << "\tOutputs potential ranges. Density is map granularity.\n";
		std::cout << "usage: " << argv[0] << " name type density lowPotential highPotential steps\n";
		std::cout << "\tOutputs file equipotential_name.xyz.\n";
		return 0;
	}
	
	///read in options
	char *name=argv[1];
	int type,steps;
	double density, lowPotential, highPotential;
	std::stringstream cmdArg;
	for(int i=1;i<argc;i++)
		cmdArg << argv[i] << ' ';
	cmdArg >> type >> density;
	if(argc==7)
		cmdArg >> lowPotential >> highPotential >> steps;
	
	///Configuration variables
	Blob<double> System;
	
	///load old variables, then initialize them
	Script<double, Blob <double> > fileIO(name,std::ios::in,&System);
	fileIO.read();
	fileIO.close();
	
	std::string newName("equipotential_");
	newName+=name;
	newName+=".xyz";
	
	position<double> *p=System.getPositions();
	int nParticles=System.readNParticles();
	threeVector<double> size=System.readSize();
	double cutoff=System.readCutoff();
	double cutoffSquared=cutoff*cutoff;
	double *constants=System.getTwoBodyUconst();
	
	//Get the parameters defining equipotential map
	std::map<unsigned int, double> mapElements;//each accessible map element
	
	fourVector<int> nCell;
	double lineDensity=pow(density, 0.3333333);
	nCell.x=size.x/lineDensity;
	nCell.y=size.y/lineDensity;
	nCell.z=size.z/lineDensity;
	nCell.t=nCell.x*nCell.y*nCell.z;
	
	threeVector<double> cellSize;
	cellSize.x=size.x/static_cast<double>(nCell.x);
	cellSize.y=size.y/static_cast<double>(nCell.y);
	cellSize.z=size.z/static_cast<double>(nCell.z);
	
	threeVector<int> octetBounds;
	octetBounds.x=cutoff/cellSize.x+1;
	octetBounds.y=cutoff/cellSize.y+1;
	octetBounds.z=cutoff/cellSize.z+1;
	
	//create map
	for(int i=0;nParticles;i++)
	{
		//get the closest coordinate
		threeVector<int> iCell;
		iCell.x=p[i].x/cellSize.x;
		iCell.y=p[i].y/cellSize.y;
		iCell.z=p[i].z/cellSize.z;
		
		//look through all nearby test points in all octets
		for(int x=iCell.x-octetBounds.x;x<=octetBounds.x;x++)
		for(int y=iCell.y-octetBounds.y;y<=octetBounds.y;y++)
		for(int z=iCell.z-octetBounds.z;z<=octetBounds.z;z++)
		{
			//physical location
			position<double> testPoint;
			testPoint.type=type;
			testPoint.x=static_cast<double>(x%nCell.x)*cellSize.x;
			testPoint.y=static_cast<double>(y%nCell.y)*cellSize.y;
			testPoint.z=static_cast<double>(z%nCell.z)*cellSize.z;
			
			//distance
			threeVector<double> d;
			d.x=p[i].x-testPoint.x;
			d.x=p[i].x-testPoint.x;
			d.x=p[i].x-testPoint.x;
			
			//don't even try if test point is out of range or zero potential
			double potential=nonBondedP(d, constants, cutoffSquared);
			if(potential!=0)
			{
				//retrieve mapped location, and test for existence
				int indexKey=x%nCells.x+(y%nCells.y)*nCells.x+(z%nCells.z)*nCells.x*nCells.y;
				
				std::map<unsigned int, double>::iterator it;
				it=mapElements.lower_bound(indexKey);
				if(!mapElements.key_comp()(it->first, indexKey) || it==mapElements.end())
					mapElements.insert(it, std::map<unsigned int, double>::value_type(indexKey,0));
				
				it->second+=potential;
			}
		}
	}
	
	//create a range map
	std::vector< position<double> > rangeMap;
	double rangeStep=(highPotential-lowPotential)/static_cast<double>(steps);
	for(std::map<unsigned int, unsigned int>::iterator current=cells.begin();
	    current!=cells.end();
	    ++current)
	{
		if(current->second>lowPotential && current->second<highPotential)
		{
			//physical location
			position<double> testPoint;
			testPoint.type=type;
			testPoint.x=static_cast<double>(x%nCell.x)*cellSize.x;
			testPoint.y=static_cast<double>(y%nCell.y)*cellSize.y;
			testPoint.z=static_cast<double>(z%nCell.z)*cellSize.z;
			rangeMap.push_back(
		
	}
	
	xyzFormat<double> xyzFile(p,nParticles);
	xyzFile.open(newName.c_str(), std::ios::out);
	xyzFile.store();
	xyzFile.close();
	return 0;
}

//this is sort of an all in one variant of a cell/verlet list
template <typename T, void FORCE(T cutoffSquared, int nT, position<T> &p1, position<T> &p2, threeVector<T> &a1, threeVector<T> &a2, T *fC)>
std::map<unsigned int, double> computeEquipotentialMap(position<T> *p, 
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
	std::map<unsigned int, double> mapElements;//each accessible map element
	
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
			
			//comapre it to previously loaded elements
			for(int j=0;j<ni;j++)
				FORCE(cutoffSquared,nTypes, tempP[ni], tempP[j], tempA[ni], tempA[j], TBFC);
		}
		
		//unhash our keys
		int x=current->first%nCells.x;
		int y=int(current->first/nCells.x)%nCells.y;
		int z=int(current->first/(nCells.x*nCells.y));
		//look at nearby keys, conservative search
		
		for(int j=0;j<27;j++)
		{
			//compute a nearby key
			int a=j%3-1;
			int b=int(j/3)%3-1;
			int c=int(j/9)-1;
			int nearbyKey=((x+a+nCells.x)%nCells.x)+\
			(((y+b+nCells.y)%nCells.y)*nCells.x)+\
			(((z+c+nCells.z)%nCells.z)*nCells.x*nCells.y);
			
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
					//reduce next acceleration element back to array
				}
			}
		}
	}
	//return input;
}
