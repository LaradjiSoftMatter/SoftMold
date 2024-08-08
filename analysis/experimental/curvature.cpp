//This is a very short and simple program that shows you how to modify a variable already present in an mpd script.
//The reason this is used is because some variables (forces, potentials, geometry, etc...) are not easily modified
// in the mpd script. Using a file that opens the mpd script, modifies, and then saves it (as below) is far easier
// than modifying the actual variable. Most of it looks like the run.cpp file, it just replaces the "run loop"
// with the modification.

//For molecular dynamics forces and potentials
#include "../include/MD.h"

//For the molecular dynamics variables
#include "../include/system.h"

struct splineType {
	std::vector<std::vector<double> > spline;
	threeVector<double> orient;
	unsigned int dim() { return 3; };
};

int main(int argc, char* argv[])
{
	if(argc!=3)
	{
		//so simple, this can't possibly mess it up
		std::cout << "usage: " << argv[0] << " name type boxSize\n";
		return 0;
	}
	
	///read in options
	char *name=argv[1];
	int type;
	double boxSize;
	std::stringstream cmdArg;
	cmdArg << argv[2] << ' ' << argv[3];
	cmdArg >> type >> boxSize;
	//int index=atoi(argv[2]);
	
	///Configuration variables
	Blob<double> System;
	
	
	///load old variables, then initialize them
	Script<double, Blob <double> > fileIO(name,std::ios::in,&System);
	fileIO.read();
	fileIO.close();
	
	threeVector<double> size=System.readSize();
	
	//get only particles we want
	std::vector<position<double> > p;
	position<double> *allP=System.getPositions();
	for(int i=0;i<System.readNParticles();i++)
		if(allP[i].type==type)
			p.push_back(allP[i]);
	
	splineType zeroSpline;
	zeroSpline.spline=std::vector<std::vector<double> > (zeroSpline.dim(), std::vector<double>(zeroSpline.dim(), 0.0));
	zeroSpline.orient=0.0;
	
	fourVector<unsigned int> nCells;//housekeeping
	//std::vector< keyVal<unsigned int,unsigned int> > hashIndices;//for setting up verlet lists
	nCells.x=size.x/boxSize;
	nCells.y=size.y/boxSize;
	nCells.z=size.z/boxSize;
	nCells.t=nCells.x*nCells.y*nCells.z;
	threeVector<T> cellSize;
	//get these a little closer to a uniform size
	cellSize.x=size.x/nCells.x;
	cellSize.y=size.y/nCells.y;
	cellSize.z=size.z/nCells.z;
	
	std::map<unsigned int, std::vector<unsigned int> > cells;//accessed by (cell,head element)
	std::map<unsigned int, splineType> splines;//accessed by (cell,spline constants)
	std::map<unsigned int, bool> flags;
	
	//get lists for all particle's locations
	for(unsigned int i=0;i<p.size();i++)
	{
		unsigned int x=int(p[i].x/cellSize.x)%nCells.x;
		unsigned int y=int(p[i].y/cellSize.y)%nCells.y;
		unsigned int z=int(p[i].z/cellSize.z)%nCells.z;
		
		unsigned int key=x+y*nCells.x+z*nCells.x*nCells.y;
		
		std::map<unsigned int, std::list<unsigned int> >::iterator it;
		//it=cells.find(hashIndices[i].key);
		it=cells.lower_bound(key);//supposedly faster (Effective STL)
		//cell linked list
		if(it!=cells.end() && !(cells.key_comp()(it->first, key)))
		{
			//already exists, just push back
			cells[key].push_back(i);
		}
		else
		{
			//initialize operation
			cells.insert(it, std::map<unsigned int, std::vector<unsigned int> >::value_type(key, i));
			
			std::map<unsigned int, bool> flagHint;
			flagHint->first=it->first;
			
			std::map<unsigned int, splineType> splineHint;
			splineHint->first=it->first;
			
			splines.insert(splineHint,std::map<unsigned int, splineType>::value_type(key,zeroSpline));
			flags.insert(flagHint,std::map<unsigned int, bool>::value_type(key,false));
		}
	}
	
	//orient the surfaces
	for(std::map<unsigned int, std::vector<unsigned int> >::iterator current=cells.begin();
	    current!=cells.end();
	    ++current)
	{
		if(!flags[current->first])
		{
			flags[current->first]=true;
			
			std::vector<unsigned int> sStack;
			sStack.push_back(current->first);
			
			//our initial guess
			std::vector<unsigned int> localBox=current->second;
			
			while(sStack.size()>0)
			{
				
			}
		}
	}
	
	//calculate spline for nearby particles, in O(ln(cells.size())+nParticles*density) time
	for(std::map<unsigned int, std::vector<unsigned int> >::iterator current=cells.begin();
	    current!=cells.end();
	    ++current)
	{
		std::vector<unsigned int> localBox=current->second;
		
		//unhash our keys
		int x=current->first%nCells.x;
		int y=int(current->first/nCells.x)%nCells.y;
		int z=int(current->first/(nCells.x*nCells.y));
		
		//look at nearby keys
		for(int a=-1;a<2;++a)
			for(int b=-1;b<2;++b)
				for(int c=-1;c<2;++c)
		{
			//compute a nearby key, k
			int nextX=x+a;
			int nextY=y+b;
			int nextZ=z+c;
			
			unsigned int nearbyKey=nextX%nCells.x+((nextY%nCells.y)*nCells.x)+((nextZ%nCells.z)*nCells.x*nCells.y);
			std::map<unsigned int, std::vector<unsigned int> >::iterator it;
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
	
	std::string newName("frames_");
	newName+=name;
	newName+=".xyz";
	
	position<double> *p=System.getPositions();
	int nParticles=System.readNParticles();
	
	xyzFormat<double> xyzFile(p,nParticles);
	xyzFile.open(newName.c_str(), std::ios::in);
	
	for(int i=0;xyzFile.load();i++)
	{
		pairInteractions.build();
		double potential=pairInteractions.computePotential(index);
		
		std::cout << System.readStoreInterval()*(double)i << '\t' << potential << '\n';
	}
	
	return 0;
}

