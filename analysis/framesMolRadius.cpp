//This is a very short and simple program that shows you how to modify a variable already present in an mpd script.
//The reason this is used is because some variables (forces, potentials, geometry, etc...) are not easily modified
// in the mpd script. Using a file that opens the mpd script, modifies, and then saves it (as below) is far easier
// than modifying the actual variable. Most of it looks like the run.cpp file, it just replaces the "run loop"
// with the modification.

//For molecular dynamics forces and potentials
#include "../include/MD.h"

//For the molecular dynamics variables
#include "../include/system.h"

struct sortByType {
	position<double> * p;
	sortByType(position<double> *particles):p(particles){};
	bool operator() (int a, int b) const noexcept {return p[a].type<p[b].type;}
};

template <typename T>
using umapfwd=typename std::unordered_map<int,std::forward_list<T>>;

int cellHash(threeVector<int> k, threeVector<int> s)
{
	return k.x+k.y*s.x+k.z*s.x*s.y;
}

threeVector<int> unCellHash(int h, threeVector<int> s)
{
	threeVector<int> k;
	k.x=h%s.x;
	k.y=(h/s.x)%s.y;
	k.z=(h/(s.x*s.y))%s.z;
	return k;
}

std::vector<std::vector<int> > getSurfaces(std::vector<std::vector<int> > mols, position<double> *p, threeVector<double> size, double cutoff, int type)
{
	umapfwd<int> cells;
	threeVector<int> s;
	s.x=std::floor(size.x/cutoff)-1;
	s.y=std::floor(size.y/cutoff)-1;
	s.z=std::floor(size.z/cutoff)-1;
	threeVector<double> cellSize;
	cellSize.x=static_cast<double>(size.x/s.x);
	cellSize.y=static_cast<double>(size.y/s.y);
	cellSize.z=static_cast<double>(size.z/s.z);
	
	for(auto mol:mols)
	{
		for(auto i:mol)
		{
			if(p[i].type==type)
			{
				threeVector<int> k;
				k.x=std::floor(p[i].x/cellSize.x);
				k.y=std::floor(p[i].y/cellSize.y);
				k.z=std::floor(p[i].z/cellSize.z);
				int h=hash(k,s);
				auto iter=cells.find(h);
				if(iter!=cells.end())
					iter->second.push_front(i);
				else
					cells[h].push_front(i);
			}
		}
	}
	
	std::vector<int> flags(mols.size(),-1);
	std::vector<std::vector<int> > surfaces;
	for(auto& cell:cells)
	{
		std::deque<int> stack;
		stack.push_back(cell.second.front());
		if(flag[cell.second.front()]==-1)
		{
			flag[cell.second.front()]=surfaces.size();
			std::vector<int> surface;
			while(stack.size()>0)
			{
				auto 
			}
			surfaces.push_back(surface);
		}
	}
	
	return surfaces;
}

int main(int argc, char* argv[])
{
	if(argc!=5)
	{
		//so simple, this can't possibly mess it up
		std::cerr << "usage: " << argv[0] << " name molIndex cutoff headType\n";
		return 0;
	}
	
	///read in options
	char *name=argv[1];
	int molIndex, headType;
	double cutoff;
	
	std::stringstream cmdArg;
	for(int i=2;i<argc;i++)
		cmdArg << argv[i] << ' ';
	cmdArg >> molIndex >> cutoff >> headType;
		
	
	///Configuration variables
	Blob<double> System;
	
	///load old variables, then initialize them
	Script<double, Blob <double> > fileIO(name,std::ios::in,&System);

	fileIO.read();

	fileIO.close();
	
	//our current size
	threeVector<double> size=System.readSize();
	
	//check if size varies
	std::string newName("size_");
	newName+=name;
	newName+=".dat";
	
	std::fstream sizeFile;
	sizeFile.open(newName.c_str(), std::ios::in);
	if(sizeFile.is_open())
	{
		double sTime=-1;
		threeVector<double> sCurrent;
		//for the first iteration, I expect this to exit the loop after one read,
		// otherwise, this is a template for subsequent reads in the loop that is
		// coming up.
		while(sTime<0 && !sizeFile.eof())
			sizeFile >> sTime >> sCurrent.x >> sCurrent.y >> sCurrent.z;
		size=sCurrent;
	}
	
	position<double> *p=System.getPositions();
	int nParticles=System.readNParticles();
	
	//Get the molecule's particle indices
	if(molIndex<0 || molIndex>=System.readNMolecules())
	{
		std::cerr << "Molecule " << molIndex << " is out of bounds!" << std::endl;
		std::cerr << "nMolecules=" << System.readNMolecules() << std::endl;
		return -1;
	}
	std::vector<std::vector<int> > molI;
	auto mol=System.getMolecules()[molIndex];
	switch(mol.readType())
	{
		case CHAIN:
			for(int j=0;j<m.readNBond();j++)
			{
				int start=mol.getBonds()[j].s[START];
				int nChains=mol.getBonds()[j].s[NCHAINS];
				int length=mol.getBonds()[j].s[CHAINLENGTH];
				for(int chain=0;chain<nChains;chain++)
				{
					std::vector<int> mI;
					for(int k=chain*chainLength+start;k<(chain+1)*chainLength+start;k++)
						mI.push_back(k);
					std::sort(mI.begin(),mI.end(),sortByType(p));
					molI.push_back(mI);
				}
			}
		break;
		//need to add more
		default:
			std::cerr << "molType not recognized!" << std::endl;
			return -1;
		break;
	}
	
	//Get our frames file
	std::string framesName("frames_");
	framesName+=name;
	framesName+=".xyz";
	
	xyzFormat<double> xyzFile(p,nParticles);
	xyzFile.open(framesName.c_str(), std::ios::in);
	
	double time=0;
	while(xyzFile.load())
	{
		//make sure the size hasn't changed
		if(sizeFile.is_open())
		{
			double sTime=0;
			threeVector<double> sCurrent=size;//just in case it is current
			//for the first iteration, I expect this to exit the loop after one read,
			// otherwise, this is a template for subsequent reads in the loop that is
			// coming up.
			while(sTime<time && !sizeFile.eof())
				sizeFile >> sTime >> sCurrent.x >> sCurrent.y >> sCurrent.z;
			size=sCurrent;
		}
		
		std::vector<std::vector<int> > surfaces;
		
		
		///Locate center of mass for centerTypes
		threeVector<double> com=0;
		int nPCom=0;
		for(int i:pI)
		{
			threeVector<double> d;
			d.x=p[i].x-p[0].x;
			d.y=p[i].y-p[0].y;
			d.z=p[i].z-p[0].z;
			
			if(d.x>size.x/2.0) p[i].x-=size.x;
			if(d.y>size.y/2.0) p[i].y-=size.y;
			if(d.z>size.z/2.0) p[i].z-=size.z;
			if(d.x<-size.x/2.0) p[i].x+=size.x;
			if(d.y<-size.y/2.0) p[i].y+=size.y;
			if(d.z<-size.z/2.0) p[i].z+=size.z;
			
			com.x+=p[i].x;
			com.y+=p[i].y;
			com.z+=p[i].z;
			nPCom++;
		}
		if(nPCom==0)
		{
			std::cerr << "Types selected are not present!\n";
			return -1;
		}
		
		com.x/=(double)nPCom;
		com.y/=(double)nPCom;
		com.z/=(double)nPCom;
		
		while(com.x>size.x) com.x-=size.x;
		while(com.y>size.y) com.y-=size.y;
		while(com.z>size.z) com.z-=size.z;
		while(com.x<0) com.x+=size.x;
		while(com.y<0) com.y+=size.y;
		while(com.z<0) com.z+=size.z;
		
		///Locate center of mass for centerTypes
		double radius=0;
		for(int i:pI)
		{
			threeVector<double> d;
			d.x=p[i].x-com.x;
			d.y=p[i].y-com.y;
			d.z=p[i].z-com.z;
			
			if(d.x>size.x/2.0) d.x-=size.x;
			if(d.y>size.y/2.0) d.y-=size.y;
			if(d.z>size.z/2.0) d.z-=size.z;
			if(d.x<-size.x/2.0) d.x+=size.x;
			if(d.y<-size.y/2.0) d.y+=size.y;
			if(d.z<-size.z/2.0) d.z+=size.z;
			
			radius+=sqrt(d.x*d.x+d.y*d.y+d.z*d.z);
		}
		radius/=(double)nPCom;
		std::cout << time << ' ' << radius << std::endl;
		std::cerr << time << std::endl;
		time+=System.readStoreInterval();
	}
	
	
	return 0;
}

