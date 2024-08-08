//Loading an MD simulation "using systemMD.h".
//Everything needed to do a simulation is already included in systemMD.h.
//Even standard libraries like iostream and fstream are inlcuded.
//But, if any extra libraries are needed you can include them here or in systemMD.h.
//Permenent libraries (included with standard library or with this package) should
// go in systemMD.h, and anything you need that is external should go in this file.

//For molecular dynamics forces and potentials
#include "../include/MD.h"

//For the molecular dynamics variables
#include "../include/system.h"

#include <deque>
#include <list>
#include <chrono>

template <typename T>
int hash(T c, T nFactor)
{
	//threeVector<int> c;
	while(c.x>nFactor.x)c.x-=nFactor.x;
	while(c.y>nFactor.y)c.y-=nFactor.y;
	while(c.z>nFactor.z)c.z-=nFactor.z;
	while(c.x<0)c.x+=nFactor.x;
	while(c.y<0)c.y+=nFactor.y;
	while(c.z<0)c.z+=nFactor.z;
	
	return c.x+c.y*nFactor.x+c.z*nFactor.x*nFactor.y;
}

template <typename T>
T unhash(int v, T nFactor)
{
	T c;
	c.x=v%nFactor.x;
	c.y=static_cast<int>(v/nFactor.x)%nFactor.y;
	c.z=static_cast<int>(v/(nFactor.x*nFactor.y))%nFactor.z;
	return c;
}



int main(int argc, char **argv)
{
	if(argc<5)
	{
		std::cerr << argv[0] << " name surfSpacing nSpacing type1 ...\n";
		return 0;
	}
	
	std::stringstream cmdArg;
	for(int i=2;i<argc;i++)
		cmdArg << argv[i] << '\t';
	
	double surfSpacing;
	int nSpacing;
	cmdArg >> surfSpacing >> nSpacing;
	double surfSpacingSqr=surfSpacing*surfSpacing;
	
	std::vector<int> types;
	
	for(int i=4;i<argc;i++)
	{
		int type;
		cmdArg >> type;
		types.push_back(type);
	}
	
	//the variables for the simulation, call with empty constructor
	Blob<double> System;
	
	//file input-output script engine
	//name part (argv[1] in this case) is used with no extension (*.mpd)
	Script<double, Blob <double> > fileIO(argv[1],std::ios::in,&System);
	
	//read the variables in
	fileIO.read();
	
	//close the file that it is reading from
	fileIO.close();
	
	//set time
	double time=0;

	//grab sizes
	threeVector<double> size=System.readSize();
	
	std::string newName("size_");
	newName+=argv[1];
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
		time=sTime;
		size=sCurrent;
		//std::cout << size.x << '\t' << size.y << '\t' << size.z << '\n';
	}
	
	position<double> *p=System.getPositions();
	int nParticles=System.readNParticles();
	
	//grab particles, by type, for surfaces
	std::vector<twoVector<int> > pIndex;
	
	for(int i=0;i<System.readNParticles();i++)
	{
		bool match=false;
		for(int tIndex=0;!match && tIndex<types.size();tIndex++)
			match=(p[i].type==types[tIndex]);
		
		if(match)
		{
			twoVector<int> buf;
			buf.x=i;
			buf.y=-1;
			pIndex.push_back(buf);
		}
	}
	
	std::cerr << "using " << pIndex.size() << " points for histogram" << std::endl;
	
	//Get our frames file
	std::string framesName("frames_");
	framesName+=argv[1];
	framesName+=".xyz";
	
	xyzFormat<double> xyzFile(p,nParticles);
	xyzFile.open(framesName.c_str(), std::ios::in);
	
	std::vector<twoVector<double> > oldHist;
	
	int frame=0;
	while(xyzFile.load())
	{
		time+=System.readStoreInterval();
		std::cerr << time << std::endl;
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
		//size=System.readSize();
		//set up a quick search using the minimum surface spacing
		threeVector<int> surfN;
		surfN.x=static_cast<int>(size.x/surfSpacing);
		surfN.y=static_cast<int>(size.y/surfSpacing);
		surfN.z=static_cast<int>(size.z/surfSpacing);
		threeVector<double> surfR;
		surfR.x=size.x/static_cast<double>(surfN.x);
		surfR.y=size.y/static_cast<double>(surfN.y);
		surfR.z=size.z/static_cast<double>(surfN.z);
		
		typedef std::vector<position<double> > cellType;//fastest is vector, deque is a little slower
		
		std::map<int,cellType > cells;
		
		for(int j=0;j<pIndex.size();j++)
		{
			int i=pIndex[j].x;
			threeVector<int> c;
			c.x=static_cast<int>(p[i].x/surfR.x)%surfN.x;
			c.y=static_cast<int>(p[i].y/surfR.y)%surfN.y;
			c.z=static_cast<int>(p[i].z/surfR.z)%surfN.z;
			int hash=c.x+c.y*surfN.x+c.z*surfN.x*surfN.y;
			//int hash=cellHash(p[i]);
			std::map<int,cellType >::iterator it=cells.lower_bound(hash);
			if(it!=cells.end() && it->first==hash)//!(cells.key_comp()(it->first, hash)))
			{
				//push operation
				it->second.push_back(p[i]);
			}
			else
			{
				//initialize operation
				cellType buf;
				buf.push_back(p[i]);
				
				cells.insert(it, std::map<int,cellType >::value_type(hash,buf));
			}
		}
		
		std::vector<twoVector<double> > histogram;
		double rSpacing=surfSpacing/static_cast<double>(nSpacing+1);
		//double rSpacingScalar=surfSpacing/M_PI;
		
		for(int i=0;i<nSpacing;i++)
		{
			twoVector<double> hist;
			hist.x=static_cast<double>(i)*rSpacing;
			hist.y=0;
			histogram.push_back(hist);
		}
		
		int nFound=0;
		//int nIndex=0;
		
		//obtaining histogram
		for(auto &cell:cells)
		{
			int currentHash=cell.first;
			cellType &cellPoints=cell.second;
			
			threeVector<int> iCoor=unhash(currentHash,surfN);
			
			for(auto &aP:cellPoints)
			{
				int neighborSum=nFound;
				//if(oldHist.size()>0)
				std::vector<double> dist;
				double minDist=0;//size.x*size.x+size.y*size.y+size.z*size.z;
				//range search
				for(int k=0;k<27;k++)
				{
					threeVector<int> nCoor=iCoor;
					nCoor.x+=k%3-1;
					nCoor.y+=int(k/3)%3-1;
					nCoor.z+=int(k/9)-1;
					int neighborHash=hash(nCoor,surfN);
					auto it=cells.lower_bound(neighborHash);
					if(it!=cells.end() && it->first==neighborHash)
					{
						cellType &neigPoints=it->second;
						for(auto &neP:neigPoints)
						{
							threeVector<double> d;
							d.x=neP.x-aP.x;
							d.y=neP.y-aP.y;
							d.z=neP.z-aP.z;
							
							while(d.x>size.x/2.0)d.x-=size.x;
							while(d.y>size.y/2.0)d.y-=size.y;
							while(d.z>size.z/2.0)d.z-=size.z;
							while(d.x<-size.x/2.0)d.x+=size.x;
							while(d.y<-size.y/2.0)d.y+=size.y;
							while(d.z<-size.z/2.0)d.z+=size.z;
							
							double dr=d.x*d.x+d.y*d.y+d.z*d.z;
							if(dr<surfSpacingSqr && dr!=0)// && dotProduct(normalSum,otherSum)>sqrt(2.0)/2.0)
							{
								dr=sqrt(dr);
								//if(dr<minDist)
									minDist+=dr;
								int histIndex=static_cast<int>(dr/rSpacing);
								histogram[histIndex].y+=1.0;
								nFound++;
								//dist.push_back(dr);
								//std::cout << dr << '\n';
							}
						}
					}
				}
				neighborSum=nFound-neighborSum;
				
				std::cout << neighborSum << ' ' << minDist/static_cast<double>(neighborSum) << '\n';
				//if(neighborSum>100)
				//	std::cout << floor((neighborSum-100)/10) << ' ' << aP.x << ' ' << aP.y << ' ' << aP.z << '\n';
			}
		}
		std::cout << std::endl;
		throw 0;
		
		double Vi=(4.0/3.0)*M_PI;
		double Vt=Vi*surfSpacing*surfSpacing*surfSpacing;
		int index=0;
		for(auto &P:histogram)
		{
			double Px=P.x+rSpacing;
			double VShell=Vi*((Px*Px*Px)-(P.x*P.x*P.x));
			P.x=P.x+rSpacing/2.0;
			P.y=P.y*static_cast<double>(nSpacing)*rSpacing*surfSpacing/(static_cast<double>(nFound)*VShell);
			if(oldHist.size()>0)
			{
				std::cout << P.x << ' ' << P.y << '\n';
			}
		}
		std::cout << std::endl;
		
		oldHist=histogram;
		
		frame++;
	}
	return 0;
}
