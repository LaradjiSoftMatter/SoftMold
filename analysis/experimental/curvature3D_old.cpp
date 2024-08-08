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

#define sThreshold 20
#define pSig 0.3

struct edgeType {
	int x,y;
};

//for sorting keys, logarithmic search
bool compByX(edgeType a, edgeType b)
{
	//if(a.x<b.x)
	//	return true;
	//return false;
	return a.x<b.x;
}

//for sorting keys, logarithmic search
bool equByXY(edgeType a, edgeType b)
{
	return (a.x==b.x) && (a.y==b.y);
}


//for simple linear search
class equByY
{
	public:
		bool operator () (edgeType a)
		{
			return a.y==y;
		};
		int x;
		int y;
};

//for simple linear search
class equByX
{
	public:
		bool operator () (edgeType a)
		{
			return a.x==x;
		};
		int x;
		int y;
};

//bool equByY(edgeType a, edgeType b)
//{
//	if(a.y==b.y)
//		true;
//	false;
//}

//to transform into unique range
//bool equByX(edgeType a, edgeType b)
//{
//	if(a.x==b.x)
//		true;
//	false;
//}

//I should double check this
template <typename T>
fourVector<double> circumscribe(T a, T b, T c)
{
	//std::cerr << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~``" << std::endl;
	//std::cerr << a.x << '\t' << a.y << '\t' << a.z << std::endl;
	//std::cerr << b.x << '\t' << b.y << '\t' << b.z << std::endl;
	
	T crossProd=crossProduct(a,b);
	
	//std::cerr << crossProd.x << '\t' << crossProd.y << '\t' << crossProd.z << std::endl;
	
	
	double twoCrossSqr=2.0*(dotProduct(crossProd,crossProd));
	//double twoCrossSqr=2.0*(magnitude(crossProd));
	
	//std::cerr << "twoCrossSqr: " << twoCrossSqr << std::endl;
	
	fourVector<double> rFact;
	
	//radius
	rFact.t=(dotProduct(a,a)*dotProduct(b,b)*dotProduct(c,c))/twoCrossSqr;
	//rFact.t=magnitude(a)*magnitude(b)*magnitude(c)/twoCrossSqr;
	
	//std::cerr << "twoCrossSqr: " << twoCrossSqr << std::endl;
	
	T aNeg=a;
	T bNeg=b;
	T cNeg=c;
	
	aNeg.x=-aNeg.x;
	aNeg.y=-aNeg.y;
	aNeg.z=-aNeg.z;
	
	bNeg.x=-bNeg.x;
	bNeg.y=-bNeg.y;
	bNeg.z=-bNeg.z;
	
	cNeg.x=-cNeg.x;
	cNeg.y=-cNeg.y;
	cNeg.z=-cNeg.z;
	
	//the offset scalers (threeVector r=rFact.x*a+rFact.y*b+rFact.z*c;)
	//rFact.x=dotProduct(b,b)*(dotProduct(a,cNeg))/twoCrossSqr;
	//rFact.y=dotProduct(cNeg,cNeg)*(dotProduct(aNeg,b))/twoCrossSqr;
	//rFact.z=dotProduct(a,a)*(dotProduct(c,bNeg))/twoCrossSqr;
	
	rFact.x=dotProduct(b,b)*(dotProduct(a,cNeg))/twoCrossSqr;
	rFact.y=dotProduct(cNeg,cNeg)*(dotProduct(aNeg,b))/twoCrossSqr;
	rFact.z=dotProduct(a,a)*(dotProduct(c,bNeg))/twoCrossSqr;
	return rFact;
	//std::cerr << "=====================================" << std::endl;
	//std::cin.get();
}



int main(int argc, char **argv)
{
	if(argc<5)
	{
		std::cerr << argv[0] << " name surfSpacing dry type1 type2 ...\n";
		std::cerr << "\tdry is for no output, if the surface spacing needs to be optimized. 0=false 1=true\n";
		return 0;
	}
	
	std::stringstream cmdArg;
	for(int i=2;i<argc;i++)
		cmdArg << argv[i] << '\t';
	
	double surfSpacing;
	int dry;
	cmdArg >> surfSpacing >> dry;
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
	std::vector<int> pIndex;
	
	
	
	for(int i=0;i<System.readNParticles();i++)
	{
		bool match=false;
		for(int tIndex=0;!match && tIndex<types.size();tIndex++)
			match=(p[i].type==types[tIndex]);
		
		if(match)
			pIndex.push_back(i);
	}
	
	//Get our frames file
	std::string framesName("frames_");
	framesName+=argv[1];
	framesName+=".xyz";
	
	//xyzFormat<double> xyzFile(p,nParticles);
	//xyzFile.open(framesName.c_str(), std::ios::in);
	
	std::vector<int> sFlagOrigin;
	int nSurfaces=0;
	
	//while(xyzFile.load())
	{
		//subsample pIndex to pFilt on a grid, weighted by center of mass
		//std::map<int, position<double> >::iterator->first is an integer for the hashed grid position
		//std::map<int, position<double> >::iterator->second is the weighted central position
		//position pos=pFilt[gridId]; would get you the actual position
		std::map<int, position<double> > pFilt;
		
		
		std::vector< std::vector< position<double> > > cellList;
		time+=System.readStoreInterval();
		std::cerr << time << '\t';
		//make sure the size hasn't changed
		if(sizeFile.is_open())
		{
			double sTime=0;
			threeVector<double> sCurrent=size;//just in case it is current
			//for the first iteration, I expect this to exit the loop after one read,
			// otherwise, this is a template for subsequent reads in the loop that is
			// coming up.
			while(sTime<time && !sizeFile.eof())
			{
				sizeFile >> sTime >> sCurrent.x >> sCurrent.y >> sCurrent.z;
				//std::cout << time << '\t' << sTime << '\t' << size.x << '\t' << size.y << '\t' << size.z << std::endl;
				
			}
			size=sCurrent;
		}
		size=System.readSize();
		//set up a quick search using the minimum surface spacing
		threeVector<int> surfN;
		surfN.x=static_cast<int>(size.x/surfSpacing);
		surfN.y=static_cast<int>(size.y/surfSpacing);
		surfN.z=static_cast<int>(size.z/surfSpacing);
		threeVector<double> surfR;
		surfR.x=size.x/static_cast<double>(surfN.x);
		surfR.y=size.y/static_cast<double>(surfN.y);
		surfR.z=size.z/static_cast<double>(surfN.z);
		
		typedef std::vector<int> cellType;//fastest is vector, deque is a little slower
		
		std::map<int,cellType > cells;
		
		for(int j=0;j<pIndex.size();j++)
		{
			int i=pIndex[j];
			threeVector<int> c;
			c.x=static_cast<int>(p[i].x/surfR.x)%surfN.x;
			c.y=static_cast<int>(p[i].y/surfR.y)%surfN.y;
			c.z=static_cast<int>(p[i].z/surfR.z)%surfN.z;
			int hash=c.x+c.y*surfN.x+c.z*surfN.x*surfN.y;
			//int hash=cellHash(p[i]);
			std::map<int,cellType >::iterator it=cells.lower_bound(hash);
			if(it!=cells.end() && it->first==hash)//!(cells.key_comp()(it->first, hash)))
			{
				//std::cerr << it->first << '\t' << hash
				
				//push operation
				it->second.push_back(j);
			}
			else
			{
				//initialize operation
				cellType buf;
				buf.push_back(j);
				//cells[hash]=buf;
				cells.insert(it, std::map<int,cellType >::value_type(hash,buf));
			}
		}
		
		std::cerr << "using " << pIndex.size() << " points for surfaces" << std::endl;
		
		//std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
		
		//get surfaces
		std::vector< std::vector<int> > surfaces;//this just grows with indices per surface
		std::vector<int> sFlag(nParticles,-1);//cooresponds to the surfaces found i.e. 0=surface 1, 1=surface 2...
		for(int i=0;i<pIndex.size();++i)
		{
			//unaccounted surface
			if(sFlag[pIndex[i]]==-1)
			{
				std::vector< edgeType > edges;
				std::vector<int> surface;//new surface
				//std::vector<int> stack;
				std::deque<int> stack;//almost the same as vector
				//std::list<int> stack;
				stack.push_back(i);
				sFlag[pIndex[i]]=surfaces.size();
				surface.push_back(pIndex[i]);
				while(stack.size()>0)
				{
					int element=stack.back();
					stack.pop_back();
					
					int j=pIndex[element];
					
					//calculate hash
					threeVector<int> iCoor;
					threeVector<double> corP;
					for(corP.x=p[j].x;corP.x>size.x;corP.x-=size.x);
					for(corP.y=p[j].y;corP.y>size.y;corP.y-=size.y);
					for(corP.z=p[j].z;corP.z>size.z;corP.z-=size.z);
					for(corP.x=corP.x;corP.x<0;corP.x+=size.x);
					for(corP.y=corP.y;corP.y<0;corP.y+=size.y);
					for(corP.z=corP.z;corP.z<0;corP.z+=size.z);
					iCoor.x=static_cast<int>(corP.x/surfR.x)%surfN.x;
					iCoor.y=static_cast<int>(corP.y/surfR.y)%surfN.y;
					iCoor.z=static_cast<int>(corP.z/surfR.z)%surfN.z;
					
					//find connected elements and push onto stack
					//check the nearby cells for traversal
					for(int k=0;k<27;k++)
					{
						int a=k%3-1;
						int b=int(k/3)%3-1;
						int c=int(k/9)-1;
						int neighborHash=
						((iCoor.x+a+surfN.x)%surfN.x)+
						(((iCoor.y+b+surfN.y)%surfN.y)*surfN.x)+
						(((iCoor.z+c+surfN.z)%surfN.z)*surfN.x*surfN.y);
						//cellHash(iCoor.x+a,iCoor.y+b,iCoor.z+c);
						//std::cerr << hash << '\t' << neighborHash << std::endl;
						
						std::map<int,cellType >::iterator it=cells.find(neighborHash);
						for(int m=0;it!=cells.end() && m<it->second.size();m++)
							//for(int m=0;m<pIndex.size();m++)
						{
							threeVector<double> d;
							
							int l=pIndex[it->second[m]];
							//int l=pIndex[m];
							d.x=p[l].x-p[j].x;
							d.y=p[l].y-p[j].y;
							d.z=p[l].z-p[j].z;
							
							while(d.x>size.x/2.0)d.x-=size.x;
							while(d.y>size.y/2.0)d.y-=size.y;
							while(d.z>size.z/2.0)d.z-=size.z;
							while(d.x<-size.x/2.0)d.x+=size.x;
							while(d.y<-size.y/2.0)d.y+=size.y;
							while(d.z<-size.z/2.0)d.z+=size.z;
							
							if(l!=j && d.x*d.x+d.y*d.y+d.z*d.z<surfSpacingSqr)
							{
								if(sFlag[pIndex[it->second[m]]]==-1)
								{
									sFlag[pIndex[it->second[m]]]=surfaces.size();
									surface.push_back(l);
									
									stack.push_back(it->second[m]);
									
									do {
										d.z=p[l].z-p[j].z;
										if(d.z>size.z/2.0)
											p[l].z-=size.z;
										d.z=p[l].z-p[j].z;
									} while(d.z>size.z/2.0);
									
									do {
										d.z=p[l].z-p[j].z;
										if(d.z<-size.z/2.0)
											p[l].z+=size.z;
										d.z=p[l].z-p[j].z;
									} while(d.z<-size.z/2.0);
								}
							}
						}
					}
				}
				surfaces.push_back(surface);
			}
		}
		
		//get a list of nearest neighbor edges
		std::vector< std::vector< edgeType > > sEdges;
		std::vector<std::vector<threeVector<int> > > sTriangles;
		
		for(int s=0;s<surfaces.size();++s)
		{
			std::vector< edgeType > edges;
			std::vector<threeVector<int> > triangles;
			for(int i=0;i<surfaces[s].size() && surfaces[s].size()>10;i++)
			{
				if(i%(surfaces[s].size()/100)==0 && i!=0)
				{
					std::cerr << (i*100/surfaces[s].size())+1 << "% of surface " << s << ".";
					std::cerr << triangles.size() << " triangles so far!" << std::endl;
				}
				int j=surfaces[s][i];
				
				std::vector<int> candidates;
				//std::vector<int> neighbors;
				//candidates.push_back(j);
				
				//calculate hash
				threeVector<int> iCoor;
				threeVector<double> corP;
				for(corP.x=p[j].x;corP.x>size.x;corP.x-=size.x);
				for(corP.y=p[j].y;corP.y>size.y;corP.y-=size.y);
				for(corP.z=p[j].z;corP.z>size.z;corP.z-=size.z);
				for(corP.x=corP.x;corP.x<0;corP.x+=size.x);
				for(corP.y=corP.y;corP.y<0;corP.y+=size.y);
				for(corP.z=corP.z;corP.z<0;corP.z+=size.z);
				iCoor.x=static_cast<int>(corP.x/surfR.x)%surfN.x;
				iCoor.y=static_cast<int>(corP.y/surfR.y)%surfN.y;
				iCoor.z=static_cast<int>(corP.z/surfR.z)%surfN.z;
				
				//check the nearby cells for traversal
				for(int k=0;k<27;k++)
				{
					int a=k%3-1;
					int b=int(k/3)%3-1;
					int c=int(k/9)-1;
					int neighborHash=
					((iCoor.x+a+surfN.x)%surfN.x)+
					(((iCoor.y+b+surfN.y)%surfN.y)*surfN.x)+
					(((iCoor.z+c+surfN.z)%surfN.z)*surfN.x*surfN.y);
					//cellHash(iCoor.x+a,iCoor.y+b,iCoor.z+c);
					//std::cerr << hash << '\t' << neighborHash << std::endl;
					
					std::map<int,cellType >::iterator it=cells.find(neighborHash);
					for(int m=0;it!=cells.end() && m<it->second.size();m++)
						//for(int m=0;m<pIndex.size();m++)
					{
						threeVector<double> d;
						
						int l=pIndex[it->second[m]];
						//int l=pIndex[m];
						d.x=p[l].x-p[j].x;
						d.y=p[l].y-p[j].y;
						d.z=p[l].z-p[j].z;
						
						//while(d.x>size.x/2.0)d.x-=size.x;
						//while(d.y>size.y/2.0)d.y-=size.y;
						//while(d.z>size.z/2.0)d.z-=size.z;
						//while(d.x<-size.x/2.0)d.x+=size.x;
						//while(d.y<-size.y/2.0)d.y+=size.y;
						//while(d.z<-size.z/2.0)d.z+=size.z;
						
						double dr=d.x*d.x+d.y*d.y+d.z*d.z;
						
						if(l>j && sFlag[pIndex[it->second[m]]]==sFlag[j] && dr<3.0*surfSpacing*surfSpacing)
						{
							edgeType edge;
							edge.x=j;
							edge.y=l;
							//if(a>l)//no duplicates
							edges.push_back(edge);
							
							//if(a!=-2 && b!=-2 && c!=-2 && a!=2 && b!=2 && c!=2)
								candidates.push_back(l);
							//neighbors.push_back(l);
						}
					}
				}
				
				//for(int j=0;j<candidates.size();j++)
				{
					for(int k=0;k<candidates.size();k++)
					{
						for(int l=k+1;l<candidates.size();l++)
						{
							int accept=0;
							int a=j;
							int b=candidates[k];
							int c=candidates[l];
							
							threeVector<double> dab,dbc,dca;
							//ab
							dab.x=p[a].x-p[b].x;
							dab.y=p[a].y-p[b].y;
							dab.z=p[a].z-p[b].z;
							//while(dab.x>size.x/2.0)dab.x-=size.x;
							//while(dab.y>size.y/2.0)dab.y-=size.y;
							//while(dab.z>size.z/2.0)dab.z-=size.z;
							//while(dab.x<-size.x/2.0)dab.x+=size.x;
							//while(dab.y<-size.y/2.0)dab.y+=size.y;
							//while(dab.z<-size.z/2.0)dab.z+=size.z;
							//bc
							dbc.x=p[b].x-p[c].x;
							dbc.y=p[b].y-p[c].y;
							dbc.z=p[b].z-p[c].z;
							//while(dbc.x>size.x/2.0)dbc.x-=size.x;
							//while(dbc.y>size.y/2.0)dbc.y-=size.y;
							//while(dbc.z>size.z/2.0)dbc.z-=size.z;
							//while(dbc.x<-size.x/2.0)dbc.x+=size.x;
							//while(dbc.y<-size.y/2.0)dbc.y+=size.y;
							//while(dbc.z<-size.z/2.0)dbc.z+=size.z;
							//ca
							dca.x=p[c].x-p[a].x;
							dca.y=p[c].y-p[a].y;
							dca.z=p[c].z-p[a].z;
							//while(dca.x>size.x/2.0)dca.x-=size.x;
							//while(dca.y>size.y/2.0)dca.y-=size.y;
							//while(dca.z>size.z/2.0)dca.z-=size.z;
							//while(dca.x<-size.x/2.0)dca.x+=size.x;
							//while(dca.y<-size.y/2.0)dca.y+=size.y;
							//while(dca.z<-size.z/2.0)dca.z+=size.z;
							
							fourVector<double> scal=circumscribe(dab,dbc,dca);
							
							position<double> sphere;
							
							threeVector<double> d;
									
							//int l=pIndex[m];
							d.x=sphere.x-p[l].x;
							d.y=sphere.y-p[l].y;
							d.z=sphere.z-p[l].z;
							
							scal.t=sqrt(d.x*d.x+d.y*d.y+d.z*d.z)-0.00001;
							
							//if(scal.t<3.0 && scal.t>1.0)//sanity check?
							{
								//a sphere bounding our 3 points, second deluancy criteria
								sphere.x=scal.x*p[a].x+scal.y*p[b].x+scal.z*p[c].x;
								sphere.y=scal.x*p[a].y+scal.y*p[b].y+scal.z*p[c].y;
								sphere.z=scal.x*p[a].z+scal.y*p[b].z+scal.z*p[c].z;
								
								
								
								//double scalTSqr=scal.t;
								//scal.t*=scal.t;
								
								//while(sphere.x>size.x)sphere.x-=size.x;
								//while(sphere.y>size.y)sphere.y-=size.y;
								//while(sphere.z>size.z)sphere.z-=size.z;
								//while(sphere.x<0)sphere.x+=size.x;
								//while(sphere.y<0)sphere.y+=size.y;
								//while(sphere.z<0)sphere.z+=size.z;
								
								for(int m=0;m<candidates.size();m++)
								{
									//threeVector<double> d;
									
									//int l=pIndex[m];
									d.x=sphere.x-p[candidates[m]].x;
									d.y=sphere.y-p[candidates[m]].y;
									d.z=sphere.z-p[candidates[m]].z;
									
									while(d.x>size.x/2.0)d.x-=size.x;
									while(d.y>size.y/2.0)d.y-=size.y;
									while(d.z>size.z/2.0)d.z-=size.z;
									while(d.x<-size.x/2.0)d.x+=size.x;
									while(d.y<-size.y/2.0)d.y+=size.y;
									while(d.z<-size.z/2.0)d.z+=size.z;
									
									double dr=d.x*d.x+d.y*d.y+d.z*d.z;
									
									if(candidates[m]!=b && candidates[m]!=c && sqrt(dr)<scal.t && sFlag[pIndex[m]]==s && candidates[m]!=a)
									{
										accept++;
										//m=pIndex.size();
									}
								}
							}
							//else
							//{
							//	accept++;
							//}
							if(accept<2)
							{
								threeVector<int> deluancy;
								/*threeVector<double> d;
								
								std::cout << sphere.x << ' ' << sphere.y << ' ' << sphere.z << '\n';
								std::cout << p[a].x << ' ' << p[a].y << ' ' << p[a].z << '\n';
								std::cout << p[b].x << ' ' << p[b].y << ' ' << p[b].z << std::endl;
								std::cout << p[c].x << ' ' << p[c].y << ' ' << p[c].z << std::endl;
								
								d.x=sphere.x-p[a].x;
								d.y=sphere.y-p[a].y;
								d.z=sphere.z-p[a].z;
								std::cout << sqrt(d.x*d.x+d.y*d.y+d.z*d.z) << std::endl;
								
								d.x=sphere.x-p[b].x;
								d.y=sphere.y-p[b].y;
								d.z=sphere.z-p[b].z;
								std::cout << sqrt(d.x*d.x+d.y*d.y+d.z*d.z) << std::endl;
								d.x=sphere.x-p[c].x;
								d.y=sphere.y-p[c].y;
								d.z=sphere.z-p[c].z;
								std::cout << sqrt(d.x*d.x+d.y*d.y+d.z*d.z) << std::endl;
								
								std::cout << scal.t << std::endl;
								
								std::cin.get();*/
								
								deluancy.x=a;
								deluancy.y=b;
								deluancy.z=c;
								triangles.push_back(deluancy);
							}
							//else
							//	std::cout << "2 " << sphere.x << ' ' << sphere.y << ' ' << sphere.z << '\n';
						}
					}
				}
				
			}
			std::cerr << "Found " << triangles.size() << " triangles!" << std::endl;
			sEdges.push_back(edges);
			sTriangles.push_back(triangles);
		}
		
		for(int i=0;i<sEdges.size();i++)
		{
			std::cerr << "Surface " << i << " has " << sEdges[i].size() << " edges." << std::endl;
		}
		
		//our first surface is assumed to be correct
		if(sFlagOrigin.size()==0)
		{
			sFlagOrigin=sFlag;
			nSurfaces=surfaces.size();
			//throw 0;
		}
		else//make new surfaces match
		{
			std::cout << std::endl;
			std::vector< std::vector<int> > cMap(nSurfaces, std::vector<int>(surfaces.size(),0) );
			for(int i=0;i<sFlag.size();i++)
			{
				cMap[sFlagOrigin[i]][sFlag[i]]++;
			}
			std::vector<int> sMap(nSurfaces, 0);
			for(int i=0;i<cMap.size();i++)
				for(int j=0;j<cMap[i].size();j++)
					sMap[i]+=cMap[i][j];
			
			std::cout << "nP\t\t";
			for(int i=0;i<surfaces.size();i++)
				std::cout << surfaces[i].size() << '\t';
			std::cout << "\n\to/n\t";
			for(int i=0;i<surfaces.size();i++)
				std::cout << i << '\t';
			std::cout << std::endl;
			for(int i=0;i<cMap.size();i++)
			{
				std::cout << sMap[i] << '\t' << i << "\t[";
				for(int j=0;j<cMap[i].size();j++)
				{
					std::cout << std::fixed << std::setprecision(4) << static_cast<double>(cMap[i][j])/static_cast<double>(sMap[i]) << '\t';
				}
				std::cout << "]" << std::endl;
			}
			//no match, swap with other
			if(cMap[0][0]<pSig*pIndex.size())
			{
				std::cerr << "\nSwapping: " << cMap[0][0] << '\t' << pSig*pIndex.size() << std::endl;
				for(int i=0;i<cMap[0].size();i++)
				{
					if(cMap[0][i]>pSig*pIndex.size())
					{
						surfaces[i].swap(surfaces[0]);
						sEdges[i].swap(sEdges[0]);
					}
				}
			}
		}
		
		//std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
		//std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
		//std::cerr << time_span.count() << std::endl;
		
		std::cerr << "number of surfaces: " << surfaces.size() << std::endl;
		
		//each coordinate is an edge in sEdges
		//std::vector<std::vector<threeVector<int> > > sTriangles;
		
		//each coordinate is a triangle in sTriangles
		//each index cooresponds to an edge in edges vector (auto &edges=sEdges[s];)
		std::vector<std::vector<std::vector<int> > > sENeighbors;

		//for(int s=0;s<sEdges.size();s++)
		for(int s=0;false;s++)
		{
			//auto &edges=sEdges[s];
			std::vector< edgeType > &edges=sEdges[s];
			
			std::cerr << "Surface " << s << " has " << edges.size() << " edges." << std::endl;
			
			typedef std::vector<edgeType >::iterator tVIter;
			
			if(edges.size()>0)
			{
				std::sort(edges.begin(), edges.end(), compByX);
				//for(tVIter edge=edges.begin();edge<edges.end();edge++)
				//	std::cout << (*edge).x << '\t' << (*edge).y << std::endl;
				//throw 0;
			}
			
			//triangles[eNeighbors[e][i]] sharing the edge.
			//index e is an index that corresponds to the edge list index (edges[e]).
			//Since this is a surface, we assume edges are only shared by 2 triangles.
			//Extremely rare cases might find more than 2 triangles, but we can leave those as orphaned vertices.
			
			std::vector<std::vector<int> > eNeighbors;
			
			if(edges.size()>0)
				eNeighbors.assign(edges.size(), std::vector<int>());
			
			std::vector<threeVector<int> > triangles;
			//enum edgeState {FREE, CROSSING, GROWTH, ENCAPSULATED, DELETED};
			
			//Edges are stored int threeVector like this:
			//x: {a,a,a,b,b,b,b,...}
			//y: {c,d,f,c,f,g,z,...}
			
			//The following converts "nearest neighbor" subgraph into deluancy triangulation graph
			//Our nearest neighbor subgraph is only valid if our bounding boxes are able to contain 4 points
			//i.e.
			//valid: -----    invalid: -----
			//      |   ..|           | .   |
			//      | . . |           |.   .|
			//       -----             -----
			//Our hash subgraph might contain only one point, but our bounding box is composed of 
			//at least 27 hash boxes. If the hash box isn't connected to another hash box, the graph is invalid.
			//With 3.11 lipids/r_min^2 and surf spacing ~ 1.75, each hash box contains about 9 lipids.
			//There are around 270 edge possibilities. This could be reduced with a smaller surf spacing, but 
			//it might not be worth it. The neighbers are cutoff a little more with a bounding sphere, 
			//surf spacing radius, leaving about 120 edges. We should reject next neighbors with distances 
			//less than the surf spacing.
			//std::vector<bool> usedEdge(edges.size(),false);
			
			//candidate triangles ({l,j,k,m}={a,b,c,d}):
			//   m
			// l<|>k
			//   j
			// mlj and jkm are the candidates
			//flipped crossing:
			//   k
			// m<|>j
			//   l
			// kml and ljk are the candidates
			//These vectors are the candidates, abcd is the normal one, bcda is the flipped one
			//They are labeled according to the diagram above
			//std::vector<fourVector<int> > ljkm;//, mljk;
			
			for(int edge=0;edge<edges.size();edge++)
			{
				//the edge is used
				//if(usedEdge[edge]==FREE || usedEdge[edge]==GROWTH)
				{
					int l=edges[edge].x;
					
					//narrow our search range
					edgeType key;
					//equByY key;
					key.x=l;
					tVIter lLow=std::lower_bound(edges.begin(),edges.end(),key, compByX);
					tVIter lHigh=std::upper_bound(lLow,edges.end(),key, compByX);
					
//std::cerr << "lHigh-lLow=" << lHigh-lLow << " lLow-edges.begin()=" << lLow-edges.begin() << std::endl;
//std::cin.get();
					
					for(int ljEdge=edge+1;ljEdge<edges.size() && edges[ljEdge].x==edges[edge].x;ljEdge++)
					{
						int j=edges[ljEdge].y;
						
						//check for a triangle in this set:
						//   m?
						// l<|   
						//   j
						for(int lmEdge=ljEdge+1;lmEdge<edges.size() && edges[lmEdge].x==edges[edge].x;lmEdge++)
						{
							int m=edges[lmEdge].y;
							
							//check for encapsulated points in this set:
							//   m
							// l<|>k?
							//   j
							
							//narrow our search range
							
							key.x=m;
							tVIter mLow=std::lower_bound(edges.begin(),edges.end(),key, compByX);
							tVIter mHigh=std::upper_bound(mLow,edges.end(),key, compByX);
//std::cerr << "mHigh-mLow=" << mHigh-mLow << " mLow-edges.begin()=" << mLow-edges.begin() << std::endl;
//std::cin.get();
							//get the potential crossing between l and k along j and m
							int jmEdge=-1;
							key.x=m;
							key.y=j;
							//tVIter jIter=std::find(mLow,mHigh,key);
							tVIter jIter=mHigh;
							for(tVIter mEdge=mLow;mEdge<mHigh;mEdge++)
							{
//std::cerr << "mEdge-mLow=" << mEdge-mLow << std::endl;
//std::cin.get();
								if((*mEdge).y==key.y)
								{
									jIter=mEdge;
								}
							}
							
							if(jIter!=mHigh)
								jmEdge=jIter-edges.begin();
							
							std::vector<int> mkEdge, jkEdge;
							
							key.x=l;
							//do j and m share another common vertex (k)? our first deluancy criteria
							for(tVIter mEdge=mLow;jmEdge!=-1 && mEdge<mHigh;mEdge++)
							{
								key.y=(*mEdge).y;
								//tVIter kIter=std::find(lLow, lHigh, key);
								tVIter kIter=lHigh;
								for(tVIter mEdge=lLow;mEdge<lHigh;mEdge++)
									if((*mEdge).y==key.y)
										kIter=mEdge;
								
								if(kIter!=lHigh && (*kIter).y!=(*jIter).y)
								{
									mkEdge.push_back(mEdge-edges.begin());
									jkEdge.push_back(kIter-edges.begin());
								}
								//this gives both mkEdge and jkEdge
							}
							
							//check our first nearest neighbor criteria
							threeVector<double> da,db,dc;
							//ab
							da.x=p[l].x-p[m].x;
							da.y=p[l].y-p[m].y;
							da.z=p[l].z-p[m].z;
							while(da.x>size.x/2.0)da.x-=size.x;
							while(da.y>size.y/2.0)da.y-=size.y;
							while(da.z>size.z/2.0)da.z-=size.z;
							while(da.x<-size.x/2.0)da.x+=size.x;
							while(da.y<-size.y/2.0)da.y+=size.y;
							while(da.z<-size.z/2.0)da.z+=size.z;
							//bc
							db.x=p[m].x-p[j].x;
							db.y=p[m].y-p[j].y;
							db.z=p[m].z-p[j].z;
							while(db.x>size.x/2.0)db.x-=size.x;
							while(db.y>size.y/2.0)db.y-=size.y;
							while(db.z>size.z/2.0)db.z-=size.z;
							while(db.x<-size.x/2.0)db.x+=size.x;
							while(db.y<-size.y/2.0)db.y+=size.y;
							while(db.z<-size.z/2.0)db.z+=size.z;
							//ca
							dc.x=p[j].x-p[l].x;
							dc.y=p[j].y-p[l].y;
							dc.z=p[j].z-p[l].z;
							while(dc.x>size.x/2.0)dc.x-=size.x;
							while(dc.y>size.y/2.0)dc.y-=size.y;
							while(dc.z>size.z/2.0)dc.z-=size.z;
							while(dc.x<-size.x/2.0)dc.x+=size.x;
							while(dc.y<-size.y/2.0)dc.y+=size.y;
							while(dc.z<-size.z/2.0)dc.z+=size.z;
							
							fourVector<double> scal=circumscribe(da,db,dc);
							
							position<double> sphere;
							
							if(mkEdge.size()>0 && scal.t<surfSpacing)
							{
								//a sphere bounding our 3 points, second deluancy criteria
								sphere.x=scal.x*p[l].x+scal.y*p[j].x+scal.z*p[m].x;
								sphere.y=scal.x*p[l].y+scal.y*p[j].y+scal.z*p[m].y;
								sphere.z=scal.x*p[l].z+scal.y*p[j].z+scal.z*p[m].z;
							}
							else//try the flipped version on another iteration, purge candidates
							{
								mkEdge.clear();
							}
							
							//does it overlap with any other points?
							if(mkEdge.size()>0)
							{
								double radSqr=scal.t*scal.t;
								//calculate hash
								threeVector<int> iCoor;
								threeVector<double> corP;
								for(corP.x=sphere.x;corP.x>size.x;corP.x-=size.x);
								for(corP.y=sphere.y;corP.y>size.y;corP.y-=size.y);
								for(corP.z=sphere.z;corP.z>size.z;corP.z-=size.z);
								for(corP.x=corP.x;corP.x<0;corP.x+=size.x);
								for(corP.y=corP.y;corP.y<0;corP.y+=size.y);
								for(corP.z=corP.z;corP.z<0;corP.z+=size.z);
								iCoor.x=static_cast<int>(corP.x/surfR.x)%surfN.x;
								iCoor.y=static_cast<int>(corP.y/surfR.y)%surfN.y;
								iCoor.z=static_cast<int>(corP.z/surfR.z)%surfN.z;
								//int hash=iCoor.x+iCoor.y*surfN.x;
								//threeVector<int> iCoor=cellHash.cellCoord(p[j]);
								//find connected elements and push onto stack
								//check the nearby cells for traversal
								for(int nB=0;nB<27;nB++)
								{
									int a=nB%3-1;
									int b=int(nB/3)%3-1;
									int c=int(nB/9)-1;
									int neighborHash=
									((iCoor.x+a+surfN.x)%surfN.x)+
									(((iCoor.y+b+surfN.y)%surfN.y)*surfN.x)+
									(((iCoor.z+c+surfN.z)%surfN.z)*surfN.x*surfN.y);
									//cellHash(iCoor.x+a,iCoor.y+b,iCoor.z+c);
									//std::cerr << hash << '\t' << neighborHash << std::endl;
									
									std::map<int,cellType >::iterator it=cells.find(neighborHash);
									for(int n=0;it!=cells.end() && n<it->second.size();n++)
										//for(int m=0;m<pIndex.size();m++)
									{
										threeVector<double> d;
										
										int q=pIndex[it->second[n]];
										//int l=pIndex[m];
										d.x=p[q].x-sphere.x;
										d.y=p[q].y-sphere.y;
										d.z=p[q].z-sphere.z;
										
										while(d.x>size.x/2.0)d.x-=size.x;
										while(d.y>size.y/2.0)d.y-=size.y;
										while(d.z>size.z/2.0)d.z-=size.z;
										while(d.x<-size.x/2.0)d.x+=size.x;
										while(d.y<-size.y/2.0)d.y+=size.y;
										while(d.z<-size.z/2.0)d.z+=size.z;
										
										if(d.x*d.x+d.y*d.y+d.z*d.z<radSqr)
										{
											//didn't work, quit!
											mkEdge.clear();
											it=cells.end();
											nB=28;
										}
									}
								}
							}
							
							//we just start on the next potential edge if no neighbors are shared
							for(int kElement=0;kElement<mkEdge.size();kElement++)
							{
								int k=edges[mkEdge[kElement]].y;
								
								bool success=true;
								
								//threeVector<double> da,db,dc;
								//ab
								da.x=p[k].x-p[m].x;
								da.y=p[k].y-p[m].y;
								da.z=p[k].z-p[m].z;
								while(da.x>size.x/2.0)da.x-=size.x;
								while(da.y>size.y/2.0)da.y-=size.y;
								while(da.z>size.z/2.0)da.z-=size.z;
								while(da.x<-size.x/2.0)da.x+=size.x;
								while(da.y<-size.y/2.0)da.y+=size.y;
								while(da.z<-size.z/2.0)da.z+=size.z;
								//bc
								db.x=p[m].x-p[j].x;
								db.y=p[m].y-p[j].y;
								db.z=p[m].z-p[j].z;
								while(db.x>size.x/2.0)db.x-=size.x;
								while(db.y>size.y/2.0)db.y-=size.y;
								while(db.z>size.z/2.0)db.z-=size.z;
								while(db.x<-size.x/2.0)db.x+=size.x;
								while(db.y<-size.y/2.0)db.y+=size.y;
								while(db.z<-size.z/2.0)db.z+=size.z;
								//ca
								dc.x=p[j].x-p[k].x;
								dc.y=p[j].y-p[k].y;
								dc.z=p[j].z-p[k].z;
								while(dc.x>size.x/2.0)dc.x-=size.x;
								while(dc.y>size.y/2.0)dc.y-=size.y;
								while(dc.z>size.z/2.0)dc.z-=size.z;
								while(dc.x<-size.x/2.0)dc.x+=size.x;
								while(dc.y<-size.y/2.0)dc.y+=size.y;
								while(dc.z<-size.z/2.0)dc.z+=size.z;
								
								scal=circumscribe(da,db,dc);
								threeVector<double> sphere2;
								//a sphere bounding our 3 points, second deluancy criteria
								sphere2.x=scal.x*p[l].x+scal.y*p[j].x+scal.z*p[m].x;
								sphere2.y=scal.x*p[l].y+scal.y*p[j].y+scal.z*p[m].y;
								sphere2.z=scal.x*p[l].z+scal.y*p[j].z+scal.z*p[m].z;
								
								double radSqr=scal.t*scal.t;
if(radSqr!=radSqr)
{
std::cerr << scal.t << " p: " << k << ' ' << m << ' ' << j << std::endl;
std::cerr << sphere.x << '\t' << sphere.y << '\t' << sphere.z << '\t' << radSqr << '\t' << scal.t << std::endl;
std::cin.get();
}
								
								//calculate hash
								threeVector<int> iCoor;
								threeVector<double> corP;
								for(corP.x=sphere2.x;corP.x>size.x;corP.x-=size.x);
								for(corP.y=sphere2.y;corP.y>size.y;corP.y-=size.y);
								for(corP.z=sphere2.z;corP.z>size.z;corP.z-=size.z);
								for(corP.x=corP.x;corP.x<0;corP.x+=size.x);
								for(corP.y=corP.y;corP.y<0;corP.y+=size.y);
								for(corP.z=corP.z;corP.z<0;corP.z+=size.z);
								iCoor.x=static_cast<int>(corP.x/surfR.x)%surfN.x;
								iCoor.y=static_cast<int>(corP.y/surfR.y)%surfN.y;
								iCoor.z=static_cast<int>(corP.z/surfR.z)%surfN.z;
								//int hash=iCoor.x+iCoor.y*surfN.x;
								//threeVector<int> iCoor=cellHash.cellCoord(p[j]);
								//find connected elements and push onto stack
								//check the nearby cells for traversal
								for(int nB=0;nB<27;nB++)
								{
									int a=nB%3-1;
									int b=int(nB/3)%3-1;
									int c=int(nB/9)-1;
									int neighborHash=
									((iCoor.x+a+surfN.x)%surfN.x)+
									(((iCoor.y+b+surfN.y)%surfN.y)*surfN.x)+
									(((iCoor.z+c+surfN.z)%surfN.z)*surfN.x*surfN.y);
									//cellHash(iCoor.x+a,iCoor.y+b,iCoor.z+c);
									//std::cerr << hash << '\t' << neighborHash << std::endl;
									
									std::map<int,cellType >::iterator it=cells.find(neighborHash);
									for(int n=0;it!=cells.end() && n<it->second.size();n++)
										//for(int m=0;m<pIndex.size();m++)
									{
										threeVector<double> d;
										
										int q=pIndex[it->second[n]];
										//int l=pIndex[m];
										d.x=p[q].x-sphere2.x;
										d.y=p[q].y-sphere2.y;
										d.z=p[q].z-sphere2.z;
										
										while(d.x>size.x/2.0)d.x-=size.x;
										while(d.y>size.y/2.0)d.y-=size.y;
										while(d.z>size.z/2.0)d.z-=size.z;
										while(d.x<-size.x/2.0)d.x+=size.x;
										while(d.y<-size.y/2.0)d.y+=size.y;
										while(d.z<-size.z/2.0)d.z+=size.z;
										
										if(d.x*d.x+d.y*d.y+d.z*d.z<radSqr)
										{
											//didn't work, quit!
											it=cells.end();
											nB=28;
											success=false;
										}
									}
								}
								
								
								
								
								if(success)
								{
									//   m
									// l<|>k
									//   j
									// need lm,lj, jm, jl, mj, and ml
									//eNeighbors[lmEdge].push_back(triangles.size());
									//eNeighbors[ljEdge].push_back(triangles.size());
									//eNeighbors[jmEdge].push_back(triangles.size());
									
									//get the flipped edges too, jl
									//key.x=j;
									//key.y=l;
									
									//tVIter jLow=std::lower_bound(edges.begin(),edges.end(),key, compByX);
									//tVIter jHigh=std::upper_bound(jLow,edges.end(),key, compByX);
									//tVIter kIter=std::find(lLow, lHigh, key);
									//tVIter kIter=lHigh;
									//for(tVIter mEdge=jLow;mEdge<jHigh;mEdge++)
									//	if((*mEdge).y==key.y)
									//		kIter=mEdge;
									
									//if(kIter==lHigh)
									//{
									//	std::cerr << "Cannot find " << key.x << " to " << key.y << " edge!" << std::endl;
									//}
									
									//eNeighbors[kIter-edges.begin()].push_back(triangles.size());
									
									//mj
									//key.x=m;
									//key.y=j;
									//tVIter mjIter=std::find(mLow,mHigh, key);
									//tVIter mjIter=lHigh;
									//for(tVIter mEdge=mLow;mEdge<mHigh;mEdge++)
									//	if((*mEdge).y==key.y)
									//		mjIter=mEdge;
									
									//if(mjIter==mHigh)
									//{
									//	std::cerr << "Cannot find " << key.x << " to " << key.y << " edge!" << std::endl;
									//}
									//eNeighbors[mjIter-edges.begin()].push_back(triangles.size());
									
									//ml
									//key.x=m;
									//key.y=l;
									//tVIter mlIter=std::find(mLow,mHigh, key);
									//tVIter mlIter=lHigh;
									//for(tVIter mEdge=mLow;mEdge<mHigh;mEdge++)
									//	if((*mEdge).y==key.y)
									//		mlIter=mEdge;
									
									//if(mlIter==mHigh)
									//{
									//	std::cerr << "Cannot find " << key.x << " to " << key.y << " edge!" << std::endl;
									//}
									//eNeighbors[mjIter-edges.begin()].push_back(triangles.size());
									
									
									//just the triangle
									threeVector<int> triangle;
									triangle.x=l;//edges[lmEdge].x;
									triangle.y=m;//edges[ljEdge].y;
									triangle.z=j;//edges[jmEdge].y;
									//if(triangle.x!=triangle.y && triangle.y!=triangle.z && triangle.x!=triangle.z)
									triangles.push_back(triangle);
									//std::cout << 1 << '\t' << sphere.x << '\t' << sphere.y << '\t' << sphere.z << std::endl;
								}
								
							}
							
						}
					}
				}
			}
			std::cerr << "Found " << triangles.size() << " triangles!" << std::endl;
			sTriangles.push_back(triangles);
			sENeighbors.push_back(eNeighbors);
		}
		
		for(int i=0;i<sEdges.size();i++)
			std::cerr << "Surface " << i << " has " << sEdges[i].size() << " edges." << std::endl;
		
		for(int s=0;s<surfaces.size() && dry==0 && surfaces.size()<10;s++)
		{
			//std::cout << time << '\t';
			if(sTriangles[s].size()>0)
			{
			
			//std::cout << gradSquared << '\t';// << hSum.x << '\t' << hSum.y << std::endl;
			//std::cerr << "Surface " << s << " gradSquared=" << gradSquared << " vectSums=" << aSum.x << ',' << aSum.y << ',' << bSum.x << ',' << bSum.y;
			
			std::string buf;
			buf="triangleMesh_N";
			std::string whatever;
			std::stringstream parseNumber;
			parseNumber << s;
			parseNumber >> whatever;
			
			buf+=whatever;
			buf+="_";
			buf+=argv[1];
			buf+=".stl";
			std::fstream dataFile;
			dataFile.open(buf.c_str(), std::ios::app | std::ios::out);
			
			dataFile << "solid " << argv[1] << '\n';
			//dump the mesh to triangle file
			for(int i=0;i<sTriangles[s].size();i++)
			{
				dataFile << "facet normal " << 1.0 << ' ' << 1.0 << ' ' << 1.0 << '\n';
				dataFile << "\touter loop " << '\n';
				for(int j=0;j<3;j++)
				{
					int k=sTriangles[s][i].s[j];
					if(k!=-1)
					{
						dataFile << "\t\tvertex " <<
							p[k].x << ' ' <<
							p[k].y << ' ' <<
							p[k].z << '\n';
					}
					else
					{
						dataFile << "\t\tnotvertex " <<
							k << ' ' << '\n';
							//k << ' ' <<
							//k << std::endl;
					}
				}
				dataFile << "\tendloop" << '\n';
				dataFile << "endfacet" << '\n';
				
			}
			dataFile << "endsolid " << argv[1] << '\n';
			dataFile.close();
			}
			/*
			buf.clear();
			buf="heightSlice_N";
			//std::string whatever;
			//std::stringstream parseNumber;
			//parseNumber << s;
			//parseNumber >> whatever;
			
			buf+=whatever;
			buf+="_";
			buf+=argv[1];
			buf+=".dat";
			//std::fstream dataFile;
			dataFile.open(buf.c_str(), std::ios::app | std::ios::out);
			
			//for(int i=0;i<heightMap.size();i++)
			{
				int i=40;
				for(int j=0;j<heightMap[i].size();j++)
				{
					
					//dataFile << "1\t" << static_cast<double>(i)*r.x+r.x/2.0 << '\t' << static_cast<double>(j)*r.y+r.y/2.0 << '\t' << heightMap[i][j] << std::endl;
					dataFile << static_cast<double>(j)*r.y+r.y/2.0 << '\t' << heightMap[i][j] << '\t' << gradSquaredMap[i][j] << std::endl;
				}
				dataFile << std::endl;
			}
			//dataFile << std::endl;
			dataFile.close();
			*/
			/*
			buf.clear();
			buf="size_N";
			
			buf+=whatever;
			buf+="_";
			buf+=argv[1];
			buf+=".xyz";
			dataFile.open(buf.c_str(), std::ios::app | std::ios::out);
			
			dataFile << time << '\t' << size.x << '\t' << size.y << std::endl;
			dataFile.close();
			*/
			/*
			buf.clear();
			buf="gradSquared_N";
			
			buf+=whatever;
			buf+="_";
			buf+=argv[1];
			buf+=".xyz";
			dataFile.open(buf.c_str(), std::ios::app | std::ios::out);
			
			dataFile << time << '\t' << gradSquared << std::endl;
			dataFile.close();
			*/
		}
		//if(surfaces.size()==2)
		//	std::cout << std::endl;
		std::cerr << std::endl;
		
	}	
	return 0;
}
