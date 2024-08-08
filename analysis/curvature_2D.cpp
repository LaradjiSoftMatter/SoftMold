//This is a very short and simple program that shows you how to modify a variable already present in an mpd script.
//The reason this is used is because some variables (forces, potentials, geometry, etc...) are not easily modified
// in the mpd script. Using a file that opens the mpd script, modifies, and then saves it (as below) is far easier
// than modifying the actual variable. Most of it looks like the run.cpp file, it just replaces the "run loop"
// with the modification.

//For molecular dynamics forces and potentials
#include "../include/MD.h"

//For the molecular dynamics variables
#include "../include/system.h"

//For cell operations
#include "cellInclude/cell.h"

//For options parsing
#include "optionsInclude/options.h"

template <typename T>
double curvature2D(const T &a, const T &b, const T &c, const int &axis)
{
	int x=(axis+1)%3;
	int y=(axis+2)%3;
	
	double x1=a.s[x];
	double y1=a.s[y];
	double x2=b.s[x];
	double y2=b.s[y];
	double x3=c.s[x];
	double y3=c.s[y];
 	double delx1=x2-x1;
	double delx2=x3-x2;
	double dely1=y2-y1;
	double dely2=y3-y2;
	double delx3=x3-x1;
	double dely3=y3-y1;
	
 	double AA=-x1*dely2+y1*delx2+x2*y3-x3*y2;
 	double BB=(x1*x1+y1*y1)*dely2-(x2*x2+y2*y2)*dely3+(x3*x3+y3*y3)*dely1;
 	double CC=-(x1*x1+y1*y1)*delx2+(x2*x2+y2*y2)*delx3-(x3*x3+y3*y3)*delx1;
 	double xc=-BB/(2.0*AA);
	double yc=-CC/(2.0*AA);
 	double rad=sqrt((xc-x1)*(xc-x1)+(yc-y1)*(yc-y1));
	double curv=-1.0/rad;
	double crossprod=-(delx1*dely2-dely1*delx2);
	if(crossprod<0.0) curv=-curv;
	//if(curv!=curv) return 0;
	return curv;
} 


//sort surface lists by size
template <typename T>
struct sortByAngle {
	T _com;
	int x,y;
	sortByAngle(T com, int axis):_com(com),x((axis+1)%3),y((axis+2)%3){}
	bool operator()(const T &a, const T &b)
	{
		T dac,dbc;
		dac.s[x]=a.s[x]-_com.s[x];
		dac.s[y]=a.s[y]-_com.s[y];
		dbc.s[x]=b.s[x]-_com.s[x];
		dbc.s[y]=b.s[y]-_com.s[y];
		return atan2(dac.s[y],dac.s[x])<atan2(dbc.s[y],dbc.s[x]);
	}
};

int main(int argc, char* argv[])
{
	//initialize options
	mpd::options opt(argc, argv, "Program to do slices along an axis and calculate curvature.");
	
	//define and pack options, similar to optarg, formatting is:
	//optionIterator(long,short,description,dependency,nArguments,[visibleBool])
	opt.emplace_back("--name","-n", "[string] Name of system. Required.","", 1);
	opt.emplace_back("--type","-t", "[integer] Particle types for surfaces. Required.","--name", 1);
	opt.emplace_back("--axisIndex","-a", "[integer] Axis along tube orientation, {0,1,2}={x,y,x}. Required.","--name", 1);
	opt.emplace_back("--cutoff","-c", "[float] Cutoff for surfaces.","--nSurfaces", 1);
	opt.emplace_back("--nSurfaces","-N", "[integer] Expected number of surfaces.","--type", 1);
	opt.emplace_back("--minSurfaceSize","-m", "[integer] Minimum particles in a surface.","--nSurfaces", 1);
	opt.emplace_back("--sliceThickness","-s", "[float] Slice thickness along axis.","--type", 1);
	opt.emplace_back("--nAngles","-A", "[integer] Number of angles to smooth around axis.","--type", 1);
	opt.emplace_back("--mergeFrames","-f", "[integer] Number of frames to merge for statistics.","--type", 1);
	opt.emplace_back("--dCurve","-d", "[float] Probability distribution step length in output.","--type", 1);
	opt.emplace_back("--orientType","-o", "[integer] Particle type for orientation.","--type", 1);
	opt.emplace_back("--maxAngle","-x", "[float] Maximum angle cutoff about orientation type.","--orientType", 1);
	
	//check for depency errors
	if(opt.dependency_error())
		return 0;
	//check for unrecognized options
	if(opt.unrecognized_option())
		return 0;
	//check if "--help" has been specified
	if(opt["--help"].active())
	{
		opt.help();
		return 0;
	}
	
	//variables
	//check if "--name" has not been specified
	if(!opt["--name"].active())
	{
		std::cerr << "--name is required!" << std::endl;
		opt.help();
		return 0;
	}
	if(!opt["--axisIndex"].active())
	{
		std::cerr << "--axisIndex is required!" << std::endl;
		opt.help();
		return 0;
	}
	if(!opt["--type"].active())
	{
		std::cerr << "--type is required!" << std::endl;
		opt.help();
		return 0;
	}
	std::string name=opt["--name"].sargs(0);
	double cutoff=2.0;
	if(opt["--cutoff"].active())
		cutoff=opt["--cutoff"].args<double>(0);
	double sliceThickness=cutoff;
	if(opt["--sliceThickness"].active())
		sliceThickness=opt["--sliceThickness"].args<double>(0);
	double dCurve=0.01;
	if(opt["--dCurve"].active())
		dCurve=opt["--dCurve"].args<double>(0);
	double maxAngle=M_PI;
	if(opt["--maxAngle"].active())
		maxAngle=opt["--maxAngle"].args<double>(0);
	int orientType=-1;
	if(opt["--orientType"].active())
		orientType=opt["--orientType"].args<int>(0);
	int type=opt["--type"].args<int>(0);
	int nSurfaces=1;
	if(opt["--nSurfaces"].active())
		nSurfaces=opt["--nSurfaces"].args<int>(0);
	int minSurfaceSize=10;
	if(opt["--minSurfaceSize"].active())
		minSurfaceSize=opt["--minSurfaceSize"].args<int>(0);
	int axisIndex=opt["--axisIndex"].args<int>(0);
	int nAngles=10;
	if(opt["--nAngles"].active())
		nAngles=opt["--nAngles"].args<int>(0);
	int mergeFrames=0;
	if(opt["--mergeFrames"].active())
		mergeFrames=opt["--mergeFrames"].args<int>(0);
	
	///Configuration variables
	Blob<double> System;
	
	///load old variables, then initialize them
	Script<double, Blob <double> > fileIO(name.c_str(),std::ios::in,&System);
	fileIO.read();
	fileIO.close();
	
	///Initialize variables and density map
	std::vector<int> indices;
	for(int i=0;i<System.readNParticles();i++)
		if(System.getPositions()[i].type==type)
			indices.push_back(i);
	
	std::vector<int> orientIndices;
	for(int i=0;i<System.readNParticles();i++)
		if(System.getPositions()[i].type==orientType)
			orientIndices.push_back(i);
	
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
	
	//Get our frames file
	std::string framesName("frames_");
	framesName+=name;
	framesName+=".xyz";
	
	position<double> *p=System.getPositions();
	int nParticles=System.readNParticles();
	
	xyzFormat<double> xyzFile(p,nParticles);
	xyzFile.open(framesName.c_str(), std::ios::in);
	
	double dAngle=2.0*M_PI/nAngles;
	std::vector<std::vector<std::vector<double>>> 
		curvatures(nSurfaces,std::vector<std::vector<double>>(nAngles,std::vector<double>()));
	std::vector<std::map<int,double>> cHistogram(nSurfaces,std::map<int,double>());
	
	double time=0;
	int frame=0;
	std::vector<std::vector<position<double>>> xyzAll;
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
		
		//make a list of positions
		std::vector<position<double>> q;
		for(auto &i:indices)
			q.push_back(p[i]);
		
		std::vector<position<double>> orientP;
		for(auto &i:orientIndices)
			orientP.push_back(p[i]);
		position<double> orientCom=centerOfMass(&(*orientP.begin()), &(*orientP.end()));
		
		//get the surfaces
		std::vector<std::vector<position<double>>> surfaces;
		if(nSurfaces>1)
		{	
			//initialize cells and cell dimensions
			threeVector<int> nCells;
			nCells.x=size.x/cutoff;
			nCells.y=size.y/cutoff;
			nCells.z=size.z/cutoff;
			threeVector<double> cSize;
			cSize.x=size.x/nCells.x;
			cSize.y=size.y/nCells.y;
			cSize.z=size.z/nCells.z;
			umapfwd<position<double>> cMap=createCMap(&(*q.begin()), &(*q.end()), nCells, cSize);
			surfaces=createSurfaces(cMap, nCells, cSize, size, cutoff);
			//sort the surfaces
			std::sort(surfaces.begin(),surfaces.end(),sortBySize);
			
			//remove surfaces smaller than or equal to minimum surface size
			while(surfaces.back().size()<=minSurfaceSize)
				surfaces.pop_back();
		}
		else
			surfaces.push_back(q);
		
		//check the constraints
		if(surfaces.size()==nSurfaces)
		{
			//std::vector<position<double>> xyzout;
			//dump the surfaces
			std::cerr << time << ' ' << surfaces.size() << ": ";
			for(auto &surface:surfaces)
				std::cerr << surface.size() << ' ';
			std::cerr << std::endl;
			
			std::vector<position<double>> xyz;
			//std::cout << time << ' ';
			int surfN=0;
			for(auto &s:surfaces)
			{
				std::vector<std::vector<position<double>>> slices=
					createAxisSlices(&(*s.begin()), &(*s.end()), axisIndex, sliceThickness, nAngles);
				//double curvature=0;
				//int sliceN=0;
				for(auto &slice:slices)
				{
					if(slice.size()>2)
					{
						position<double> com=centerOfMass(&(*slice.begin()),&(*slice.end()));
						double orientAngle=0;
						if(orientType!=-1)
							orientAngle=getAngle(orientCom,com,axisIndex);
						std::sort(slice.begin(),slice.end(),sortByAngle<position<double>>(com,axisIndex));
						for(auto pos=slice.begin();pos!=slice.end();++pos)
						{
							int offset=pos-slice.begin();
							position<double> aa=*pos;
							//xyzout.push_back(aa);
							//xyzout.back().type=sliceN;
							position<double> bb=*((offset+1)%slice.size()+slice.begin());
							position<double> cc=*((offset+2)%slice.size()+slice.begin());
							//curvature+=curvature2D(aa, bb, cc, axisIndex);
							double curv=curvature2D(aa, bb, cc, axisIndex);
							double angle=getAngle(bb,com,axisIndex);
							angle-=orientAngle;
							angle-=angle>M_PI?2.0*M_PI:0.0;
							angle+=angle<-M_PI?2.0*M_PI:0.0;
							
							if(mergeFrames!=0)
							{
								int i=(angle+M_PI)/dAngle;
								if(i==nAngles)
									i=0;
								//if(bb.x<orientCom.x+5 && bb.x>orientCom.x-5)
									curvatures[surfN][i].push_back(curv);
							}
							else if(fabs(angle)<maxAngle)
							{
								int hIndex=curv/dCurve;
								auto bin=cHistogram[surfN].find(hIndex);
								if(bin==cHistogram[surfN].end())
									cHistogram[surfN][hIndex]=1;
								else
									bin->second++;
							}
						}
					}
				}
				surfN++;
			}
		}
		else//ignore surfaces that don't match constraints
		{
			std::cerr << time << " Continuous surface!" << std::endl;
		}
		if(mergeFrames!=0)
		{
			if(frame%mergeFrames==0 && frame!=0)
			{
				for(int s=0;s<nSurfaces;s++)
				{
					for(int i=0;i<curvatures[s].size();i++)
					{
						double angle=dAngle*i-M_PI+dAngle*0.5;
						double curv=0;
						for(auto &j:curvatures[s][i])
							curv+=j;
						if(curvatures[s][i].size()>0)
							curv/=curvatures[s][i].size();
						std::cout << angle << ' ' << curv << std::endl;
						curvatures[s][i].clear();
					}
					std::cout << std::endl;
				}
			}
		}
		else
		{
			for(auto &cHist:cHistogram)
			{
				double total=0;
				for(auto &hist:cHist)
					total+=hist.second;
				for(auto &hist:cHist)
				{
					double curv=hist.first*dCurve;
					double count=hist.second;
					std::cout << curv << ' ' << count/(total*dCurve) << std::endl;
				}
				std::cout << std::endl;
				cHist.clear();
			}
		}
		time+=System.readStoreInterval();
		frame++;
	}
	
	return 0;
}

