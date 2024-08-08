#include "include.h"

#define BB 0
#define BM 1
#define JANUSTYPE 5

int main(int argc, char **argv)
{
	if(argc!=14)
	{
		//14, brush + NP
		std::cout << "Usage: " << argv[0] << " name seed nBrushes arealDensity UminBrushNP nTessellations ratio nNanoparticles nanoRadius UminBrushBrush UmaxBrushBrush nMonomers bondLength temp" << std::endl;
		return 0;
	}
	char *name=argv[1];
	std::stringstream cmdArg;
	
	for(int i=2;i<argc;i++)
		cmdArg << argv[i] << ' ';
	
	int seed,nBrushes,nMonomers,nNanoparticles,nTessellations;
	double arealDensity, UminBrushNP, UminBrushBrush,UmaxBrushBrush,bondLength,temp,nanoRadius,ratio;
	
	cmdArg >> seed >> nBrushes >> arealDensity >> UminBrushNP >> nTessellations >> ratio >> nNanoparticles >> nanoRadius >> UminBrushBrush >> UmaxBrushBrush >> nMonomers >> bondLength >> temp;
	
	//the variables for the simulation
	Blob<double> System;
	
	System.setGamma(1.0);
	System.setNTypes(8);
	System.setSeed(seed);
	
	threeVector<bool> wrap;//periodic boundaries flag, be forwarned, this doesn't work yet, it is always on
	wrap.x=true;
	wrap.y=true;
	wrap.z=true;
	System.setPeriodic(wrap);
	System.setCutoff(2.0);
	
	System.setInitialTime(0);
	System.setFinalTime(100);
	System.setDeltaT(0.02);
	System.setStoreInterval(10);
	System.setMeasureInterval(1);
	System.setInitialTemp(temp);
	System.setFinalTemp(temp);
	
	//initialize constants (force and potential constants)
	//generate force constants
	std::vector<double> Umax(System.readNTypes()*System.readNTypes());
	std::vector<double> Umin(System.readNTypes()*System.readNTypes());
	
	//could be one loop, but it's good for an example
	//of how to set constants
	for(int i=0;i<System.readNTypes();i++)
	{
		for(int j=0;j<System.readNTypes();j++)
		{
			//i=first type, j=second type, indexed grid
			int k=i+j*System.readNTypes();
			Umin[k]=0;
			Umax[k]=100;
		}
	}
	
	//Brush monomers
	Umin[BM+BM*System.readNTypes()]=UminBrushBrush;
	Umax[BM+BM*System.readNTypes()]=UmaxBrushBrush;
	
	Umin[BM+JANUSTYPE*System.readNTypes()]=UminBrushNP;
	Umin[JANUSTYPE+BM*System.readNTypes()]=UminBrushNP;
	
	//two body constants
	std::vector<double> twoBodyFconst=mpd::laradjiRevaleeFC(Umax,Umin,CUTOFF,RMIN);
	for(auto &v:twoBodyFconst)
		System.addTwoBodyFconst(v);
	
	std::vector<double> twoBodyUconst=mpd::laradjiRevaleePC(Umax,Umin,CUTOFF,RMIN);
	for(auto &v:twoBodyUconst)
		System.addTwoBodyUconst(v);
	
	//double radius=sqrt(((double)nBrushes/arealDensity)/(4.0*M_PI));//density adjusted
	threeVector<double>size;
	size=0;
	size.z=static_cast<double>(nMonomers)*bondLength*2.0;
	System.setSize(size);
	
	threeVector<double> pos;
	pos.x=0;
	pos.y=0;
	pos.z=1.0;
	
	double *constants=new double[4];
	constants[0]=bondLength;
	constants[1]=100;
	constants[2]=1;
	constants[3]=0;
	
	//simple brush
	brush<double>(System, nBrushes, nMonomers, pos, bondLength, arealDensity, constants, 4, 1.0, BB, BM);
	
	//Add a nanoparticle
	MTRand randNum(System.readSeed());
	
	std::vector< threeVector<double> > nanoPos;
	
	int nanoOffset=System.readNParticles();
	
	for(int i=0;i<nNanoparticles;i++)
	{
		std::cerr << "Placing nanoparticle " << i << "!" << std::endl;
		threeVector<double> toPlace;
		bool overlap;
		do
		{
			overlap=false;
			toPlace.x=(System.readSize().x-2.0*nanoRadius)*randNum.rand53()+nanoRadius;
			toPlace.y=(System.readSize().y-2.0*nanoRadius)*randNum.rand53()+nanoRadius;
			toPlace.z=static_cast<double>(nMonomers+2)*bondLength+nanoRadius+1.0;
			for(int j=0;j<nanoPos.size();j++)
			{
				threeVector<double> d;
				d.x=nanoPos[j].x-toPlace.x;
				d.y=nanoPos[j].y-toPlace.y;
				d.z=nanoPos[j].z-toPlace.z;
				d.x-=(d.x>System.readSize().x/2.0)?System.readSize().x:0;
				d.x+=(d.x<-System.readSize().x/2.0)?System.readSize().x:0;
				
				d.y-=(d.y>System.readSize().y/2.0)?System.readSize().y:0;
				d.y+=(d.y<-System.readSize().y/2.0)?System.readSize().y:0;
				
				d.z-=(d.z>System.readSize().z/2.0)?System.readSize().z:0;
				d.z+=(d.z<-System.readSize().z/2.0)?System.readSize().z:0;
				
				if(sqrt(d.x*d.x+d.y*d.y+d.z*d.z)<2.0*nanoRadius+0.5)
					overlap=true;
			}
		} while(overlap);
		std::cerr << i << '\t' << toPlace.x << '\t' << toPlace.y << '\t' << toPlace.z << std::endl;
		nanoPos.push_back(toPlace);
	}
	
	for(auto &nPos:nanoPos)
		janus(System,1200.0,45.0,250.0,0.5,ratio,nPos,nanoRadius,nTessellations);
	
	molecule< double, fourVector<int> > boundaryZ;
	molecule< double, fourVector<int> > boundaryZ2;
	boundaryZ.setType(BOUNDARY);
	boundaryZ2.setType(OFFSET_BOUNDARY);
	for(int i=0;i<System.readNParticles();i++)
	{
		fourVector<int> buf;
		buf.s[0]=i;
		if(System.getPositions()[i].type==JANUSTYPE)
			boundaryZ2.addBond(buf);
		else
			boundaryZ.addBond(buf);
	}
	boundaryZ.addConstant(2);//along Z, direction
	boundaryZ.addConstant(0.0);//at 0=z position, z_0
	boundaryZ.addConstant(4);//new:Doesn't matter old:half width of z, cutoff
	boundaryZ.addConstant(50);//k, epsilon
	
	boundaryZ2.addConstant(2);//along Z, direction
	boundaryZ2.addConstant(0.0);//at 0=z position, z_0
	boundaryZ2.addConstant(nanoRadius);//new:size of nanoparticle
	boundaryZ2.addConstant(50);//k, epsilon
	
	if(boundaryZ.readNBond()>0)
		System.addMolecule(boundaryZ);
	if(boundaryZ2.readNBond()>0)
		System.addMolecule(boundaryZ2);
	
	size=System.readSize();
	size.z*=2.0;
	System.setSize(size);
	
	std::cout << "Storing configuration...";
	
	//Write configuration
	Script<double,Blob <double> > output(name,std::ios::out,&System);
	output.write();
	
	//This stores the xyz file to check new configuration
	position<double> *p=System.getPositions();
	int nParticles=System.readNParticles();
	
	std::string newName("");
	newName+=name;
	newName+=".xyz";
	
	xyzFormat<double> xyzFile(p, nParticles);
	xyzFile.open(newName.c_str(), std::ios::out);
	xyzFile.store();
	xyzFile.close();
	
	std::cout << "Done.\nExiting...\n";
	
	delete constants;
	return 0;
}


