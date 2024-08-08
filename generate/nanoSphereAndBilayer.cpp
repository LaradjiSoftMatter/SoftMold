/**
 * \brief Creates a bilayer membrane with nanoparticles.
 * Nanoparticle chain placement.
 */

#include "include.h"

int main(int argc, char **argv)
{
	if(argc!=10)
	{
		//7, bilayer, nanoparticle
		std::cout << "Usage: " << argv[0] << " name seed anchorHeadUmin nLipids lipidArealDensity ";
		std::cout << "nanoRadius nanoHeadUmin nNanoparticles nanoNanoUmin\n";
		return 0;
	}

	char *name=argv[1];
	double anchorHeadUmin=atof(argv[3]);
	int nLipids=atoi(argv[4]);
	double arealDensity=atof(argv[5]);
	double nanoRadius=atof(argv[6]);
	double UminNanoHead=atof(argv[7]);
	int nNanoparticles=atoi(argv[8]);
	double nanoNanoUmin=atof(argv[9]);
	
	//the variables for the simulation
	Blob<double> System;
	
	System.setGamma(1.0);
	System.setNTypes(6+nNanoparticles);
	System.setSeed(atoi(argv[2]));
	
	threeVector<bool> wrap;//periodic boundaries flag, be forwarned, this doesn't work yet, it is always on
	wrap.x=true;
	wrap.y=true;
	wrap.z=true;
	System.setPeriodic(wrap);
	System.setCutoff(2.0);
	
	System.setInitialTime(0);
	System.setFinalTime(50000);
	System.setDeltaT(0.01);
	System.setStoreInterval(100);
	System.setMeasureInterval(10);
	System.setDeltaLXY(0.01);//if you don't set it, it doesn't activate
	System.setInitialTemp(3.0);
	System.setFinalTemp(3.0);
	
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
	
	//tail types hydrophobicity
	Umin[TAIL+TAIL*System.readNTypes()]=-6;
	Umax[TAIL+TAIL*System.readNTypes()]=200;
	
	//Nanoparticle-nanoparticle attraction
	for(int i=4;i<4+nNanoparticles;i++)
	{
		for(int j=4;j<4+nNanoparticles;j++)
		{
			if(i!=j)
			{
				Umin[i+j*System.readNTypes()]=nanoNanoUmin;
				Umax[i+j*System.readNTypes()]=100;
			}
		}
	
		Umin[i+HEAD*System.readNTypes()]=UminNanoHead;
		Umax[i+HEAD*System.readNTypes()]=200;
		
		Umin[HEAD+i*System.readNTypes()]=Umin[i+HEAD*System.readNTypes()];
		Umax[HEAD+i*System.readNTypes()]=Umax[i+HEAD*System.readNTypes()];
	}	
	
	
	//two body constants
	std::vector<double> twoBodyFconst=mpd::laradjiRevaleeFC(Umax,Umin,CUTOFF,RMIN);
	for(auto &v:twoBodyFconst)
		System.addTwoBodyFconst(v);
	
	std::vector<double> twoBodyUconst=mpd::laradjiRevaleePC(Umax,Umin,CUTOFF,RMIN);
	for(auto &v:twoBodyUconst)
		System.addTwoBodyUconst(v);
	
	//double nanoRadius=sqrt(((double)nLipids/arealDensity)/(4.0*M_PI));//density adjusted
	threeVector<double>size;
	size=0;
	size.z=80;
	System.setSize(size);
	
	threeVector<double> pos;
	pos.x=0;//System.readSize().x/2.0;
	pos.y=0;//System.readSize().y/2.0;
	pos.z=System.readSize().z/2.0;
	double bondLength=0.7;
	
	double *constants=new double[4];
	constants[0]=bondLength;
	constants[1]=100;
	constants[2]=1;
	constants[3]=100;
	
	//bilayer with cytoskeleton, regular hexagonal cytoskeleton has an aspect ratio sqrt(3/2)=y/x
	bilayer<double>(System, nLipids, 3, pos, bondLength, arealDensity, constants, 4, 1.0);
	
	//Add a nanoparticle
	int nanoType=4;//could use different nanotypes at different lattice positions
	
	position<double> *unitCell=new position<double>[4];
	position<double> *bonded=new position<double>[12];
	
	constants[0]=sqrt(2.0)/4.0;//bondlength=latticeLength*sqrt(2.0)/2.0 where latticeLength=0.5
	constants[1]=2800;//kbond
	
	std::vector< threeVector<double> > nanoPos;
	MTRand randNum(System.readSeed());
	//1.31345933134 is an expansion factor for bond relaxation at 2800 kBond and 0.5 bondLength
	for(int i=0;i<nNanoparticles;i++)
	{
		unitCell[0].x=0;
		unitCell[0].y=0;
		unitCell[0].z=0;
		unitCell[0].type=nanoType+i;
		unitCell[1].x=0.5;
		unitCell[1].y=0.5;
		unitCell[1].z=0;
		unitCell[1].type=nanoType+i;
		unitCell[2].x=0.5;
		unitCell[2].y=0;
		unitCell[2].z=0.5;
		unitCell[2].type=nanoType+i;
		unitCell[3].x=0;
		unitCell[3].y=0.5;
		unitCell[3].z=0.5;
		unitCell[3].type=nanoType+i;
		bonded[0].x=0.5;
		bonded[0].y=0.5;
		bonded[0].z=0;
		bonded[0].type=nanoType+i;
		bonded[1].x=0.5;
		bonded[1].y=0;
		bonded[1].z=0.5;
		bonded[1].type=nanoType+i;
		bonded[2].x=0;
		bonded[2].y=0.5;
		bonded[2].z=0.5;
		bonded[2].type=nanoType+i;
		bonded[3].x=-0.5;
		bonded[3].y=-0.5;
		bonded[3].z=0;
		bonded[3].type=nanoType+i;
		bonded[4].x=-0.5;
		bonded[4].y=0;
		bonded[4].z=-0.5;
		bonded[4].type=nanoType+i;
		bonded[5].x=0;
		bonded[5].y=-0.5;
		bonded[5].z=-0.5;
		bonded[5].type=nanoType+i;
		bonded[6].x=0.5;
		bonded[6].y=-0.5;
		bonded[6].z=0;
		bonded[6].type=nanoType+i;
		bonded[7].x=0.5;
		bonded[7].y=0;
		bonded[7].z=-0.5;
		bonded[7].type=nanoType+i;
		bonded[8].x=0;
		bonded[8].y=0.5;
		bonded[8].z=-0.5;
		bonded[8].type=nanoType+i;
		bonded[9].x=-0.5;
		bonded[9].y=0.5;
		bonded[9].z=0;
		bonded[9].type=nanoType+i;
		bonded[10].x=-0.5;
		bonded[10].y=0;
		bonded[10].z=0.5;
		bonded[10].type=nanoType+i;
		bonded[11].x=0;
		bonded[11].y=-0.5;
		bonded[11].z=0.5;
		bonded[11].type=nanoType+i;
		std::cout << "Placing nanoparticle " << i << "!" << std::endl;
		threeVector<double> toPlace;
		
		
		bool overlap;
		do
		{
			overlap=false;
			double dNanoparticles=1.0;
			double nanoDensity=0;
			if(dNanoparticles<=0)//no chain
			{
				toPlace.x=(System.readSize().x-2.0*nanoRadius)*randNum.rand53()+nanoRadius;
				toPlace.y=(System.readSize().y-2.0*nanoRadius)*randNum.rand53()+nanoRadius;
				if(nanoDensity<=0)//bilayer
					toPlace.z=pos.z+nanoRadius+1.0+3.0*2.0*0.7;
				else//no bilayer
					toPlace.z=(System.readSize().y-2.0*nanoRadius)*randNum.rand53()+nanoRadius;
			}
			else//chain placement
			{
				double xPos, yPos;
				bool outOfBounds;
				if(i!=0)
				{
					do
					{
						outOfBounds=false;
						double theta=2.0*M_PI*randNum.rand53();
						xPos=sin(theta)*(dNanoparticles+2.0*nanoRadius);
						yPos=cos(theta)*(dNanoparticles+2.0*nanoRadius);
						if(nanoPos[i-1].x+xPos>System.readSize().x-nanoRadius)
							outOfBounds=true;
						if(nanoPos[i-1].y+yPos>System.readSize().y-nanoRadius)
							outOfBounds=true;
						if(nanoPos[i-1].x+xPos<nanoRadius)
							outOfBounds=true;
						if(nanoPos[i-1].y+yPos<nanoRadius)
							outOfBounds=true;
					} while(outOfBounds);
					
					toPlace.x=xPos+nanoPos[i-1].x;
					toPlace.y=yPos+nanoPos[i-1].y;
					toPlace.z=pos.z+nanoRadius+1.0+3.0*2.0*0.7;
				}
				else
				{
					toPlace.x=(System.readSize().x-2.0*nanoRadius)*randNum.rand53()+nanoRadius;
					toPlace.y=(System.readSize().y-2.0*nanoRadius)*randNum.rand53()+nanoRadius;
					toPlace.z=pos.z+nanoRadius+1.0+3.0*2.0*0.7;
				}
			}

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
				
				if(sqrt(d.x*d.x+d.y*d.y+d.z*d.z)<2.0*nanoRadius+1.0)
					overlap=true;
			}
		} while(overlap);
		nanoPos.push_back(toPlace);
		nanoSphere<double>(System, nanoPos[i], nanoRadius, 0.5*1.31345933134, unitCell, 4, bonded, 12, constants, 2);
	}
	
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
	delete bonded,unitCell;
	return 0;
}


