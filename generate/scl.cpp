#include "include.h"

int main(int argc, char **argv)
{
	if(argc<3)
	{
		std::cout << "Usage: " << argv[0] << " name seed\n";
		return 0;
	}

	char *name=argv[1];
	
	//the variables for the simulation
	Blob<double> System;
	
	System.setGamma(1.0);
	System.setNTypes(2);
	System.setSeed(atoi(argv[2]));
	
	threeVector<bool> wrap;//periodic boundaries flag, be forwarned, this doesn't work yet, it is always on
	wrap.x=true;
	wrap.y=true;
	wrap.z=true;
	System.setPeriodic(wrap);
	System.setCutoff(2.0);
	
	System.setInitialTime(0);
	System.setFinalTime(1000);
	System.setDeltaT(0.02);
	System.setStoreInterval(0.02);
	System.setMeasureInterval(0.02);
	//System.setDeltaLXY(0.01);//if you don't set it, it doesn't activate
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
	
	//two body constants
	std::vector<double> twoBodyFconst=mpd::laradjiRevaleeFC(Umax,Umin,CUTOFF,RMIN);
	for(auto &v:twoBodyFconst)
		System.addTwoBodyFconst(v);
	
	std::vector<double> twoBodyUconst=mpd::laradjiRevaleePC(Umax,Umin,CUTOFF,RMIN);
	for(auto &v:twoBodyUconst)
		System.addTwoBodyUconst(v);
	
	threeVector<double>size;
	size.x=50;
	size.y=50;
	size.z=50;
	System.setSize(size);
	
	//Adding a simple cubic lattice, 
	double latticeLength=0.5;
	double tolerance=latticeLength*0.1;//10 percent tolerance for connections
	
	molecule<double,fourVector<int> > m;
	m.addConstant(latticeLength);
	m.addConstant(100);
	m.setType(BOND);

	position<double> unitCell[3];//3 vectors from the base particle
	unitCell[0].x=latticeLength*1.0;
	unitCell[0].y=latticeLength*0;
	unitCell[0].z=latticeLength*0;

	unitCell[1].x=latticeLength*0;
	unitCell[1].y=latticeLength*1.0;
	unitCell[1].z=latticeLength*0;
	
	unitCell[2].x=latticeLength*0;
	unitCell[2].y=latticeLength*0;
	unitCell[2].z=latticeLength*1.0;
	
	MTRand randNum(System.readSeed());
	double Vrms=sqrt(3.0*System.readInitialTemp());
	fourVector<int> bond;
	position<double> vertex;
	vertex.type=1;
	threeVector<double> v,a;
	a=0;

	//particle offset
	int offset=System.readNParticles();
	System.allocParticle(1000);//should be enough
	
	for(int i=0;i<10;i++)
	{
		for(int j=0;j<10;j++)
		{
			for(int k=0;k<10;k++)
			{
				//To check other positions
				position<double> *p=System.getPositions();
				
				//base initially assumed to be a particle we are placing
				int baseIndex=System.readNParticles();

				vertex.x=double(i)*latticeLength+5.0;
				vertex.y=double(j)*latticeLength+5.0;
				vertex.z=double(k)*latticeLength+5.0;
				
				//for velocity
				double theta=M_PI*randNum.rand53();
				double phi=M_PI*2.0*randNum.rand53();
				v.x=Vrms*cos(phi)*sin(theta);
				v.y=Vrms*sin(phi)*sin(theta);
				v.z=Vrms*cos(theta);
				
				//check if already placed, distance should be within tolerance
				//This method is extremely crude and inefficient
				bool withinTolerance=false;
				
				for(int n=offset;n<System.readNParticles();n++)
				{
					//our distance
					threeVector<double> d;
					d.x=p[n].x-vertex.x;
					d.y=p[n].y-vertex.y;
					d.z=p[n].z-vertex.z;
					if(sqrt(d.x*d.x+d.y*d.y+d.z*d.z)<tolerance)
					{
						withinTolerance=true;
						baseIndex=n;//reset base index because one is already present
					}
				}

				if(!withinTolerance)
					System.addParticle(vertex,v,a);
				
				//in case this changes when you add a particle
				p=System.getPositions();
				
				//for the arms
				for(int arm=0;arm<3;arm++)
				{
					//nearby
					vertex.x=double(i)*latticeLength+unitCell[arm].x+5.0;
					vertex.y=double(j)*latticeLength+unitCell[arm].y+5.0;
					vertex.z=double(k)*latticeLength+unitCell[arm].z+5.0;
					
					//for velocity
					theta=M_PI*randNum.rand53();
					phi=M_PI*2.0*randNum.rand53();
					v.x=Vrms*cos(phi)*sin(theta);
					v.y=Vrms*sin(phi)*sin(theta);
					v.z=Vrms*cos(theta);
					

					//only add an arm if it isn't at a terminating face
					if(double(i)+unitCell[arm].x<=9.0 && 
						double(j)+unitCell[arm].y<=9.0 &&
						double(k)+unitCell[arm].z<=9.0)
					{
						//check if already placed, distance should be within tolerance
						//This method is extremely crude and inefficient
						withinTolerance=false;
		
						int armIndex=System.readNParticles();

						//check nearby
						for(int n=offset;n<System.readNParticles();n++)
						{
							//our distance
							threeVector<double> d;
							d.x=p[n].x-vertex.x;
							d.y=p[n].y-vertex.y;
							d.z=p[n].z-vertex.z;
							if(sqrt(d.x*d.x+d.y*d.y+d.z*d.z)<tolerance)
							{
								withinTolerance=true;
								armIndex=n;//reset arm index
							}
						}
		
						if(!withinTolerance)
							System.addParticle(vertex,v,a);

						bond.s[0]=baseIndex;
						bond.s[1]=armIndex;
						m.addBond(bond);
					}					
					//in case this changes when you add a particle
					p=System.getPositions();
				}
				
			}
		}
	}
	
	System.addMolecule(m);
	
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
	
	return 0;
}
