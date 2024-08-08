#include <list>
#include <deque>
#include <cstdlib>

/**
 * \brief Data extraction while the program runs.
 * This is the data extraction header. It is currently used for molecular dynamics
 * (see MD.cpp), but it is probably general enough for other usage, like DPD. This
 * is a user modifiable header that can extract quantities like kinetic or potential
 * energy, center of mass, tension, etc... Some objects are defined twice, once here
 * and once in MD.cpp; that shouldn't be a problem unless you run out of memory. The
 * basic usage of this is as follows:
 * 	1. Build your system and this object.
 * 	2. Initialize some values.
 *	3. Run MD until measure interval is complete.
 * 	4. Execute the compute member function of this object.
 * 	5. repeat step 2 unless system is done.
 * 	6. Destroy objects and exit.
 * 
 * What you are modifying here is three parts:
 * 	1. The constructor member is run during build phase (step 1 above).
 * 	2. The initialize member is run during initialization (step 2 above).
 * 	2. The compute member is executed every measure interval (step 4 above).\n
 * 	3. The destructor member is run prior to exiting (step 6 above).
 * 
 * Also, note that you are passing the blob object into the object's constructor; so
 * at an interval, you should have all the information about the system. 
 * 
 * This is also a template. T is the numeric type (float, double, int, etc...) and
 * U is the system blob type (see systemMD.h and systemDPD.h for an example of the blob).
 * 
 * If you plan on modifying the object, then you should probably put quickly computed 
 * variables in the compute member, and persistent variables, which require information
 * from the previous steps, as private members of this object. 
 */

template <typename T, typename U>
class dataExtraction {
	public:
		/**
		 * The constructor, where blob is a pointer (don't forget the reference operator (&)!)
		 * and name is a pointer to a '\0' terminated char string. Execute this after you build
		 * the system (after reading in an mpd file, for example). 
		 */
		dataExtraction(U *blob, char *name);
		
		/**
		 * Alternative constructor, where blob is a pointer (don't forget the reference operator (&)!)
		 * and name is a pointer to a '\0' terminated char string. Execute this after you build
		 * the system (after reading in an mpd file, for example). Also accepts the absolute position
		 * of a particle list, so that you can track diffusion.
		 */
		dataExtraction(U *blob, char *name, position<T> *absolutePosition);
		
		/**
		 * Alternative constructor, where blob is a pointer (don't forget the reference operator (&)!)
		 * and name is a pointer to a '\0' terminated char string. Execute this after you build
		 * the system (after reading in an mpd file, for example). Also accepts a cutoff parameter.
		 */
		dataExtraction(U *blob, char *name, double cutoffParameter);
		
		/**
		 * The destructor, where old memory goes to die. This should be executed after it falls
		 * out of scope, like when a function (e.g. main()) exits.
		 */
		~dataExtraction();
		
		/**
		 * An initialize member. Call this after you build it. It automatically calls compute too.
		 * Mainly, it just does some extra setup prior to your molecular dynamics loop.
		 */
		void initialize();
		
		/**
		 * This does some computations faster. Data extraction is still put off until compute is called.
		 */
		void computeFast();
		
		/**
		 * A compute member. Call this every measure interval to update the data extraction. When
		 * you modify it, don't forget to close files when you are done using them. Typical extraction
		 * is as follows:
		 * 	1. Set values you want.
		 * 	2. Compute against data in the system blob.
		 * 	3. Open a file.
		 * 	4. Write computed data to the file.
		 * 	5. Close the file.
		 * 
		 * Although, you could swap steps 2 and 3.
		 */
		void compute();
		
		/**
		 * Starts calculation of diffusion.
		 */
		void startDiffusion();
		
		
	private:
		//These variables are typically set during construction:
			//The system blob.
			U *System;
			
			//The name of the system.
			char *name;
			
			//list of beads
			std::vector<int> beads;
			std::vector<T> beadRadius;
			
			std::deque< std::vector< std::vector<T> > > profilesToAvg;
			
			//for profiles
			double cutoffParameter;
};

template <typename T, typename U>
dataExtraction<T,U>::dataExtraction(U *blob, char *name, position<T> *absolutePosition)
{
	System=blob;
	this->name=name;
	this->cutoffParameter=1.2;
}

template <typename T, typename U>
dataExtraction<T,U>::dataExtraction(U *blob, char *name, double cutoffParameter)
{
	System=blob;
	this->name=name;
	this->cutoffParameter=cutoffParameter;
	
}

template <typename T, typename U>
dataExtraction<T,U>::dataExtraction(U *blob, char *name)
{
	System=blob;
	this->name=name;
	this->cutoffParameter=1.2;
}

template <typename T, typename U>
dataExtraction<T,U>::~dataExtraction()
{
}

template <typename T, typename U>
void dataExtraction<T,U>::initialize()
{
	//useful references
	position<T> *p=(*System).getPositions();
	threeVector<T> *v=(*System).getVelocities();
	threeVector<T> *a=(*System).getAccelerations();
	threeVector<T> s=(*System).readSize();
	molecule<T, fourVector<int> > *m=(*System).getMolecule();
	int nMolecules=(*System).readNMolecules();
	T cutoffSqr=(*System).readCutoff();
	cutoffSqr*=cutoffSqr;
	/*
	std::fstream dataFile;
	std::string buf("beadNeighbors_");
	buf+=name;
	buf+=".dat";
	dataFile.open(buf.c_str(), std::ios::app | std::ios::out);
	dataFile << "#time nBeadBeadNeighbors dBeadBeadNeighbors nBeadParticleNeighbors dBeadParticleNeighbors" << std::endl;
	dataFile.close();
	*/
	//gather large bead particles
	for(int i=0;i<nMolecules;i++)
	{
		if(m[i].readType()==BEAD)
		{
			for(int j=0;j<m[i].readNBond();j++)
			{
				beads.push_back(m[i].getBonds()[j].x);
				beadRadius.push_back(m[i].getConstants()[BEADRADIUS]);
				//std::cerr << *(beadRadius.end()-1) << std::endl;
				//std::cin.get();
			}
		}
	}
}

template <typename T, typename U>
void dataExtraction<T,U>::compute()
{
	//useful references
	int nP=(*System).readNParticles();
	position<T> *p=(*System).getPositions();
	threeVector<T> *v=(*System).getVelocities();
	threeVector<T> *a=(*System).getAccelerations();
	threeVector<T> s=(*System).readSize();
	T cutoffSqr=(*System).readCutoff();
	cutoffSqr*=cutoffSqr;
	std::vector< threeVector<T> > beadTriplets;
	std::vector< std::vector<T> > chainDistributions;
	std::vector< fourVector<T> > beadNeighbors;
	
	//the dereferenced values from both of these is and index of beads[value]
	std::vector< std::vector<int> > clusters;
	std::vector< std::list<int> > chainLists;
	//std::vector< std::vector<int> > clusterNeighbors((*System).readNTypes(),std::vector<int>());
	//std::vector< std::vector<int> > clusterDuplicates((*System).readNTypes(),std::vector<int>());
	
	std::vector< std::vector<int> > clusterNeighbors(beads.size(),std::vector<int>());
	std::vector< std::vector<int> > clusterDuplicates(beads.size(),std::vector<int>());
	
	//profiles[i] and beads[j], where i=j
	std::vector< std::vector<T> > profiles;
	T beadPotential=0;
	
	//Go through all beads
	for(int i=0;i<beads.size();i++)
	{
		clusters.push_back(std::vector<int>());
		//against all other beads
		for(int j=0;j<beads.size();j++)
		{
			if(i!=j)
			{
				threeVector<T> d;
				d.x=p[beads[i]].x-p[beads[j]].x;
				if(d.x>s.x/2.0) d.x-=s.x;
				if(d.x<-s.x/2.0) d.x+=s.x;
				d.y=p[beads[i]].y-p[beads[j]].y;
				if(d.y>s.y/2.0) d.y-=s.y;
				if(d.y<-s.y/2.0) d.y+=s.y;
				d.z=p[beads[i]].z-p[beads[j]].z;
				if(d.z>s.z/2.0) d.z-=s.z;
				if(d.z<-s.z/2.0) d.z+=s.z;
				
				T cutoffSquared2=beadRadius[i]+beadRadius[j]+2.0;
				cutoffSquared2*=cutoffSquared2;
				
				if(d.x*d.x+d.y*d.y+d.z*d.z<cutoffSquared2)
				{
					clusters[i].push_back(j);
					std::cerr << j << '\t';
				}
			}
		}
		std::cerr << std::endl;
	}
	
	//find all particles that neighbor within cutoff
	for(int j=0;j<nP;j++)
	{
		for(int i=0;i<beads.size();i++)
		{
			if(beads[i]!=j && p[j].type==HEAD)
			{
				threeVector<T> d;
				d.x=p[beads[i]].x-p[j].x;
				if(d.x>s.x/2.0) d.x-=s.x;
				if(d.x<-s.x/2.0) d.x+=s.x;
				d.y=p[beads[i]].y-p[j].y;
				if(d.y>s.y/2.0) d.y-=s.y;
				if(d.y<-s.y/2.0) d.y+=s.y;
				d.z=p[beads[i]].z-p[j].z;
				if(d.z>s.z/2.0) d.z-=s.z;
				if(d.z<-s.z/2.0) d.z+=s.z;
				//std::cerr << d.x*d.x+d.y*d.y+d.z*d.z << ' ' << beadRadius[i]+cutoffParameter << std::endl;
				T cutoffSquared2=beadRadius[i]+cutoffParameter;
				cutoffSquared2*=cutoffSquared2;
				
				if(d.x*d.x+d.y*d.y+d.z*d.z<cutoffSquared2)
				{
					/*
					clusterNeighbors[p[j].type].push_back(j);
					if(clusterNeighbors[p[j].type].size()>1)
					{
						if(*(clusterNeighbors[p[j].type].end()-1)==*(clusterNeighbors[p[j].type].end()-2))
						{
							clusterNeighbors[p[j].type].pop_back();
							clusterDuplicates[p[j].type].push_back(j);
						}
					}
					*/
					
					clusterNeighbors[i].push_back(j);
					if(clusterNeighbors[i].size()>1)
					{
						if(*(clusterNeighbors[i].end()-1)==*(clusterNeighbors[i].end()-2))
						{
							clusterNeighbors[i].pop_back();
							clusterDuplicates[i].push_back(j);
								//std::cerr << "ClusterDup:" << clusterNeighbors[i].size() << ' ' << clusterDuplicates[i].size() << std::endl;
						}
					}
					
				}
			}
		}
	}

	//keeps track of used chain elements
	std::vector<int> chainIndex(clusters.size(), -1);
	
	//go through all particles' lists, N
	for(int i=0;i<clusters.size();i++)
	{
		if(chainIndex[i]==-1 && clusters[i].size()<3 && clusters[i].size()>0)
		{
			//mark it
			chainIndex[i]=chainLists.size();
			
			//start a new list
			chainLists.push_back(std::list<int>());
			
			//with the ith particle's index
			chainLists.back().push_front(i);
			std::cerr << i << '\t';
			
			//recursive traversal requires a stack, although this really isn't needed
			// for linear chains, but it could be extended to branched chains
			std::vector<int> pStack;
			
			//"left"
			pStack.push_back(clusters[i][0]);
			std::cerr << "left " << clusters[i][0] << '\t';
			chainLists.back().push_front(clusters[i][0]);
			
			//traverse, ln(N)->1 in linear case
			while(pStack.size()>0)
			{
				int currentCluster=pStack.back();
				pStack.pop_back();
				
				chainIndex[currentCluster]=chainLists.size()-1;
				
				//grab all ends or non branching chains
				for(int j=0;j<clusters[currentCluster].size() 
					&& clusters[currentCluster].size()<3;j++)
				{
					if(chainIndex[clusters[currentCluster][j]]==-1)
					{
						//chainIndex[clusters[currentCluster][j]]=chainLists.size()-1;
						pStack.push_back(clusters[currentCluster][j]);
						std::cerr << clusters[currentCluster][j] << "!\t";
						chainLists.back().push_front(clusters[currentCluster][j]);
					}
				}
				//chainIndex[currentCluster]=chainLists.size()-1;
			}
			
			//"right"
			if(clusters[i].size()>1)
			{
				pStack.push_back(clusters[i][1]);
				std::cerr << "right " << clusters[i][1] << '\t';
				chainLists.back().push_back(clusters[i][1]);
			}
			
			//traverse, ln(N)->1 in linear case
			while(pStack.size()>0)
			{
				int currentCluster=pStack.back();
				pStack.pop_back();
				
				chainIndex[currentCluster]=chainLists.size()-1;
				
				//grab all ends or non-branching chains
				for(int j=0;j<clusters[currentCluster].size() 
					&& clusters[currentCluster].size()<3;j++)
				{
					if(chainIndex[clusters[currentCluster][j]]==-1)
					{
						//chainIndex[clusters[currentCluster][j]]=chainLists.size()-1;
						pStack.push_back(clusters[currentCluster][j]);
						std::cerr << clusters[currentCluster][j] << "!\t";
						chainLists.back().push_back(clusters[currentCluster][j]);
					}
				}
				//chainIndex[currentCluster]=chainLists.size()-1;
			}
		}
		std::cerr << std::endl;
	}
	
	std::vector<threeVector<T> > membraneNormals;
	
	for(int i=0;i<beads.size();i++)
	{
		T cutoffSquared2=beadRadius[i]+2.0;
		cutoffSquared2*=cutoffSquared2;
		
		std::vector<threeVector<T> > nearbyHeads;
		
		// and all the neighboring head vectors
		threeVector<T> membraneNormal(0);
		T vectorCount=0;
		
		for(int j=0;j<clusterNeighbors[i].size();j++)
		{
			int k=clusterNeighbors[i][j];
			if(p[k].type==HEAD)
			{
				membraneNormal.x+=p[k].x;
				membraneNormal.y+=p[k].y;
				membraneNormal.z+=p[k].z;
				vectorCount++;
			}
		}
		
		if(vectorCount>0)
		{
			membraneNormal.x/=vectorCount;
			membraneNormal.y/=vectorCount;
			membraneNormal.z/=vectorCount;
		}
		
		//get the local membrane normal vector relative to bead<----membrane
		membraneNormal.x=p[beads[i]].x-membraneNormal.x;
		membraneNormal.y=p[beads[i]].y-membraneNormal.y;
		membraneNormal.z=p[beads[i]].z-membraneNormal.z;
		
		//normalize it
		if(magnitude(membraneNormal)>0)
			membraneNormal=unitVector(membraneNormal);
		else
			membraneNormal=0;
		membraneNormals.push_back(membraneNormal);
		
		
		//assume our periodicity is omega=N*M_PI
		T omega=static_cast<T>(clusters[i].size())*M_PI;
		
		//profile is regularly spaced by delta phi
		T dTheta=sin(1.0/beadRadius[i]);
		std::vector<T> profile(static_cast<int>(2.0*M_PI/dTheta)+1,0);
		
//std::cout << i+beads.size()+1 << '\t' << membraneNormal.x << '\t' << membraneNormal.y << '\t' << membraneNormal.z << std::endl;
		
		for(int j=0;j<clusterNeighbors[i].size();j++)
		{
			int k=clusterNeighbors[i][j];
			if(p[k].type==HEAD)
			{
				//eliminate boundaries
				threeVector<T> translated;
				translated.x=p[k].x-p[beads[i]].x;
				if(translated.x>s.x/2.0) translated.x-=s.x;
				if(translated.x<-s.x/2.0) translated.x+=s.x;
				translated.y=p[k].y-p[beads[i]].y;
				if(translated.y>s.y/2.0) translated.y-=s.y;
				if(translated.y<-s.y/2.0) translated.y+=s.y;
				translated.z=p[k].z-p[beads[i]].z;
				if(translated.z>s.z/2.0) translated.z-=s.z;
				if(translated.z<-s.z/2.0) translated.z+=s.z;
				
				translated=unitVector(translated);
//std::cout << i << '\t' << translated.x << '\t' << translated.y << '\t' << translated.z << std::endl;
				
				T phi=0;
				T theta=0;
				
				//they aren't just zero (vertical)
				//if(magnitude(membraneNormal+translated)!=0)
				{
					//translated=unitVector(membraneNormal+translated);
					
					theta=atan2(translated.y,translated.x)+M_PI;
					
					//phi=asin(translated.z);//(z/r) where r=1
					phi=dotProduct(translated,membraneNormal);
					
//std::cout << theta << '\t' << phi << std::endl;
					//fix quadrant issue
//					if(translated.y>0 && translated.x<0)
//						theta+=M_PI/2.0;
//					if(translated.y<0 && translated.x<0)
//						theta+=M_PI;
//					if(translated.y<0 && translated.x>0)
//						theta+=3.0*M_PI/2.0;
//std::cout << theta << '\t' << phi << std::endl;
				}
				
				//the highest line around nanoparticle
				int pPos=static_cast<int>(theta/dTheta);
				if(pPos>=0)
					profile[pPos]=(phi>profile[pPos])?phi:profile[pPos];
				else
					std::cerr << "pPos: " << pPos << std::endl;
			}
		}
		profiles.push_back(profile);
//std::cout << std::endl;
	}
//throw 0;
	
	std::fstream dataFile, dataFile2;
	std::string buf;
	/*
	
	for(int i=0;i<chainLists.size();i++)
	{
		std::cerr << "Chain: ";
		for(std::list<int>::iterator it=chainLists[i].begin();
			it!=chainLists[i].end();it++)
				std::cerr << *it << '\t';
		std::cerr << std::endl;
		if(chainLists[i].size()>2)
		{
			buf.clear();
			buf="beadChains_N";
			std::string whatever;
			std::stringstream parseNumber;
			parseNumber << chainLists[i].size();
			parseNumber >> whatever;
			
			buf+=whatever;
			buf+="_";
			buf+=name;
			buf+=".dat";
			dataFile.open(buf.c_str(), std::ios::app | std::ios::out);
			
			buf.clear();
			buf="beadChainsDistance_N";
			buf+=whatever;
			buf+="_";
			buf+=name;
			buf+=".dat";
			dataFile2.open(buf.c_str(), std::ios::app | std::ios::out);
			
			//if(oldNChains!=chainLists.size())
			//	dataFile << std::endl;
			
			dataFile << (*System).readInitialTime();
			dataFile2 << (*System).readInitialTime();
			
			for(std::list<int>::iterator it=chainLists[i].begin(), last=chainLists[i].begin();
				last!=chainLists[i].end();it++)
			{
				last=it;
				int first=beads[*last];
				last++;
				int second=beads[*last];
				last++;
				int third=beads[*last];
				last++;
				
				threeVector<T> da,db;
				
				da.x=p[first].x-p[second].x;
				if(da.x>s.x/2.0) da.x-=s.x;
				if(da.x<-s.x/2.0) da.x+=s.x;
				da.y=p[first].y-p[second].y;
				if(da.y>s.y/2.0) da.y-=s.y;
				if(da.y<-s.y/2.0) da.y+=s.y;
				da.z=p[first].z-p[second].z;
				if(da.z>s.z/2.0) da.z-=s.z;
				if(da.z<-s.z/2.0) da.z+=s.z;
				
				db.x=p[second].x-p[third].x;
				if(db.x>s.x/2.0) db.x-=s.x;
				if(db.x<-s.x/2.0) db.x+=s.x;
				db.y=p[second].y-p[third].y;
				if(db.y>s.y/2.0) db.y-=s.y;
				if(db.y<-s.y/2.0) db.y+=s.y;
				db.z=p[second].z-p[third].z;
				if(db.z>s.z/2.0) db.z-=s.z;
				if(db.z<-s.z/2.0) db.z+=s.z;
				T avgD=(magnitude(da)+magnitude(db))/2.0;
				da=unitVector(da);
				db=unitVector(db);
				
				//corrected angle
				T angle=acos(-dotProduct(da,db));
				
				dataFile << '\t' << angle;
				dataFile2 << '\t' << avgD;
				
			}
			dataFile << std::endl;
			dataFile.close();
			dataFile2 << std::endl;
			dataFile2.close();
		}
	}
	*/
	/*
	//all possible three body angles
	for(int i=0;i<clusters.size();i++)
	{
		if(clusters[i].size()>1)
		{
			buf.clear();
			buf="beadBendClusters_N";
			std::string whatever;
			std::stringstream parseNumber;
			parseNumber << clusters[i].size();
			parseNumber >> whatever;
			buf+=whatever;
			buf+="_";
			buf+=name;
			buf+=".dat";
			dataFile.open(buf.c_str(), std::ios::app | std::ios::out);
			
			buf.clear();
			buf="beadBendDistances_N";
			buf+=whatever;
			buf+="_";
			buf+=name;
			buf+=".dat";
			dataFile2.open(buf.c_str(), std::ios::app | std::ios::out);
			
			dataFile2 << (*System).readInitialTime();
			dataFile << (*System).readInitialTime();
			
			for(int j=0;j<clusters[i].size();j++)
			{
				for(int k=j+1;k<clusters[i].size();k++)
				{
					int second=beads[i];
					int first=beads[clusters[i][j]];
					int third=beads[clusters[i][k]];
					
					threeVector<T> da,db;
					
					da.x=p[first].x-p[second].x;
					if(da.x>s.x/2.0) da.x-=s.x;
					if(da.x<-s.x/2.0) da.x+=s.x;
					da.y=p[first].y-p[second].y;
					if(da.y>s.y/2.0) da.y-=s.y;
					if(da.y<-s.y/2.0) da.y+=s.y;
					da.z=p[first].z-p[second].z;
					if(da.z>s.z/2.0) da.z-=s.z;
					if(da.z<-s.z/2.0) da.z+=s.z;
					
					db.x=p[second].x-p[third].x;
					if(db.x>s.x/2.0) db.x-=s.x;
					if(db.x<-s.x/2.0) db.x+=s.x;
					db.y=p[second].y-p[third].y;
					if(db.y>s.y/2.0) db.y-=s.y;
					if(db.y<-s.y/2.0) db.y+=s.y;
					db.z=p[second].z-p[third].z;
					if(db.z>s.z/2.0) db.z-=s.z;
					if(db.z<-s.z/2.0) db.z+=s.z;
					
					T avgD=(magnitude(da)+magnitude(db))/2.0;
					
					da=unitVector(da);
					db=unitVector(db);
					
					//corrected angle
					T angle=acos(-dotProduct(da,db));
					
					dataFile << '\t' << angle;
					dataFile2 << '\t' << avgD;
				}
			}
			dataFile << std::endl;
			dataFile.close();
			dataFile2 << std::endl;
			dataFile2.close();
		}
		if(clusters[i].size()==1)
		{
			buf.clear();
			buf="beadBendDistances_N1_";
			buf+=name;
			buf+=".dat";
			dataFile.open(buf.c_str(), std::ios::app | std::ios::out);
			
			dataFile << (*System).readInitialTime();
			
			int second=beads[i];
			int first=beads[clusters[i][0]];
			
			threeVector<T> da;
					
			da.x=p[first].x-p[second].x;
			if(da.x>s.x/2.0) da.x-=s.x;
			if(da.x<-s.x/2.0) da.x+=s.x;
			da.y=p[first].y-p[second].y;
			if(da.y>s.y/2.0) da.y-=s.y;
			if(da.y<-s.y/2.0) da.y+=s.y;
			da.z=p[first].z-p[second].z;
			if(da.z>s.z/2.0) da.z-=s.z;
			if(da.z<-s.z/2.0) da.z+=s.z;
			
			dataFile << (*System).readInitialTime() << '\t' << magnitude(da) << std::endl;
			dataFile.close();
		}
	}
	
	
	buf.clear();
	buf="simpleBeadBeadCluster_";
	buf+=name;
	buf+=".dat";
	dataFile.open(buf.c_str(), std::ios::app | std::ios::out);
	
	double avgClusterSize=0;
	for(int i=0;i<clusters.size();i++)
		avgClusterSize+=static_cast<double>(clusters[i].size());
	if(clusters.size()>0)
		avgClusterSize/=clusters.size();
	dataFile << (*System).readInitialTime() << '\t' << avgClusterSize << std::endl;

	dataFile.close();
	*/
	/*
	position<T> zero;
	zero.type=1;
	zero.x=0;
	zero.y=0;
	zero.z=0;
	std::vector< position<T> > contact(nP*3, zero);
	for(int i=nP;i<nP*2;i++)
		contact[i].type=2;
	for(int i=nP*2;i<nP*3;i++)
		contact[i].type=3;
	*/

	
	//std::cerr << clusterNeighbors[HEAD].size() << std::endl;
	//for(int j=0;j<clusterNeighbors.size();j++)
	if(clusterNeighbors.size()>0)
	{
		std::cerr << "HEADS: " << clusterNeighbors.size() << '\t' << clusterNeighbors[0].size() << std::endl;
		buf.clear();
			buf="headNeighborsPerBead_";
		buf+=name;
		buf+=".dat";
		dataFile.open(buf.c_str(), std::ios::app | std::ios::out);
		dataFile << (*System).readInitialTime() << '\t';
		for(int i=0;i<clusterNeighbors.size();i++)
			dataFile << clusterNeighbors[i].size() << '\t';
		dataFile << std::endl;
		//dataFile << (*System).readInitialTime() << '\t' << clusterNeighbors[HEAD].size();
		/*
		for(std::vector<int>::iterator i=clusterNeighbors[HEAD].begin();i<clusterNeighbors[HEAD].end();i++)
		{
			//contact.push_back(p[(*i)]);
			contact[(*i)+nP]=p[(*i)];
			(*(contact.end()-1)).type=2;
		}
		for(std::vector<int>::iterator i=clusterDuplicates[HEAD].begin();i<clusterDuplicates[HEAD].end();i++)
		{
			//contact.push_back(p[(*i)]);
			contact[(*i)+nP*2]=p[(*i)];
			(*(contact.end()-1)).type=3;
		}
		*/
		
			
		//dataFile << '\t' << clusterDuplicates[HEAD].size() << std::endl;
		
		dataFile.close();
	}
	/*
	std::fstream outXYZ;
	
	buf.clear();
	buf="contactParticles_";
	buf+=name;
	buf+=".xyz";
	outXYZ.open(buf.c_str(), std::ios::out | std::ios::app);
	outXYZ << contact.size() << "\nabc\n";
	//for(std::vector<int>::iterator i=clusterNeighbors[HEAD].begin();i<it;i++)
	//{
		//outXYZ << "before " << (*i) << '\t' << p[(*i)].type << '\t' << p[(*i)].x << '\t' << p[(*i)].y << '\t' << p[(*i)].z << std::endl;
	//	outXYZ << (*i) << '\t' << p[(*i)].type << '\t' << p[(*i)].x << '\t' << p[(*i)].y << '\t' << p[(*i)].z << std::endl;
	//}
	//for(std::vector<int>::iterator i=it;i<clusterNeighbors[HEAD].end();i++)
	//{
	//	outXYZ << "after " << (*i) << '\t' << p[(*i)].type << '\t' << p[(*i)].x << '\t' << p[(*i)].y << '\t' << p[(*i)].z << std::endl;
	//}
	
	for(int i=0;i<contact.size();i++)
		outXYZ << contact[i].type << '\t' << contact[i].x << '\t' << contact[i].y << '\t' << contact[i].z << std::endl;
	outXYZ.close();
	*/
	//throw 1;
	
	profilesToAvg.push_back(profiles);
	double runningAvgLength=10;
	if(profilesToAvg.size()==static_cast<int>(runningAvgLength))
	{
		//profiles of phi as a function of theta around particles
		for(int i=0;i<profiles.size();i++)
		{
			buf.clear();
			buf="beadProfiles_N";
			std::string whatever;
			std::stringstream parseNumber;
			parseNumber << i;
			parseNumber >> whatever;
			buf+=whatever;
			buf+="_";
			buf+=name;
			buf+=".dat";
			dataFile.open(buf.c_str(), std::ios::app | std::ios::out);
			
			//profile is regularly spaced by delta theta
			T dTheta=sin(1.0/beadRadius[i]);
			
			//first perform running average over set
			for(int j=0;j<profiles[i].size();j++)
			{
				profiles[i][j]/=runningAvgLength;
				for(int k=0;k<profilesToAvg.size()-1;k++)
					profiles[i][j]+=profilesToAvg[k][i][j]/runningAvgLength;
			}
			
			//now output our
			for(int j=0;j<profiles[i].size();j++)
			{
				dataFile << dTheta*static_cast<T>(j) << '\t' << profiles[i][j] << std::endl;
			}
			dataFile << std::endl;
			dataFile.close();
		}
		profilesToAvg.pop_front();
	}
	
	/*
	//Go through all molecule structures
	for(int k=0;k<(*System).readNMolecules();k++)
	{
		//pick a structure by type
		switch((*System).getMolecule()[k].readType())
		{
			case BOND:
			{
				break;
			}
			case BEND:
			{
				break;
			}
			case CHAIN:
			{
				//potential+=(*System).doChainPotential(k);
				break;
			}
			case BEAD:
			{
				T tempBeadPotential=(*System).doBeadPotential(k);
				beadNeighbors.push_back((*System).doBeadNeighbors(k));
				beadPotential+=tempBeadPotential;
				break;
			}
			default:
			{
				//does nothing
				break;
			}
		}
	}
	
	buf.clear();
	buf+="beadPotential_";
	buf+=name;
	buf+=".dat";
	dataFile.open(buf.c_str(), std::ios::app | std::ios::out);
	dataFile << (*System).readInitialTime() << '\t' << beadPotential << std::endl;
	dataFile.close();
	
	buf.clear();
	buf="beadNeighbors_";
	buf+=name;
	buf+=".dat";
	dataFile.open(buf.c_str(), std::ios::app | std::ios::out);
	for(int i=0;i<beadNeighbors.size();i++)
	{
		dataFile << (*System).readInitialTime() << '\t';
		dataFile << beadNeighbors[i].x << '\t';
		dataFile << beadNeighbors[i].y << '\t';
		dataFile << beadNeighbors[i].z << '\t';
		dataFile << beadNeighbors[i].t;
		if(i<beadNeighbors.size()-1)
			dataFile << '\t';
	}
	dataFile << std::endl;
	dataFile.close();
	*/
}
