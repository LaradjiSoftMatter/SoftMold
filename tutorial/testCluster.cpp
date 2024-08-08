/** \brief Tests clustering algorithms and shows how they are used.
 * There are many different clustering algorithms. The two tested here are exact
 * clustering algorithms based on nearby elements. They will return the groups of
 * clusters to be analyzed. The first algorithm tests all elements against all
 * other elements. The second tests only nearby elements.
 */
#include "../include/algorithms/functions.h"
#include "../include/algorithms/cell.h"

/** \brief An object for comparing elements of an object
 *
 *
 */
class nearbyCompare {
	public:
		nearbyCompare(std::vector<position<double> > &elements, double cutoff, threeVector<double> size)
		{
			e=elements;
			cutoffSqr=cutoff*cutoff;
			pairInteractions.initialize(&e[0], e.size(), cutoff, size);
			pairInteractions.build();
		};
		
		inline int nElements()
		{
			return e.size();
		};
		
		inline bool compare(int &i, int &j)
		{
			bool result=false;
			threeVector<double> d;
			d.x=e[i].x-e[j].x;
			d.y=e[i].y-e[j].y;
			d.z=e[i].z-e[j].z;
			if(d.x*d.x+d.y*d.y+d.z*d.z<cutoffSqr)
			{
				result=true;
			}
			return result;
		};
		
		inline int next(int &i)
		{
			return pairInteractions.query(i);
		}
		
	private:
		
		std::vector<position<double> > e;
		double cutoffSqr;
		Cell<double> pairInteractions;
};


int main()
{
	std::vector<position<double> > p;
	double cutoff=2;
	MTRand randNum(49586);
	int boxLength=10;
	threeVector<double> size;
	size.x=boxLength*5+2.0;
	size.y=boxLength*5+2.0;
	size.z=boxLength*5+2.0;
	std::cout << "Creating " << boxLength*boxLength*boxLength << " clusters..." << std::endl;
	
	for(int a=0;a<boxLength;a++)
	for(int b=0;b<boxLength;b++)
	for(int c=0;c<boxLength;c++)
	{
		for(int i=0;i<50;i++)
		{
			position<double> buf;
			buf.x=static_cast<double>(a)*5.0+randNum.rand53();
			buf.y=static_cast<double>(b)*5.0+randNum.rand53();
			buf.z=static_cast<double>(c)*5.0+randNum.rand53();
			buf.type=1;
			p.push_back(buf);
			//std::cout << buf.type << '\t' << buf.x << '\t' << buf.y << '\t' << buf.z << std::endl;
		}
	}
	
	nearbyCompare toCompare(p, cutoff, size);
	
	std::cout << "Searching for clusters with slow algorithm..." << std::endl;
	
	int begin=time(NULL);
	std::vector< std::vector<int> > clusters=cluster(toCompare);
	int end=time(NULL);
	
	std::cout << "Found " << clusters.size() << " clusters in " << end-begin << " seconds!" << std::endl;
	
	std::cout << "Searching for clusters with fast algorithm..." << std::endl;
	
	begin=time(NULL);
	std::vector< std::vector<int> > fastClusters=clusterNext(toCompare);
	end=time(NULL);
	
	std::cout << "Found " << fastClusters.size() << " clusters in " << end-begin << " seconds!" << std::endl;
	
	/*
	std::fstream outputClusters;
	outputClusters.open("clusters.xyz", std::ios::out);
	if(outputClusters.is_open())
	{
		outputClusters << p.size() << "\ntest\n";
		for(int i=0;i<clusters.size();i++)
		{
			for(int j=0;j<clusters[i].size();j++)
			{
				int index=clusters[i][j];
				outputClusters << i+1 << '\t' 
					<< p[index].x << '\t' 
					<< p[index].y << '\t' 
					<< p[index].z << '\n';
			}
		}
		outputClusters.close();
	}
	*/
	return 0;
}
