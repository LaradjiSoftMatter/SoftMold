#include "functions.h"
#include "dataTypes.h"
#include <vector>
#include <math.h>
#include <fstream>
#ifndef TESSELASPHERE_MPD
#define TESSELASPHERE_MPD

template<typename T>
struct Edge
    {
        position<T> head, tail;
        
	/*union
        {
            struct {position<T> head, tail;};
		position<T> s[2];
        };*/

   };

template <typename T>
struct Triangle
    {
        position<T> Fvertex, Svertex, Tvertex;
        /*union
        {
            struct{position<T> Fvertex, Svertex, Tvertex;};
            position<T> s[3];
        };*/
    };


template <typename T>
class tesselaSphere
{
private:


    position<T> initialPoint; // bottom vertex
    position<T> finalPoint; // top vertex
    std::vector<Triangle<T>> SphereTriangle;
    std::vector<position<T>> SphereVetices;
    position<T> CenterPoint;
    int type1 = 1;//set type name of particle
    int type2 = 2;
    std::vector<std::vector<position<T>>> neighbors;
    std::vector<T> EdgeLengths;//this two vector classfiy the Lengths 
    std::vector<std::vector<fourVector<int>>> ClassfiedEdges; //save the Classfied Edges with the vertex index 
    std::vector<fourVector<int>> TriangleIndexList;
    std::vector<fourVector<int>> angleIndexList;
    std::vector<T> cosList;
    std::vector<T> cos_Os; // save the different angle
    std::vector<std::vector<fourVector<int>>> ClassfiedangleIndexList; //classfiy the angle with a two D array
    std::vector<T> vertexArea; //average area for each vertex


    T JanusRatio;
    T SphereRadius;
    T IcoEdgeLength;
    T IcoVolume;

    int TesselateTimes;
    std::vector<std::vector<bool>> EdgeMap;
    std::vector<std::vector<T>> distanceMap;
    T o, m, n; //basic parameter for coordinate
    std::vector<position<T>> IcosahedronVectics;
    std::vector<Triangle<T>> IcosahedronTriangle;
    std::vector<Edge<T>> IcosahedronEdge;

public:
    tesselaSphere(position<T> Center, T radius, int ReduceTimes);
    void InitialTesselaSphere();
    inline Edge<T> buildEdge(position<T> head, position<T> tail);
    inline position<T> midPoint(position<T> firstPoint, position<T> secondPoint);
    inline Triangle<T> buildTriangel(position<T> Fvertex, position<T> Svertex, position<T> Tvertex);
    inline void PrintVertex();
    inline std::vector<Triangle<T>> divideTriangle(std::vector<Triangle<T>> LastTriangle, int dTimes, T Radius);
    inline std::vector<position<T>> arrangeVertex(std::vector<Triangle<T>> CurrentTriangle);
    inline std::vector<std::vector<bool>> arrangeEdge(std::vector<Triangle<T>> CurrentTriangle, std::vector<position<T>> CurrentVertics);
    std::vector<std::vector<position<T>>> PrintNeighbors(std::vector<std::vector<bool>> CurrentEdgeMap,std::vector<position<T>> CurrentVertex);
    std::vector<std::vector<position<T>>> PrintNeighbors1();
    
    void arrangeAngleList();
    void arrangeCosList();
    void PrintIcosahedronVectics();
    void ClassifyCos_O();
    std::vector<position<T>> PrintSphereVertex();
    inline position<T> normalize(position<T> a);
    void classfiedge();
    void JanusRatioSet(T ratio);
    void CalneighborDistance();
    void CoutEdgeLengths();
    void vertOut();
    void EdgeOut();
    std::vector<position<T>> ReturnVertex();
    std::vector<std::vector<T>> ReturndistanceMap();
    std::vector<T> ReturnEdgeLengths();
    std::vector<std::vector<fourVector<int>>> ReturnClassfiedEdges();
    std::vector<fourVector<int>> ReturnangleIndexList();
    std::vector<T> ReturnCosList();//The cos_0 of each angle
    std::vector<T> Returncos_Os(); // Return different cos_Os
    std::vector<std::vector<fourVector<int>>> ReturnClassfiedangleIndexList(); //Return classied index
    void RelcateSphere(std::vector<position<T>> newPosition); //relocate the position of particles, includes the vertices and the center position
    void settype1(int t);
    void settype2(int t);
    void setcenterpointtype(int t);
    void CalVertexArea();
    std::vector<T> ReturnVertexArea();
    position<T> ReturnCenterPoint();
    T ReturnRadius();
    ~tesselaSphere();
};

template <typename T>
tesselaSphere<T>::tesselaSphere (position<T> Center, T radius, int ReduceTimes)
{

    this->IcosahedronVectics.resize(12);
    this->IcosahedronTriangle.resize(20);
    this->IcosahedronEdge.resize(30);
    this->TesselateTimes = ReduceTimes;
    this->SphereRadius = radius;
    this->IcoEdgeLength = 4/std::sqrt(10+2*std::sqrt(5));
    this->o = 0;
    this->m = std::sqrt(50-10*std::sqrt(5))/10*radius;
    this->n = std::sqrt(50+10*std::sqrt(5))/10*radius;
    IcosahedronVectics[0].x = -m;
    IcosahedronVectics[0].y = o;
    IcosahedronVectics[0].z = n;

    IcosahedronVectics[1].x = m;
    IcosahedronVectics[1].y = o;
    IcosahedronVectics[1].z = n;

    IcosahedronVectics[2].x = -m;
    IcosahedronVectics[2].y = o;
    IcosahedronVectics[2].z = -n;

    IcosahedronVectics[3].x = m;
    IcosahedronVectics[3].y = o;
    IcosahedronVectics[3].z = -n;

    IcosahedronVectics[4].x = o;
    IcosahedronVectics[4].y = n;
    IcosahedronVectics[4].z = m;

    IcosahedronVectics[5].x = o;
    IcosahedronVectics[5].y = n;
    IcosahedronVectics[5].z = -m;

    IcosahedronVectics[6].x = o;
    IcosahedronVectics[6].y = -n;
    IcosahedronVectics[6].z = m;

    IcosahedronVectics[7].x = o;
    IcosahedronVectics[7].y = -n;
    IcosahedronVectics[7].z = -m;

    IcosahedronVectics[8].x = n;
    IcosahedronVectics[8].y = m;
    IcosahedronVectics[8].z = o;

    IcosahedronVectics[9].x = -n;
    IcosahedronVectics[9].y = m;
    IcosahedronVectics[9].z = o;

    IcosahedronVectics[10].x = n;
    IcosahedronVectics[10].y = -m;
    IcosahedronVectics[10].z = o;

    IcosahedronVectics[11].x = -n;
    IcosahedronVectics[11].y = -m;
    IcosahedronVectics[11].z = o;

    for(int i = 0; i < IcosahedronVectics.size(); i++)
    {
        IcosahedronVectics[i].type = ReduceTimes+1;
    }
    //PrintIcosahedronVectics();

    //build triangle
    IcosahedronTriangle[0] = buildTriangel(IcosahedronVectics[1], IcosahedronVectics[4], IcosahedronVectics[0]);
    IcosahedronTriangle[1] = buildTriangel(IcosahedronVectics[4], IcosahedronVectics[9], IcosahedronVectics[0]);
    IcosahedronTriangle[2] = buildTriangel(IcosahedronVectics[4], IcosahedronVectics[5], IcosahedronVectics[9]);
    IcosahedronTriangle[3] = buildTriangel(IcosahedronVectics[8], IcosahedronVectics[5], IcosahedronVectics[4]);
    IcosahedronTriangle[4] = buildTriangel(IcosahedronVectics[1], IcosahedronVectics[8], IcosahedronVectics[4]);
    IcosahedronTriangle[5] = buildTriangel(IcosahedronVectics[1], IcosahedronVectics[10], IcosahedronVectics[8]);
    IcosahedronTriangle[6] = buildTriangel(IcosahedronVectics[10], IcosahedronVectics[3], IcosahedronVectics[8]);
    IcosahedronTriangle[7] = buildTriangel(IcosahedronVectics[8], IcosahedronVectics[3], IcosahedronVectics[5]);
    IcosahedronTriangle[8] = buildTriangel(IcosahedronVectics[3], IcosahedronVectics[2], IcosahedronVectics[5]);
    IcosahedronTriangle[9] = buildTriangel(IcosahedronVectics[3], IcosahedronVectics[7], IcosahedronVectics[2]);
    IcosahedronTriangle[10] = buildTriangel(IcosahedronVectics[3], IcosahedronVectics[10], IcosahedronVectics[7]);
    IcosahedronTriangle[11] = buildTriangel(IcosahedronVectics[10], IcosahedronVectics[6], IcosahedronVectics[7]);
    IcosahedronTriangle[12] = buildTriangel(IcosahedronVectics[6], IcosahedronVectics[11], IcosahedronVectics[7]);
    IcosahedronTriangle[13] = buildTriangel(IcosahedronVectics[6], IcosahedronVectics[0], IcosahedronVectics[11]);
    IcosahedronTriangle[14] = buildTriangel(IcosahedronVectics[6], IcosahedronVectics[1], IcosahedronVectics[0]);
    IcosahedronTriangle[15] = buildTriangel(IcosahedronVectics[10], IcosahedronVectics[1], IcosahedronVectics[6]);
    IcosahedronTriangle[16] = buildTriangel(IcosahedronVectics[11], IcosahedronVectics[0], IcosahedronVectics[9]);
    IcosahedronTriangle[17] = buildTriangel(IcosahedronVectics[2], IcosahedronVectics[11], IcosahedronVectics[9]);
    IcosahedronTriangle[18] = buildTriangel(IcosahedronVectics[5], IcosahedronVectics[2], IcosahedronVectics[9]);
    IcosahedronTriangle[19] = buildTriangel(IcosahedronVectics[11], IcosahedronVectics[2], IcosahedronVectics[7]);

    //build edge

    //edge from 0 first layer of Pentagon, vertex from anti-clock direction (1,4,9,11,6)
    IcosahedronEdge[0] = buildEdge(IcosahedronVectics[0], IcosahedronVectics[1]);
    IcosahedronEdge[1] = buildEdge(IcosahedronVectics[0], IcosahedronVectics[4]);
    IcosahedronEdge[2] = buildEdge(IcosahedronVectics[0], IcosahedronVectics[9]);
    IcosahedronEdge[3] = buildEdge(IcosahedronVectics[0], IcosahedronVectics[11]);
    IcosahedronEdge[4] = buildEdge(IcosahedronVectics[0], IcosahedronVectics[6]);
    //edges of first Pentagon, vertex from anti-clock direction 1->4->9->11->6->1
    IcosahedronEdge[5] = buildEdge(IcosahedronVectics[1], IcosahedronVectics[4]);
    IcosahedronEdge[6] = buildEdge(IcosahedronVectics[4], IcosahedronVectics[9]);
    IcosahedronEdge[7] = buildEdge(IcosahedronVectics[9], IcosahedronVectics[11]);
    IcosahedronEdge[8] = buildEdge(IcosahedronVectics[11], IcosahedronVectics[6]);
    IcosahedronEdge[9] = buildEdge(IcosahedronVectics[6], IcosahedronVectics[1]);
    //edges from first Pentagon to second Pentagon, anti-clock direction (1,4,9,11,6) -> (10,8.5,2,7)
    IcosahedronEdge[10] = buildEdge(IcosahedronVectics[1], IcosahedronVectics[10]);
    IcosahedronEdge[11] = buildEdge(IcosahedronVectics[1], IcosahedronVectics[8]);
    IcosahedronEdge[12] = buildEdge(IcosahedronVectics[4], IcosahedronVectics[8]);
    IcosahedronEdge[13] = buildEdge(IcosahedronVectics[4], IcosahedronVectics[5]);
    IcosahedronEdge[14] = buildEdge(IcosahedronVectics[9], IcosahedronVectics[5]);
    IcosahedronEdge[15] = buildEdge(IcosahedronVectics[9], IcosahedronVectics[2]);
    IcosahedronEdge[16] = buildEdge(IcosahedronVectics[11], IcosahedronVectics[2]);
    IcosahedronEdge[17] = buildEdge(IcosahedronVectics[11], IcosahedronVectics[7]);
    IcosahedronEdge[18] = buildEdge(IcosahedronVectics[6], IcosahedronVectics[7]);
    IcosahedronEdge[19] = buildEdge(IcosahedronVectics[6], IcosahedronVectics[10]);
    //edges of Second Pentagon, vertex from anti-clock direction 10->8->5->2->7->10
    IcosahedronEdge[20] = buildEdge(IcosahedronVectics[10], IcosahedronVectics[8]);
    IcosahedronEdge[21] = buildEdge(IcosahedronVectics[8], IcosahedronVectics[5]);
    IcosahedronEdge[22] = buildEdge(IcosahedronVectics[5], IcosahedronVectics[2]);
    IcosahedronEdge[23] = buildEdge(IcosahedronVectics[2], IcosahedronVectics[7]);
    IcosahedronEdge[24] = buildEdge(IcosahedronVectics[7], IcosahedronVectics[10]);
    //edges from second layer to the bottom point 3, vertex from anti-direction 10->8->5->2->7->10
    IcosahedronEdge[25] = buildEdge(IcosahedronVectics[10], IcosahedronVectics[3]);
    IcosahedronEdge[26] = buildEdge(IcosahedronVectics[8], IcosahedronVectics[3]);
    IcosahedronEdge[27] = buildEdge(IcosahedronVectics[5], IcosahedronVectics[3]);
    IcosahedronEdge[28] = buildEdge(IcosahedronVectics[2], IcosahedronVectics[3]);
    IcosahedronEdge[29] = buildEdge(IcosahedronVectics[7], IcosahedronVectics[3]);

    this->CenterPoint = Center;
    this->CenterPoint.type = type1;
    	//std::cout<< "in1" << std::endl;

  /*
  for(int i = 0; i < IcosahedronTriangle.size();i++)
        {
            std::cout << IcosahedronTriangle[i].Fvertex.x << " " <<IcosahedronTriangle[i].Fvertex.y << " " << IcosahedronTriangle[i].Fvertex.z << std::endl;
            std::cout << IcosahedronTriangle[i].Svertex.x << " " <<IcosahedronTriangle[i].Svertex.y << " " << IcosahedronTriangle[i].Svertex.z << std::endl;
            std::cout << IcosahedronTriangle[i].Tvertex.x << " " <<IcosahedronTriangle[i].Tvertex.y << " " << IcosahedronTriangle[i].Tvertex.z << std::endl;
            std::cout << IcosahedronTriangle[i].Fvertex.x << " " <<IcosahedronTriangle[i].Fvertex.y << " " << IcosahedronTriangle[i].Fvertex.z << std::endl;
            std::cout << std::endl;
        }
   */     
SphereTriangle = divideTriangle(IcosahedronTriangle, TesselateTimes, SphereRadius);
	//std::cout<< "in2" << std::endl;
      
        /*
        std::cout << SphereTriangle.size() << std::endl;
        for(int i = 0; i < SphereTriangle.size();i++)
        {
            std::cout << SphereTriangle[i].Fvertex.x << " " <<SphereTriangle[i].Fvertex.y << " " << SphereTriangle[i].Fvertex.z << std::endl;
            std::cout << SphereTriangle[i].Svertex.x << " " <<SphereTriangle[i].Svertex.y << " " << SphereTriangle[i].Svertex.z << std::endl;
            std::cout << SphereTriangle[i].Tvertex.x << " " <<SphereTriangle[i].Tvertex.y << " " << SphereTriangle[i].Tvertex.z << std::endl;
            std::cout << SphereTriangle[i].Fvertex.x << " " <<SphereTriangle[i].Fvertex.y << " " << SphereTriangle[i].Fvertex.z << std::endl;
            std::cout << std::endl;
        }
        */

SphereVetices = arrangeVertex(SphereTriangle);
	//std::cout<< "in3" << std::endl;
        
      /*  
    std::cout << SphereVetices.size() << std::endl;
     std::cout << "txt" << std::endl; 
        for(int i = 0; i < SphereVetices.size(); i++)
        {
            std::cout << SphereVetices[i].type << " " << SphereVetices[i].x << " " <<SphereVetices[i].y << " " << SphereVetices[i].z << std::endl;
        }
        */

EdgeMap = arrangeEdge(SphereTriangle, SphereVetices);

neighbors = PrintNeighbors(EdgeMap, SphereVetices);
	//std::cout<< "in4" << std::endl;
CalneighborDistance();
CoutEdgeLengths();
/*
for(int i = 0; i < EdgeLengths.size();i++)
{
    std::cout << i << " " << EdgeLengths[i]<< std::endl;
}
*/
//	std::cout<< "in5" << std::endl;
classfiedge();

	//std::cout<< "in6" << std::endl;

arrangeAngleList();
	//std::cout<< "in7" << std::endl;
arrangeCosList();
    //std::cout<< "in8" << std::endl;
ClassifyCos_O();
    //std::cout<< "in9" << std::endl;

}

template <typename T>
inline Edge<T> tesselaSphere<T>::buildEdge(position<T> head, position<T> tail)
{
    Edge<T> tempEdge;
    tempEdge.head = head;
    tempEdge.tail = tail;
    return tempEdge;
}


template <typename T>
inline position<T> tesselaSphere<T>::midPoint(position<T> firstPoint, position<T> secondPoint)
{
    position<T> mid;
    mid.x = (firstPoint.x + secondPoint.x)/2;
    mid.y = (firstPoint.y + secondPoint.y)/2;
    mid.z = (firstPoint.z + secondPoint.z)/2;
    return mid;
}

template <typename T>
inline Triangle<T> tesselaSphere<T>::buildTriangel(position<T> Fvertex, position<T> Svertex, position<T> Tvertex)
{
        Triangle<T> tempTriangle;
        tempTriangle.Fvertex = Fvertex;
        tempTriangle.Svertex = Svertex;
        tempTriangle.Tvertex = Tvertex;
        return tempTriangle;
}

template <typename T>
void tesselaSphere<T>::PrintVertex()
{
    std::cout << SphereVetices.size() << std::endl;
    std::cout << "txt" << std::endl;

    for(int i = 0; i < SphereVetices.size(); i++)
    {
        std::cout << SphereVetices[i].type << " " << SphereVetices[i].x << " " << SphereVetices[i].y << " " << SphereVetices[i].z << std::endl;
    }
}

template <typename T>
inline std::vector<Triangle<T>> tesselaSphere<T>::divideTriangle(std::vector<Triangle<T> > LastTriangle, int dTimes, T Radius)
{
    std::vector<Triangle<T>> temT;
    T R = Radius;
    //std::cout << "Radius" << Radius << std::endl;
    if(dTimes > 0)
    {
    for(int i = 0; i < LastTriangle.size(); i++)
    {
        position<T> v1, v2, v3, m12, m13, m23;
        v1 = LastTriangle[i].Fvertex;
        v2 = LastTriangle[i].Svertex;
        v3 = LastTriangle[i].Tvertex;
        m12.x = (v1.x+ v2.x)/2;
        m12.y = (v1.y+ v2.y)/2;
        m12.z = (v1.z+ v2.z)/2;
        m12 = normalize(m12);
        //if(m12.x < 0.00001) m12.x = 0;
        //if(m12.y < 0.00001) m12.y = 0;
        //if(m12.z < 0.00001) m12.z = 0; 
        //m12 = cartesianToSpherical<position<T>>(m12); 
        //m12.x = Radius;
        //std::cout << "Radius" << Radius << " " << m12.x << std::endl;
        //m12 = sphericalToCartesian<position<T>>(m12);
        m13.x = (v1.x+ v3.x)/2;
        m13.y = (v1.y+ v3.y)/2;
        m13.z = (v1.z+ v3.z)/2;
        m13 = normalize(m13);
        //if(m13.x < 0.00001) m13.x = 0;
        //if(m13.y < 0.00001) m13.y = 0;
        //if(m13.z < 0.00001) m13.z = 0;
        //std::cout << "m13b" << " "<< m13.x << " " << m13.y << " " << m13.z << std::endl; 
        //m13 = cartesianToSpherical<position<T>>(m13);
        //std::cout << "m13b1" << " "<< m13.x << " " << m13.y << " " << m13.z << std::endl; 
        //m13.x = Radius;
        //std::cout << "m13b2" << " "<< m13.x << " " << m13.y << " " << m13.z << std::endl; 
        //m13 = sphericalToCartesian<position<T>>(m13);
        //std::cout << "m13b3"<< " " << m12.x << " " << m13.y << " " << m13.z << std::endl; 
            //std::cout << "end" << std::endl; 
        m23.x = (v2.x+ v3.x)/2;
        m23.y = (v2.y+ v3.y)/2;
        m23.z = (v2.z+ v3.z)/2;
        m23 = normalize(m23);
        //if(m23.x < 0.00001) m23.x = 0;
        //if(m23.y < 0.00001) m23.y = 0;
        //if(m23.z < 0.00001) m23.z = 0; 
        //m23 = cartesianToSpherical<position<T>>(m23);
        //m23.x = Radius;
        //m23 = sphericalToCartesian<position<T>>(m23);
        m12.type = dTimes;
        m13.type = dTimes;
        m23.type = dTimes;
        temT.push_back(buildTriangel(v1,m12,m13));
        temT.push_back(buildTriangel(v2,m12,m23));
        temT.push_back(buildTriangel(v3,m13,m23));
        temT.push_back(buildTriangel(m12,m13,m23));
    }
    }
    
    dTimes--;
    if(dTimes > 0)
    {
        return divideTriangle (temT, dTimes, R);
    }
            //std::cout << temT.size() << std::endl;
    return temT;
}

template<typename T>
inline std::vector<position<T>> tesselaSphere<T>::arrangeVertex(std::vector<Triangle<T>> CurrentTriangle)
{
    std::vector<position<T>> Vertices;
    for(int i = 0; i < CurrentTriangle.size(); i++)
    {
        bool Flag = true;
        for(int j = 0; j < Vertices.size(); j++)
        {

            if(Vertices[j].x == CurrentTriangle[i].Fvertex.x && 
            Vertices[j].y == CurrentTriangle[i].Fvertex.y &&
            Vertices[j].z == CurrentTriangle[i].Fvertex.z)
            {
                Flag = false;
                break;
            }
        }
        if(Flag == true)
        {
          Vertices.push_back(CurrentTriangle[i].Fvertex);          
        }

        Flag = true;
        for(int j = 0; j < Vertices.size(); j++)
        {

            if(Vertices[j].x == CurrentTriangle[i].Svertex.x && 
            Vertices[j].y == CurrentTriangle[i].Svertex.y &&
            Vertices[j].z == CurrentTriangle[i].Svertex.z)
            {
                Flag = false;
                break;
            }
        }
        if(Flag == true)
        {
          Vertices.push_back(CurrentTriangle[i].Svertex);          
        }
        Flag = true;
        for(int j = 0; j < Vertices.size(); j++)
        {

            if(Vertices[j].x == CurrentTriangle[i].Tvertex.x && 
            Vertices[j].y == CurrentTriangle[i].Tvertex.y &&
            Vertices[j].z == CurrentTriangle[i].Tvertex.z)
            {
                Flag = false;
                break;
            }
        }
        if(Flag == true)
        {
          Vertices.push_back(CurrentTriangle[i].Tvertex);          
        }

    }

    return Vertices;
}

template<typename T>
inline std::vector<std::vector<bool>> tesselaSphere<T>::arrangeEdge(std::vector<Triangle<T>> CurrentTriangle, std::vector<position<T>> CurrentVertics)
{
    std::vector<std::vector<bool>> EMap;
    EMap.resize(CurrentVertics.size());
    for(int i = 0; i < EMap.size(); i++)
    {
        EMap[i].resize(CurrentVertics.size());
        for(int j = 0; j < CurrentVertics.size();  j++)
        {
            EMap[i][j] = false;
        }
    }
    for(int i =0; i < CurrentTriangle.size(); i++)
    {
        int m, n, l;
        for(int j = 0; j < CurrentVertics.size(); j++)
        {
            if(std::abs(CurrentVertics[j].x - CurrentTriangle[i].Fvertex.x) < 0.001 &&
            std::abs(CurrentVertics[j].y - CurrentTriangle[i].Fvertex.y) < 0.001 &&
            std::abs(CurrentVertics[j].z - CurrentTriangle[i].Fvertex.z) < 0.001)
            {
                m = j;
            }
            else if(std::abs(CurrentVertics[j].x - CurrentTriangle[i].Svertex.x) < 0.001 &&
            std::abs(CurrentVertics[j].y - CurrentTriangle[i].Svertex.y) < 0.001 &&
            std::abs(CurrentVertics[j].z - CurrentTriangle[i].Svertex.z) < 0.001)
            {
                n = j;
            }
            else if(std::abs(CurrentVertics[j].x - CurrentTriangle[i].Tvertex.x) < 0.001 &&
            std::abs(CurrentVertics[j].y - CurrentTriangle[i].Tvertex.y) < 0.001 &&
            std::abs(CurrentVertics[j].z - CurrentTriangle[i].Tvertex.z) < 0.001)
            {
                l = j;
            }
        }
        EMap[m][n] = true;
        EMap[n][m] = true;
        EMap[m][l] = true;
        EMap[l][m] = true;
        EMap[n][l] = true;
        EMap[l][n] = true;
    }
    return EMap;
}

template<typename T>
std::vector<std::vector<position<T>>> tesselaSphere<T>::PrintNeighbors(std::vector<std::vector<bool>> CurrentEdgeMap, std::vector<position<T>> CurrentVertex)
{
    std::vector<std::vector<position<T>>> Neighbors;
    Neighbors.resize(CurrentVertex.size());
    for(int i = 0; i < CurrentEdgeMap.size(); i++)
    {
        Neighbors[i].push_back(CurrentVertex[i]);
        for(int j = 0; j < CurrentEdgeMap.size();j++)
        {
            if(CurrentEdgeMap[i][j] == true)
            {
                Neighbors[i].push_back(CurrentVertex[j]);
            }
        }
    }
    return Neighbors;
}

template<typename T>
inline position<T> tesselaSphere<T>::normalize(position<T> a)
{
    T r;
    r = std::sqrt(a.x*a.x + a.y*a.y + a.z*a.z);
    a.x = a.x/r*this->SphereRadius;
    a.y = a.y/r*this->SphereRadius;
    a.z = a.z/r*this->SphereRadius;
    return a;
}

template<typename T>
void tesselaSphere<T>::JanusRatioSet(T ratio)
{
    T bottomX = 0, topX = 0.0;
    for(int i = 0; i < SphereVetices.size(); i++)
    {
        if(SphereVetices[i].x < bottomX) bottomX = SphereVetices[i].x;
        if(SphereVetices[i].x > topX) topX = SphereVetices[i].x;
    }
    T d = topX - bottomX;
    T interface = bottomX + d*ratio;
    for(int i = 0; i < SphereVetices.size(); i++)
    {
        if(SphereVetices[i].x <= interface) SphereVetices[i].type =type1;
        if(SphereVetices[i].x > interface) SphereVetices[i].type =type2;
    }
}

template<typename T>
void tesselaSphere<T>::CalneighborDistance()
{
    this->distanceMap.resize(EdgeMap.size());
    for(int i = 0; i < this->EdgeMap.size(); i++)
    {
        this->distanceMap[i].resize(EdgeMap[i].size());
        for(int j = 0; j < this->EdgeMap[i].size(); j++)
        {
            distanceMap[i][j] =  0;
            if(EdgeMap[i][j]==true)
            {
                T a = SphereVetices[i].x - SphereVetices[j].x;
                T b = SphereVetices[i].y - SphereVetices[j].y;
                T c = SphereVetices[i].z - SphereVetices[j].z;
                distanceMap[i][j] = std::sqrt(a*a + b*b + c*c);
            }
        }
    }

    //return this->distanceMap;
}



template<typename T>
void tesselaSphere<T>::PrintIcosahedronVectics()
{
    
    std::cout << "12" << std::endl;
     std::cout << "txt" << std::endl;   
    for(int i = 0; i < this->IcosahedronVectics.size(); i++)
    {
        std::cout << "1" << " " << this->IcosahedronVectics[i].x << " " << this->IcosahedronVectics[i].y << " " << this->IcosahedronVectics[i].z << std::endl;
    }
}

template<typename T>
std::vector<std::vector<position<T>>> tesselaSphere<T>::PrintNeighbors1()
{
    return neighbors;

}

template<typename T>
std::vector<position<T>> tesselaSphere<T>::PrintSphereVertex()
{
    return SphereVetices;
}

template<typename T>
void tesselaSphere<T>::vertOut()
{
    std::fstream verO("SphereVertices.dat", std::ios::out);
    for(int i = 0; i < SphereVetices.size();i++)
    {
        verO << "4" << " " << SphereVetices[i].x << " " << SphereVetices[i].y << " " << SphereVetices[i].z << std::endl;
    }
    verO << "4" << " " << 0.0 << " " << 0.0 << " " << 0.0 << std::endl;
    verO.close();

}

template<typename T>
void tesselaSphere<T>::EdgeOut()
{
    //print edge on the sphere
    int EdgeCounter = 0;
    std::fstream EdgeO("SphereEdges.dat",std::ios::out);
    for(int i = 0; i < EdgeMap.size(); i++)
    {
        for(int j = 0; j < i; j++)
        {
            if(EdgeMap[i][j] == true)
            {
                EdgeO << i << " " << j << std::endl;
                EdgeCounter++;
            }
        }
    }
    std::cout << "2nd Part Start Index: " << EdgeCounter << std::endl; 
    //print the edge between index and center point
    for(int i = 0; i < EdgeMap.size(); i++)
    {
        EdgeO << i << " " << EdgeMap.size() << std::endl;
    }
    EdgeO.close();

}

template<typename T>
void tesselaSphere<T>::CoutEdgeLengths()
{
    for(int i = 0; i < distanceMap.size(); i++)
    {
        for(int j = 0; j < i; j++)
        {
            if(distanceMap[i][j] != 0)
            {
            bool EFlag = true;
            for(int k =0; k < EdgeLengths.size(); k++)
            {
                if(abs(distanceMap[i][j] - EdgeLengths[k]) < 0.0000001)
                {
                    EFlag = false;
                }
            }
            if(EFlag == true)
            {
                EdgeLengths.push_back(distanceMap[i][j]);
            }
            }
            
        }
    }
}

template<typename T>
void tesselaSphere<T>::classfiedge()
{
    	ClassfiedEdges.resize(EdgeLengths.size());
        for(int i = 0; i < distanceMap.size(); i++)
	{
			//std::cout<< "c1" << std::endl;
        for(int j = 0; j < i; j++)
		{
				//std::cout<< "c2" << std::endl;
            if(distanceMap[i][j] != 0)
			{	
                	//std::cout<< "c3" << std::endl;
            bool errorFlag = false;//check if some edges didnot record corre
			for(int k = 0; k < EdgeLengths.size(); k++)
			{
					//std::cout<< "c4" << std::endl;
                if(abs(distanceMap[i][j] - EdgeLengths[k])<0.0000001)
                {
                    	//std::cout << i << " " << j << " " << distanceMap[i][j] << std::endl;
                        //std::cout<< "c5" << std::endl;
                    fourVector<int> tempedge;
                    tempedge.s[0] = i;
                    tempedge.s[1] = j;
                    	//std::cout<< "c6" << std::endl;
                    ClassfiedEdges[k].push_back(tempedge);
                    errorFlag = true;
                    	//std::cout<< "c7" << std::endl;
                    break;
                }
			}
            if(errorFlag == false)
            {
                std::cerr << "Something happenend when classfing edges!" << std::endl;
            }
			}

		}
	}

}


template<typename T>
void tesselaSphere<T>::arrangeAngleList()
{
    //Creat triangle index list
    for(int i = 0; i < SphereTriangle.size(); i++)
    {
        fourVector<int> tempangle;
        bool flag1 = false, flag2 = false, flag3 = false;
        for(int j = 0; j < SphereVetices.size(); j++)
        {
            if(std::abs(SphereVetices[j].x - SphereTriangle[i].Fvertex.x) < 0.001 &&
            std::abs(SphereVetices[j].y - SphereTriangle[i].Fvertex.y) < 0.001 &&
            std::abs(SphereVetices[j].z - SphereTriangle[i].Fvertex.z) < 0.001)
            {
                tempangle.s[0] = j;
                flag1 = true;
                continue;
            }
            else if (std::abs(SphereVetices[j].x - SphereTriangle[i].Svertex.x) < 0.001 &&
            std::abs(SphereVetices[j].y - SphereTriangle[i].Svertex.y) < 0.001 &&
            std::abs(SphereVetices[j].z - SphereTriangle[i].Svertex.z) < 0.001)
            {
                tempangle.s[1] = j;
                flag2 = true;
                continue;
            }
            else if (std::abs(SphereVetices[j].x - SphereTriangle[i].Tvertex.x) < 0.001 &&
            std::abs(SphereVetices[j].y - SphereTriangle[i].Tvertex.y) < 0.001 &&
            std::abs(SphereVetices[j].z - SphereTriangle[i].Tvertex.z) < 0.001)
            {
                tempangle.s[2] = j;
                flag3 = true;
                continue;
            }
            
            if(flag1 == true && flag2 == true && flag3 == true)
            {
                break;
            }
            
            if(j == SphereVetices.size() && (flag1 == false || flag2 == false || flag3 == false))
            {
                std::cerr << "A triangle " << i << " did not match the vertices "  << std::endl;
            }
        
        }
        TriangleIndexList.push_back(tempangle);
    }

    //Creat angle list
    for(int i = 0; i < TriangleIndexList.size(); i++)
    {
        int a,b,c;
        fourVector<int> tempangle;
        a = TriangleIndexList[i].s[0];
        b = TriangleIndexList[i].s[1];
        c = TriangleIndexList[i].s[2];
        tempangle.s[0] = a;
        tempangle.s[1] = b;
        tempangle.s[2] = c;
        angleIndexList.push_back(tempangle);
        tempangle.s[0] = b;
        tempangle.s[1] = c;
        tempangle.s[2] = a;
        angleIndexList.push_back(tempangle);
        tempangle.s[0] = c;
        tempangle.s[1] = a;
        tempangle.s[2] = b;
        angleIndexList.push_back(tempangle);
    }


}

template<typename T>
void tesselaSphere<T>::arrangeCosList()
{
    for(int i = 0; i < angleIndexList.size();i++)
    {
        int a = angleIndexList[i].s[0];
        int b = angleIndexList[i].s[1];
        int c = angleIndexList[i].s[2];
        position<T> V1,V2;
        V1.x = SphereVetices[a].x - SphereVetices[b].x;
        V1.y = SphereVetices[a].y - SphereVetices[b].y;
        V1.z = SphereVetices[a].z - SphereVetices[b].z;
        V2.x = SphereVetices[c].x - SphereVetices[b].x;
        V2.y = SphereVetices[c].y - SphereVetices[b].y;
        V2.z = SphereVetices[c].z - SphereVetices[b].z;
        T tempcos = dotProduct(V1, V2)/(magnitude<position<T>>(V1)*magnitude<position<T>>(V2));
        cosList.push_back(tempcos);
    }
}

template<typename T>
void tesselaSphere<T>::ClassifyCos_O()
{
    for(int i = 0; i < cosList.size(); i++)
    {
        bool cosFlag = true;
        for(int j = 0; j < cos_Os.size(); j++)
        {
            if(abs(cosList[i] - cos_Os[j]) < 0.000001)
            {
                cosFlag = false;
                break;
            }
        }
        if(cosFlag == true)
        {
            cos_Os.push_back(cosList[i]);
        }
    }
    ClassfiedangleIndexList.resize(cos_Os.size());
    for(int i = 0; i < cosList.size(); i++)
    {
        bool cosflag = false;
        for(int j = 0; j < cos_Os.size();j++)
        {
            if(abs(cosList[i] - cos_Os[j]) < 0.000001)
            {
                ClassfiedangleIndexList[j].push_back(angleIndexList[i]);
                cosflag = true;
                break;
            }
        }
        if(cosflag == false)
        {
            std::cerr << "cos_O " << i << " " << cosList[i] << "did not match!" << std::endl;
        }
    }
    
}


template<typename T>
void tesselaSphere<T>::RelcateSphere(std::vector<position<T>> newPosition)
{
    if((newPosition.size()-1) != SphereVetices.size())
    {
        std::cout << "This array did not match sphere" << std::endl;
    }
    //relocate the position of vertices
    for(int i = 0; i < SphereVetices.size(); i++)
    {
        SphereVetices[i] = newPosition[i];
    }
    //relocate the position of center point
    CenterPoint = newPosition.back();


    CalneighborDistance();


} 

template<typename T>
void tesselaSphere<T>::CalVertexArea()
{
    vertexArea.resize(SphereVetices.size());
    std::vector<std::vector<double>> tempAreas;
    tempAreas.resize(SphereVetices.size());
    for(int i = 0; i < TriangleIndexList.size(); i++)
    {
        position<T> a = SphereVetices[TriangleIndexList[i].s[0]];
        position<T> b = SphereVetices[TriangleIndexList[i].s[1]];
        position<T> c = SphereVetices[TriangleIndexList[i].s[2]];
        threeVector<T> V1, V2, V3;
        V1.x = b.x - a.x;
        V1.y = b.y - a.y;
        V1.z = b.z - a.z;
        V2.x = c.x - a.x;
        V2.y = c.y - a.y;
        V2.z = c.z - a.z;
        V3 = crossProduct<T>(V1, V2);
        //std::cout << V3.x << " " <<V3.y << " " << V3.z << std::endl; 
        double tempArea = 0.5*magnitude(V3);// the area of the current triangle
        tempAreas[TriangleIndexList[i].s[0]].push_back(tempArea);
        tempAreas[TriangleIndexList[i].s[1]].push_back(tempArea);
        tempAreas[TriangleIndexList[i].s[2]].push_back(tempArea);
    }

    //calculate the average area for each vertex
    for(int i = 0; i < vertexArea.size(); i++)
    {
        vertexArea[i] = 0;
        for(int j = 0; j < tempAreas[i].size(); j++)
        {
            vertexArea[i] += tempAreas[i][j];
        }
        vertexArea[i] /= 3;
        //std::cout << i << " " << vertexArea[i] << std::endl;
    }
    
    //calculate the total area
    double TotalArea = 0;
    for(int i = 0; i < vertexArea.size(); i++)
    {
        TotalArea += vertexArea[i];
    }
    //std::cout << TotalArea << std::endl;

}


template<typename T>    
std::vector<std::vector<T>> tesselaSphere<T>::ReturndistanceMap()
{
    return distanceMap;
}

template<typename T>
   std::vector<T> tesselaSphere<T>::ReturnEdgeLengths()
{
    return EdgeLengths;
}

template<typename T>
std::vector<position<T>> tesselaSphere<T>::ReturnVertex()
{
    return SphereVetices;
}

template<typename T>
std::vector<std::vector<fourVector<int>>> tesselaSphere<T>::ReturnClassfiedEdges()
{
    return ClassfiedEdges;
}


template<typename T>
void tesselaSphere<T>::settype1(int t)
{
    type1 = t;
}

template<typename T>
void tesselaSphere<T>::settype2(int t)
{
    type2 = t;
}

template<typename T>
position<T> tesselaSphere<T>::ReturnCenterPoint()
{
    return CenterPoint;
}

template<typename T>
T tesselaSphere<T>::ReturnRadius()
{
    return SphereRadius;
}

template<typename T>
void tesselaSphere<T>::setcenterpointtype(int t)
{
    CenterPoint.type = t;
}
    
template<typename T>
std::vector<fourVector<int>> tesselaSphere<T>::ReturnangleIndexList()
{
    return angleIndexList;
}

template<typename T>
std::vector<T> tesselaSphere<T>::ReturnCosList()
{
    return cosList;
}

template<typename T>
 std::vector<T> tesselaSphere<T>::Returncos_Os()
 {
     return cos_Os;
 }
    
template<typename T>
std::vector<std::vector<fourVector<int>>> tesselaSphere<T>::ReturnClassfiedangleIndexList()
{
    return ClassfiedangleIndexList;
}


template<typename T>
std::vector<T> tesselaSphere<T>::ReturnVertexArea()
{
    return vertexArea;
}


template <typename T>
tesselaSphere<T>::~tesselaSphere()
{

}
#endif
