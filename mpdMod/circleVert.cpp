#include <fstream>
#include <cmath>

int main()
{
	int n=100;
	//std::cout << n << "\ntest\n";
	double dTheta=2.0*M_PI/(n);
	double radius=10;
	std::fstream vertOut("circleVert.dat", std::ios::out);
	for(double theta=0;theta<2.0*M_PI-0.00001;theta+=dTheta)
		vertOut << "4 " << cos(theta)*radius << ' ' << sin(theta)*radius << " 0" << std::endl;
	
	std::fstream edgeOut("circleEdge.dat", std::ios::out);
	for(int i=0;i<n;i++)
		edgeOut << i << ' ' << (i+1)%n << std::endl;
	
	return 0;
}
