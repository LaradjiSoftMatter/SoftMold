#include <iostream>
#include <sstream>
#include <cmath>

/*
 *      * <----- Angle from here
 *     /_\
 *leg /   \ leg
 *   /     \
 *  *-------*
 *   stride
 */
double isoApexAngle(double strideLength, double legLength)//aka opposite and radius
{
	double val=(M_PI-2.0*acos(strideLength/(2.0*legLength)));
	if(val!=val)
		return M_PI;
	return val;
}

int main(int argc, char **argv)
{
	if(argc!=3)
	{
		std::cerr << "Usage: " << argv[0] << " strideLength legLength" << std::endl;
		std::cerr << "Returns apex angle of isosceles triangle. " << std::endl;
		return 0;
	}
	
	std::stringstream cmdArg;
	cmdArg << argv[1] << ' ' << argv[2];
	double strideLength,legLength;
	cmdArg >> strideLength >> legLength;
	std::cout << isoApexAngle(strideLength,legLength) << std::endl;
	return 0;
}
