#include <iostream>
#include <stringstream>
#include <fstream>

int main(int argc, char **argv)
{
	if(argc!=2)
	{
		std::cerr << "Usage: " << argv[0] << " trajectory.dat" << std::endl;
		std::cerr << "\ntrajectory.dat is a file with columns {time,x,y,z}
	}

}
