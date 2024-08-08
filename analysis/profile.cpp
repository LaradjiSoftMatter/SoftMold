//Loading an MD simulation "using systemMD.h".
//Everything needed to do a simulation is already included in systemMD.h.
//Even standard libraries like iostream and fstream are inlcuded.
//But, if any extra libraries are needed you can include them here or in systemMD.h.
//Permenent libraries (included with standard library or with this package) should
// go in systemMD.h, and anything you need that is external should go in this file.

#include "../include/MD.h"
#include "../include/system.h"

int main(int argc, char **argv)
{
	if(argc!=4)
	{
		std::cerr << argv[0] << " name direction deltaL\n";
		return 0;
	}
	
	std::stringstream cmdArg;
	double deltaL;
	int direction;
	cmdArg << argv[2] << '\t' << argv[3];
	cmdArg >> direction >> deltaL;
	if(direction>2 || direction<0)
	{
		std::cerr << "direction must be 0, 1, or 2." << std::endl;
		return 0;
	}
	
	//the variables for the simulation, call with empty constructor
	Blob<double> System;
	
	//file input-output script engine
	//name part (argv[1] in this case) is used with no extension (*.mpd)
	Script<double, Blob <double> > fileIO(argv[1],std::ios::in,&System);
	
	//read the variables in
	fileIO.read();
	
	//close the file that it is reading from
	fileIO.close();
	
	position<double> *p=System.getPositions();
	std::vector<double> density;
	double comD=0;
	for(int i=0;i<System.readNParticles();i++)
	{
		int index=floor(p[i].s[direction]/deltaL);
		if(density.size()<=index)
			density.resize(index+1,0.0);
		density[index]++;
		comD+=p[i].s[direction];
	}
	std::cout << comD/static_cast<double>(System.readNParticles()) << std::endl << std::endl;
	for(int i=0;i<density.size();i++)
		std::cout << static_cast<double>(i)*deltaL << '\t' << density[i] << std::endl;

	return 0;
}
