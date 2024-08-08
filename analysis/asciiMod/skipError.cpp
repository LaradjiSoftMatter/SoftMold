#include <iostream>
#include <sstream>
#include <string>

int main(int argc,char **argv)
{
	int skip=5;
	if(argc!=2)
	{
		std::cout << "Usage: " << argv[0] << " nSkip < inFile > outFile " << std::endl;
		return 0;
	}
	std::stringstream cmdArg;
	cmdArg << argv[1];
	cmdArg >> skip;
	
	std::string buf;
	for(long line=0;std::getline(std::cin,buf);line++)
	{
		std::stringstream s;
		s << buf;
		double a,b,c;
		s >> a >> b >> c;
		if(line%skip==0)
			std::cout << a << '\t' << b << '\t' << c << std::endl;
		else
			std::cout << a << '\t' << b << '\t' << 0.0 << std::endl;
	}
	return 0;
}
