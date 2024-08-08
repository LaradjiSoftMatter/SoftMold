//This is a very short and simple program that shows you how to modify a variable already present in an mpd script.
//The reason this is used is because some variables (forces, potentials, geometry, etc...) are not easily modified
// in the mpd script. Using a file that opens the mpd script, modifies, and then saves it (as below) is far easier
// than modifying the actual variable. Most of it looks like the run.cpp file, it just replaces the "run loop"
// with the modification.

//For molecular dynamics forces and potentials
#include "../include/MD.h"

//For the molecular dynamics variables
#include "../include/system.h"

char decodeDim(int dim)
{
	char out='X';
	if(dim%3==1)
		out='Y';
	if(dim%3==2)
		out='Z';
	return out;
}

int main(int argc, char* argv[])
{
	if(argc<2)
	{
		//so simple, this can't possibly mess it up
		std::cerr << "usage: " << argv[0] << " name [options]\n";
		std::cerr << "\t[options]:\n";
		std::cerr << "\t\t-t [int] [char]\n";
		std::cerr << "\t\t\tUse type [int] with mask [char].\n";
		std::cerr << "\t\t-d [char]\n";
		std::cerr << "\t\t\tUse dimension [char] for projection. Can be x, y, z, X, Y, or Z. Default is z.\n";
		std::cerr << "\t\t-r [double]\n";
		std::cerr << "\t\t\tUse displacement [double] for each 'pixel'. Default is 2.0.\n";
		std::cerr << "\t\t-s [char]\n";
		std::cerr << "\t\t\tUse sign [char] for each 'pixel' direction. Can be '+' (looking up) or '-' (looking down). Default is '-'.\n";
		std::cerr << "\t\t-b\n";
		std::cerr << "\t\t\tUse pbm format.\n";
		std::cerr << "\t\t-e\n";
		std::cerr << "\t\t\tDo not drop empty lines.\n";
		return 0;
	}
	
	std::stringstream cmdArg;
	for(int i=2;i<argc;i++)
		cmdArg << argv[i] << ' ';
	std::string option;
	std::vector<std::pair<unsigned,char> > typeMasks;
	int dim=2;
	double deltaD=2;
	double sign=-1.0;
	bool pbm=false, dropEmpty=true;
	///read in options
	while(cmdArg >> option)
	{
		if(option=="-t" || option=="--typeMask")
		{
			std::pair<int,char> tm;
			if(cmdArg >> tm.first >> tm.second)
				typeMasks.push_back(tm);
			else
			{
				std::cerr << "Bad or incomplete type mask argument!" << std::endl;
				return 1;
			}
		} 
		else if(option=="-b" || option=="--pbm")
		{
			pbm=true;
		}
		else if(option=="-e" || option=="--dropEmpty")
		{
			dropEmpty=false;
		}
		else if(option=="-r" || option=="--displacement")
		{
			if(cmdArg >> deltaD)
				continue;
			else
			{
				std::cerr << "Bad displacement!" << std::endl;
				return 1;
			}
		}
		else if(option=="-s" || option=="--sign")
		{
			char s;
			if(cmdArg >> s)
			{
				if(s=='+')
					sign=1.0;
				else if(s=='-')
					sign=-1.0;
				else
				{
					std::cerr << "Bad sign!" << std::endl;
					return 1;
				}
			}
			else
			{
				std::cerr << "Bad sign argument!" << std::endl;
				return 1;
			}
		}
		else if(option=="-d" || option=="--dimension")
		{
			char d;
			if(cmdArg >> d)
			{
				switch(d)
				{
					case 'x':
					case 'X':
						dim=0;
						break;
					case 'y':
					case 'Y':
						dim=1;
						break;
					case 'z':
					case 'Z':
						dim=2;
						break;
					default:
						std::cerr << "Dimension '" << d << "' not recognized!" << std::endl;
						return 1;
				}
			}
			else
			{
				std::cerr << "Bad dimension argument!" << std::endl;
				return 1;
			}
		}
		else
		{
			std::cerr << "Unrecognized option '" << option << "'." << std::endl;
			return 1;
		}
	}
	
	char *name=argv[1];
	
	///Configuration variables
	Blob<double> System;
	
	///load old variables, then initialize them
	Script<double, Blob <double> > fileIO(name,std::ios::in,&System);
	fileIO.read();
	fileIO.close();
	
	//std::string newName("frames_");
	//newName+=name;
	//newName+=".xyz";
	
	position<double> *p=System.getPositions();
	int nParticles=System.readNParticles();
	
	int nTypes=System.readNTypes();
	for(auto& tm:typeMasks)
	{
		if(tm.first>nTypes)
		{
			std::cerr << "Type '" << tm.first << "' out of range!" << std::endl;
			return 1;
		}
	}
	
	threeVector<double> size=System.readSize();
	twoVector<unsigned> plane;
	plane.x=(dim+1)%3;
	plane.y=(dim+2)%3;
	if(deltaD<0 || size.s[plane.x]<deltaD || size.s[plane.y]<deltaD)
	{
		std::cerr << "Displacement is outside system dimensions!" << std::endl;
		return 1;
	}
	
	molecule<double,fourVector<int>> *m=System.getMolecule();
	int nMol=System.readNMolecules();
	std::vector<int> beads;
	std::vector<double> radii;
	for(int i=0;i<nMol;i++)
	{
		if(m[i].readType()==BEAD)
		{
			fourVector<int> *bond=m[i].getBonds();
			
			for(int j=0;j<m[i].readNBond();j++)
			{
				radii.push_back(m[i].getConstants()[BEADRADIUS]);
				beads.push_back(m[i].getBonds()[j].x);
			}
		}
	}
	//xyzFormat<double> xyzFile(p,nParticles);
	//xyzFile.open(newName.c_str(), std::ios::in);
	
	//while(xyzFile.load())
	{
		twoVector<unsigned> n;
		twoVector<double> c;
		n.x=static_cast<unsigned>(size.s[plane.x]/deltaD);
		n.y=static_cast<unsigned>(size.s[plane.y]/deltaD);
		c.x=size.s[plane.x]/static_cast<double>(n.x);
		c.y=size.s[plane.y]/static_cast<double>(n.y);
		double zeroDepth=0;
		if(sign<0)
			zeroDepth=sign*size.s[dim];
		std::vector<std::vector<int>> image(n.x,std::vector<int>(n.y,0));
		std::vector<std::vector<double>> depth(n.x,std::vector<double>(n.y,zeroDepth));
		for(int i=0;i<nParticles;i++)
		{
			twoVector<unsigned> cell;
			cell.x=static_cast<unsigned>(p[i].s[plane.x]/c.x);
			cell.y=static_cast<unsigned>(p[i].s[plane.y]/c.y);
			if(depth[cell.x][cell.y]<sign*p[i].s[dim])
			{
				depth[cell.x][cell.y]=sign*p[i].s[dim];
				image[cell.x][cell.y]=p[i].type;
			}
			for(double radius=0;radius<0.5;radius+=deltaD)
			{
				double deltaTheta=deltaD/(2.0*radius);
				for(double theta=0;theta<2.0*M_PI;theta+=deltaTheta)
				{
					twoVector<double> r;
					r.x=p[i].s[plane.x]+radius*cos(theta);
					r.y=p[i].s[plane.y]+radius*sin(theta);
					r.x-=(r.x>=size.s[plane.x])?size.s[plane.x]:0.0;
					r.y-=(r.y>=size.s[plane.y])?size.s[plane.y]:0.0;
					r.x+=(r.x<0.0)?size.s[plane.x]:0.0;
					r.y+=(r.y<0.0)?size.s[plane.y]:0.0;
					twoVector<unsigned> rCell;
					rCell.x=static_cast<unsigned>(r.x/c.x);
					rCell.y=static_cast<unsigned>(r.y/c.y);
					double height=sqrt(0.25-radius*radius);
					if(depth[rCell.x][rCell.y]<sign*(p[i].s[dim]+height))
					{
						depth[rCell.x][rCell.y]=sign*(p[i].s[dim]+height);
						image[rCell.x][rCell.y]=p[i].type;
					}
				}
			}
			for(int j=0;j<beads.size();j++)
			{
				if(i==beads[j])
				{
					cell.x=static_cast<unsigned>(p[i].s[plane.x]/c.x);
					cell.y=static_cast<unsigned>(p[i].s[plane.y]/c.y);
					if(depth[cell.x][cell.y]<sign*(p[i].s[dim]+radii[j]))
					{
						depth[cell.x][cell.y]=sign*(p[i].s[dim]+radii[j]);
						image[cell.x][cell.y]=p[i].type;
					}
					for(double radius=0;radius<radii[j];radius+=deltaD)
					{
						double deltaTheta=deltaD/(2.0*radius);
						for(double theta=0;theta<2.0*M_PI;theta+=deltaTheta)
						{
							twoVector<double> r;
							r.x=p[i].s[plane.x]+radius*cos(theta);
							r.y=p[i].s[plane.y]+radius*sin(theta);
							r.x-=(r.x>=size.s[plane.x])?size.s[plane.x]:0.0;
							r.y-=(r.y>=size.s[plane.y])?size.s[plane.y]:0.0;
							r.x+=(r.x<0.0)?size.s[plane.x]:0.0;
							r.y+=(r.y<0.0)?size.s[plane.y]:0.0;
							twoVector<unsigned> rCell;
							rCell.x=static_cast<unsigned>(r.x/c.x);
							rCell.y=static_cast<unsigned>(r.y/c.y);
							double height=sqrt(radii[j]*radii[j]-radius*radius);
							if(depth[rCell.x][rCell.y]<sign*(p[i].s[dim]+height))
							{
								depth[rCell.x][rCell.y]=sign*(p[i].s[dim]+height);
								image[rCell.x][rCell.y]=p[i].type;
							}
						}
					}
				}
			}
		}
		double minDepth=zeroDepth+size.s[dim];
		double maxDepth=zeroDepth;
		for(auto& rDepth:depth)
		{
			for(auto& dDepth:rDepth)
			{
				if(dDepth<minDepth)
					minDepth=dDepth;
				if(dDepth>maxDepth)
					maxDepth=dDepth;
			}
		}
		double midDepth=(maxDepth+minDepth)/2.0;
		double nPBM=4096;
		double pbmSteps=(maxDepth-minDepth)/(nPBM-1);
		//std::cout << minDepth << ' ' << maxDepth << ' ' << pbmSteps << std::endl;
		//return 0;
		if(pbm)
		{
			std::cout << "P2" << std::endl;
			std::cout << n.x << ' ' << n.y << std::endl;
			std::cout << int(nPBM) << std::endl;
			std::cout << "# This is a raw depth image for" << std::endl;
			std::cout << "# " << decodeDim(dim+1) << "-" << decodeDim(dim+2);
			std::cout << " plane of " << argv[1] << " with sign " << sign << std::endl;
			for(int i=0;i<n.x;i++)
			{
				for(int j=0;j<n.y;j++)
					if(depth[i][j]!=zeroDepth)
						std::cout << static_cast<int>((depth[i][j]-minDepth)/pbmSteps) << ' ';
						//std::cout << static_cast<int>(depth[i][j]) << ' ';
					else
						std::cout << nPBM << ' ';
				std::cout << std::endl;
			}
		}
		else
		{
			std::cout << decodeDim(dim+2) << "*--->" << std::endl;
			for(int i=0;i<n.y;i++)
			{
				if(i%10==0)
				{
					std::cout << "|";
					std::cout.width(9);
					std::cout << std::left << static_cast<double>(i)*c.y;
				}
			}
			std::cout << '|' << size.s[(dim+2)%3] << std::endl;
			for(int i=0;i<n.x;i++)
			{
				std::vector<char> line;
				for(int j=0;j<n.y;j++)
				{
					char out=' ';
					if(image[i][j]!=0)
					{
						if(depth[i][j]<midDepth)
							out='a'+image[i][j]%26;
						else
							out='A'+image[i][j]%26;
						for(auto& tm:typeMasks)
							if(tm.first==image[i][j])
								out=tm.second;
					}
					line.push_back(out);
				}
				if(i==0) line.push_back(decodeDim(dim+1));
				if(i==1) line.push_back('*');
				if(i==2) line.push_back('|'); 
				if(i==3) line.push_back('|'); 
				if(i==4) line.push_back('V');
				if(i%10==0)
				{
					if(i!=0)
						line.push_back('|');
					std::stringstream number;
					number << static_cast<double>(i)*c.x;
					char nC;
					while(number >> nC)
						line.push_back(nC);
				}
				bool empty=true;
				for(auto& out:line)
					if(out!=' ')
						empty=false;
				if(!empty || dropEmpty)
				{
					for(auto& out:line)
						std::cout << out;
					std::cout << std::endl;
				}
			}
		}
	}
	
	return 0;
}

