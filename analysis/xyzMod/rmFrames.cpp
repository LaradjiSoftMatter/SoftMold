#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <sstream>

struct pos {
	int t;
	double x,y,z;
};

struct frame {
	std::string header;
	std::vector<pos> p;
	void clear()
	{
		this->p.clear();
		this->header="";
	}
	pos operator- (frame f)
	{
		if(f.p.size()!=this->p.size())
		{
			throw 1;
		}
		pos d;
		for(int i=0;i<this->p.size();i++)
		{
			double buf=this->p[i].t-f.p[i].t;
			d.t+=buf*buf;
			buf=this->p[i].x-f.p[i].x;
			d.x+=buf*buf;
			buf=this->p[i].y-f.p[i].y;
			d.y+=buf*buf;
			buf=this->p[i].z-f.p[i].z;
			d.z+=buf*buf;
		}
		d.t/=this->p.size();
		d.x/=this->p.size();
		d.y/=this->p.size();
		d.z/=this->p.size();
		d.t=sqrt(d.t);
		d.x=sqrt(d.x);
		d.y=sqrt(d.y);
		d.z=sqrt(d.z);
		return d;
	}	
};

void frameErrors(int n)
{
	switch(n)
	{
		case 1:
			std::cerr << "Wrong frame size!" << std::endl;
			break;
	}
}

template <typename STREAM>
STREAM& operator << (STREAM& stream, pos p)
{
	stream << p.t << '\t' << p.x << '\t' << p.y << '\t' << p.z;
	return stream;
}

template <typename STREAM>
STREAM& operator >> (STREAM& stream, pos &p)
{
	stream >> p.t;
	if(stream.bad() || stream.eof())
	{
		throw 1;
	}
	stream >> p.x;
	if(stream.bad() || stream.eof())
	{
		throw 2;
	}
	stream >> p.y;
	if(stream.bad() || stream.eof())
	{
		throw 3;
	}
	stream >> p.z;
	if(stream.bad())
	{
		throw 4;
	}
	return stream;
}

void streamErrors(int n)
{
	switch(n)
	{
		case 1:
			std::cerr << "Error reading type!" << std::endl;
			break;
		case 2:
			std::cerr << "Error reading x!" << std::endl;
			break;
		case 3:
			std::cerr << "Error reading y!" << std::endl;
			break;
		case 4:
			std::cerr << "Error reading z!" << std::endl;
			break;
		default:
			std::cerr << "Unknown error." << std::endl;
	}
}

template <typename STREAM>
STREAM& operator >> (STREAM& stream, frame &f)
{
	int n;
	stream >> n;
	if(stream.eof())
		return stream;
	stream >> f.header;
	for(int i=0;i<n;i++)
	{
		pos p;
		try {
			stream >> p;
		} catch(int e) {
			std::cerr << "Error at particle " << i << std::endl;
			streamErrors(e);
			throw 0;
		}
		f.p.push_back(p);
	}
	return stream;
}

template <typename STREAM>
STREAM& operator << (STREAM& stream, frame f)
{
	stream << f.p.size() << '\n';
	stream << f.header << '\n';
	for(auto& p:f.p)
		stream << p << '\n';
	return stream;
}

int main(int argc, char **argv)
{
	if(argc!=2 && argc!=3 && argc!=4 && argc!=5)
	{
		std::cerr << "Usage: " << argv[0] << " file.xyz nFrames offset threshold" << std::endl;
		return 0;
	}
	
	std::fstream fin(argv[1],std::ios::in);
	std::stringstream cmdArg;
	for(int i=2;i<argc;i++)
		cmdArg << argv[i] << ' ';
	int nFrames=0, offset=0;
	double threshold=-1.0;
	bool testFlag=false;
	if(argc>=3)
		cmdArg >> nFrames;
	if(argc>=4)
		cmdArg >> offset;
	if(argc>=5)
		cmdArg >> threshold;
	
	std::vector<frame> f;
	try {
		frame buf;
		while(fin >> buf)
		{
			if(f.size()>0)
				std::cerr << "Average rms displacement: " << buf-f.back() << std::endl;
			f.push_back(buf);
			buf.clear();
		}
	} catch (int e)
	{
		frameErrors(e);
	}
	return 0;
}
