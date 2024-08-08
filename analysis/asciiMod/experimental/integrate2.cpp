#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <string>
#include <limits>

struct twoVector {
	typedef int size_type;
	typedef double value_type;
	double& x() {return s[0];}
	double& y() {return s[1];}
	double x() const {return s[0];}
	double y() const {return s[1];}
	double s[2];
	twoVector(double x, double y){s[0]=x;s[1]=y;}
	twoVector(const twoVector &init){s[0]=init.s[0];s[1]=init.s[1];}
	twoVector()=default;
};

//template <typename T, typename V>
//class polygonArea {
//	public:
	
		template <typename V>
		//T operator() (const V v1)
		typename V::value_type polygonArea(const V v1)
		{
			return 0;
		}
		template <typename V, typename ...Args>
		//T operator() (const V v1, const V v2, Args... args)
		typename V::value_type polygonArea(const V v1, const V v2, Args... args)
		{
			return polygonArea<V>(v2, args...)+(v1.y()+v2.y())*(v1.x()-v2.x())*0.5;
		}
		template <typename V>
		typename V::value_type polygonArea(const std::vector<V> v)
		{
			typename V::value_type sum=0;
			for(int i=0;i<v.size()-1;i++)
				sum+=(v[i].y()+v[i+1].y())*(v[i].x()-v[i+1].x());
			return sum*0.5;
		}
//};

/*
template <typename T>
class trapazoidArea {
	trapazoidArea(T dx, T dy)
	{
		return (b.y()+a.y())*0.5*(b.x()-a.x());
	}
};

template <typename T, V=twoVector<T> >
T rectangleArea(V a, V b)
{
	return (b.x()-a.x())*(b.y()-a.y());
}

template <typename T=double, typename A=trapazoidArea>
struct twoVectorIntegrator {
	T x() {return xC;}
	T y() {return ySum;}
	T xC, xS, xE, ySum;
	twoVectorIntegrator(const twoVector &init, const T &xStart, const T &xEnd):
		xS(xStart),xE(xEnd)
	{}
	twoVectorIntegrator(const twoVector &init, const T &xEnd):
		xS(init.x()),xC(init.x()),xE(xEnd),ySum
	{
		xC=init.s[0];
		s[1]=init.s[1];
	}
	twoVectorIntegrator(const twoVectorIntegrator &init) //copy constructor
	{
		s[0]=init.s[0];
		s[1]=init.s[1];
		xS=init.xS;
		xE=init.xE;
	}
};
*/
bool sortByX(twoVector a, twoVector b)
{
	return a.x()<b.x();	
}

/*
template<typename T>
class iter
{
	protected:
	T value;
		
		using tr = iterator_traits<T>;
		
	public:
		using iterator_type	 = T;
		using iterator_category = typename tr::iterator_category;
		using value_type		= typename tr::value_type;
		using difference_type   = typename tr::difference_type;
		using reference		 = typename tr::reference;
		using pointer		   = typename tr::pointer;
		
		iter() : value(T()) { }
		
		explicit iter(const T& value) : value(value) { }
		
		// Forward iterator requirements
		reference operator*() const { return *value; }
		
		pointer operator->() const { return value; }
		
		iter& operator++() { ++value; return *this; }
		
		iter operator++(int) { return iter(value++); }
		
		// Bidirectional iterator requirements
		iter& operator--() { --value; return *this; }
		
		iter operator--(int) { return iter(value--); }
		
		// Random access iterator requirements
		reference operator[](const difference_type& n) const { return value[n]; }
		
		iter& operator+=(const difference_type& n) { value += n; return *this; }
		
		iter operator+(const difference_type& n) const { return iter(value + n); }
		
		iter& operator-=(const difference_type& n) { value -= n; return *this; }
		
		iter operator-(const difference_type& n) const { return iter(value - n); }
		
		const T& base() const { return value; }
};
*/
/*
template<typename T>
class integrateVector
{
	protected:
		T value;
		
		using tr = iterator_traits<T>;
		
	public:
		using iterator_type	 = T;
		using iterator_category = typename tr::iterator_category;
		using value_type		= typename tr::value_type;
		using difference_type   = typename tr::difference_type;
		using reference		 = typename tr::reference;
		using pointer		   = typename tr::pointer;
		
		iter() : value(T()) { }
		
		explicit iter(const T& value) : value(value) { }
		
		// Forward iterator requirements
		reference operator*() const { return *value; }
		
		pointer operator->() const { return value; }
		
		iter& operator++() { ++value; return *this; }
		
		iter operator++(int) { return iter(value++); }
		
		// Bidirectional iterator requirements
		iter& operator--() { --value; return *this; }
		
		iter operator--(int) { return iter(value--); }
		
		// Random access iterator requirements
		reference operator[](const difference_type& n) const { return value[n]; }
		
		iter& operator+=(const difference_type& n) { value += n; return *this; }
		
		iter operator+(const difference_type& n) const { return iter(value + n); }
		
		iter& operator-=(const difference_type& n) { value -= n; return *this; }
		
		iter operator-(const difference_type& n) const { return iter(value - n); }
		
		const T& base() const { return value; }
};

template <typename V, typename I, typename R>
class integrate {
	public:
		integrate(V::iterator first, V::iterator last, int fOffset=0, int lOffset=0):
			a(first),b(last),c(first),sum(I(a)),aO(fOffset),bO(lOffset)
		{
			if(fOffset<0)
				throw ERROR("integrate: First offset is less than zero!");
			if(lOffset>0)
				throw ERROR("integrate: Last offset is greater than zero!");
		}
		V::iterator begin()
		{
			return a+aO;
		}
		V::iterator end()
		{
			return b+bO;
		}
		V::iterator operator++ ()
		{
			sum+=R(c);
			return c++;
		}
		V::iterator operator-- ()
		{
			if(c!=a)
				sum-=R(c);
			else
				sum=I(c);
			return c--;
		}
		V::value_type operator() ()
		{
			return sum;
		}
		
	private:
		V::value_type sum;
		int count, aO, bO,total;
		V::iterator a,b,c;
}

template <typename V, typename I, typename M, typename L>
V integrate(V::iterator begin, V::iterator end)
{
	V vSum;
	vSum.push_back(*begin)
	for(auto vIt=begin+1;vIt!=end-1;++vIt)
	{
		auto vS=*vIt;
		vS.x()=v[i].x();
		vS.y()+=(v[i].y()+v[i-1].y())*0.5*(v[i].x()-v[i-1].x());
		vSum.push_back(vS);
		if(sumOnlyPrefix.size()==0)
			std::cout << vSum.x() << ' ' << vSum.y()+constant << std::endl;
	}
	return vSum;
}
*/
void showSwitches(int argc, char **argv)
{
	std::cerr << "Integrates a set of values from stdin to stdout. Default is cumulative sum." << std::endl;
	std::cerr << "Usage: " << argv[0] << " [switches] < inFile > outFile" << std::endl;
	std::cerr << "Usage: inStream | " << argv[0] << " [switches] | outStream" << std::endl;
	std::cerr << "Usage: inStream | " << argv[0] << " [switches] > outFile" << std::endl;
	std::cerr << "\nSwitches:" << std::endl;
	
	std::cerr << "--deltaX [float]" << std::endl;
	std::cerr << "\tX values are not given, use deltaX instead." << std::endl;
	std::cerr << "\tLimitation: Can only use --firstX or --lastX." << std::endl;
	std::cerr << "\tLimitation: Assumes point in range [firstX, ...] or [..., lastX]." << std::endl;
	
	std::cerr << "--noX" << std::endl;
	std::cerr << "\tX values are not given." << std::endl;
	std::cerr << "\tLimitation: Must use both --firstX and --lastX." << std::endl;
	std::cerr << "\tLimitation: Assumes points in range [firstX, lastX] are equidistant." << std::endl;
	
	std::cerr << "--firstX [float]" << std::endl;
	std::cerr << "\tX values are given, but start at [float] for partial integration." << std::endl;
	
	std::cerr << "--lastX [float]" << std::endl;
	std::cerr << "\tX values are given, but end at [float] for partial integration." << std::endl;
	
	std::cerr << "--sumOnly [string]" << std::endl;
	std::cerr << "\tOnly output the total sum over the range using [string] as a prefix." << std::endl;
	
	std::cerr << "--constant [float]" << std::endl;
	std::cerr << "\tConstant added to integration." << std::endl;
	
	std::cerr << "--help" << std::endl;
	std::cerr << "\tShow command usage (this)." << std::endl;
	
	twoVector v0(10,10), v1(10,17), v2(17,17), v3(17,10);
	/*************************************************************
		(10,17)*--*(17,17)
		       |  |
		(10,10)*--*(17,10)
	*************************************************************/
	std::cerr << polygonArea(v3,v2,v1,v0,v3) << std::endl;
	std::vector<twoVector> v;
	v.push_back(v3);
	v.push_back(v2);
	v.push_back(v1);
	v.push_back(v0);
	v.push_back(v3);
	std::cerr << polygonArea(v) << std::endl;
}

int main(int argc, char **argv)
{
	std::cout.precision(std::numeric_limits<double>::max_digits10);
	//Process Command Line Arguments:	
	std::stringstream cmdArg;
	for(int i=1;i<argc;i++)
		cmdArg << argv[i] << ' ';
	enum SWITCHESS {DELTAX_S, NOX_S, FIRSTX_S, LASTX_S, SUMONLY_S, CONSTANT_S, HELP_S};
	
	std::map<std::string, SWITCHESS> switches;
	switches[std::string("--deltaX")]=DELTAX_S;
	switches[std::string("--noX")]=NOX_S;
	switches[std::string("--firstX")]=FIRSTX_S;
	switches[std::string("--lastX")]=LASTX_S;
	switches[std::string("--sumOnly")]=SUMONLY_S;
	switches[std::string("--constant")]=CONSTANT_S;
	switches[std::string("--help")]=HELP_S;
	
	double deltaX=0, firstX=0, lastX=0, constant=0;
	std::string sumOnlyPrefix("");
	bool noX=false, fX=false, lX=false;
	
	//std::string option;
	//while(cmdArg >> option)
	for(std::string option; cmdArg >> option;)
	{
		if(switches.find(option)==switches.end())
		{
			std::cerr << "Error: " << option << ": Unknown switch!" << std::endl;
			showSwitches(argc,argv);
			return -1;
		}
		switch(switches[option])
		{
			case DELTAX_S:
				noX=true;
				cmdArg >> deltaX;
				if(cmdArg.eof())
				{
					std::cerr << "Error: " << option << " needs [float] value!" << std::endl;
					showSwitches(argc,argv);
					return -1;
				}
				break;
			case NOX_S:
				noX=true;
				break;
			case FIRSTX_S:
				cmdArg >> firstX;
				fX=true;
				if(cmdArg.eof())
				{
					std::cerr << "Error: " << option << " needs [float] value!" << std::endl;
					showSwitches(argc,argv);
					return -1;
				}
				break;
			case LASTX_S:
				cmdArg >> lastX;
				lX=true;
				if(cmdArg.eof())
				{
					std::cerr << "Error: " << option << " needs [float] value!" << std::endl;
					showSwitches(argc,argv);
					return -1;
				}
				break;
			case SUMONLY_S:
				std::cerr << cmdArg.eof() << ' ' << cmdArg.bad() << ' ' << cmdArg.fail() << std::endl;
				cmdArg >> sumOnlyPrefix;
				std::cerr << "sumOnlyPrefix is " << sumOnlyPrefix << std::endl;
				std::cerr << cmdArg.eof() << ' ' << cmdArg.bad() << ' ' << cmdArg.fail() << std::endl;
				if(cmdArg.eof())
				{
					std::cerr << "Error: " << option << " needs [string] value!" << std::endl;
					showSwitches(argc,argv);
					return -1;
				}
				break;
			case CONSTANT_S:
				cmdArg >> constant;
				if(cmdArg.eof())
				{
					std::cerr << "Error: " << option << " needs [float] value!" << std::endl;
					showSwitches(argc,argv);
					return -1;
				}
				break;
			case HELP_S:
				showSwitches(argc,argv);
				return 0;
				break;
		}
	}
	
	//Error Checking:
	if(deltaX!=0 && fX && lX)
	{
		std::cerr << "Error: Cannot use --deltaX, --firstX, and --lastX together!" << std::endl;
		showSwitches(argc, argv);
		return -1;
	}
	
	if(deltaX==0 && noX && !fX && !lX)
	{
		std::cerr << "Error: --noX requires either --deltaX or --firstX or --lastX!" << std::endl;
		showSwitches(argc, argv);
		return -1;
	}
	
	if(fX && lX && firstX>=lastX)
	{
		std::cerr << "Error: --firstX must be less than lastX!" << std::endl;
		showSwitches(argc, argv);
		return -1;
	}
	
	//Collect Values From stdin:
	std::vector<twoVector> v;
	if(noX)
	{
		for(double value;std::cin >> value;)
		{
			twoVector newV;
			if(v.size()==0)
				newV.x()=firstX;
			else
				newV.x()=v.back().x()+deltaX;
			newV.y()=value;
			v.push_back(newV);
		}
		
		if(v.size()>0)
		{
			if(deltaX==0)
			{
				deltaX=(lastX-firstX)/static_cast<double>(v.size()-1);
				double currentX=firstX;
				for(auto& thisV:v)
				{
					thisV.x()=currentX;
					currentX+=deltaX;
				}
			}
			else
			{	
				if(!lX)
				{
					lastX=v.back().x();
				}
				else if(!fX)
				{
					v.back().x()=lastX;
					for(int i=v.size()-2;i>=0;--i)
						v[i].x()=v[i+1].x()-deltaX;
					firstX=v.front().x();
				}
			}
		}
	}
	else
	{
		int iter=0;
		for(double value;std::cin >> value;)
		{
			if(iter%2==0)
			{
				twoVector newV;
				newV.x()=value;
				newV.y()=0;
				if(
					(fX && newV.x()>=firstX && !lX) || 
					(lX && newV.x()<=lastX && !fX) || 
					(!lX && !fX) ||
					(fX && newV.x()>=firstX && lX && newV.x()<=lastX)  
				) 
					v.push_back(newV);
			}
			if(iter%2==1 && v.size()>0)
				v.back().y()=value;
			iter++;
		}
		if(iter%2!=0)
		{
			std::cerr << "Wrong offset or missing column!" << std::endl;
			return -1;
		}
		if(v.size()>0)
		{
			std::sort(v.begin(), v.end(), sortByX);
			firstX=v.front().x();
			lastX=v.back().x();
		}
	}
	
	//std::vector<twoVector> vSum;
	twoVector vSum;
	vSum.x()=0;
	vSum.y()=0;
	for(int i=0;i<v.size();i++)
	{
		vSum.x()=v[i].x();
		vSum.y()+=(v[i].y()+v[i-1].y())*0.5*(v[i].x()-v[i-1].x());
		//vSum.push_back(vS);
		if(sumOnlyPrefix.length()==0)
			std::cout << vSum.x() << ' ' << vSum.y()+constant << std::endl;
	}
	if(sumOnlyPrefix.length()>0)
		std::cout << sumOnlyPrefix << ' ' << vSum.y()+constant << std::endl;
	
	return 0;
}
