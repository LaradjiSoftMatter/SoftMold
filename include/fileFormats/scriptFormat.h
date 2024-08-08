//This is a template for a flexible format.
//It doesn't require a specific number of parameters, types of parameters, or anything.
//Just give it the blob and it will extract elements from the file until the blob decides otherwise.
#include "../algorithms/functions.h"

#define END_INPUT "end"

//this extracts elements from a file
template <typename T, typename U>
class Script {
	public:
		//blob is special class
		Script(const char *name, std::ios_base::openmode mode, U *blob);
		//these are just the current requirements
		~Script();
		bool open(const char *name, std::ios_base::openmode mode);
		void close();
		void write();//write to the file
		void read();//read from the file
	private:
		U *b;//the blob, a group of elements with special member functions
		std::fstream file;//the file
		std::stringstream current;//the current line
};

template <typename T, typename U>
Script<T,U>::Script(const char *name, std::ios_base::openmode mode, U *blob)
{
	b=blob;
	this->open(name,mode);
}

template <typename T, typename U>
Script<T,U>::~Script()
{
	this->close();
}

template <typename T, typename U>
bool Script<T,U>::open(const char *name, std::ios_base::openmode mode)
{
	if(file.is_open())
		file.close();//screw you! I'm going to close this when you didn't expect it!
		
	std::stringstream filename;
	filename << name << ".mpd";
	file.open(filename.str().c_str(),mode);
	if(!file.is_open())
	{
		std::cout << "Could not open " << filename.str() << " in Script class!\n";
		throw 0;
	}
	return true;
}

template <typename T, typename U>
void Script<T,U>::close()
{
	if(file.is_open())
		file.close();
}

template <typename T, typename U>
void Script<T,U>::write()
{
	std::string out("");
	(*b).reset();
	while(out!=END_INPUT)
	{
		out=(*b).output();
		//interestingly, if you want your script to be readable by a human being
		//, you should include newlines in your output
		if(out!=END_INPUT)
			file << out << ' ';
	}
	file << '\n';
}

template <typename T, typename U>
void Script<T,U>::read()
{
	std::string in("");
	file.seekg(std::ios::beg);
	file.clear();
	(*b).reset();
	while(!file.eof() && in!=END_INPUT)
	{
		file >> in;
		if(!file.eof())
			(*b).input(in);
	}
	//this actually checks for failure, otherwise, you would just guess
	(*b).input(END_INPUT);
}
