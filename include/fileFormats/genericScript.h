//This is a more flexible version of scriptFormat.h
//It just takes a variable list.
//See dataTypes.h for variable list and anonymous data types.
#include "../algorithms/dataTypes.h"

#ifndef genericScript
#define genericScript

#define END_INPUT "end"

//this extracts elements from a file
template <typename T, typename U>
class genericScript {
	public:
		//blob is special class
		genericScript(const char *name, std::ios_base::openmode mode, U *blob);
		//these are just the current requirements
		~genericScript();
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
genericScript<T,U>::genericScript(const char *name, std::ios_base::openmode mode, U *blob)
{
	b=blob;
	this->open(name,mode);
}

template <typename T, typename U>
genericScript<T,U>::~genericScript()
{
	this->close();
}

template <typename T, typename U>
bool genericScript<T,U>::open(const char *name, std::ios_base::openmode mode)
{
	if(file.is_open())
		file.close();//screw you! I'm going to close this when you didn't expect it!
		
	std::stringstream filename;
	filename << name << ".mpd";
	file.open(filename.str().c_str(),mode);
	if(!file.is_open())
	{
		std::cout << "Could not open " << filename.str() << " in genericScript class!\n";
		throw 0;
	}
	return true;
}

template <typename T, typename U>
void genericScript<T,U>::close()
{
	if(file.is_open())
		file.close();
}

template <typename T, typename U>
void genericScript<T,U>::write()
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
void genericScript<T,U>::read()
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

#endif