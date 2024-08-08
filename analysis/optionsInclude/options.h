#ifndef OPTIONS_MPD
#define OPTIONS_MPD
#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>
#include <unordered_map>
#include <list>
#include <vector>
#include <algorithm>

namespace mpd {
/** \b A data type for reading and checking in options from command line.
 *  Important bits are emplace_back for adding options and iterators over options.
 *
**/
struct options {
	/** \b The data element of an option in the command line.
	 *  This is an element of options.
	**/
	struct option {
		/// The value this structure works with is std::string.
		using value_type= std::string;
		
		/// The reference to std::string.
		using reference = value_type&;
		
		/// The const reference for std::string.
		using const_reference= const value_type&;
		
		/// The iterator to the arguments vector.
		using iterator  = typename std::vector<value_type>::iterator;
		
		/// The const iterator to the arguments vector.
		using const_iterator  = typename std::vector<value_type>::const_iterator;
		
		//constructors
		/// Construct an option with dependencies (requirements).
		option(value_type longOpt, value_type shortOpt, value_type desc, 
			value_type req, size_t max):
			_longOption(longOpt),_shortOption(shortOpt), 
			_required(req), _description(desc), _maxArgs(max), 
			_visible(true), _eof(false), _active(0), _good(true),
			_bad(false), _fail(false)
		{}
		
		/// Construct an option without dependencies (requirements).
		option(value_type longOpt, value_type shortOpt, 
			value_type desc, size_t max):
			_longOption(longOpt),_shortOption(shortOpt), 
			_description(desc), _maxArgs(max), 
			_visible(true), _eof(false), _active(0), _good(true),
			_bad(false), _fail(false)
		{}
		
		/// Construct an option with visibility.
		option(value_type longOpt, value_type shortOpt, value_type desc, 
			value_type req, size_t max, bool visible):
			_longOption(longOpt),_shortOption(shortOpt), 
			_required(req), _description(desc), _maxArgs(max), 
			_visible(visible), _eof(false), _active(0), _good(true),
			_bad(false), _fail(false)
		{}
		
		//access data
		/// Description of the option, for use with --help option.
		value_type description() const {return _description;}
		
		/// Long option name.
		value_type longOption() const {return _longOption;}
		
		/// Short option name.
		value_type shortOption() const {return _shortOption;}
		
		/// Required (dependent) option, for use with dependency tree.
		value_type required() const {return _required;}
		
		//range access
		/// Iterator to the begining of arguments. Use with for(arg:args()).
		iterator begin(){return _args.begin();}
		
		/// Iterator to the begining of arguments. Use with for(arg:args()).
		const_iterator cbegin(){return _args.cbegin();}
		
		/// Iterator to the ending of arguments. Use with for(arg:args()).
		iterator end(){return _args.end();}
		
		/// Iterator to the ending of arguments. Use with for(arg:args()).
		const_iterator cend(){return _args.cend();}
		
		//modifiers
		/// Push a read argument onto option.
		void push_back(const_reference in){_args.push_back(in);}
		
		//capacity access
		/// Return the number of arguments read.
		size_t size() const {return _args.size();}
		
		/// Return the number of arguments read.
		size_t max_size() const {return _maxArgs;}
		
		//flags
		/// Check if option is visible.
		bool visible() const {return _visible;}
		
		/// Count active options.
		size_t active() const {return _active;}
		
		/// Check if option is at end of input.
		bool eof() const {return _eof;}
		
		/// Check if option is still good.
		bool good() const {return _good;}
		
		/// Check if option is bad (std::stringstream).
		bool bad() const {return _bad;}
		
		/// Check if option is failed (std::stringstream).
		bool fail() const {return _fail;}
		
		/// Set end of file flag.
		void set_eof() {_eof=true; _good=false;}
		
		/// Increment active flag.
		void set_active() {++_active;}
		
		//convert data
		/// Convert all arguments to some type.
		template <typename T>
		std::vector<T> args()
		{
			std::vector<T> out;
			for(auto &arg:_args)
			{
				std::stringstream convertArg(arg);
				T buf;
				convertArg >> buf;
				if(convertArg.bad())
				{
					std::cerr << "Warning: Bad conversion for argument \"" 
						<< arg << "\" of option " << _longOption 
						<< " or " << _shortOption << "!" << std::endl;
					_good=false;
					_bad=true;
				}
				if(convertArg.fail())
				{
					std::cerr << "Warning: Failed conversion for argument \"" 
						<< arg << "\" of option " << _longOption 
						<< " or " << _shortOption << "!" << std::endl;
					_good=false;
					_fail=true;
				}
				out.push_back(buf);
			}
			return out;
		}
		
		/// Convert one indexed argument to some type.
		template <typename T>
		T args(size_t i)
		{
			T out;
			std::stringstream convertArg(_args[i]);
			convertArg >> out;
			if(convertArg.bad())
			{
				std::cerr << "Warning: Bad conversion for argument \"" 
					<< _args[i] << "\" of option " << _longOption 
					<< " or " << _shortOption << "!" << std::endl;
				_good=false;
				_bad=true;
			}
			if(convertArg.fail())
			{
				std::cerr << "Warning: Failed conversion for argument \"" 
					<< _args[i] << "\" of option " << _longOption 
					<< " or " << _shortOption << "!" << std::endl;
				_good=false;
				_fail=true;
			}
			return out;
		}
		
		/// Return all arguments as strings.
		std::vector<value_type> sargs(){return _args;}
		
		/// Return one argument as a string.
		value_type sargs(size_t i){return _args[i];}
		
		private:
			std::vector<value_type> _args;
			value_type _longOption, _shortOption, _required, _description;
			size_t _maxArgs;
			bool _visible, _eof, _good, _fail, _bad;
			size_t _active;
	};
	
	/// String type is std::string.
	using str = std::string;
	
	/// This structure utilizes the internally defined option type.
	using value_type= option;
	
	/// The reference to the option type.
	using reference= option&;
	
	/// The iterator to a list of option.
	using iterator  = typename std::list<value_type>::iterator;
	
	/// The const iterator to a list of option.
	using const_iterator  = typename std::list<value_type>::const_iterator;
	
	/// Map of options to easily find a particular option.
	using optionMap = typename std::unordered_map<str,iterator>;
	
	//constructors
	/// Empty constructor. Use push_back to add elements.
	options(){};
	
	/// Constructor taking in c/c++ arguments.
	options(int argc, char ** argv, str description):_argc(argc)
	{
		std::stringstream cmdArgs;
		for(int i=0;i<argc;i++)
			cmdArgs << argv[i] << ' ';
		cmdArgs >> _callProgramName;
		str buf;
		while(cmdArgs >> std::quoted(buf))
		{
			_argv.push_back(buf);
			_argvUsed.push_back(false);
		}
		
		//default option for a globally required value
		this->emplace_back(_callProgramName, "", description, "", 0, false);
		this->emplace_back("--help", "-h", "Show this help option", "", 0);
	}
	
	/// Constructor with a list of options. Probably not needed.
	options(const std::list<value_type> &opts):_opts(opts){};
	
	//access data
	/// Locate an option by std::string (aka option::value_type).
	iterator find(option::value_type in)
	{
		auto out=_oMap.find(in);
		if(out==_oMap.end())
			return _opts.end();
		return out->second;
	}
	
	/// Locate an option by std::string (aka option::value_type).
	const_iterator find(option::value_type in) const
	{
		auto out=_oMap.find(in);
		if(out==_oMap.end())
			return _opts.cend();
		return out->second;
	}
	
	/// Retrieve option state by [] operator
	reference operator [] (option::value_type in) {return *_oMap[in];}
	
	/// Return std::string arguments from an option. Returns an empty vector if option doesn't exist.
	std::vector<option::value_type> get(option::value_type in)
	{
		auto buf=this->find(in);
		if(buf==_opts.end())
			return std::vector<option::value_type>();
		return buf->args<option::value_type>();
	}
	
	//range access
	/// An iterator to the beginning of the options, for use with for(option:options).
	iterator begin(){return _opts.begin();}
	
	/// A const iterator to the beginning of the options, for use with for(option:options).
	const_iterator cbegin(){return _opts.cbegin();}
	
	/// An iterator to the beginning of the options, for use with for(option:options).
	iterator end(){return _opts.end();}
	
	/// A const iterator to the beginning of the options, for use with for(option:options).
	const_iterator cend(){return _opts.cend();}
	
	//capacity
	/// Number of options.
	size_t size() const {return _opts.size();}
	
	/// Maximum number of options.
	size_t max_size() const {return _opts.max_size();}
	
	//modifiers
	/// Add an option. Not sure if needed.
	void push_back(const value_type &in) { _opts.push_back(in);}
	
	/// Construct an option.
	template <typename ...Args>
	void emplace_back(Args &&...args) 
	{ 
		_opts.emplace_back(args...);
		this->_process(--_opts.end());
	}
	
	/// Check input for correctness and contiguity.
	bool dependency_error()
	{
		for(auto &o:_opts)
		{
			if(o.active() && o.visible())
			{
				auto buf=_oMap.find(o.required());
				if(buf!=_oMap.end())
				{
					if(!buf->second->active() && buf->second->visible())
					{
						std::cerr << "Warning: Option " << o.longOption() << " or "
							<< o.shortOption() << " requires inactive option "
							<< buf->second->longOption();
						if(buf->second->shortOption()!="")
							std::cerr << " or " 
							<< buf->second->shortOption();
						std::cerr << "." << std::endl;
						return true;
					}
				}
			}
		}
		return false;
	}
	
	/// Check if option string hasn't been added.
	bool unrecognized_option()
	{
		bool test=false;
		for(auto arg:this->other())
		{
			std::cerr << "Warning: String \"" << arg << "\" not recognized." << std::endl;
			test=true;
		}
		return test;
	}
	
	/// Get unrecognized options.
	std::vector<str> other()
	{
		std::vector<str> out;
		for(auto buf=_argvUsed.begin();buf!=_argvUsed.end();++buf)
			if(!*buf) out.push_back(_argv[buf-_argvUsed.begin()]);
		return out;
	}
	
	//output
	/// Show help.
	void help()
	{
		for(auto &o:_opts)
		{
			if(o.longOption()!=_callProgramName)
				std::cerr << '\t';
			std::cerr << o.longOption();
			for(int i=0;i<o.max_size();i++)
				std::cerr << " [" << i << "]";
			std::cerr << ": " << o.description() << std::endl;
		}
	}
	
	private:
		std::list<value_type> _opts;
		optionMap _oMap;
		int _argc;
		std::vector<str> _argv;
		std::vector<bool> _argvUsed;
		str _callProgramName;
		void _process(std::list<value_type>::iterator in)
		{
			bool checkOnce=false;
			if(in->longOption().size()>0) _oMap[in->longOption()]=in;
			
			for(auto buf=std::find(_argv.begin(),_argv.end(),in->longOption());
			buf!=_argv.end();buf=std::find(buf++,_argv.end(),in->longOption()))
			{
				_argvUsed[buf-_argv.begin()]=true;
				in->set_active();
				++buf;
				for(int i=0;i<in->max_size();++i)
					if(buf!=_argv.end()) 
					{
						_argvUsed[buf-_argv.begin()]=true;
						in->push_back(*buf);
						++buf;
					} else in->set_eof();
			}
			if(!in->good())
			{
				checkOnce=true;
				std::cerr << "Warning: ";
				if(in->eof())
					std::cerr << "End of line reached for option " << in->longOption();
				if(in->bad())
					std::cerr << "Bad input reached for option " << in->longOption();
				if(in->fail())
					std::cerr << "Failed to reach input for option " << in->longOption();
				std::cerr << "!" << std::endl;
			}
			
			if(in->shortOption().size()>0) _oMap[in->shortOption()]=in;
			
			for(auto buf=std::find(_argv.begin(),_argv.end(),in->shortOption());
			buf!=_argv.end();buf=std::find(buf++,_argv.end(),in->shortOption()))
			{
				_argvUsed[buf-_argv.begin()]=true;
				in->set_active();
				++buf;
				for(int i=0;i<in->max_size();++i)
					if(buf!=_argv.end()) 
					{
						_argvUsed[buf-_argv.begin()]=true;
						in->push_back(*buf);
						++buf;
					} else in->set_eof();
			}
			
			if(!in->good() && !checkOnce)
			{
				std::cerr << "Warning: ";
				if(in->eof())
					std::cerr << "End of line reached for option " << in->shortOption();
				if(in->bad())
					std::cerr << "Bad input reached for option " << in->shortOption();
				if(in->fail())
					std::cerr << "Failed to reach input for option " << in->shortOption();
				std::cerr << "!" << std::endl;
			}
			
		}
	};
}
#endif
