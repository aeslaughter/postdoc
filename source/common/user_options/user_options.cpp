//! \file user_options.cpp Source code for the UserOptions class

// Use my common library and namespace
#include "common/include.h"
using namespace SlaughterCommon;

// Class constructor
UserOptions :: UserOptions(string message) : opt_list(message), vis_list(message){    
	is_slave = false;
	hidden = false;
}

// An option for showing hidden options
void UserOptions :: show_hidden(){
	std::cout << title;
	std::cout << opt_list << "\n";
	exit(0);
}

// Adds an option without an associated value (e.g., --help)
void UserOptions :: add_flag(string flag, string message){
	
	// Adds options to the complete, public list
    opt_list.add_options()(flag.c_str(), message.c_str());
    
    // Adds options to the visable list
    vis_list.add_options()(flag.c_str(), message.c_str());
}
 
// Implements the command-line arguments 
void UserOptions :: apply_options(int argc, char* argv[]){

	// Restricts this function if the class is a slave of another
	slave_test("apply_options(...)");

	// Create positional descriptions from the storage vectors
	opt::positional_options_description pos_desc;
	for (unsigned int i = 0; i < pos_str.size(); ++i){
		pos_desc.add(pos_str[i].c_str(), pos_int[i]);
	}

	// Parse the command line and store the options
	opt::store(opt::command_line_parser(argc, argv).options(opt_list).positional(pos_desc).run(), opt_map);

	// Parse the configuration file and store the options
	if (opt_map.count("config")){
		string config = get<string>("config");
		opt::store(opt::parse_config_file<char>(config.c_str(), opt_list, false), opt_map);
	}

	// Prepare the options for use
	opt::notify(opt_map);
	   
	// Print the options information if --help is used   
    if (opt_map.count("help")) {
		std::cout << title;
		std::cout << vis_list << "\n";
		exit(0);
	}	
}    

// Adds options from other instances of the UserOptions class
UserOptions & UserOptions :: add(UserOptions & new_opt){
	
	// Define the added class as a slave
	new_opt.is_slave = true;
	
	// Insert the position options into this class
	pos_str.insert(pos_str.end(), new_opt.pos_str.begin(), new_opt.pos_str.end());
	pos_int.insert(pos_int.end(), new_opt.pos_int.begin(), new_opt.pos_int.end());
	
	// Add the class to the complete list
	opt_list.add(new_opt.opt_list);	
	
	// If the list is visable add it, otherwise it is skipped
	if (!new_opt.hidden){
		vis_list.add(new_opt.vis_list);	
	}

	// Return the class itself
	return *this;
}

// Adds positional information to the storage vectors
void UserOptions :: add_positional(string handle, int pos){
	
	// Remove the ",x" portion of the handle input
	size_t idx = handle.rfind(",");
	if (idx > 0 && idx < handle.size()){
		handle.erase(idx, handle.size());
	}
				
	// Append the vectors			
	pos_str.push_back(handle);
	pos_int.push_back(pos);
}

// Checks for existance of an option
bool UserOptions :: exist(string handle){
	return opt_map.count(handle.c_str());	
}

// Retuns a true falue if was present in the command line
bool UserOptions :: get_flag(string flag){
	
    if (opt_map.count(flag)) {
		return true;
	} else {
		return false;
	}	
}

// Tests if the current class is a slave, it produces an error if it is
void UserOptions :: slave_test(const char* func_name){
	if (is_slave){
		printf("ERROR: member %s is not available.\nThis instance of UserOptions was included in an call to add() from another instance, thus it has been idenfied as a slave class. The %s member is not available for slave instances.\n", func_name, func_name);		
		exit(2);
	}
}

// Adds header text to be printed prior to the list with --help
void UserOptions :: add_title(string str){
	slave_test("add_title(...)");
	title = str;
}

