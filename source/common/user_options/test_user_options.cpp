/** \example test_user_options.cpp 
 * \anchor test_user_options
 * A test function for the UserOptions class.
 * 
 * see \ref test_user_options for more details
 */
 
//![include]

// Standard include
#include <vector>

// Include my common library, which contains the UserOptions class
#include "common/include.h"
using namespace SlaughterCommon;

// The main function for test_user_options.cpp
int main(int argc, char* argv[]){
//![include]

//![main]
// Define the master class for storing and accessing options
UserOptions main("General Options");

// Add an extra title to the help output
main.add_title("\nThis is a test program, the options that are accessible\nfrom the command line are listed below.\n\nFor example:\n./test_user_options --value=5\n\n");

// Add a help option
main.add_flag("help,h", "If you need help");

// Add a verbose options that will all the results when used
main.add_flag("verbose", "Display's the values of all the variables");

// Assign a configuration file, including a default file
main.add_option<string>("config", "../data/test_user_options/test_user_options.cfg", "Specify a configuration file");

// Some value that you may want to used
main.add_option<double>("value,v", 0, "The value of some importart variable");
//![main]

//![options]
// Define another instance of UserOptions that includes io info
UserOptions io("Input/Ouput Options");
io.add_option<string>("input,i", "input.txt", "The input filename", 1);
io.add_option<string>("output,o", "default.txt", "The output file for widget document", 1);
io.add_option<std::vector<string> >("file-list","A list of files", 3);

// Define some advanced options (the configuration file sets some of these)
UserOptions adv("Advanced Options");
main.add_option<double>("advanced-value", 0, "The value of some really importart variable");
adv.add_option<string>("input.path","The path to the input file");
adv.add_option<string>("input.name","The name of the input file, w/o the extension");
adv.add_option<string>("input.ext","The extension to the input file");
adv.add_option<std::vector<int> >("many,m", "Many of these are allowed");
adv.add_option<std::vector<double> >("multi","This input can contain multiple values");

std::vector<int> def;
def.push_back(5);
def.push_back(6);
adv.add_option<std::vector<int> >("many-multi", def, "You can combine behavior and list default(s)", "[5,6]");
//![options]

//![hidden]
// Add some hidden options
UserOptions hide("Hidden Options");
hide.add_option<int>("big-red-button", 0, "Don't change this this value!");
hide.add_flag("show-hidden","If you know about it then you know what it does");
hide.hidden = true; // this hides the hide instance
//![hidden]

//![apply]
// Group the various classes and apply the command line inputs
main.add(io).add(adv).add(hide);
main.apply_options(argc, argv);
//![apply]

//![use]
// If --show-hidden is used display the available hidden options
if (main.get_flag("show-hidden")){
	std::cout << hide.opt_list << std::endl;
	return 0;
}
// If --verbose is used, show all the values
bool display = main.get_flag("verbose");
if (display){
	// Show the main options
	printf("config: %s\n", (main.get<string>("config")).c_str());
	printf("value: %f\n", main.get<double>("value"));
	
	// Show the input/output options
	printf("input: %s\n", (main.get<string>("input")).c_str());
	printf("output: %s\n", (main.get<string>("output")).c_str());
	
	if(main.exist("file-list")){
		std::vector<string> vec0 = main.get<std::vector<string> >("file-list");
		for(int i = 0; i < vec0.size(); i++){
			printf("file-list[%d] = %s\n", i, vec0[i].c_str());
		}	
	}
	
	// Show the advanced options
	printf("advanced-value: %f\n", main.get<double>("advanced-value"));
	printf("input.path: %s\n", (main.get<string>("input.path")).c_str());
	printf("input.name: %s\n", (main.get<string>("input.name")).c_str());
	printf("input.ext: %s\n", (main.get<string>("input.ext")).c_str());
	printf("bit-red-button: %d\n", main.get<int>("big-red-button"));
	
	if (main.exist("many")){
		std::vector<int> vec1 = main.get<std::vector<int> >("many");
		for(int i = 0; i < vec1.size(); i++){
			printf("many[%d] = %d\n", i, vec1[i]);
		}
	}
	
	if (main.exist("multi")){
		std::vector<double> vec2 = main.get<std::vector<double> >("multi");
		for(int i = 0; i < vec2.size(); i++){
			printf("multi[%d] = %f\n", i, vec2[i]);
		}
	}
	
	if (main.exist("many-multi")){
		std::vector<int> vec3 = main.get<std::vector<int> >("many-multi");
		for(int i = 0; i < vec3.size(); i++){
			printf("many-multi[%d] = %d\n", i, vec3[i]);
		}	
	}
	
	
} else {
	printf("The program worked great!\n");
}

return 0;
}
//![use]
