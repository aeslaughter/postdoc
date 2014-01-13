//! \file file_parts.cpp Source code for FileParts class.

// Standard library includes
#include <iostream>     // input/output class (printf)
#include <sstream>		// allows for stream string creation
#include <stdio.h>		// input/output to string (sprintf) 
#include <fstream>      // file reading class (fopen)

// FileParts include statements
#include "common/file_parts/file_parts.h"  
using namespace SlaughterCommon;

// Overloaded version for handling char* input
FileParts :: FileParts(const char* c_str){
    std::string str(c_str);
    init(str);
}

// FileParts Class constructor
FileParts :: FileParts(std::string str){
    init(str);
}

// String version of assign	
void FileParts :: assign(std::string str){
	init(str);
}
	
// Char* version of assign
void FileParts :: assign(const char* c_str){
	std::string str(c_str);
    init(str);
}

// Return a string with a time step
std::string FileParts :: add_tstep(int tstep, int pad, const char* prfx){
	
	// Build a string stream containing the number
	std::stringstream tss;
	tss.width(pad);
	tss.fill('0');
	tss << tstep;

	// Build the string for outpout
	std::string str;
	str.append(path);
	str.append(name);	
	str.append(prfx);
	str.append(tss.str());
	str.append(ext);
	
	// Return the string
	return str;
}

// Intialization function
int FileParts :: init(std::string str){
	
    // Define the complete file path
    full = str;

    // Define variables
    size_t dot, slash;

    // Locate the final slash
    slash = str.rfind("/");
    
    // Locate the last period
    dot = str.find(".", slash);
    
    // Extract the extension. The length parameter eliminates a problem with finding a value for "found" outside of the string.
    if (dot > 0 && dot < str.length()){	// an extensions exists
        ext = str.substr(dot);
    }

    // Extract the path and the file name witout the extension. 
    if (slash > 0 && slash < str.length()){	
        name = str.substr(slash+1, dot - slash - 1);
        path = str.substr(0, slash+1);
    }

    // Check if the file exists
    FILE* fp = fopen(str.c_str(), "r");
    if (fp){
        exist = true;
        fclose(fp);
    } else {
        exist = false;
    }
  
    return 0;
}

// A function for display the file parts
void FileParts :: display(){
    // Print all of the parts to the screen
    printf("ext = %s\n", ext.c_str());
    printf("name = %s\n", name.c_str());
    printf("path = %s\n", path.c_str()); 
    printf("full = %s\n", full.c_str());
    printf("exist = %s\n", (exist)?"true":"false"); 
}

// A function for updating the full file path
void FileParts :: update(){
	full.assign(path);
	full.append(name);
	full.append(ext);	
}
