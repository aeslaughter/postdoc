/** \example test_file_parts.cpp 
 * A test function for the FileParts class.
 * This is the source code for a simple test function that tests the 
 * functionality of the FileParts class. To run this example simply
 * execute the \c test_file_parts executable. It will dispay the the 
 * components of the test file and path for both types of input to the 
 * class: \c std::string and \c char*.
 */
 
// General includes
#include <string>                   // use the standard string class
#include <stdio.h>                  // allows use of printf for output

// FileParts related includes and namespace
#include "common/file_parts/file_parts.h"  // header for FileParts class
using namespace SlaughterCommon;

//! Main program for testing the FileParts class functionality
/** This function is the main program executed to test that the FileParts
 * class is operating as expected. It tests the two input variants of the 
 * class.
 */ 
int main(){
    
    // This tests the std::string input
    std::string str("./this/is/a/test/dir/and/file/test.tar.gz");
    FileParts str_file(str);
  
    // Test the char* input
    FileParts char_file(str.c_str());
    
    // Output the results from the string input
    printf("\n\n");
    printf("FileParts with a std::string input:\n");
    str_file.display();
    printf("\n\n");
    
    // Output the results from the char* input
    printf("FileParts with a char* input:\n");
    char_file.display();
    printf("\n\n");
    
    // Alter the filename
    char_file.name.append("(test)");
    char_file.update();
    char_file.display();
    printf("\n\n");
    
    // Test = operator
    FileParts eq_file;
    eq_file.assign("a/path/to/file/assigned/with/equal/sign.txt");
    eq_file.display();
    
    // Test the time stamp behavior
    FileParts time_file;
    time_file.assign("this/is/a/time/series/filename.dat");
    std::string newfile = time_file.add_tstep(234, 5, "_");
    printf("\n\n%s\n", newfile.c_str());
    
    newfile.assign(time_file.add_tstep(234, 5));
    printf("%s\n", newfile.c_str());
    
    newfile.assign(time_file.add_tstep(234));
    printf("%s\n\n", newfile.c_str());
    
    
    return 0;
}
