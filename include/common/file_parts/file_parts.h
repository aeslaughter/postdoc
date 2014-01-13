//! \file file_parts.h Header file for FileParts class source code

// Avoid multiple includes
#ifndef file_parts_h
#define file_parts_h

// Required includes
#include <string>	// the stardard library string

// Add to my common namespace  
namespace SlaughterCommon{

/*! \class FileParts file_parts.h "common/file_parts/file_parts.h"
 * \brief A class for simple file name handling.
 * 
 * On creation it seperates the file into components as well as tests 
 * for it existance.
 * 
 * Example usage:
 * \code
 * FileParts filename(/my/path/and/file.txt);
 * FILE* fid = fopen(filename.full,'r');
 * 	 ... gather some data ...
 * fclose(fid);
 * \endcode
 * 
 */
class FileParts{
public:
    std::string full;   //!< The complete file name and path as input
    std::string path;   //!< The directory of the file name
    std::string name;   //!< The file name without the path or extension
    std::string ext;    //!< The file name extension, including the period  
    bool exist;         //!< A boolean flag indicating if the file exists
    
    //! Default constructor
    FileParts(){};
    
    //! Cconstructor, string input version
    /*! Creates a FileParts object using a std::string input
     * \param str A standard string containing the complete file name and path
     */
	FileParts(std::string str);
    
    //! Constructor, char input version
    /*! Creates a FileParts object with a const char* input
     * \param c_str A char* containing the complete file name and path
     */ 
    FileParts(const char* c_str);    
    
    //! A function that displays the various parts of the file
    void display();
    
    //! A function for updating the full file path
    /*! This allows for the user to alter the components and then
     * create a full file path from these new complnents. For example:
     * \code
     * FileParts filename("path/to/a/file.txt");
     * filename.name.append("2");
     * filename.update();
	 * \endcode
     */ 
    void update();
    
	//! Assign member allows user to initilize class after decleration
	/*! In some instances the class must be declared before the 
     * complete file path is known. As such, it should be possible to
     * do the following, which this operator enables.
     * \param str A standard string containing the complete file name and path
     * \code
     * FileParts filename;
     * filname.assign("path/to/a/file.txt");
     * \endcode
     */
	void assign(std::string str);
	
	//! Assign member using char* as input.
	/*! \param c_str A char* containing the complete file name and path
	 */
	void assign(const char* c_str); 
	
	//! Special function for inserting a time series stamp
	/*! \param tstep An interger value containing the time step
	 *  \param pad An interger value containing the zero padding level,
	 * the default is 3.
	 *  \param prfx A \c char that is inserted between the filename
	 * and the numeric time step. The default is an empty value.
	 * 
	 * The following code retuns the following strings with a time stamp 
	 * inserted, it does not edit the object itself: \n
	 * \c path/to/a/file_0222.txt 
	 * 
	 * \code
	 * FileParts filename;
	 * filename.assign("path/to/a/file.txt");
	 * std::string str = filename.add_tstep(222, 4, "_");
	 * printf("%s\n", str.c_str());
	 * \endcode
	 */ 
	std::string add_tstep(int tstep, int pad = 3, const char* prfx = "");
	
private:
    //! A function that initializes the public attributes
    /*! 
     * \param str A standard string containing the complete file name and path
     */
    int init(std::string str);
}; // FileParts class
}  // slaughter namespace
#endif
