/*! \file user_options.h
\brief Header file for UserOptions class, with source file UserOptions.cpp.
*/

// Avoid multiple includes
#ifndef user_options_h
#define user_options_h

// Boost includes related to program_options  
#include <boost/program_options.hpp>	
#include <boost/program_options/positional_options.hpp>
#include <boost/program_options/parsers.hpp>

// Boost includes related to custom functionality of this class
#include <boost/any.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

// Boost exception handling
#include <boost/exception/all.hpp>

// Standard library includes
#include <stdio.h>						
#include <iostream>		
#include <string>		
#include <vector>
#include <typeinfo>

// Definitions and namespaces uses
namespace opt = boost::program_options;	// short-hand for namespace
using std::string;

 // UserOptions is part of the SlaughterCommon namespace; this is documented in File_parts.h
namespace SlaughterCommon{

/*! \class UserOptions user_options.h "common/user_options/user_options.h"
 * \brief A class for handling command line specified options.
 * 
 * This class is basically a wrapper for the Boost.Program_options 
 * library. As such, the syntax for this class was designed to be somewhat
 * similar to the Boost libraries. 
 * 
 * The example program detailed in \ref test_user_options has extensive
 * documentation for this class including. The resulting executable
 * allows the user to test the various features.
 * 
 * The functionality of the Boost.Program_options is not restrained in 
 * anyway. It is possible to work directly with the Boost classes and 
 * this class. The members \c opt_list and \c opt_map may be used as
 * shown in the Boost documentation they represent the
 * \c boost::program_options::options_description and 
 * \c boost::program_options::variables_map classes, respectively.
 * 
 * \n<b>Special Behavior: </b>\n
 * The option "config", if added, is automatically considered to contain
 * the name of a configuration file that contains program options. This
 * options must use the string type. The options in the coinfiguration 
 * file are always overwritten by those on the command line. The 
 * following code demonstrates this behavior.
 * 
 * \code
 * UserOptions opt("General Options");
 * opt.add_option<string>("config", "default.cfg", "Configuration file name");
 * \endcode
 * 
 * The options "help", if added with the \c add_flag function, will
 * cause the program to automatically exit after displaying all of the
 * available user options. See the documentation for \c add_flag for
 * more details.
 * 
 * Grouping is accomplished using the \c add member function. Example
 * code for creating groups are provided in the class documentation. It
 * is possible to hide groups, just set the \c hidden attribute to \c
 * false. Note, the master class (the one to which all the \c add
 * commands are linked) will always display its options, only classes
 * that are attached can be hidden.
 */
class UserOptions{
    public:  
		//! A flag for hidding the options associated with this class.
		bool hidden;
    
		//! A method for showing all options
		/*!
		 * This function exists to allow for a \c - -advanced type 
		 * flag that will display all the options. This must be called
		 * after the apply_option member is called.
		 * 
		 * \see v2m.cpp
		 */ 
		void show_hidden();
    
		//! Boost::program_options::options description class associated with this UserOptions instance.
        opt::options_description opt_list; // contains everything

        //! Boost::program_options::variables_map class associated with this UserOptions instance.
		opt::variables_map opt_map;
			
		//! Class constructor
		/*! This constrcutor requires a single argument, the message
		 * that will be displayed when the command-line options are
		 * displayed on the screen.
		 * \param message A std::string containing the desired message.
		 * 
		 * <b>&nbsp; Example Syntax:</b>
		 * \code{.cpp}
		 * UserOptions opt("Available options");
		 * \endcode
		 */
        UserOptions(string message);
        
		//! Function for adding a new command-line option, without a default.
		/*! \tparam Type The type of variable desired, e.g., \b double 
		 *
		 * \param handle The text reference to the added option, uses
		 * the same syntax as Boost
		 * \param message A description of the option
		 * 
		 * <b>&nbsp; Example Syntax:</b>
		 * \code
		 * UserOptions opt("Available options");
		 * opt.add_option<double>("angle,a", 35, "Interface angle");
		 * \endcode
		 */
		template<typename Type> void add_option(string handle, string message){  

			// Add options to complete set of program options
			opt_list.add_options()(handle.c_str(), opt::value<Type>()->multitoken(), message.c_str());
		
			// Add the option to the visable list
			vis_list.add_options()(handle.c_str(), opt::value<Type>()->multitoken(), message.c_str());
		}
		
		//! Function for adding a new command-line option, with a default
		/*! \tparam Type The type of variable desired, e.g., \b double 
		 *
		 * \param handle The text reference to the added option, uses
		 * the same syntax as Boost
		 * \param dvalue The default value of the same Type
		 * \param message A description of the option
		 * 
		 * <b>&nbsp; Example Syntax:</b>
		 * \code{.cpp}
		 * UserOptions opt("Available options");
		 * opt.add_option<double>("angle,a", 35, "Interface angle");
		 * \endcode
		 */
        template<typename Type> void add_option(string handle, Type dvalue, string message){
			
			// Add options to complete set of program options
			opt_list.add_options()(handle.c_str(), opt::value<Type>()->default_value(dvalue)->multitoken(), message.c_str());
		
			// Add the option to the visable list
			vis_list.add_options()(handle.c_str(), opt::value<Type>()->default_value(dvalue)->multitoken(), message.c_str());	
		}
		
		//! Function for adding a new command-line option, with a default and default text
		/*! \tparam Type The type of variable desired, e.g., \b double 
		 *
		 * \param handle The text reference to the added option, uses
		 * the same syntax as Boost
		 * \param dvalue The default value of the same Type
		 * \param message A description of the option
		 * \param dtext A string containing the default text to display; useful when ouputing a vector
		 * 
		 * <b>&nbsp; Example Syntax:</b>
		 * \code{.cpp}
		 * UserOptions opt("Available options");
		 * std::vector vec;
		 * vec.push_back(35);
		 * vec.push_back(50);
		 * opt.add_option<std::<double> >("angle,a", vec, "Interface angle", "[35,50]");
		 * \endcode
		 */
		template<typename Type> void add_option(string handle, Type dvalue, string message, string dtext){
		
			// Add options to complete set of program options
			opt_list.add_options()(handle.c_str(), opt::value<Type>()->default_value(dvalue, dtext.c_str())->multitoken(), message.c_str());
		
			// Add the option to the visable list
			vis_list.add_options()(handle.c_str(), opt::value<Type>()->default_value(dvalue, dtext.c_str())->multitoken(), message.c_str());
		
		}
		
		//! Function for adding a new command-line option, without a default but with a positional argument
		/*! \tparam Type The type of variable desired, e.g., \b double 
		 *
		 * \param handle The text reference to the added option, uses
		 * the same syntax as Boost
		 * \param message A description of the option
		 * \param pos Positional specification, refer to 
		 * Boost.Program_options documentation for more information.
		 * 
		 * <b>&nbsp; Example Syntax:</b>
		 * \code{.cpp}
		 * UserOptions opt("Available options");
		 * opt.add_option<string>("input,i", "Input file name", 1);
		 * \endcode
		 */
		template<typename Type> void add_option(string handle, string message, int pos){  
			
			// Add options to complete set of program options
			opt_list.add_options()(handle.c_str(), opt::value<Type>(), message.c_str());
			
			// Add the option to the visable list
			vis_list.add_options()(handle.c_str(), opt::value<Type>(), message.c_str());
			
			// Append the position vectors with the new entry
			add_positional(handle, pos);
		}
			
		//! Function for adding a new command-line option, with a default and with a positional argument
		/*! \tparam Type The type of variable desired, e.g., \b double 
		 *
		 * \param handle The text reference to the added option, uses
		 * the same syntax as Boost
		 * \param dvalue The default value of the same Type
		 * \param message A description of the option
		 * \param pos Positional specification, refer to 
		 * Boost.Program_options documentation for more information.
		 * 
		 * <b>&nbsp; Example Syntax:</b>
		 * \code
		 * UserOptions opt("Available options");
		 * opt.add_option<string>("input,i", "default.txt", "Input file name", 1);
		 * \endcode
		 */
        template<typename Type> void add_option(string handle, Type dvalue, string message, int pos){  
			
			// Add options to complete set of program options
			opt_list.add_options()(handle.c_str(), opt::value<Type>()->default_value(dvalue), message.c_str());
			
			// Add the option to the visable list
			vis_list.add_options()(handle.c_str(), opt::value<Type>()->default_value(dvalue), message.c_str());

			// Append the position vectors with the new entry
			add_positional(handle, pos);
		}
		
		//! Tests for the existance of an option
		/*! Allows the user to test if an option contains a value,
		 * which may not be the case for an option that is specified
		 * without a default value.
		 */ 
		bool exist(string handle);

		//! Function for returning a command-line specified option
		/*! This allows the users to extract the default value or the
		 * user specified value from the command line. Note, that the
		 * \c apply_options member function must be called prior to using
		 * this function. 
		 * \tparam Type The type of variable, e.g., \b double 
		 * 
		 * \param handle The text reference to the options, do not include
		 * the short-hand value.
		 * 
		 * <b>&nbsp; Example Syntax:</b>
		 * \code{.cpp}
		 * UserOptions opt("Available options");
		 * opt.add_option<string>("input,i", "default.txt", "Input file name", 1);
		 * opt.apply_options(argc, argv);
		 * std::string filename = opt.get<string>("input");
		 * \endcode
		 */
        template<typename Type> Type get(string handle){  
			if(opt_map.count(handle.c_str())){
				return boost::any_cast<Type>(opt_map[handle.c_str()].as<Type>());
			} else {
				printf("ERROR: The option, --%s, was not found, it either does not exist or a value has not been assigned.\n", handle.c_str());
				exit(1);
			}
		}

		//! Function for adding a flag (no associated value) options
		/*! \param flag Text to use for accessing help, e.g., "help,h"
		 * \param message associated with the flag option
		 * 
		 * <b>&nbsp; Example Syntax:</b>
		 * \code{.cpp}
		 * UserOptions opt("Available options");
		 * opt.add_help("help,h", "Show all available options.");
		 * \endcode
		 * 
		 * A flag does not contain any other input, just the command 
		 * itself (e.g., - -help or - -verbose).
		 * 
		 * A true|false value is returned from the \c get_flag() function
		 * if the options is specified.
		 * 
		 * Note, "help" is a special case. When this options is given
		 * it automatically lists the various command-line options
		 * and stops execution of the program.
		 */
        void add_flag(string flag, string message);
        
        //! Function for testing if an flag was supplied on the command-line
		/*! \param flag Text to use for accessing flag value, e.g., "help"
         */ 
        bool get_flag(string flag);
        
        //! Function for collecting and applying command-line inputs
        /*! This function should be the last to be called, as shown in
         * the example.
         * \param argc Number of command-line inputs, from \c main
         * \param argv Character array storing the inputs, from \c main
         * 
		 * <b>&nbsp; Example Syntax:</b>
		 * \code{.cpp}
		 * UserOptions opt("Available options");
		 * opt.add_option<string>("input,i", "default.txt", "Input file name", 1);
		 * opt.apply_options(argc, argv);
		 * std::string filename = opt.get<string>("input");
		 * \endcode
		 */
        void apply_options(int argc, char* argv[]);
        
		//! A function for grouping UserOptions together
		/*! It is possible to seperate user options into groups, as 
		 * shown in the Boost.Program_options documentation. This can
		 * also be done with the UserOptions class in a similar fashion
		 * by creating several instances of the UserOptions class.
		 * 
		 * \param new_opt The instance of a UserOptions that should be added
		 * to the calling class, the class inputed here is considered then
		 * to be a slave class to the calling class. Slave classes have 
		 * the apply_options member disabled.
		 * 
		 * \code{.cpp}
		 * UserOptions opt("General Options");
		 * opt.add_help("help,h", "Display the available options");
		 * opt.add_option<double>("angle,a", 35, "Interface angle");
		 * 
		 * UserOptions opt2("Advanced Options");
		 * opt2.add_options<int>("advanced",0, "Use advanced options");
		 * opt.add(opt2);
		 * opt.apply_options(argc, argv);
		 * \endcode
		 */
        UserOptions& add(UserOptions & new_opt);
        
        //! Allows user to add title text that prints with - -help
        /*! \param str A string contiaing the text to display.
         * 
		 * \code{.cpp}
         * UserOptions opt("Available Options");
         * opt.add_title("\nThis is my wonderful program.\n\nCopyright 2012\n");
         * opt.add_flag("help","Display the available options");
         * opt.add_option("input","The input file name");
         * opt.apply_options(argc, argv);
         * \endcode
         */
        void add_title(string str);
            
    private:
    	//! Options description class containing visible items only
        opt::options_description vis_list;
    
		//! Vector for storing positional option names
		std::vector<string> pos_str;
		
		//! Vector for storing positional option item counts
		std::vector<int> pos_int;
		
		//! Storage location for title string
		string title;
		
		//! Flag indicating if class instance is a slave to another
		bool is_slave;
		
		//! Function for handling slave classes
		void slave_test(const char* func_name);

		//! Function for adding a positional option
		void add_positional(std::string handle, int pos);
		
}; // UserOptions

} // Common namespace

#endif
