/*! \file v2m.cpp
 * Source code for a program, \c v2m, to create a mesh from 2D
 * images slices. The program is utilzes the Vol2mesh library.
 * The basis for the development comes from the CGAL example: Section 
 * 50.3.3 Domains From Segmented 3D Images. 
 * 
 * This file includes the main function and uses the Vol2mesh class, the
 * class handles most of the behavior a majority of this code is for 
 * handling the command-line inputs. 
 * 
 * This use of executable associated with this source code is detailed 
 * here: \ref v2m.
 */
 
 /*! \page v2m Image Mesh Generation: v2m 
 * \tableofcontents 
 * 
 * \section v2m_intro Introduction
 * The program, v2m, converts 3-D "binary" image files such as 
 * those gathered from a CT-scanner and converts them into a 3-D 
 * tetrahedral mesh for the use in finite element programs. The images 
 * must be comprised of 0's and 1's, where the value of 1 is the solid 
 * that is meshed.
 * 
 * This program relies on the CGAL (<a href="http://www.cgal.org/">www.cgal.org</a>) 
 * and VTK (<a href="http://www.vtk.org/">www.vtk.org</a>) libaries. 
 * The main meshing capabilities are demonstrated in the CGAL 
 * documentation, in particluar <a href="http://www.cgal.org/Manual/latest/doc_html/cgal_manual/Mesh_3/Chapter_main.html#Subsection_50.3.3"> Section 50.3.3, "Domains From Segmented 3D Images." </a>
 * 
 * \section v2m_syntax Syntax for Executable
 * The main executable, v2m, is controlled via command-line inputs,
 * which allow the user to specify the input and output files as well as
 * various meshing criteria parameters. The executable is in the \c exec
 * directory. To view the available parameters and the associated 
 * default values use the help command:
 * 
 * \code
 * v2m --help
 * \endcode
 *
 * It is also to specify a configuration file with the \c --config option
 * followed by the configuration filename. The configuration file should 
 * follow the format of the Boost Program Optoins library: http://www.boost.org/doc/libs/1_49_0/doc/html/program_options.html .
 * 
 * \subsection v2m_io Input and Output
 * In the simplest form, only the input file name must be specified, as 
 * follows. Currently, *.inr, *.inr.gz, *.tif, and *.tiff input files 
 * are supported. Sample files are available in the \c data/v2m 
 * directory. In this case the output file is placed in the same 
 * directory in the *.vtu format.
 * 
 * \code
 * v2m ../data/v2m/randvol_30.inr
 * \endcode
 * 
 * The output format can be specified by using the \c --output-file-format
 * flag (\c -f). The period in the file type must be included. The 
 * following four file types are supported: 
 * - .mesh Medit format 
 * - .vtu VTK Unstructured XML format
 * - .vtk legacy VTK unstructured format
 * - .ex2 ExodusII format
 * 
 * \code
 * v2m ../data/v2m/randvol_30.inr -f .vtk
 * \endcode
 * 
 * The complete path and file for the output may also be specified with
 * the \c --output-file (\c -o) flag as such or without any flag. The
 * following examples are all valid methods for specifing the input and
 * output files. The flags may be used, but by default the first input
 * without a flag is considered the input file and the second without 
 * a flag the output file.
 * \code
 * v2m -i ../data/v2m/randvol_30.inr -o ./randvol.vtk
 * v2m ../data/v2m/randvol_30.inr -o ./randvol.vtk
 * v2m ../data/v2m/randvol_30.inr ./randvol.vtk
 * v2m --normalize ../data/v2m/randvol_30.inr --enable-file ./randvol.vtk
 * \endcode
 * 
 * \subsection v2m_criteria Mesh Criteria
 * The meshing behavior may be modified with the meshing criteria list
 * below. Details of the various critiera may be found in the CGAL 
 * documentation, <a href="http://www.cgal.org/Manual/latest/doc_html/cgal_manual/Mesh_3/Chapter_main.html"> Chapter 50, 3D Mesh Generation.</a> 
 * - \c --facet-angle
 * - \c --facet-size
 * - \c --facet-distance
 * - \c --cell-size
 * - \c --cell-radius-edge-ratio 
 * 
 * For example, the following changes the facet angle to 20. 
 * \code
 * v2m ../data/v2m/input/randvol_30.inr --facet-angle=20
 * \endcode
 * 
 * By default \c v2m prints the input settings and meshing results
 * to the screen. These results are briefly detailed in \ref v2m_results.
 * This ouptut may be disable by using the \c --disable-screen
 * flag. It is also possible to print this information to a file
 * with the \c --enable-file flag. By setting this flag a file is
 * created with the same name as the output with a .info extension 
 * added.
 * 
 * \code
 * v2m ../data/input/v2m/randvol_30.inr --disable-screen
 * v2m ../data/input/v2m/randvol_30.inr --enable-file
 * \endcode
 * 
 * An additional option, \c --normalize, exists for setting the parameters
 * to match the dimensional shifts of the image (see \ref v2m_dim).
 * For example, the following two commands will produce very different
 * meshes because of the change in dimension. The first image is given
 * dimensions of cube of size 30 by default, the second the same image
 * is set to a unit cube. 
 * 
 * \code
 * v2m ../data/input/v2m/randvol_30.inr
 * v2m ../data/input/v2m/randvol_30.inr --xdim=1 --ydim=1 --zdim=1
 * \endcode 
 * 
 * The \c --normalize flag removes this difference
 * by normalizing the mesh criteria by the voxel dimensions (except the
 * \c --cell-radius-edge-ratio). This option allows for the criteria to 
 * adjust to the image and allows the same mesh criteria values to apply
 * to images of various sizes. For example, the following commands will 
 * produce nearly identical meshes.
 * 
 * \code
 * v2m ../data/input/v2m/randvol_30.inr
 * v2m ../data/input/v2m/randvol_30.inr --normalize --xdim=1 --ydim=1 --zdim=1
 * \endcode 
 * 
 * \subsection v2m_dim Mesh Dimensions
 * The \c v2m program has two methods for specifing the image 
 * dimensions. By default \c v2m assumes that the dimensions are
 * equivalent to the number of pixels (e.g., a 30 pixel cube would be 
 * dimensioned as 30x30x30).
 * 
 * The first method to specify the image dimensions is by the \c xdim,
 * \c ydim, and \c zdim command-line arguments. For example, the command
 * below sets the dimensions of the image to a unit cube.
 * 
 * \code
 * v2m ../data/input/v2m/randvol_30.inr --xdim=1 --ydim=1 --zdim=1
 * \endcode
 * 
 * The second method defines the image size based on the size of a pixel. 
 * This method supercedes the aforemention method, thus if both are usec
 * the voxel (pixel) dimension is used. The following command sets the
 * dimensions of the complete image to 15 x 15 x 15, by setting the voxel
 * size to a 0.5 cube.
 * 
 * \code
 * v2m ../data/input/v2m/randvol_30.inr --vx=0.5 --vy=0.5 --vz0.5
 * \endcode
 * 
 * \subsection v2m_optimize Mesh Optimization
 * GCAL offers four optimization routines for generating meshes, which
 * are detailed in the <a href="http://www.cgal.org/Manual/latest/doc_html/cgal_manual/Mesh_3/Chapter_main.html#Subsection_50.2"> Section 50.2, "Interface" </a> of the CGAL user manual. \ref v2m_example2 "Example 2"
 * demonstrates the use and importance of mesh optimization.
 * 
 * By default the perturb and exude options are used, these can be disabled
 * with the \c --disable-exude and \c --disable-perturb flags.
 *
 * The Lloyd and Odt optimization are not used by default becuase they add 
 * signficant computation time to the meshing. These options may be enabled
 * with the \c --enable-lloyd and \c --enable-odt flags.
 *
 * Also, the \c --enable-all and --disable-all flags allow for the user to toggle
 * on or off all of the optimization routines. 
 * 
 * Each of the methods has a number of options associated, these values are discussed
 * in the CGAL documentation in detail and the available values may be 
 * viewed from the command line by using the \c --advanced flag.
 *
 * \subsection v2m_subdomain Meshing Subdomains
 * It is also possible to mesh subdomains and set varying mesh criteria
 * on each subdomain. \ref v2m_example1 demonstrates the usage of this 
 * feature by mimicing the example in the CGAL documentation.
 * 
 * A subdomain is defined in the image file by the values associated 
 * with the pixels, which should be between 0 through 255. 
 * Three examples of adding a subdomain with the \c --subdomain command 
 * as follow. The first indicates that all pixels with a value of 20 
 * will have a \c cell-size of 2, all others will use the default value.
 * 
 * The second command sets both the default value and the subdomain value;
 * the default \c cell-size is set to 4 and the subdomain value is again
 * set to 2.
 * 
 * The third command demonstates that there are two subdomains, pixel
 * values of 20 and 40. The default \c cell-size is set to 4 and set to
 * 2 and 1 for the two subdomains defined, respectively.
 * 
 * \code
 * v2m ./input.inr --subdomain 20 --cell-size 2
 * v2m ./input.inr --subdomain 20 --cell-size 4 2
 * v2m ./input.inr --subdomain 20 40 --cell-size 4 2 1
 * \endcode
 * 
 * Currently only the \c cell-size criteria may be different
 * on each subdomain, all other values will use the default. So, if 
 * you are adding the subdomain feature and want to change the criteria
 * for the \c --facet-size you must specify to values, one to set the 
 * default (2) and the other for the subdomain, the latter of which will not 
 * be used (9 in the command below). It is planned to have adjustable 
 * parameters for each criteria, but it may not be possible with the 
 * CGAL library.
 * 
 * \code
 * v2m ./input.inr --subdomain 20 --cell-size 2
 * v2m ./input.inr --subdomain 20 --cell-size 4 2 --facet-size 1 9
 * \endcode
 *
 * The \c --subdomain option is only need to change the meshing criteria
 * for the various domains. The program automatically defines the 
 * subdomains, unless the \c --disable-subdomains flag is used. In this
 * case the ids are eliminated from the file. 
 * 
 * Finally, only the *.medit and *.ex2 output format support the ability 
 * to write the subdomains to a file.   
 * 
 * \section v2m_results Mesh Results
 * As mentioned above \c v2m outputs various parameters that
 * describe the resulting mesh. An example output is given below. 
 * 
 * \verbinclude v2m_example_info_file.txt
 * 
 * Notice that the output includes a table of mesh quality results. 
 * These are computed using the VTK mesh quality class
 * (http://www.vtk.org/doc/nightly/html/classvtkMeshQuality.html). The 
 * name quality refers to the quality measure member function used to
 * compute the values in that row. For example, the "Edge ratio" uses
 * the \c SetTetQualityMeasureToEdgeRatio() member function of the
 * \c vtkMeshQuality class. 
 * 
 * The quality parameters are calculated for 
 * each element, the data presented in the table summarizes the results
 * and include the lower and upper bounds of the computed values, the 
 * average value for the entire mesh, the standard deviation, and 
 * coeffient of variation.
 * 
 * \section v2m_inr INR Format
 * The *.inr format is the main format used by the CGAL library.
 * Currently, \c v2m only supports direct input of *.inr files and 
 * *.tif files (8-bit grey scale only). The *.tif support uses the ITK (www.ITK.org)
 * libraries to covert the *.tif into a raw file format (*.inr without 
 * the header) that CGAL can read.
 * 
 * For convenience a set of MATLAB functions was developed from
 * for converting tiff files to inr format (e.g., /matlab/v2m/tif2inr.m). 
 * For example, the following MATLAB commands convert the randvol_30.tif 
 * file in the data directory into a *.inr file.
 * 
 * \code
 * >> tif2inr('../../data/v2m/randvol_30.tif');
 * \endcode
 * 
 * Additional functions exists within the \c /matlab/v2m directory
 * for generating random volumes, viewing volumes, and save and converting
 * between *.tif and *.inr files. The conversion functions were developed
 * from the ISO2Mesh library (http://iso2mesh.sourceforge.net/cgi-bin/index.cgi).
 * 
 * \section v2m_examples Example Usage
 * The following examples utilize various example input data files, these
 * files are contained in the \c /data/v2m/input directory. In all
 * cases the commands export the meshes to the \c /data/v2m/output
 * directory. The lloyd and odt optimization is were not used 
 * in all cases except in \ref v2m_example2 "Example 2"; enabling these 
 * options causes the  * mesh time required to increase significantly.
 * 
 * \subsection v2m_example1 Example 1: CGAL Liver
 * The CGAL library includes various meshing examples, this section 
 * demonstrates how these examples may be reproduced using v2m. The
 * first example produces a 3D mesh of a liver in *.vtu format. This 
 * code should be executed from the \c /bin directory. The resulting
 * 3D mesh is displayed below the code snippet.
 * 
 * \code
 * ./v2m ../data/v2m/input/liver.inr.gz ../data/v2m/output/liver.vtu 
 * \endcode
 * \image html v2m_liver.png "Mesh of CGAL Liver"
 * 
 * The resulting mesh has 16,702 elements and 4,017 external faces. The
 * resulting mesh had a mean condition number of 1.27 with a standard
 * deviation of 0.075 and coefficient of variation of 5.9%. 
 * 
 * The liver image contains sub-domains, the following command will
 * change the \c cell-size to 2 for the kidney portion of the image as
 * shown below.
 * 
 * \code
 * ./v2m ../data/v2m/input/liver.inr.gz ../data/v2m/output/liver.ex2 --subdomain 127 --cell-size 2
 * \endcode
 * \image html v2m_liver_subdomain.png "Mesh of CGAL Liver with subdomain" 
 * 
 * \subsection v2m_example2 Example 2: CGAL Brain
 * In similar example to \ref v2m_example1 "Example 1", the CGAL
 * brain.inr image may be meshed using the following commands. In this 
 * case, the output format is changed to *.vtk format and the \c --facet-size 
 * is reduced to 3 from the default of 6.
 * 
 * \code
 * ./v2m ../data/v2m/input/brain.inr ../data/v2m/output/brain.vtk --facet-size 4 
 * \endcode
 * \image html v2m_brain.png "Mesh of CGAL Brain"
 *
 * The resulting mesh has 21,019 elements and 6,716 external faces. The
 * resulting mesh had a mean condition number of 1.44 with a standard
 * deviation of 0.21 and coefficient of variation of 14.5%. This mesh
 * took 15.8 sec. to generate.
 * 
 * While little difference visually is detectable, the importance of 
 * using optimization techniques is apparent when 
 * compared. The brain image was executed with various optimization 
 * settings. The mesh criteria remained the same as above. The results 
 * of these runs are shown in the following table.
 * 
 * \htmlinclude v2m_optimization_table.txt
 * 
 * \subsection v2m_example3 Example 3: Porous Volume
 * This example demonstrates the usage of \c v2m to mesh a 3D image
 * of a random porous volume. It provides examples of how the various
 * mesh parameters can effect the mesh. Each mesh created in this section
 * is labeled with a reference (e.g., 3a and 3b). This reference is used
 * to report various mesh statistics and quality measures in the table 
 * at the end of the section.
 * 
 * The first step is to develop the volume, a 30 pixel cube is used
 * here. This volume is generated using the MATLAB
 * function \c randvol.m is uses (this function is located in the 
 *	\c ~/matlab/v2m directory). The following MATLAB commands, 
 * executed from this directory created the base image used here: 
 * \c randvol_30.inr, which is plotted using MATLAB following the code.
 * 
 * \code
 * >> v = randvol(30);
 * >> showvol(v);
 * >> saveinr(v, '../../data/v2m/input/randvol_30.inr');
 * \endcode
 * 
 * Using the default settings for \c v2m, as follows, creates a mesh
 * (Mesh 3a) that obviously does not capture the detail of the image. 
 * 
 * \code
 * ./v2m ../data/v2m/input/randvol_30.inr ../data/v2m/output/randvol_30.vtu
 * \endcode
 * 
 * One method to improve the detail is to adjust the facet size as shown
 * in the following images (Mesh 3b and 3c, respectively). It is also 
 * possible to modifiy the cell-size to achieve more detailed mesh
 * (Mesh 3d).
 * 
 * \code
 * ./v2m ../data/v2m/input/randvol_30.inr ../data/v2m/output/randvol_30.vtu --facet-size=2 
 * ./v2m ../data/v2m/input/randvol_30.inr ../data/v2m/output/randvol_30.vtu --facet-size=0.25
 * ./v2m ../data/v2m/input/randvol_30.inr ../data/v2m/output/randvol_30.vtu --facet-size=2 --cell-size=0.5
 * \endcode
 * 
 * \htmlonly
 * <table align=center>
 * <tr>
 * <td style="width:320px"> <img src="v2m_ex3_matlab.png" style="width:300px; height:300px;" /> </td>
 * <td style="width:320px"><img src="v2m_ex3_1.png" style="width:300px; height:300px;" /> </td>
 * </tr><tr>
 * <td align=center><b> Matlab random volume </b></td>
 * <td align=center><b> Default v2m mesh (Mesh 3a) </b></td>
 * </tr>
 * </table>
 * \endhtmlonly
 * 
 * \htmlonly
 * <table align=center>
 * <tr>
 * <td style="width:320px"> <img src="v2m_ex3_2.png" style="width:300px; height:300px; float:center" /> </td>
 * <td style="width:320px"> <img src="v2m_ex3_3.png" style="width:300px; height:300px; float:center" /> </td>
 * </tr><tr>
 * <td align=center><b> facet-size = 2 (Mesh 3b) </b></td>
 * <td align=center><b> facet-size = 0.25 (Mesh 3c) </b></td>
 * </tr>
 * </table>
 * \endhtmlonly 
 * \n\n
 * \image html v2m_ex3_4.png "facet-size = 2; cell-size = 0.5 (Mesh 3d)"
 * 
 * Changing the facet-distance allows for a mesh (Mesh 3e) that applies detail 
 * where it is needed. For example:
 * 
 * \code
 * ./v2m ../data/v2m/input/randvol_30.inr ../data/v2m/output/randvol_30.vtu --facet-size=2  --facet-distance=0.1
 * \endcode
 * 
 * \image html v2m_ex3_5.png "facet-size = 2; facet-distance = 0.1 (Mesh 3e)"
 * 
 * The table below summarizes a portion of the mesh results as reported
 * by \c v2m. Obviously, the possible settings are endless and 
 * dependant on the problem, but hopefully this example demonstrated 
 * how to use \c v2m and modify the command-line to fit your specific needs.
 * 
 * \htmlinclude v2m_example3_table.txt
 * 
 * \n\n\n
 */
  
// Header files for the Vol2mesh class
#include "vol2mesh/include.h" 
using namespace SlaughterVol2mesh;

//! A subfunction for defining and gathering command-line options
/*! \param argc Number of commandline parameters
 *  \param argv Character array containing the command line text
 * 
 * This function implements the UserOptions class for defining 
 * command line inputs, see \ref test_user_options for details on 
 * using this class.
 */ 
UserOptions v2m_command_line_options(int argc, char** argv);

//! Function for getting the output file
/*!
 * The \c v2m program will automatically generate a file for outputting
 * the mesh using the \c --input-file and \c --output-format by simply
 * changing the extension of the input file to that specified by
 * the format.
 * 
 * \param opt UserOptions class containing the user specified
 */ 
string v2m_get_output_file(UserOptions& opt);

//! Function for setting up the image based on the user input
/*! The meshing behavior may be defined by the user from the command
 * line, this functions uses the supplied (or default) values and 
 * inserts the values into the Vol2mesh object.
 *
 * \param v2m Vol2mesh class initilized with the image file
 * \param opt UserOptions class containing the user specified
 * 	values from the command-line
 */
void v2m_initialize_image(Vol2mesh& v2m, UserOptions& opt); 

//! Function for setting the image dimensions
/*! 
 * \param v2m Vol2mesh class initilized with the image file
 * \param opt UserOptions class containing the user specified
 * 	values from the command-line
 */ 
void v2m_set_dimensions(Vol2mesh& v2m, UserOptions& opt);

//! Function for setting the mesh criteria from command line options
/*! 
 * \param v2m Vol2mesh class initilized with the image file
 * \param opt UserOptions class containing the user specified
 * 	values from the command-line
 */ 
void v2m_set_mesh_criteria(Vol2mesh& v2m, UserOptions& opt); 
 
//! Main program for generating meshes from 3D image data
/*! 
 * This function is the main program executed when using the v2m
 * executable, details are provided here: \ref v2m.
 * 
 * The command-line arguments are parsed into the two input variables 
 * of this program:
 * \param argc The command-line argument count
 * \param argv An array of characters that contain the command-line text, 
 * this text is parsed with the UserOptions class in the v2m_command_line_options
 * sub-function.
 */ 
int main(int argc, char* argv[]){
    // Use CGALs timer to keep track of execution time
	CGAL::Timer t;
	t.start();    
	
	// Get the user options from the command line
	UserOptions opt = v2m_command_line_options(argc, argv);
	
	// Get the output file name
	string output_file;
	output_file = v2m_get_output_file(opt);

	// Vol2mesh class
	Vol2mesh vol2mesh(opt.get<string>("input-file"));

	// Setup the image based on the command-line options
	v2m_initialize_image(vol2mesh, opt); 

	// Generate the mesh
	vol2mesh.generate_mesh();

    // Output the mesh to a file
    vol2mesh.write(output_file);

	// Stop the timer and display the results  
	t.stop();
    vol2mesh.print_results(t.time(), opt.get_flag("disable-screen"), opt.get_flag("enable-file"));
    
	// Program complete
    return 0;
}

// A sub-function for defining and gathering command-line options
UserOptions v2m_command_line_options(int argc, char** argv){
	
	// The general options and title
	UserOptions gen("General Options");
	gen.add_title("\nThis function builds a mesh from 2D images slices.The basis for the\ndevelopment comes from the following CGAL library example:\nSection 50.3.3 Domains From Segmented 3D Images\n(http://www.cgal.org/Manual/latest/doc_html/cgal_manual/Mesh_3/Chapter_main.html)\n\n");
	
	gen.add_flag("help,h", "List the available options");
	gen.add_flag("advanced", "Show the complete list of options");
	gen.add_option<string>("config", "Specify a configuration file");
	gen.add_flag("enable-file", "Create a *.*.info file of mesh results");
	gen.add_flag("disable-screen", "Disable printing the mesh results to the screen");
	
	// Input and output related items
	UserOptions io("Input/Output Options");
	io.add_option<string>("input-file,i", "Name of the input file, it must be a *.tiff, *.inr, or *.inr.gz file", 1);
	io.add_option<string>("output-file,o", "Output filename", 1);
	io.add_option<string>("output-format,f", ".ex2", "Output file format");
	
	// Define default vectors and text for meshing critiera
	vector<double> a_vec(1,30);	string a_txt("[30]");
	vector<double> s_vec(1,6);	string s_txt("[6]");
	vector<double> d_vec(1,4);	string d_txt("[4]");
	vector<double> r_vec(1,3);	string r_txt("[3]");
	vector<double> c_vec(1,8);	string c_txt("[8]");
	
	// Meshing related options
	UserOptions crit("Meshing Criteria Options");
	crit.add_option<vector<double> >("facet-angle,a", a_vec, "Facet angle mesh criteria for CGAL", a_txt);
	crit.add_option<vector<double> >("facet-size,s", s_vec, "Facet size mesh criteria for CGAL", s_txt);
	crit.add_option<vector<double> >("facet-distance,d", d_vec, "Facet distance mesh criteria for CGAL", d_txt);
	crit.add_option<vector<double> >("cell-radius-edge-ratio,r", r_vec, "Cell radius to edge ratio mesh criteria for CGAL", r_txt);
	crit.add_option<vector<double> >("cell-size,c", c_vec, "Cell size mesh criteria for CGAL", c_txt);
	crit.add_flag("normalize","Normalize mesh criteria based on the of pixels in each direction");
	
	// Mesh dimensions
	UserOptions dim("Mesh Dimension Options");
	dim.add_option<double>("xdim", "x-dimension of the image (superseded by vx)");
	dim.add_option<double>("ydim", "y-dimension of the image (superseded by vy)");
	dim.add_option<double>("zdim", "z-dimension of the image (superseded by vz)");
	dim.add_option<double>("vx", 1, "voxel (pixel) size in x-direction");
	dim.add_option<double>("vy", 1, "voxel (pixel) size in y-direction");
	dim.add_option<double>("vz", 1, "voxel (pixel) size in z-direction");

	// Subdomain options
	UserOptions dom("Subdomain Options");
	dom.add_option<vector<int> >("subdomain","List of subdomain ids");
	dom.add_flag("disable-subdomains", "Exclude subdomain index information from output");
	
	// Mesh optimization options
	UserOptions opt("Meshing Optimization Options");
	opt.add_flag("enable-lloyd", "Enable Lloyd optimization");
	opt.add_flag("enable-odt", "Enable Odt optimization");
	opt.add_flag("disable-perturb", "Disable Perturb optimization");
	opt.add_flag("disable-exude", "Disable Exude optimization");
	opt.add_flag("disable-all", "Disable all optimization routines");
	opt.add_flag("enable-all", "Enable all optimization routines");

	// Lloyd mesh optimization settings
	UserOptions set0("Lloyd Optimization Settings");
	set0.add_option<double>("lloyd.time-limit", 0, "GCAL Lloyd optimization time limit (0 = no limit)");
	set0.add_option<std::size_t>("lloyd.max-iteration", 0, "GCAL Lloyd optimization max iterations allowed (0 = no limit)");
	set0.add_option<double>("lloyd.convergence", 0.02, "GCAL Lloyd optimization convergence limit");
	set0.add_option<double>("lloyd.freeze-bound", 0, "GCAL Lloyd optimization freeze bound limit");
	set0.hidden = true;
	
	// Odt mesh optimization settings
	UserOptions set1("Odt Optimization Settings");
	set1.add_option<double>("odt.time-limit", 0, "GCAL Odt optimization time limit (0 = no limit)");
	set1.add_option<std::size_t>("odt.max-iteration", 0, "GCAL Odt optimization max iterations allowed (0 = no limit)");
	set1.add_option<double>("odt.convergence", 0.02, "GCAL Odt optimization convergence limit");
	set1.add_option<double>("odt.freeze-bound", 0, "GCAL Odt optimization freeze bound limit");	
	set1.hidden = true;	
		
	// Perturb mesh optimization settings
	UserOptions set2("Perturb Optimization Settings");
	set2.add_option<double>("perturb.time-limit", 0, "GCAL Perturb optimization time limit (0 = no limit)");
	set2.add_option<double>("perturb.sliver-bound", 0, "GCAL Perturb sliver lower bounds");
	set2.hidden = true;
	
	// Perturb mesh optimization settings
	UserOptions set3("Exude Optimization Settings");
	set3.add_option<double>("exude.time-limit", 0, "GCAL Exude optimization time limit (0 = no limit)");
	set3.add_option<double>("exude.sliver-bound", 0, "GCAL Exude sliver lower bounds");
	set3.hidden = true;
	
	// Link the groups together
	gen.add(io).add(crit).add(dim).add(dom).add(opt).add(set0).add(set1).add(set2).add(set3);

	// Apply the command-line options
	gen.apply_options(argc, argv);

	// If --advanced flag is used display the all of the options
	if (gen.get_flag("advanced")){
		gen.show_hidden();
	}

	// Return the main class
	return gen;
}

// Return the output filename
string v2m_get_output_file(UserOptions& opt){
	
	// String to contain output file
	string output;

	// Case when the file is specified
	if (opt.exist("output-file")){
		output.assign(opt.get<string>("output-file"));
	
	// Case when the file is NOT specified
	} else {
		FileParts in(opt.get<string>("input-file"));
		in.ext.assign(opt.get<string>("output-format"));
		in.update();
		output.assign(in.full);
	}

	// Return the string
	return output;
}

// A sub-function for initilizing the Vol2mesh class based on the user options
void v2m_initialize_image(Vol2mesh& v2m, UserOptions& opt){
    
    // Get the number of pixels in the image
	vector<int> pix = v2m.pixels();
	
	// Set the optimization behavior
	v2m.set_optimization("lloyd", opt.get_flag("enable-lloyd"));
	v2m.set_optimization("odt", opt.get_flag("enable-odt"));
	v2m.set_optimization("perturb", !opt.get_flag("disable-perturb"));
	v2m.set_optimization("exude", !opt.get_flag("disable-exude"));

	// Override individual settings with --disable-all
	if (opt.get_flag("disable-all")){
		v2m.set_optimization("lloyd", false);
		v2m.set_optimization("odt", false);
		v2m.set_optimization("perturb", false);
		v2m.set_optimization("exude", false);
	}

	// Override individual settings with --enable-all
	if (opt.get_flag("disable-all")){
		v2m.set_optimization("lloyd", true);
		v2m.set_optimization("odt", true);
		v2m.set_optimization("perturb", true);
		v2m.set_optimization("exude", true);
	}

	// Apply the user-defined optimizaiton settings
	// Set the Lloyd settings
	v2m.lloyd.time_limit = opt.get<double>("lloyd.time-limit");
	v2m.lloyd.max_iteration = opt.get<std::size_t>("lloyd.max-iteration");
	v2m.lloyd.convergence = opt.get<double>("lloyd.convergence");
	v2m.lloyd.freeze_bound = opt.get<double>("lloyd.freeze-bound");
	
	// Set the Odt settings
	v2m.odt.time_limit = opt.get<double>("odt.time-limit");
	v2m.odt.max_iteration = opt.get<std::size_t>("odt.max-iteration");
	v2m.odt.convergence = opt.get<double>("odt.convergence");
	v2m.odt.freeze_bound = opt.get<double>("odt.freeze-bound");

	// Set the perturb settings
	v2m.perturb.time_limit = opt.get<double>("perturb.time-limit");
	v2m.perturb.sliver_bound = opt.get<double>("perturb.sliver-bound");
	
	// Set the exude settings
	v2m.exude.time_limit = opt.get<double>("exude.time-limit");
	v2m.exude.sliver_bound = opt.get<double>("exude.sliver-bound");

	// Set the disble-subdomain behavior
	v2m.disable_subdomains(opt.get_flag("disable-subdomains"));

	// Set the image dimensions
	v2m_set_dimensions(v2m, opt);

	// Set the meshing criteria
	v2m_set_mesh_criteria(v2m, opt); 

}

// A sub-function for setting the image dimensions
void v2m_set_dimensions(Vol2mesh& v2m, UserOptions& opt){
	
	// Define a storage location for the image dimensions and voxel sizes
	vector<double> dim;
	vector<double> vox;
	
	// Vector of option handles for the user supplied image dimensions
	vector<string> dtxt;
	dtxt.push_back("xdim"); dtxt.push_back("ydim"); dtxt.push_back("zdim");	
	
	// Vector of option handles for the user supplied voxel dimensions
	vector<string> vtxt;
	vtxt.push_back("vx"); vtxt.push_back("vy"); vtxt.push_back("vz");	
	
	// Loop through each dimesion
	for (int i = 0; i < dtxt.size(); ++i){
	
		// If the option exists, set the voxel using the image dimension
		if (opt.exist(dtxt[i].c_str())){
			v2m.set_dimension(i, opt.get<double>(dtxt[i].c_str()));
		}
			
		// If the option exists, set this voxel with the user defined voxel size
		// (this overwrites the above if defined for the same direction)
		if (opt.exist(vtxt[i].c_str())){
			v2m.set_voxel(i,opt.get<double>(vtxt[i].c_str())); 
		}	
	}	
} 

// A sub-function for applying the command-line mesh criteria
void v2m_set_mesh_criteria(Vol2mesh& vol2mesh, UserOptions& opt){
	
	// Trigger normalization if desired
	vol2mesh.normalize(opt.get_flag("normalize"));

	// Extract user supplied the mesh criteria into a vector of vectors
	vector<vector<double> > mat;
	mat.reserve(5);
	mat.push_back(opt.get<vector<double> >("facet-angle")); 			// a	
	mat.push_back(opt.get<vector<double> >("facet-size"));  			// b
	mat.push_back(opt.get<vector<double> >("facet-distance")); 			// c
	mat.push_back(opt.get<vector<double> >("cell-radius-edge-ratio")); 	// d
	mat.push_back(opt.get<vector<double> >("cell-size")); 				// e

	// Determine the number of subdomains
	int n = 0;
	if (opt.exist("subdomain")){
		vector<int> id = opt.get<vector<int> >("subdomain");
		n = id.size();
	}
	
	// Get the default values, as an object and vector
	Vol2meshMeshCriteria dmc = vol2mesh.get_default_criteria();
	vector<double> dmc_vec = dmc.get_vector();

	// Loop through the user supplied vectors
	for (int i = 0; i < mat.size(); i++){
		// If the supplied vector exceeds the number of subdomains replace the
		// default value with the first value of the user supplied values
		if (mat[i].size() > n){
			dmc_vec.at(i) = mat[i][0];
			mat[i].erase(mat[i].begin());
		}
	}

	// Update the default mesh criteria
	dmc.init(dmc_vec);
	vol2mesh.set_default_criteria(dmc);

	// Add subdomain critiera, if they do not exist we are done
	if (opt.exist("subdomain")){
		
		// Extract a vector of ids
		vector<int> id = opt.get<vector<int> >("subdomain");
		
		// Pad the user vectors with the defaults
		for (int i = 0; i < mat.size(); i++){
			int n = id.size() - mat[i].size();			// num. of values to add
			mat[i].insert(mat[i].end(), n, dmc_vec[i]); // insert the defaults at the end
		}
		
		// Add the criteria to the vol2mesh object
		for (int i = 0; i < id.size(); i++){
			Vol2meshMeshCriteria smc(mat[0][i], mat[1][i], mat[2][i], mat[3][i], mat[4][i]);
			vol2mesh.add_subdomain(id[i], smc);
		}
	}		 
}
