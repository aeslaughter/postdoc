/** \file vol2mesh.h
 * Header file for source code, vol2mesh.cpp, which has the necessary
 * includes and function prototypes for the source code.
 */
 
// Includes for 3D mesh generation (CGAL)
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Labeled_image_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/Image_3.h>

// Includes for exporting the file as a *.vtu or *.vtk file
#include <CGAL/IO/Complex_3_in_triangulation_3_to_vtk.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkXMLDataSetWriter.h>
#include <vtkExodusIIWriter.h>
#include <vtkAlgorithm.h>
#include <vtkUnstructuredGridAlgorithm.h>
 
// Mesh optimization
#include <CGAL/refine_mesh_3.h>
#include <CGAL/lloyd_optimize_mesh_3.h>
#include <CGAL/odt_optimize_mesh_3.h>
#include <CGAL/perturb_mesh_3.h>
#include <CGAL/exude_mesh_3.h>

// VTK includes for mesh quality calculation
#include <vtkGeometryFilter.h>
#include <vtkMeshQuality.h>
#include <vtkFieldData.h>

// ITK includes
#include <itkImageFileReader.h>
#include <itkImage.h> 

// Include CGAL timer
#include <CGAL/Timer.h>

// Boost libraries
#include <boost/shared_ptr.hpp>

// Include the string and io stardard libraries
#include <string>
#include <stdio.h>
using std::string; 

// Use the standard library vector class
#include <vector>
using std::vector;

// Include the subdomain capable CGAL to VTK conversion function
#include "vol2mesh/Complex_3_subdomain_in_triangulation_3_to_vtk.h"

// Include the Vol2meshMeshCriteria class
#include "vol2mesh/vol2mesh_mesh_criteria.h"

// Include the my common library
#include "common/include.h"
using namespace SlaughterCommon;

// Add to my vol2mesh namespace
namespace SlaughterVol2mesh{

// Define type definitions for varius CGAL objects
//! Short-hand for CGAL kernel library class
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

//! Short-hand for 3D image domain CGAL library class 
typedef CGAL::Labeled_image_mesh_domain_3<CGAL::Image_3,K> Mesh_domain;

//! Short-hand for CGAL 3D triangulation library class
typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;

//! Short-hand for CGAL 3D complex triangulation library class
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;

//! Short-hand for CGAL 3D mesh critiera library class 
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

//! Short-hand for creating subdomain mech critiera variable
typedef CGAL::Mesh_constant_domain_field_3<Mesh_domain::R, Mesh_domain::Index> Sizing_field;

//! Short-hand notation for a vector of vectors (used in print_results)
typedef vector< vector<double> > matrix;

//! A container for storing CGAL Lloyd and Odt optimization settings
/*!
 * The Lloyd and Odt optimization methods with respect to mesh
 * generation are detailed here:
 * http://www.cgal.org/Manual/latest/doc_html/cgal_manual/Mesh_3_ref/Function_make_mesh_3.html
 * 
 * Details of the Lloyd optimization are given here:
 * http://www.cgal.org/Manual/latest/doc_html/cgal_manual/Mesh_3_ref/Function_lloyd_optimize_mesh_3.html
 * 
 * Details of the Odt optimization are given here:
 * http://www.cgal.org/Manual/latest/doc_html/cgal_manual/Mesh_3_ref/Function_odt_optimize_mesh_3.html
 */ 
struct odt_lloyd_settings{
	//! Time limit for optimization process (CGAL default 0; no limit)
	double time_limit;
	
	//! Maximum number of optimization iterations (CGAL default 0; no limit)
	std::size_t max_iteration;
	
	//! Convergence criteria (CGAL default 0.02)
	double convergence;
	
	//! Method for helping reduce run time (CGAL default 0.01)
	double freeze_bound;
};

//! A class for containing CGAL Perturb and Exude optimization settings
/*!
 * The Perturb and Exude optimization methods with respect to mesh
 * generation are detailed here:
 * http://www.cgal.org/Manual/latest/doc_html/cgal_manual/Mesh_3_ref/Function_make_mesh_3.html
 * 
 * Details of the Perturb optimization are given here:
 * http://www.cgal.org/Manual/latest/doc_html/cgal_manual/Mesh_3_ref/Function_perturb_mesh_3.html
 * 
 * Details of the Odt optimization are given here:
 * http://www.cgal.org/Manual/latest/doc_html/cgal_manual/Mesh_3_ref/Function_exude_mesh_3.html
 */ 
struct perturb_exude_settings{
	//! Time limit for optimization process (CGAL default 0; no limit)
	double time_limit;
	
	//! Lower bound on dihedral angles of mesh cells (CGAL default 0; no limit)
	double sliver_bound;
};


/*! \class Vol2mesh vol2mesh.h "vol2mesh/vol2mesh.h"
 * \brief A class for generating a 3D mesh from 2D image slices. This 
 * class acts as a wrapper to the CGAL, VTK, and ImageMagick++ libraries
 * to read, build, and export tetrahedral meshes from a pixel image. The
 * class is used by the \c v2m executable (\ref v2m) and supports
 * many of the features demonstrated in the CGAL documentation.
 *
 * CGAL Documentation:\n
 * http://www.cgal.org/Manual/latest/doc_html/cgal_manual/Mesh_3/Chapter_main.html
 * 
 * \see v2m.cpp
 */ 
class Vol2mesh{
	public:
		//! Class constructor
		/*!
		 * This constructor reads the image generating the CGAL::Image_3
		 * object, which is available through the \c image() member
		 * of this class. It also sets the pixel dimensions which is
		 * available from the \c pixels() member of this class. It sets
		 * the default mesh criteria (see Vol2meshMeshCriteria) values
		 * to the values in the CGAL examples:
		 * -# facet-angle = 30
		 * -# facet-size = 6
		 * -# facet-distance = 4
		 * -# cell-radius-edge-ratio = 3
		 * -# cell-size = 8
		 * 
		 * The normalization flag is set to false and the 
		 * optimization flags are set to the default values of the CGAL
		 * \c Make_mesh_3 function (see the set_optimization member).
		 * 
		 * The output file name is also set to "output.ex2."
		 * 
		 * The default optimization parameters are set to that of the CGAL
		 * library, the optimization methods to be used may be set 
		 * with the set_optimization function. The parameters for each
		 * method may be set with the corresponding data structure. Links
		 * to the settings are provided in the documentation for
		 * the storage structures: perturb_exude_settings and 
		 * odt_lloyd_settings.
		 * 
		 * The optimization routines are implemented in order as done 
		 * in the CGAl library: llyod, odt, perturb, exude. 
		 * 
		 * \param infile A std::string containg the name of the image 
		 * file to open; *.tiff, *.inr, and *.inr.gz files are supported
		 * 
		 * Note, this class contains no public attributes, access to any 
		 * attributes that are to be changed by the user is done
		 * with member functions designed to control the behavior 
		 * correctly. All private attributes are labeled with "_" 
		 * proceeding the name.
		 * 
		 * \see v2m.cpp
		 * \see odt_lloyd_settings
		 * \see perturb_exude_settings
		 * \see set_optimization
		 * \see Vol2meshMeshCriteria
		 */ 
		Vol2mesh(string infile);
		
		//! Sets the correct voxel size based on supplied value
		/*! The mesh dimensions are based on the dimensions of the 
		 * image voxels (pixel). This function sets the image voxel size
		 * for the desired dimension.
		 * 
		 * \param i The direction of interest: 0 = x; 1 = y; 2 = z
		 * \param vox The dimension of the pixel in specified direction
		 * 
		 * \see set_voxel(vector<double> vox)
		 * \see voxels()
		 */ 
		void set_voxel(int i, double vox);
		
		//! Sets the correct voxel size based on supplied vector
		/*! The mesh dimensions are based on the dimensions of the 
		 * image voxels (pixel). This function sets the image voxel size
		 * for based on the supplied vector of dimensions
		 * 
		 * \param vox The dimensions of the pixel in vector format
		 * 
		 * \see set_voxel(int i, double vox)
		 */ 
		void set_voxel(vector<double> vox);
		
		//! Sets the correct voxel size based on the overall image dimension
		/*! 
		 * The mesh dimensions are based on the dimensions of the 
		 * image voxels (pixel). This function sets the image voxel size
		 * for based on the supplied image dimension
		 * 
		 * \param i The direction of interest: 0 = x; 1 = y; 2 = z
		 * \param dim The dimension of the image in specified direction
		 * 
		 * \see set_dimension(vector<double> vox)
		 * \see set_voxel(vector<double> vox)
		 * \see set_voxel(int i, double vox)
		 */ 	
		void set_dimension(int i, double dim);
		
		//! Sets the correct voxel size based on vector of overall image dimensions
		/*! 
		 * The mesh dimensions are based on the dimensions of the 
		 * image voxels (pixel). This function sets the image voxel size
		 * for based on the supplied vector of image dimensions
		 * 
		 * \param dim The dimension of the image in vector format
		 * 
		 * \see set_dimension(int i, double dim)
		 * \see set_voxel(vector<double> vox)
		 * \see set_voxel(int i, double vox)
		 */ 
		void set_dimension(vector<double> dim);		
		
		//! Sets the output file
		/*!
		 * Allows the user to specify the output file seperate
		 * from the write function. Note, that the name given in the 
		 * write function overwrites this value.
		 * 
		 * \param output_file the output filename
		 * 
		 * \see write()
		 * \see write(string output_file)
		 */ 
		void set_output_file(string output_file);
		
		//! Returns a vector of the image size in pixels
		vector<int> pixels();
		
		//! Returns a vector of the voxel dimensions
		/*! 
		 * The voxel dimension dictate the dimensions of the image, this
		 * provides access to the current values.
		 * 
		 * \return A vector containg the dimensions of the image in pixels
		 * \see set_voxel(vector<double> vox)
		 * \see set_voxel(int i, double vox)
		 */ 
		vector<double> voxels();		
		
		//! Set default CGAL meshing criteria parameters
		/*! 
		 * Replaces the existing default meshing criteria with those
		 * defined in the instance of the Vol2meshMeshCriteria class
		 * supplied.
		 * 
		 * This is useful when using different criteria on sub domains,
		 * if the criteria are not defined on a subdomain with the 
		 * add_subdomain member function then the default values are 
		 * used, this allows the user to change them from what was 
		 * defined in the constructor.
		 * 
		 * If no subdomains this function is used to set the desired
		 * mesh criteria desired.
		 * 
		 * \param c Vol2meshMeshCriteria instance to replace the default 
		 * mesh criteria
		 * 
		 * \see get_default_criteria
		 * \see Vol2meshMeshCriteria
		 */ 
		void set_default_criteria(Vol2meshMeshCriteria c);
		
		//! Return a reference to the default meshing criteria
		/*!
		 * \return Reference to the default meshing criteria, a Vol2meshMeshCriteria class
		 * \see set_default_criteria
		 * \see Vol2MeshCriteria
		 */ 
		Vol2meshMeshCriteria& get_default_criteria();		
		
		//! Adds meshing criteria for a subdomain
		/*! 
		 * The actual values for each pixel (0 to 255) may be used
		 * to define subdomains within the image. This function 
		 * allows the user to specify different meshing criteria for the 
		 * subdomain.
		 * 
		 * \param id A value identifing the subdomain (0 to 255)
		 * \param c The mesh criteria to use for this subdomain
		 * 
		 * \see Vol2meshMeshCriteria
		 * \see disable_subdomains
		 */ 
		void add_subdomain(int id, Vol2meshMeshCriteria c);
		
		//! Returns reference to vector of Vol2meshMeshCriteria
		/*!
		 * The mesh criteria for the various subdomains that were added
		 * to this class using the add_subdomain function are stored in
		 * a vector, a reference is returned here to provide access
		 * to these criteria.
		 * 
		 * \return A vector of the subdomain mesh criteria
		 * 
		 * \see add_subdomain
		 */
		vector<Vol2meshMeshCriteria>& get_subdomain_criteria();
		
		//! Sets the mesh criterial normalization behavior
		/*! 
		 * The mesh criteria are defined with units that are proportional
		 * to the voxel size. This function will normalize the 
		 * supplied mesh criteria base on the voxel size so that the 
		 * criteria may be specified on a unit square basis, as such
		 * the criteria will behave similarily across images regardless
		 * of the pixel dimensions.
		 * 
		 * \param value A true/false value indicating to normalize or not
		 * 
		 * \see Vol2meshMeshCriteria
		 */ 
		void normalize(bool value);
		
		//! Disable all subdomain behavior
		/*!
		 * Setting this value to true will eliminate all subdomain
		 * options and generate the mesh with all the same criteria
		 * as defined by the class defaults. The mesh will also be exported
		 * without defining subdomain parameters.
		 *
		 * \param value A true/false value for disable/enabling subdomains
		 * 
		 * \see add_subdomian
		 * \see set_default_criteria
		 * \see get_default_criteria
		 */ 
		void disable_subdomains(bool value);
		
		//! Function for changing the mesh optimization settings
		/*!
		 * CGAL defines four methods for performing mesh optimization,
		 * this function allows the useage of these tools to be toggled
		 * on or off. 
		 * 
		 * Details of the optimization may be found in the CGAL
		 * documentation: 
		 * http://www.cgal.org/Manual/latest/doc_html/cgal_manual/Mesh_3/Chapter_main.html#Section_50.2
		 * 
		 * \param type A std::string that may be one of four values:
		 * 	"lloyd", "odt", "perturb", or "exude"
		 * \param value A boolean indicating whether the specified 
		 * optimization should be used.
		 */ 
		void set_optimization(string type, bool value);
		
		//! Returns a reference to the CGAL image object
		/*!
		 * The CGAL library coverts the image file into a Image_3 object,
		 * this function returns a reference to that object. See the 
		 * CGAL documentation for more details regarind working with
		 * this object.
		 * 
		 * \return Reference to CGAL::Image_3 object
		 * 
		 * \see Vol2mesh()
		 */ 
		CGAL::Image_3& image();

		//! Generates the 3D mesh
		/*!
		 * Until this function is called no mesh exists, as such
		 * the c3t3() function will return an empty reference. This
		 * function applies the defined mesh criteria and optimization
		 * to generate a mesh.
		 * 
		 * This function must also be called before the write() and 
		 * print_results() members.
		 * 
		 * \see set_optimization
		 * \see set_default_criteria
		 * \see get_default_criteria
		 * \see add_subdomain
		 * \see write
		 * \see print_results
		 */ 
		void generate_mesh();
			
		//! Returns a reference to the CGAL mesh
		/*! 
		 * The generate_mesh() function must be called prior to accessing
		 * the mesh.
		 *  
		 * \return A reference to the CGAL::C3t3 object containing the
		 * generated mesh
		 * 
		 * \see generate_mesh()
		 */ 	
		C3t3& c3t3();
				
		//! Exports the CGAL mesh to a file
		/*!
		 * This generate_mesh function must be called prior to writing
		 * the mesh to a file. Three file types are supported:
		 * -# *.medit - This is the only file format natively supported
		 * by the CGAL library, the following three types are all generated
		 * by convering the CGAL mesh to a VTK mesh.
		 * -# *.vtk - The legacy version of the VTK format
		 * -# *.vtu - The XML format of the VTK format
		 * -# *.ex2 - An Exodus II file.
		 * 
		 * \param output_file A string containg the filename to which the
		 * mesh will be exported, the extension dictates the format
		 * 
		 * \see write()
		 * \see set_input_file
		 * \see generate_mesh
		 */ 		
		void write(string output_file);

		//! Exports the CGAL mesh to a file
		/*!
		 * This generate_mesh function must be called prior to writing
		 * the mesh to a file. Three file types are supported:
		 * -# *.medit - This is the only file format natively supported
		 * by the CGAL library, the following three types are all generated
		 * by convering the CGAL mesh to a VTK mesh.
		 * -# *.vtk - The legacy version of the VTK format
		 * -# *.vtu - The XML format of the VTK format
		 * -# *.ex2 - An Exodus II file.
		 * 
		 * \see set_input_file
		 * \see write(string output_file)
		 * \see generate_mesh
		 */ 		
		void write();

		//! Exports the meshing results to a file and/or screen
		/*!
		 * The settings specified in the Vol2mesh class such as the 
		 * filenames and mesh criteria may be exported to the screen
		 * and or file using this function.
		 * 
		 * Also, the quality of mesh statistics are reported by this 
		 * function. These are computed using the VTK library
		 * 
		 * \param t The execution time
		 * \param disable_screen Set this to true to exclude screen display
		 * \param enable_file Set this to true to print a information
		 * file, it uses the *.info extension appeneded to the output filename
		 * ("output.ex2" by default)
		 * 
		 * \see write
		 */ 
		void print_results(double t, bool disable_screen = false, bool enable_file = false);
				
		//! Storage structure for Odt optimization settings
		/*!
		 * \see odt_lloyd_options
		 */ 
		odt_lloyd_settings odt;
		
		//! Storage structure for Lloyd optimization settings
		/*!
		 * \see odt_lloyd_options
		 */ 
		odt_lloyd_settings lloyd;
		
		//! Storage structure for Perturb optimization methods
		/*!
		 * \see perturb_exude_options
		 */ 
		perturb_exude_settings perturb;
				
		//! Storage structure for Exude optimization methods
		/*!
		 * \see perturb_exude_options
		 */ 
		perturb_exude_settings exude;	
		
	private:
		//! Default Vol2MeshCriteria class
		Vol2meshMeshCriteria default_criteria_;
		
		//! A vector of the meshing criteria for each subdomain added with add_subdomain
		std::vector<Vol2meshMeshCriteria> subdomain_criteria_;
		
		//! A flag for normalizing the mesh criteria
		bool normalize_;
		
		//! A flag for disabling the subdomain behavior
		bool disable_subdomains_;
		
		//! The input filename
		string input_file_;
		
		//! The mesh output filename
		string output_file_;

		//! Dimensions of pixels in the image
		vector<double> voxels_; 
		
		//! Number of pixels in each direction of entire image
		vector<int> pixels_; 	
		
		//! CGAL mesh object; created by generate_mesh
		C3t3 c3t3_;	
		
		//! CGAL image object; created by Vol2mesh()
		CGAL::Image_3 image_;

		//! Flag for toggling Lloyd optimization
		bool use_lloyd_;
		
		//! Flag for toggling Odt optimization
		bool use_odt_;
		
		//! Flag for toggling Perturb optimization
		bool use_perturb_;
		
		//! Flag for toggling Exude optimization
		bool use_exude_;

		//! Initilization function
		void init();

		//! Functoin for reading the 3D image into CGAL
		void read_image();		
		
		//! A function for reading a .tiff file (not working)
		void read_tiff();
		
		//! Returns the normalization value
		double get_normalize_value();

		//! Function for collecting the current VTK mesh quality statitics
		/*! 
		 * The VTK mesh quality calculation relies on setting the desired
		 * mesh quality measure (see \c get_all_quality_stats). After the 
		 * desired measure is set, this function returns the various statistics
		 * available. These values are returned in a vector object.
		 * 
		 * \param q A pointer to vtkMeshQuality class, this is created in the
		 * \c print_results sub-function
		 * \param name The name of the quality statitics to extract. By default
		 * it is set to "Mesh Tetrahedron Quality", which is the only valid 
		 * value for \c vol2mesh as it only generates tetrahedral meshes. This
		 * functionality was included for future expansion that may include
		 * mixed or different element types.
		 * \return A std::vector object containing the following mesh
		 * quality statistics:
		 * - v[0] = lower bounds of range
		 * - v[1] = uper bounds of range 
		 * - v[2] = average 
		 * - v[3] = standard deviation
		 * 
		 * \see get_all_quality_stats
		 */
		vector<double> get_current_stats(vtkMeshQuality* q, const char *name = "Mesh Tetrahedron Quality");
		
		//! Function for gathering all of the available mesh quality statistics
		/*! 
		 * The mesh quality of the CGAL generated mesh is computed using the 
		 * VTK library (www.vtk.org). This function extracts all of the available
		 * mesh quality measures in the vtkMeshQuality class
		 *  (www.vtk.org/doc/nightly/html/classvtkMeshQuality.html)
		 * and stores the names and values in the supplied vectors. The code for
		 * this function was adapted from the following VTK test function:
		 * VTK/Graphics/Testing/Cxx/MeshQuality.cxx.
		 * 
		 * \param q A pointer to vtkMeshQuality class, this is created in the
		 * \c print_results sub-function
		 * \param q_mat A vector of vectors that contains the various mesh 
		 * statistics for each of the available methods
		 * \param q_name A vector that contains the names of each of the mesh
		 * quality measures used
		 * 
		 * \see get_current_stats
		*/
		void get_all_quality_stats(vtkMeshQuality* q, matrix& q_mat, vector<string>& q_name);

};

} // SlaughterVol2Mesh namespace 
