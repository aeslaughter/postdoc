/*! \file vol2mesh_mesh_criteria.h
 * Header file for the Vol2meshMeshCriteria class
 */

// Standard includes
#include <vector>
using std::vector;

// Add this to the SlaughterVol2mesh namespace
namespace SlaughterVol2mesh{
	
/*! \class Vol2meshMeshCriteria vol2mesh_mesh_criteria.h "vol2mesh/vol2mesh_mesh_criteria.h"
 * \brief A class for storing and accessing mesh criteria data for the Vol2mesh class.
 * The Vol2mesh class relies on instances of this class to store
 * and retrieve the CGAL meshing criteria. Details of these critiera are provided in the 
 * CGAL documentation (see link below). There are five criteria, these terms are provided
 * here in order:
 * -# facet_angle
 * -# facet_size
 * -# facet_distance
 * -# cell_radius_edge_ratio
 * -# cell_size
 *
 * CGAL Documentation:\n
 * http://www.cgal.org/Manual/latest/doc_html/cgal_manual/Mesh_3/Chapter_main.html
 * 
 * \see Vol2mesh
 */  
class Vol2meshMeshCriteria{
	public:
		//! Default constructor
		/*! Sets all meshing criteria to zero.
		 * 
		 * \see Vol2meshMeshCriteria(vector<double>)
		 * \see Vol2meshMeshCriteria(double a, double b, double c, double d, double e)
		 */ 
		Vol2meshMeshCriteria();
		
		//! Vector input constructor
		/*!
		 * Constructor that accepts a vector, the values from the 
		 * vector are inserted in to the mesh criteria attributes.
		 * 
		 * \param vec A vector of mesh criteria; values must be in order
		 * shown in the main class documentation.
		 * 
		 * \see Vol2meshMeshCriteria()
		 * \see Vol2meshMeshCriteria(double a, double b, double c, double d, double e)
		 */ 
		Vol2meshMeshCriteria(vector<double> vec);
		
		//! Inidividual input constructor
		/*!
		 * Constructor that accepts a five doubles, the values are 
		 * inserted in to the mesh criteria attributes
		 * in order as listed in order:
		 * 
		 * \param a facet_angle
		 * \param b facet_size
		 * \param c facet_distance
		 * \param d cell_radius_edge_ratio
		 * \param e cell_size
		 * 
		 * \see Vol2meshMeshCriteria()
		 * \see Vol2meshMeshCriteria(vector<double>)
		 */
		Vol2meshMeshCriteria(double a, double b, double c, double d, double e);
		
		//! A method for extracting the meshing criteria as vector
		/*!
		 * \return A vector containing the five meshing criteria in
		 * order as specified in the main class documentation.
		 */ 
		vector<double> get_vector();
		
		//! A method for extracing a single mesh criteria by index
		/*!
		 * \param i Index of the desired criteria (0 through 4); criteria
		 * ordered as in the main class documentation
		 * 
		 * \return A double of the desired parameter
		 */ 
		double get_value(int i);
		
		//! Initilization function
		/*!
		 * Sets the mesh criteria to the values given in the vector
		 * in the same fashion as the vector version of the constructor.
		 * 
		 * \see Vol2meshMeshCriteria(vector<double>)
		 */
		void init(vector<double>);
		
		//! The CGAL facet angle mesh criteria
		double facet_angle;
		
		//! The CGAL facet size mesh criteria
		double facet_size;
		
		//! The CGAL facet distance mesh criteria
		double facet_distance;
		
		//! The CGAL cell radius to edge ratio mesh criteria
		double cell_radius_edge_ratio;
		
		//! The CGAL cell size meshing criteria
		double cell_size;
		
		//! The subdomain id
		int id;		
}; // end class
}  // end namespace
