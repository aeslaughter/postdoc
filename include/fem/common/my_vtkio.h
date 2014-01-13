/*! \file my_vtkio.h
  * \brief Header  for MyVTKIO class.
  */
  
// Avoid multiple includes
#ifndef my_vtkio_h
#define my_vtkio_h

// Standard C++ includes
#include <stdio.h>
#include <vector>
#include <string>
using std::vector;
using std::string;

// libMesh includes
#include <libmesh.h>
#include <vtk_io.h>
#include <mesh.h>
#include <mesh_data.h>
#include <equation_systems.h>
using namespace libMesh;

// My includes
#include "common/include.h"
using namespace SlaughterCommon;

// Add this to my FEM namespace
namespace SlaughterFEM{
/*! \class MyVTKIO my_vtkio.h "fem/common/my_vtkio.h"
 * \brief An extension of the libMesh::VTKIO class.
 */  
class MyVTKIO : public VTKIO{
	public:
	
		//! Class constructor
		/*!
		  * A class for enhancing the behavior of libMesh::VTKIO to create temporal *.vtu files, with a 
		  * single *.pvd file.
		  *	
		  * 
		  * \param filename String containing the output filename, the extension should be *.vtu and it 
		  * should not contain any information regarding the time step, this is added automatically
		  * \param es A reference to the EquationSystems object that is being output
		  */
		MyVTKIO(string filename, EquationSystems& es);
	
		//! Writes the data to a file
		/*!
		 * Automatically outputs all of the data in the EquationSystem to a new *.vtu file. It also creates a new *.pvd file
		 * with the correct links to all the files in the series.
		 * \param t The current time
		 */
		void write(double t);
		
		//! Sets the number of zeros to append onto the *.vtu files
		/*!
		 * \param pad The number of digits to include in the *.vtu filename (defaults is 3)
		  */
		void set_padding(unsigned int pad);
		
	private:
	
		//! Storage location for the times supplied by calls to the write() function
		vector<double> _time;
		
		//! The *.vtu filename (this will be what is appended with time)
		FileParts _vtu;
		
		//! The *.pvd filename (this will reference the *.vtu files)
		FileParts _pvd;
		
		//! Reference to the libMesh::EquationSystems that contains the data
		const EquationSystems& _es;
		
		//! The number of digits to append to the filename of the *.vtu files
		unsigned int _pad;
		
		//! Creates the *.pvd file containing the list of *.vtu files and the correct times
		void write_pvd();
	
}; // class
} //namespace
#endif
