/*! \file my_vtkio.cpp
  * \brief Source code for MyVTKIO class.
  */
  
// Include the header
#include "fem/common/my_vtkio.h"
using namespace SlaughterFEM;
using namespace SlaughterCommon;
  
// Class constructor
MyVTKIO :: MyVTKIO(string vtu, EquationSystems& es) : VTKIO(es.get_mesh()), _vtu(vtu), _pvd(vtu), _es(es){
	
	// Check for the correct extension
	if (_vtu.ext.compare(".vtu") != 0){
		printf("ERROR: The supplied file must have a *.vtu extension.\n");
		libmesh_error();
	}
	
	// Set the default level of padding
	_pad = 3;
	
	// Create the correct file extensions for the reference file
	_pvd.ext = ".pvd";
	_pvd.update();
}

// Change how many numbers are appended to the name
void MyVTKIO :: set_padding(unsigned int pad){
	_pad = pad;
}

// Write the data to a file
void MyVTKIO :: write(double t){

	/* libMesh::VTKIO creates the following file structure:
	 * 		filename.vtu
	 *      filename_0.vtu
	 * filename.vtu referes to filename_0.vtu. 
	 * 
	 * Currently, this function renames filename_0.vtu to filename.vtu.
	 */ 

	// Append the time vector
	_time.push_back(t);
		
	// The desired filename
	FileParts newname(_vtu.add_tstep((int)_time.size()-1, _pad, "_"));
	
	// The filename with the libMesh added _0 extension
	FileParts oldname(newname);
	oldname.ext = "_0.vtu";
	oldname.update();
	
	// Write the *.vtu for the current time
	VTKIO::write_equation_systems(newname.full.c_str(), _es);

	// Replace the libMesh created file with the actual file
	rename(oldname.full.c_str(), newname.full.c_str());

	// Write the *.pvd that lists all of the files
	write_pvd();
}

void MyVTKIO :: write_pvd(){

	// Open the file
	FILE* fp = fopen(_pvd.full.c_str(), "w");	

	// Write the collection of files
	fprintf(fp, "<?xml version=\"1.0\"?>\n"); 	
	fprintf(fp, "\t<VTKFile type=\"Collection\" version=\"0.1\">\n"); 	
	fprintf(fp, "\t\t<Collection>\n");
	for (int i = 0; i < _time.size(); i++){
		FileParts fname(_vtu.add_tstep(i, _pad, "_"));
		fprintf(fp, "\t\t\t<DataSet timestep=\"%f\" part=\"0\" file=\"%s\" />\n", _time[i], 
			fname.name.append(fname.ext).c_str());
	}
	fprintf(fp, " \t\t</Collection>\n");
	fprintf(fp, "\t</VTKFile>\n");
	fclose(fp);
}
