/*! \file vol2mesh.cpp
 * Source code for a class, \c Vol2mesh, for creating a 3D mesh from 2D
 * images slices. The program is utilzes the CGAL, VTK, and ImageMagick++
 * libraries. The basis for the development comes from the CGAL example: 
 * Section 50.3.3 Domains From Segmented 3D Images. 
 * 
 * This use of executable associated with this source code is detailed 
 * here: \ref v2m.
 */
 
// Include the necessary headers
#include "vol2mesh/vol2mesh.h"

// Include this code in SlaughterVol2Mesh namespace
namespace SlaughterVol2mesh{ 
 
// Class constructor for string input										   
Vol2mesh :: Vol2mesh(string infile) : input_file_(infile){
	init();
}

// Initilization function
void Vol2mesh :: init(){

	// Set the default voxel dimensions
	voxels_.assign(3,1);

	// Read the image
	read_image();
	
	// Set the image pixel dimensions
	pixels_.push_back(image_.xdim());
    pixels_.push_back(image_.ydim());
	pixels_.push_back(image_.zdim());
	
	// Set the default mesh criteria
	Vol2meshMeshCriteria c(30, 6, 4, 3, 8);
	default_criteria_ = c;	
	
	// Set the normalize value to false
	normalize_ = false;
	
	// Set the default optimization behavior
	use_lloyd_ = false;
	use_odt_ = false;
	use_perturb_ = true;
	use_exude_ = true;
	
	// Set the default optimization settings (CGAL defaults)
	lloyd 	= {0, 0, 0.02, 0.01};
	odt 	= {0, 0, 0.02, 0.01};
	perturb = {0, 0};
	exude   = {0, 0};
	
	// A default output file
	output_file_ = "output.ex2";	
	
};

// Function for reading an image file and storing as a CGAL::Image_3 object (image_)
void Vol2mesh :: read_image(){

    // Create File_parts object for file names    
    FileParts infile(input_file_);

    // Loads a *.inr or *.inf.gz image
    if (infile.ext.compare(".inr") == 0 || infile.ext.compare(".inr.gz") == 0){
        image_.read(infile.full.c_str());
        
    // Load a *.tiff or *.tif image
    } else if (infile.ext.compare(".tif") == 0 ||  infile.ext.compare(".tiff") == 0){
        read_tiff();
       // printf("TIFF files are not yet support, future suport is planned.\n");
       // exit(100);
        
    // Returns an error if the file type is not understood    
    } else {
        printf("ERROR: The input file type (*%s) is not supported (see vol2mesh.cpp)\n", infile.ext.c_str());
        exit(101);
    }
}						   

// Sets the output file
void Vol2mesh :: set_output_file(string output_file){
	output_file_.assign(output_file);
};	
	
// Return a vector of the pixel dimensions	
vector<int> Vol2mesh :: pixels(){
    return pixels_; 
}

// Return a vector of the voxel dimensions
vector<double> Vol2mesh :: voxels(){
	return voxels_;
}

// Set the normalize flag
void Vol2mesh :: normalize(bool value){
	normalize_ = value;
}

// Set the disable_subdomain flag
void Vol2mesh :: disable_subdomains(bool value){
	disable_subdomains_ = value;
}

// Set the optimization flags
void Vol2mesh :: set_optimization(string type, bool value){
	
	if(type.compare("lloyd") == 0){
		use_lloyd_ = value;
	
	} else if (type.compare("odt") == 0) {
		use_odt_ = value;
		
	} else if (type.compare("perturb") == 0) {
		use_perturb_ = value;
		
	} else if (type.compare("exude") == 0) {
		use_exude_ = value;
	
	} else {
		printf("ERROR: Optimiatation type %s was not recongnized", type.c_str());
	}
}

// Retun the value for normalizing the mesh criteria
double Vol2mesh :: get_normalize_value(){
	double N = 1;
	if (normalize_){
		vector<double> v = voxels();
		N = (v[0] + v[1] + v[2]) / 3;
	}
	return N;	
}

// Set the voxel dimension (single value)
void Vol2mesh :: set_voxel(int i, double vox){
	
	// Update the storage vector
	voxels_.at(i) = vox;
	
    // Update the image
    set_voxel(voxels_);
}

// Set the voxel dimensions (vector)
void Vol2mesh :: set_voxel(vector<double> vox){    
    
    // Update the storage vector
    voxels_ = vox;
    
    // Update the actual CGAL image, note the _image type from CGAL
    // is not related to the private _image value used in this class
    _image* img = image_.image();
    img->vx = vox[0];
    img->vy = vox[1];
    img->vz = vox[2];	
}

// Set the voxel dimension using overall image dimension (single value)
void Vol2mesh :: set_dimension(int i, double dim){
	
	// Extract the number of pixels in the image
	vector<int> pix = pixels();
	
	// Update the image with computed voxel size
	set_voxel(i, dim / pix[i]);
}

// Set the voxel dimensions using a vector of the overall image dimensions
void Vol2mesh :: set_dimension(vector<double> dim){
	
	// Define storage location for calculated voxel size
	vector<double> vox;	
	
	// Extract the number of pixels in the image
	vector<int> pix = pixels();
	
	// Computed voxel size
	for (int i = 0; i < 3; i++){
		vox.push_back(dim[i] / pix[i]);
	}
	
	// Update the image
	set_voxel(vox);
}

// Set the default meshing criteria
void Vol2mesh :: set_default_criteria(Vol2meshMeshCriteria c){
	default_criteria_ = c;
}

// Return the default mesh criteria	
Vol2meshMeshCriteria& Vol2mesh :: get_default_criteria(){
	return default_criteria_;
}	

// Add a subdomain id and associate mesh criteria
void Vol2mesh :: add_subdomain(int id, Vol2meshMeshCriteria c){
	
	// Test for existance of domain
	for (int i = 0; i < subdomain_criteria_.size(); i++){
		if( id == subdomain_criteria_[i].id){
			printf("ERROR: The subdomain id, %d, was already specified!",id);
			exit(1);
		}
	}
	
	// Set the id of object and append the storage vector
	c.id = id;
	subdomain_criteria_.push_back(c);
}

// Return a reference to the subdomain criteria
vector<Vol2meshMeshCriteria>& Vol2mesh :: get_subdomain_criteria(){
	return subdomain_criteria_;
}

// Generate the CGAL mesh
void Vol2mesh :: generate_mesh(){

	// Create CGAL domain object
	Mesh_domain domain(image_); 

	// Get a vector of the default mesh critiera
	vector<double> dmc_vec = default_criteria_.get_vector();

	// Compute the normalizing value
	double N = get_normalize_value();

	// Create a smart pointer to a Mesh_criteria
	boost::shared_ptr<Mesh_criteria> p_criteria;

	// Create CGAL Mesh_critiera object, ignoring subdomains
	if (disable_subdomains_ || subdomain_criteria_.empty()){	

		// Set the mesh criteria to the defaults
		boost::shared_ptr<Mesh_criteria> p(new Mesh_criteria(
		CGAL::parameters::facet_angle = dmc_vec[0] / N, 
        CGAL::parameters::facet_size = dmc_vec[1] / N, 
        CGAL::parameters::facet_distance = dmc_vec[2] / N,
        CGAL::parameters::cell_radius_edge_ratio = dmc_vec[3], 
        CGAL::parameters::cell_size = dmc_vec[4] / N));
        
        // Copy the local pointer
        p_criteria = p;

	 // Create CGAL Mesh_critiera object, with subdomains
     } else {
		// Create CGAL sizing fields, specifing the default at creation
		vector<Sizing_field> fields;
		for (int i = 0; i < dmc_vec.size(); i++){

			// Create sizing field
			Sizing_field f(dmc_vec[i] / N);
			
			// Apply the subdomain criteria to the field
			for (int j = 0; j < subdomain_criteria_.size(); j++){
				f.set_size(subdomain_criteria_[j].get_value(i), 3, 
					domain.index_from_subdomain_index(subdomain_criteria_[j].id));
			}
			
			// Store the field in a vector
			fields.push_back(f);
		}
		
		//Set the meshing criteria, based on command-line input
		boost::shared_ptr<Mesh_criteria> p(new Mesh_criteria(
			CGAL::parameters::facet_angle = dmc_vec[0], 
			CGAL::parameters::facet_size = dmc_vec[1], 
			CGAL::parameters::facet_distance = dmc_vec[2],
			CGAL::parameters::cell_radius_edge_ratio = dmc_vec[3],
			CGAL::parameters::cell_size = fields[4]));
			
        // Copy the local pointer
        p_criteria = p;		
	}

	// Create the C3t3 mesh object w/o optimization
	c3t3_ = CGAL::make_mesh_3<C3t3>(domain, *p_criteria, 
			CGAL::parameters::no_lloyd(), 
			CGAL::parameters::no_odt(), 
			CGAL::parameters::no_perturb(), 
			CGAL::parameters::no_exude());
		
	// Apply the Lloyd optimization, if desired
	if (use_lloyd_){
		CGAL::refine_mesh_3(c3t3_, domain, *p_criteria, 
			CGAL::parameters::no_odt(), 
			CGAL::parameters::no_perturb(), 
			CGAL::parameters::no_exude(),
			CGAL::parameters::lloyd(
				CGAL::parameters::time_limit = lloyd.time_limit,
				CGAL::parameters::max_iteration_number = lloyd.max_iteration,
				CGAL::parameters::convergence = lloyd.convergence,
				CGAL::parameters::freeze_bound = lloyd.freeze_bound));
	}

	// Apply the Odt optimization, if desired
	if (use_odt_){
		CGAL::refine_mesh_3(c3t3_, domain, *p_criteria, 
			CGAL::parameters::no_lloyd(), 
			CGAL::parameters::no_perturb(), 
			CGAL::parameters::no_exude(),
			CGAL::parameters::odt(
				CGAL::parameters::time_limit = odt.time_limit,
				CGAL::parameters::max_iteration_number = odt.max_iteration,
				CGAL::parameters::convergence = odt.convergence,
				CGAL::parameters::freeze_bound = odt.freeze_bound));
	}

	// Apply the perturb, if desired
	if (use_perturb_){
		CGAL::refine_mesh_3(c3t3_, domain, *p_criteria, 
			CGAL::parameters::no_odt(), 
			CGAL::parameters::no_lloyd(), 
			CGAL::parameters::no_exude(),
			CGAL::parameters::perturb(
				CGAL::parameters::time_limit = perturb.time_limit,
				CGAL::parameters::sliver_bound = perturb.sliver_bound));
	}

	// Apply the exude, if desired
	if (use_exude_){
		CGAL::refine_mesh_3(c3t3_, domain, *p_criteria, 
			CGAL::parameters::no_odt(), 
			CGAL::parameters::no_lloyd(), 
			CGAL::parameters::no_perturb(),
			CGAL::parameters::exude(
				CGAL::parameters::time_limit = exude.time_limit,
				CGAL::parameters::sliver_bound = exude.sliver_bound));
	}
}

// Return a reference to the CGAL mesh
C3t3& Vol2mesh :: c3t3(){
	return c3t3_;
}

// Return a reference to the CGAL image
CGAL::Image_3& Vol2mesh :: image(){
	return image_;
}

// Write the CGAL mesh to a file specified
void Vol2mesh :: write(string output_file){
	
	// Create file parts object of the output file
    FileParts outfile(output_file);
    
    // Store the output file in the private member
    output_file_.assign(output_file);

	// Ouptut the file in *.mesh format
    if (outfile.ext.compare(".mesh") == 0){		
        std::ofstream medit_file(output_file.c_str());
        c3t3_.output_to_medit(medit_file);
    
    // Output the file in *.vtu, *.vtk, or *.ex2 format (uses VTK library)
    } else {
        
        // Convert CGAL object to a vtkUnstructuredGrid (http://cgal-discuss.949826.n4.nabble.com/mesh-to-vtk-output-td3586974.html)
        vtkUnstructuredGrid *output; //vtk object to convert into
        
        // Do not include sub-domains if disabled
        if (disable_subdomains_){
			output = CGAL::output_c3t3_to_vtk_unstructured_grid(c3t3_); 
		} else {
			output = output_c3t3_subdomain_to_vtk_unstructured_grid(c3t3_);  
		}
        output->Squeeze();

		// Ouput in the *.vtu format
        if (outfile.ext.compare(".vtu") == 0){
            vtkXMLDataSetWriter *vtu = vtkXMLDataSetWriter::New(); 
            vtu->SetInput(output);
            vtu->SetFileName(output_file.c_str());
            vtu->Write(); 
            
        // Output in the *.vtk format    
        } else if (outfile.ext.compare(".vtk") == 0){
            vtkUnstructuredGridWriter *vtk = vtkUnstructuredGridWriter::New();
            vtk->SetInput(output);
            vtk->SetFileName(output_file.c_str());
            vtk->Write();
        
        // Output in the ExodusII (*.ex2) format    
        } else if (outfile.ext.compare(".ex2") == 0){  
			vtkExodusIIWriter *exII = vtkExodusIIWriter::New();	
			exII->SetFileName(output_file.c_str());
            exII->SetInput(output);        
            exII->Write();

        // Produce an error if the file extension is not understood
        } else {
            printf("ERROR: The desired output file type (*.%s) is not supported (see vol2mesh.cpp)\n", outfile.ext.c_str());
            exit(102);
        }
    }
}

// Write the CGAL mesh with previously defined filename
void Vol2mesh :: write(){
	if(!output_file_.empty()){
		write(output_file_);
	} else {
		printf("ERROR: An output file was not specified\n");
	}
}	

// A function for printing the results of the mesh to a file and/or the screen
void Vol2mesh :: print_results(double t, bool disable_screen, bool enable_file){
       
	// Convert CGAL object to a vtkUnstructuredGrid 
	// (http://cgal-discuss.949826.n4.nabble.com/mesh-to-vtk-output-td3586974.html)
	vtkUnstructuredGrid *uGrid;
	uGrid = CGAL::output_c3t3_to_vtk_unstructured_grid(c3t3_); 
	uGrid->Squeeze();
              
    // Compute mesh quality information 
    vtkMeshQuality *q = vtkMeshQuality::New();
    q->SetInput(uGrid);
    
    // Variables for storing quality statistics
    matrix q_mat;
    vector<string> q_name;
    
    // Gather statistics for the mesh
    get_all_quality_stats(q, q_mat, q_name);
    
    // Define variables for building the vector or strings for display    
    int w = 85; // must be greater than 85
    char c[w];    
    vector<string> s;
    string tmp;
    
    // File information header
    tmp.assign("FILE INFORMATION ");
    tmp.append(w - tmp.size() - 1, '-');
    tmp.append("\n");
    s.push_back(tmp);
    
    // File input and output names
    sprintf(c, " %12s: %s\n", "input-file", input_file_.c_str());         
		s.push_back(c);
    sprintf(c, " %12s: %s\n\n", "output-file", output_file_.c_str());     
		s.push_back(c);
   
    // Input paramaters header
    tmp.assign("DEFAULT MESH CRITERIA ");
    tmp.append(w - tmp.size() - 1, '-');
    tmp.append("\n");
    s.push_back(tmp);
 
    // User suplied options
    sprintf(c, " %23s: %6.3f\n", "facet-angle", default_criteria_.facet_angle);            
		s.push_back(c);
    sprintf(c, " %23s: %6.3f\n", "facet-size", default_criteria_.facet_size);              
		s.push_back(c);
    sprintf(c, " %23s: %6.3f\n", "facet-distance", default_criteria_.facet_distance);      
		s.push_back(c);
    sprintf(c, " %23s: %6.3f\n", "cell-radius-edge-ratio", 
		default_criteria_.cell_radius_edge_ratio);      
        s.push_back(c);
    sprintf(c, " %23s: %6.3f\n\n", "cell-size", default_criteria_.cell_size); 
		s.push_back(c);

	// Mesh results header
    tmp.assign("MESH RESULTS ");
    tmp.append(w - tmp.size() - 1, '-');
    tmp.append("\n");
    s.push_back(tmp);
     
    // Mesh results
    vector<int> pix = pixels();
    vector<double> vox = voxels();

    sprintf(c, " %23s: %6.3f\n", "execution time (sec.)", t);            
		s.push_back(c);
	sprintf(c, " %23s: %d, %d, %d\n", "num. of pixels (x,y,z)",
		pix[0], pix[1], pix[2]);
		s.push_back(c);	
	sprintf(c, " %23s: %6.3f, %6.3f, %6.3f\n", "pixel dim. (x,y,z)",
		vox[0], vox[1], vox[2]);
		s.push_back(c);	
	sprintf(c, " %23s: %6.3f, %6.3f, %6.3f\n", "image dim. (x,y,z)",
		pix[0]*vox[0], pix[1]*vox[1], pix[2]*vox[2]);
		s.push_back(c);				
    sprintf(c, " %23s: %d\n", "num. of elements", (int)c3t3_.number_of_cells()); 
		s.push_back(c);
	sprintf(c, " %23s: %d\n\n", "num. of faces", (int)c3t3_.number_of_facets()); 
		s.push_back(c);
		
	// Mesh quality header
	tmp.assign("TETRAHEDRAL QUALITY ");
    tmp.append(w - tmp.size() - 1, '-');
    tmp.append("\n");
    s.push_back(tmp);			
		
	// Print the mesh quality table labels
	sprintf(c,"%24s%10s%10s%10s%10s%10s\n", 
		"Name", "Lower", "Upper", "Average", "Std. dev.", "COV (%)");
		s.push_back(c);
		
	// Print each of the mesh quality results			
	for(int i = 0; i < q_name.size(); ++i){
		sprintf(c,"%24s%10.3f%10.3f%10.3f%10.3f%10.3f\n", q_name[i].c_str(), 
			q_mat[i][0], q_mat[i][1], q_mat[i][2], q_mat[i][3],
			q_mat[i][3] / q_mat[i][2] * 100);
		s.push_back(c);		
	}

	// Add sub-domain mesh criteria
	if(!subdomain_criteria_.empty()){
		for(int i = 0; i < subdomain_criteria_.size(); i++){
			// Input paramaters header
			sprintf(c, "SUBDOMAIN %d: MESH CRITERIA ", subdomain_criteria_[i].id);
			tmp.assign(c);
			tmp.append(w - tmp.size() - 1, '-');
			tmp.append("\n");
			s.push_back(tmp);

			// User suplied options
			sprintf(c, " %23s: %6.3f\n", "facet-angle", subdomain_criteria_[i].facet_angle);            
				s.push_back(c);
			sprintf(c, " %23s: %6.3f\n", "facet-size", subdomain_criteria_[i].facet_size);              
				s.push_back(c);
			sprintf(c, " %23s: %6.3f\n", "facet-distance", subdomain_criteria_[i].facet_distance);      
				s.push_back(c);
			sprintf(c, " %23s: %6.3f\n", "cell-radius-edge-ratio", 
				subdomain_criteria_[i].cell_radius_edge_ratio);      
				s.push_back(c);
			sprintf(c, " %23s: %6.3f\n\n", "cell-size", subdomain_criteria_[i].cell_size); 
				s.push_back(c);	
		}
	}

	// Output the message to the screen
    if (!disable_screen){
        printf("\n\n");
        for (int i = 0; i < s.size(); ++i){
            printf("%s", s[i].c_str());
        }
        printf("\n\n");
    }

	// Output the message to a file
    if (enable_file){
        string hdr_file;
        hdr_file = output_file_ + (string)".info";
        FILE* fid = fopen(hdr_file.c_str(), "w");
        for (int i = 0; i < s.size(); ++i){
            fprintf(fid, "%s", s[i].c_str());
        }
    }
}

// Returns a vector containing the vtk mesh quality statistics
vector<double> Vol2mesh :: get_current_stats(vtkMeshQuality* q, const char *name){
	vector<double> v;
	v.push_back(q->GetOutput()->GetFieldData()->GetArray(name)->GetComponent( 0, 0 ));
	v.push_back(q->GetOutput()->GetFieldData()->GetArray(name)->GetComponent( 0, 2 ));
	v.push_back(q->GetOutput()->GetFieldData()->GetArray(name)->GetComponent( 0, 1 ));
	v.push_back(q->GetOutput()->GetFieldData()->GetArray(name)->GetComponent( 0, 3 ));
	return v;
}

// Computes the quality statistics for the tetrahedral mesh
void Vol2mesh :: get_all_quality_stats(vtkMeshQuality* q, matrix & q_mat, vector<string> & q_name){

    // Edge ratio
    q->SetTetQualityMeasureToEdgeRatio();
    q->Update();
    q_name.push_back("Edge ratio");
    q_mat.push_back(get_current_stats(q));

	// Aspect ratio
    q->SetTetQualityMeasureToAspectRatio();
    q->Update();
    q_name.push_back("Aspect ratio");
    q_mat.push_back(get_current_stats(q));

    // Radius ratio
    q->SetTetQualityMeasureToRadiusRatio();
    q->Update();
    q_name.push_back("Radius ratio");
    q_mat.push_back(get_current_stats(q));
    
    // Aspect Frobenius
    q->SetTetQualityMeasureToAspectFrobenius();
    q->Update();
    q_name.push_back("Aspect Frobenius");
    q_mat.push_back(get_current_stats(q));
    
    // Minimal dihedral angle
    q->SetTetQualityMeasureToMinAngle();
    q->Update();
    q_name.push_back("Minimal dihedral angle");
    q_mat.push_back(get_current_stats(q));
    
    // Collapse ratio
    q->SetTetQualityMeasureToCollapseRatio();
    q->Update();
    q_name.push_back("Collapse ratio");
    q_mat.push_back(get_current_stats(q));
    
    // Aspect beta
    q->SetTetQualityMeasureToAspectBeta();
    q->Update();
    q_name.push_back("Aspect beta");
    q_mat.push_back(get_current_stats(q));

    // Volume
    q->SetTetQualityMeasureToVolume();
    q->Update();
    q_name.push_back("Volume");
    q_mat.push_back(get_current_stats(q));
    
    // Condition
    q->SetTetQualityMeasureToCondition();
    q->Update();
    q_name.push_back("Condition");
    q_mat.push_back(get_current_stats(q));
    
    // Jacobian
    q->SetTetQualityMeasureToJacobian();
    q->Update();
    q_name.push_back("Jacobian");
    q_mat.push_back(get_current_stats(q));
    
    // Scaled jacobian
    q->SetTetQualityMeasureToScaledJacobian();
    q->Update();
    q_name.push_back("Scaled jacobian");
    q_mat.push_back(get_current_stats(q));
    
    // Shape
    q->SetTetQualityMeasureToShape();
    q->Update();
    q_name.push_back("Shape");
    q_mat.push_back(get_current_stats(q));

    // Relative size squared
    q->SetTetQualityMeasureToRelativeSizeSquared();
    q->Update();
    q_name.push_back("Relative size squared");
    q_mat.push_back(get_current_stats(q));
    
    // Shape and size
    q->SetTetQualityMeasureToShapeAndSize();
    q->Update();
    q_name.push_back("Shape and size");
    q_mat.push_back(get_current_stats(q));
    
    // Distortion
    q->SetTetQualityMeasureToDistortion();
    q->Update();
    q_name.push_back("Distortion");
    q_mat.push_back(get_current_stats(q));
}

// This function is not working correctly
void Vol2mesh :: read_tiff(){

	// Typedefs for ITK image and associated reader
	typedef unsigned char PixelType;
	typedef itk::Image< PixelType, 3 > ImageType;
	typedef itk::ImageFileReader< ImageType > ReaderType;
	ImageType::IndexType pixelIndex;

	// Create a new reader to be accessed with a pointer
	ReaderType::Pointer reader = ReaderType::New();
	
	// Set the image input file to the string stored in this class
	reader->SetFileName(input_file_.c_str());

	// Set the ITK image pointer to the output of the reader
	ImageType::Pointer itk_image = reader->GetOutput();
	
	// Update the ITK image
	itk_image->Update();
	itk_image->UpdateOutputData();
	itk_image->UpdateOutputInformation();

	// Extract the number of dimensions and compents
	int ndim = itk_image->GetImageDimension();
	int nv = itk_image->GetNumberOfComponentsPerPixel();
	
	// Extract the dimensions in each direction
	vector<int> dim;
	for (int i = 0; i < ndim; i++){
		dim.push_back(itk_image->GetLargestPossibleRegion().GetSize()[i]);
	}

	// Open my output file
	string tmp("vol2mesh.tmp");
	FILE* fid;
	fid = fopen(tmp.c_str(),"w");

	// Loop through each pixel, extract the value and write to the file
	for (int z = 0; z < dim[2]; z++){
		for (int y = 0; y < dim[1]; y++){
			for (int x = 0; x < dim[0]; x++){
				
				// Extracts the pixel at the x,y,z location
				pixelIndex[0] = x;
				pixelIndex[1] = y;
				pixelIndex[2] = z;
				PixelType pixelValue = itk_image->GetPixel(pixelIndex);				
			
				// Writes the value to the file
				fprintf(fid, "%c", pixelValue);
			}
		}
	}

	// Closes my output file
	fclose(fid);

	// Loads image into the CGAL::Image_3 variable of this class
	image_.read_raw(tmp.c_str(), dim[0], dim[1], dim[2], 1, 1, 1);
	
	// Remove the temporary file
	remove(tmp.c_str());
}

} //namespace SlaughterVol2Mesh

