/*! \file vol2mesh_mesh_criteria.cpp
 * Source code for the Vol2meshMeshCriteria class
 */

// Standard includes
#include <stdio.h>

// Header file for this class
#include "vol2mesh/vol2mesh_mesh_criteria.h" 

// Include the this class source code in namespace
namespace SlaughterVol2mesh{
				
// Default constructor that sets all values to zero				
Vol2meshMeshCriteria :: Vol2meshMeshCriteria(){
	vector<double> vec(5,0);
	init(vec);			
}					

// Constructor that sets the values equal to the supplied vector
Vol2meshMeshCriteria :: Vol2meshMeshCriteria(vector<double> vec){
	init(vec);
}	
	
// Constructor with individually specified values	
Vol2meshMeshCriteria :: Vol2meshMeshCriteria(double a, double b, double c, double d, double e){
	vector<double> vec;
	vec.push_back(a);
	vec.push_back(b);
	vec.push_back(c);
	vec.push_back(d);
	vec.push_back(e);
	init(vec);
}		
				
// Initialization function the places vector components into the values				
void Vol2meshMeshCriteria :: init(vector<double> vec){

	// Report error if vector is not the correct size
	if (vec.size() != 5){
		printf("ERROR: The supplied vector must contain 5 values, but %d were supplied.", (int)vec.size());
		return;
	}

	// Insert the values
	facet_angle = vec[0];
	facet_size = vec[1];
	facet_distance = vec[2];
	cell_radius_edge_ratio = vec[3];
	cell_size = vec[4];
}
	
// Return the mesh criteria in a vector format
vector<double> Vol2meshMeshCriteria :: get_vector(){
	vector<double> vec;
	vec.push_back(facet_angle);
	vec.push_back(facet_size);
	vec.push_back(facet_distance);
	vec.push_back(cell_radius_edge_ratio);
	vec.push_back(cell_size);
	return vec;
}

// Return the i-th value (values ordered as in init)
double Vol2meshMeshCriteria :: get_value(int i){
	
	double value;
	if (i == 0){
		value = facet_angle; 
	
	} else if (i == 1){ 
		value = facet_size; 
	
	} else if (i == 2){
		value = facet_distance;
	
	} else if (i == 3){
		value = cell_radius_edge_ratio;

	} else if (i == 4){
		value = cell_size;
		
	} else {
		printf("ERROR: The index %d is out of range, only 0 through 4 accepted.", i);
	}
	
	return value;
}

} // end namespace
