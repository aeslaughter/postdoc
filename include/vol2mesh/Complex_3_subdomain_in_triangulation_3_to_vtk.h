/*! \file Complex_3_subdomain_in_triangulation_3_to_vtk.h
/*! 
 * This is a modified version of the CGAL fuction \c output_c3t3_to_vtk_unstructured_grid.h
 * that includes the subdomain index in the conversion to a vtk
 * unstructured grid.
 * 
 */

// Original Header:
//----------------------------------------------------------------------
// Copyright (c) 2008-2009  GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: svn+ssh://scm.gforge.inria.fr/svn/cgal/branches/next/Mesh_3/include/CGAL/IO/Complex_3_in_triangulation_3_to_vtk.h $
// $Id: Complex_3_in_triangulation_3_to_vtk.h 67117 2012-01-13 18:14:48Z lrineau $
//
//
// Author(s)     : Laurent Rineau
//----------------------------------------------------------------------

// Define the function if it has not been defined
#ifndef CGAL_COMPLEX_3_SUBDOMAIN_IN_TRIANGULATION_3_TO_VTK
#define CGAL_COMPLEX_3_SUBDOMAIN_IN_TRIANGULATION_3_TO_VTK

// Standard includes
#include <map>

// VTK includes
#include <vtkPoints.h>
#include <vtkUnstructuredGrid.h>
#include <vtkCellArray.h>
#include <vtkType.h>

#include <vtkCellData.h>
#include <vtkIntArray.h>

// Add this function to my namespace
namespace SlaughterVol2mesh {

//! A function for converting a CGAL mesh to VTK including the subdomain
/*! 
 * A template function that modifies the CGAL function (\c output_c3t3_to_vtk_unstructured_grid)
 * to include the exporting of the subdomain_index, which in the case
 * of the input files is the pixel value.
 * 
 * The header for the original CGAL file is included in the source code
 * for this file.
 * 
 * \param c3t3 CGAL::C3t3 mesh object
 * \param grid A pointer to the vtk grid; be default the grid is generated
 * 
 * \return A pointer to the vtkUnstructuredGrid object
 * 
 * \see v2m.cpp
 * \see Vol2mesh
 * \see Vol2mesh::generate_mesh
 */ 
template <typename C3T3>
vtkUnstructuredGrid* 
output_c3t3_subdomain_to_vtk_unstructured_grid(const C3T3& c3t3, 
                                     vtkUnstructuredGrid* grid = 0)
{
	typedef typename C3T3::Triangulation Triangulation;
	typedef typename Triangulation::Vertex_handle Vertex_handle;

	const Triangulation& tr = c3t3.triangulation();

	vtkPoints* const vtk_points = vtkPoints::New();
	vtkCellArray* const vtk_facets = vtkCellArray::New();
	vtkCellArray* const vtk_cells = vtkCellArray::New();

	// Create an array for storing the subdomain index
	vtkSmartPointer<vtkIntArray> blockIds = vtkSmartPointer<vtkIntArray>::New();
	blockIds->SetName("ElementBlockIds");

	vtk_points->Allocate(c3t3.triangulation().number_of_vertices());
	vtk_facets->Allocate(c3t3.number_of_facets_in_complex());
	vtk_cells->Allocate(c3t3.number_of_cells_in_complex());

	std::map<Vertex_handle, vtkIdType> V;
	vtkIdType inum = 0;

	for(typename Triangulation::Finite_vertices_iterator 
		vit = tr.finite_vertices_begin(),
		end = tr.finite_vertices_end();
	  vit != end;
	  ++vit)
	{
	typedef typename Triangulation::Point Point;
	const Point& p = vit->point();
	vtk_points->InsertNextPoint(CGAL::to_double(p.x()),
								CGAL::to_double(p.y()),
								CGAL::to_double(p.z()));
	V[vit] = inum++;
	}
	for(typename C3T3::Facets_in_complex_iterator 
		fit = c3t3.facets_in_complex_begin(),
		end = c3t3.facets_in_complex_end();
	  fit != end; ++fit) 
	{
	vtkIdType cell[3];
	int j=0;
	for (int i = 0; i < 4; ++i)
	  if (i != fit->second)
		cell[j++] =  V[(*fit).first->vertex(i)];
	CGAL_assertion(j==3);
	vtk_facets->InsertNextCell(3, cell);
	}
	for(typename C3T3::Cells_in_complex_iterator 
		cit = c3t3.cells_in_complex_begin(),
		end = c3t3.cells_in_complex_end();
	  cit != end; ++cit) 
	{
	  
	// Add the subdomain index to the vtk array
	blockIds->InsertNextValue( (int)c3t3.subdomain_index(cit) );	  
	  
	vtkIdType cell[4];
	for (int i = 0; i < 4; ++i)
	  cell[i] =  V[cit->vertex(i)];
	vtk_cells->InsertNextCell(4, cell);
	}

	if(!grid) {
	grid = vtkUnstructuredGrid::New();
	}

	grid->SetPoints(vtk_points);
	vtk_points->Delete();

	grid->SetCells(VTK_TRIANGLE, vtk_facets);
	grid->SetCells(VTK_TETRA, vtk_cells);
	
	// Add the subdomain index to the vtk grid
	grid->GetCellData()->AddArray(blockIds);

	vtk_facets->Delete();
	vtk_cells->Delete();

	return grid;
}

} // end namespace 
#endif // CGAL_COMPLEX_3_SUBDOMAIN_IN_TRIANGULATION_3_TO_VTK
