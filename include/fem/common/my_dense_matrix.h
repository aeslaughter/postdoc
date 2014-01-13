 /*! \file my_dense_matrix.h
  * \brief Header file for MyDenseMatrix class.
  */
  
// Avoid multiple includes
#ifndef my_dense_matrix_h
#define my_dense_matrix_h

// Boost includes
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>
using namespace boost::numeric;

// libMesh includes
#include <libmesh.h>
#include <libmesh_common.h>
#include <dense_matrix.h>
using namespace libMesh;

// Add to my FEM namespace
namespace SlaughterFEM{

/*! \class MyDenseMatrix my_dense_matrix.h "fem/common/my_dense_matrix.h"
 * \brief A class that adds the capability of computing the inverse of 
 * a DenseMatrix using Boost::uBlas 
 * \see LevelSetSystem
 * \tparam The variable type, e.g. double
 * \ingroup FEMcommon
 */ 
template <typename Type> class MyDenseMatrix : public DenseMatrix<Type>{
	
	public:
	
		//! Constructor
		/*!
		 * Calls the libMesh::DenseMatrix constructor
		 */ 
		MyDenseMatrix(const unsigned int m = 0, const unsigned int n = 0) : DenseMatrix<Type>(m, n){};
	
		//! Compute the inverse of the matrix
		/*!
		 * Utilize the boost::numeric::ublas functions to invert the matrix,
		 * the inverted matrix replaces the existing matrix.
		 * 
		 * This function creates a copy of the original matrix and another 
		 * matrix of the solution, so it should not be used for large 
		 * matrices where memory may be an issue.
		 */ 
		void inverse(){
		
			// Check that the matrix is square
			libmesh_assert(this->m() == this->n());
			
			// Create a copy of the DenseMatrix and a matrix for the solution
			ublas::matrix<Type> A(this->m(), this->n());
			ublas::matrix<Type> B(this->m(), this->n());
			
			// Copy the dense matrix to the ublas::matrix
			for (unsigned int i = 0; i < this->m(); i++){
				for (unsigned int j = 0; j < this->n(); j++){
					A(i,j) = this->el(i,j);
				}
			}

			// Create a permutation matrix for the LU-factorization
			ublas::permutation_matrix<std::size_t> pm(A.size1());
			
			// Perform LU-factorization
			int res = lu_factorize(A, pm);
			if (res != 0){
				libmesh_error();
			}
			
			// Create identity matrix of "inverse"
			B.assign(ublas::identity_matrix<Type> (A.size1()));

			// Backsubstitute to get the inverse
			ublas::lu_substitute(A, pm, B);
			
			// Copy the solution to the current DenseMatrix
			for (unsigned int i = 0; i < this->m(); i++){
				for (unsigned int j = 0; j < this->n(); j++){
					this->el(i,j) = B(i,j);
				}
			}
		};
	
}; // class
} // namespace
#endif
