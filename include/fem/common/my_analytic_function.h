// An libMesh::AnalyticFunction using boost::functions as input. The 
// code presented here is modeled after the libMesh library AnalyticFunction
// class in the analytic_function.h file.

// Copyright (C) 2012 Andrew E Slaughter

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

#ifndef my_analytic_function_h
#define my_analytic_function_h

#include <boost/function.hpp>
#include <boost/bind.hpp>

#include "libmesh.h"
#include "dense_vector.h"
#include "point.h"
#include "function_base.h"
using namespace libMesh;

namespace SlaughterFEM{

//! A class for using boost::functions with libmesh.
/*!
 * This class mimics the behavior of libMesh::AnalyticFunction but 
 * utilizes the boost::function behavior, as such both class members 
 * and functions can be used.
 * 
 * \ingroup FEMcommon
 */
template <typename Output = Number>
class MyAnalyticFunction : public FunctionBase<Output>{
	public:

		//! Class constructor for functions with scalar return values
		/*! 
		 * \param fptr A boost::function object that acts like a function pointer, this
		 * version outputs the type defined by the template parameter and 
		 * takes a libMesh Point and Real (time) as input.
		 */
		MyAnalyticFunction(boost::function< Output (const Point&, const Real) > fptr) :
			FunctionBase<Output> (),
			_number_fptr (fptr),
			_vector_fptr (NULL){
				
				libmesh_assert (fptr != NULL);
				this->_initialized = true;
		}

		//! Class constructor for functions with vector output, void return
		/*! 
		 * \param fptr A boost::function object that acts like a function pointer,
		 * this version output values through the DenseVector reference
		 * input (the type is defined as the template parameter) that is
		 * computed from the libMesh Point and Real (time) inputs.
		 */   		       
		MyAnalyticFunction(boost::function< void (DenseVector<Output>&, 
			const Point&, const Real)> fptr) :
			FunctionBase<Output> (),
			_number_fptr (NULL),
			_vector_fptr (fptr){
		
				libmesh_assert (fptr != NULL);
				this->_initialized = true;
		}			      
   		      
		//! Class destructor
		~MyAnalyticFunction (){};

		/*! Boost::function that points to a user provided class or function
		 * that computes the boundary values and has a scalar output.
		 */
		boost::function< Output (const Point&, const Real) > _number_fptr;

		/*! Boost::function that points to a user provided class or function
		 * that computes the boundary values and outputs a vector via 
		 * the DenseVector input reference.
		 */
		boost::function< void (DenseVector<Output>&, const Point&, const Real) > _vector_fptr;

		//! An initialization function
		/*! 
		 * This checks if either the _number_fptr or the _vector_fptr
		 * is valid, if so the class is marked as intilized.
		 */
		void init (){
			// dumb double-test
			libmesh_assert ((_number_fptr != NULL) || (_vector_fptr != NULL));

			// definitely ready
			this->_initialized = true;
		}

		//! Clears the function.
		/*! 
		 * This sets the two boost::function objects to NULL, thus 
		 * the function is set to be un-initilized.
		 */
		void clear (){
			// We probably need a method to reset these later...
			_number_fptr = NULL;
			_vector_fptr = NULL;

			// definitely not ready
			this->_initialized = false;
		}

		//! Returns a new deep copy of the function.
		/*! 
		 * \return A new copy of the object, the new object is 
		 * created with either the _number_fptr or _vector_ptr depending
		 * on which is initilized.
		 */
		virtual AutoPtr<FunctionBase<Output> > clone () const{
		
			return AutoPtr<FunctionBase<Output> > ( _number_fptr ?
				new MyAnalyticFunction<Output>(_number_fptr) :
				new MyAnalyticFunction<Output>(_vector_fptr) );
		}
		
		//! Allows the class to behave like a function
		/*! 
		 * /param p A libMesh point
		 * /param t The time
		 * /return The value of the function to which this class points
		 * at point \p p and time \p t (defaults to zero).
		 */ 
		Output operator() (const Point& p, const Real t = 0.){
			
			libmesh_assert (this->initialized());
			
			return (this->_number_fptr(p, t));
		}

		//! Allows the class to behave like a function
		/* 
		 * /param p A libMesh point
		 * /param t The time
		 * /param output A vector reference, the value of the function 
		 * to which this class points at point \p p and time \p t is 
		 * returned using this term.
		 */ 
		void operator() (const Point& p, const Real t, DenseVector<Output>& output){
			
			libmesh_assert (this->initialized());
			
			this->_vector_fptr(output, p, t);
			
			return;
		}
		
		
		
		
};

} // namespace SlaughterFEM

#endif

