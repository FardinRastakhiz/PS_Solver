#pragma once


#ifndef SES_UTILITIES

#ifndef VIENNACL_WITH_OPENCL
#define VIENNACL_WITH_OPENCL
#endif

#include "viennacl/sliced_ell_matrix.hpp"
#include "viennacl/scalar.hpp"
#include "viennacl/vector.hpp"
#include "viennacl/coordinate_matrix.hpp"
#include "viennacl/compressed_matrix.hpp"
#include "petscksp.h"

#define CHECK_DOUBLE_TYPES(T, U) (numeric_limits<T>::is_specialized && numeric_limits<U>::is_specialized && numeric_limits<T>::is_double && numeric_limits<U>::is_integer && (numeric_limits<T>::is_signed == numeric_limits<U>::is_signed) && (numeric_limits<T>::digits == numeric_limits<U>::digits))


namespace ses {
	
	/*				Shared Types			*/
	// Other Types
#define USE_FLOAT 0
#if USE_FLOAT
	typedef float LocalType;
#else
	typedef double LocalType;
#endif

	// Available Matrix Types
	typedef viennacl::sliced_ell_matrix<LocalType>  VI_SELL_MAT;
	typedef viennacl::compressed_matrix<LocalType>  VI_COMP_Mat;
	typedef viennacl::coordinate_matrix<LocalType>  VI_COO_Mat;
	typedef Mat PETSC_MAT;
	// Available Vector Types
	typedef viennacl::vector<LocalType> VI_VEC;
	typedef Vec PETSC_VEC;


	/*				Save To Files			*/
	// Saving Vector!!
	template<typename T>
	void save_vector(viennacl::vector<T> x, const char* file_name);

	void save_vector_pointer(double* x, size_t size, const char* file_name);
	void save_vector_pointer(float* x, size_t size, const char* file_name);

	// implementing equivalent of "is" is c# for type checking
	template < class T, class U >
	bool is_inst_of(U u);


	/*				Converters				*/
	LocalType* cast_to_local(double* arr, size_t size);

	template <typename out_T>
	out_T* cast_to(LocalType* arr, size_t size);

}
#define SES_UTILITIES
#endif