//#pragma once
//
//#ifndef VIENNACL_WITH_OPENCL
//#define VIENNACL_WITH_OPENCL
//#endif
//
//#include "viennacl/sliced_ell_matrix.hpp"
//#include "viennacl/vector.hpp"
//#include "viennacl/coordinate_matrix.hpp"
//#include "viennacl/compressed_matrix.hpp"
//
//#define CHECK_DOUBLE_TYPES(T, U) (numeric_limits<T>::is_specialized && numeric_limits<U>::is_specialized && numeric_limits<T>::is_double && numeric_limits<U>::is_integer && (numeric_limits<T>::is_signed == numeric_limits<U>::is_signed) && (numeric_limits<T>::digits == numeric_limits<U>::digits))
//
//namespace ses {
//	// Other Types
//	typedef double LocalType;
//
//	// Available Matrix Types
//	typedef viennacl::sliced_ell_matrix<LocalType>  VI_SELL_MAT;
//	typedef viennacl::compressed_matrix<LocalType>  VI_COMP_Mat;
//	typedef viennacl::coordinate_matrix<LocalType>  VI_COO_Mat;
//
//	// Available Vector Types
//	typedef viennacl::vector<LocalType> VI_VEC;
//
//
//}