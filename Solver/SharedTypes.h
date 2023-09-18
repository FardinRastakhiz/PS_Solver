#pragma once

#include "viennacl/sliced_ell_matrix.hpp"
#include "viennacl/vector.hpp"
#include "viennacl/coordinate_matrix.hpp"
#include "viennacl/compressed_matrix.hpp"

namespace ses {
	// Other Types
	typedef float ScalarType;

	// Available Matrix Types
	typedef viennacl::sliced_ell_matrix<ScalarType>  VI_SELL_MAT;
	typedef viennacl::compressed_matrix<ScalarType>  VI_COMP_Mat;
	typedef viennacl::coordinate_matrix<ScalarType>  VI_COO_Mat;

	// Available Vector Types
	typedef viennacl::vector<ScalarType> VI_VEC;
}