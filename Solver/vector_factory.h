#pragma once

#include "utilities.h"
#ifndef SES_VEC_FACTORY


//#include "viennacl/forwards.h"
//#include "viennacl/vector.hpp"
//#include "viennacl/tools/tools.hpp"
//#include "viennacl/linalg/sparse_matrix_operations.hpp"
#include "viennacl/scalar.hpp"
#include "viennacl/vector.hpp"
#include <vector>
#include "petscksp.h"

namespace ses {
	void create_vector(int size, LocalType* values, VI_VEC& out_vec) {
		out_vec = VI_VEC(size);
		//std::vector<int, PreAllocator<int>> cpu_vec(0, PreAllocator<int>(&values2[0], size));
		std::vector<int> cpu_vec(size);
		for (size_t i = 0; i < size; i++)
		{
			cpu_vec[i] = values[i];
		}
		viennacl::copy(cpu_vec.begin(), cpu_vec.end(), out_vec.begin());
	}
	//void fill_petsc_vector(int size, LocalType* values, PETSC_VEC out_vec) {
	//	// Set rhs values
	//	for (int i = 0; i < size; i++) { VecSetValue(out_vec, i, values[i], INSERT_VALUES); }
	//	VecAssemblyBegin(out_vec);
	//	VecAssemblyEnd(out_vec);
	//}
	//void create_petsc_vector(int size, PETSC_VEC out_vec) {
	//	// Create vectors
	//	VecCreate(PETSC_COMM_WORLD, &out_vec);
	//	//VecSetType(b, VECVIENNACL);
	//	VecSetSizes(out_vec, PETSC_DECIDE, size);
	//	VecSetFromOptions(out_vec);
	//}
}

#define SES_VEC_FACTORY
#endif