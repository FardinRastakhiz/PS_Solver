#include "pch.h"
#include "vector_factory.h"

//#include "viennacl/forwards.h"
//#include "viennacl/vector.hpp"
//#include "viennacl/tools/tools.hpp"
//#include "viennacl/linalg/sparse_matrix_operations.hpp"

#include "petscksp.h"
#include "viennacl/scalar.hpp"
#include "viennacl/vector.hpp"
#include <vector>

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
	void fill_petsc_vector(int size, LocalType* values, PETSC_VEC& out_vec) {
		// Set rhs values
		for (int i = 0; i < size; i++) { VecSetValue(out_vec, i, values[i], INSERT_VALUES); }
		VecAssemblyBegin(out_vec);
		VecAssemblyEnd(out_vec);
	}
	void create_petsc_vector(int size, PETSC_VEC& out_vec) {
		// Create vectors
		VecCreate(PETSC_COMM_WORLD, &out_vec);
		VecSetSizes(out_vec, PETSC_DECIDE, size);
		VecSetFromOptions(out_vec);
	}
	//void create_petsc_vectors(int iterations ,int size, std::vector<PETSC_VEC> out_vecs) {
	//	for (int i = 0; i < iterations; i++) {
	//		VecCreate(PETSC_COMM_WORLD, &(out_vecs[i]));
	//		VecSetSizes(out_vecs[i], PETSC_DECIDE, size);
	//		VecSetFromOptions(out_vecs[i]);
	//	}
	//}
	//void fill_petsc_vectors(int iterations, int size, LocalType** values, std::vector<PETSC_VEC> out_vecs) {
	//	for (int i = 0; i < iterations; i++) {
	//		for (int j = 0; j < size; j++) { VecSetValue(out_vecs[i], j, values[i][j], INSERT_VALUES); }
	//		VecAssemblyBegin(out_vecs[i]);
	//		VecAssemblyEnd(out_vecs[i]);
	//	}

	//}

	
}