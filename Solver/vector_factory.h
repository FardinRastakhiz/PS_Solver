#pragma once

#ifndef SES_VEC_FACTORY
#include "utilities.h"

namespace ses {
	void create_vector(int size, LocalType* values, VI_VEC& out_vec);
	void fill_petsc_vector(int size, LocalType* values, PETSC_VEC& out_vec);
	void create_petsc_vector(int size, PETSC_VEC& out_vec);
	void create_petsc_vectors(int iterations, int size, std::vector<PETSC_VEC> out_vecs);
	void fill_petsc_vectors(int iterations, int size, LocalType** values, std::vector<PETSC_VEC> out_vecs);
}

#define SES_VEC_FACTORY
#endif