#pragma once

#ifndef SES_VEC_FACTORY
#include "utilities.h"

namespace ses {
	void create_vector(int size, LocalType* values, VI_VEC& out_vec);
	void fill_petsc_vector(int size, LocalType* values, PETSC_VEC& out_vec);
	void create_petsc_vector(int size, PETSC_VEC& out_vec);
}

#define SES_VEC_FACTORY
#endif