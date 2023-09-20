#pragma once

#ifndef SES_MAT_FACTORY
#include "utilities.h"

namespace ses {

	void create_matrix(int num_rows, int num_cols, int nnz,
		int* rows, int* cols, LocalType* values, VI_SELL_MAT& matrix);

	void create_matrix(int num_rows, int num_cols, int nnz,
		int* rows, int* cols, LocalType* values, VI_COMP_Mat& matrix);


	void create_matrix(int num_rows, int num_cols, int nnz,
		int* rows, int* cols, LocalType* values, VI_COO_Mat& matrix);
}

#define SES_MAT_FACTORY
#endif