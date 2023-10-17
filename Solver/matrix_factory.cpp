#include "pch.h"
#include "matrix_factory.h"

//
// ublas includes
//
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/operation_sparse.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <numeric>
#include "petscksp.h"

// IMPORTANT: Must be set prior to any ViennaCL includes if you want to use ViennaCL algorithms on Eigen objects
#define VIENNACL_WITH_UBLAS 1

using namespace boost::numeric;

namespace ses {

	typedef ublas::mapped_matrix<LocalType> CPUMatType;


	ublas::mapped_matrix<LocalType> create_ublas_matrix(int num_rows, int num_cols, int nnz, int* rows, int* cols, LocalType* values) {

		ublas::mapped_matrix<LocalType> matrix(num_rows, num_cols, nnz);
		//int min_row = 100000;

		for (int i = 0; i < nnz; i++) {
			if (i % 1000000 == 0)
				std::cout << i << std::endl;
			matrix(rows[i] - 1, cols[i] - 1) = values[i];
			/*if (rows[i] < min_row)
				min_row = rows[i];
			if (rows[i] != cols[i])
			{
				matrix(cols[i] - 1, rows[i] - 1) = values[i];
			}*/
			//else
			//{
			//	matrix(cols[i] - 1, rows[i] - 1) += values[i];
			//}
		}
		return matrix;
	}


	// The main function to create CPU Matrix
	CPUMatType create_cpu_matrix(int num_rows, int num_cols, int nnz, int* rows, int* cols, LocalType* values) {
		return create_ublas_matrix(num_rows, num_cols, nnz, rows, cols, values);
	}

	void create_matrix(int num_rows, int num_cols, int nnz, int* rows, int* cols, LocalType* values, VI_SELL_MAT& matrix) {
		CPUMatType ublas_mat = create_cpu_matrix(num_rows, num_cols, nnz, rows, cols, values);
		std::cout <<"ublas_mat1: " << ublas_mat(0, 0) << std::endl;
		//matrix = VI_SELL_MAT(num_rows, num_cols, 64);
		viennacl::copy(ublas_mat, matrix);
		std::cout << "matrix row per block: " << matrix.rows_per_block() << std::endl;
	}

	void create_matrix(int num_rows, int num_cols, int nnz, int* rows, int* cols, LocalType* values, VI_COMP_Mat& matrix) {
		
		viennacl::copy(create_cpu_matrix(num_rows, num_cols, nnz, rows, cols, values), matrix);
	}


	void create_matrix(int num_rows, int num_cols, int nnz, int* rows, int* cols, LocalType* values, VI_COO_Mat& matrix) {
		
		viennacl::copy(create_cpu_matrix(num_rows, num_cols, nnz, rows, cols, values), matrix);
	}

	void create_matrix(int num_rows, int num_cols, int nnz, int* rows, int* cols, LocalType* values, PETSC_MAT& matrix) {
		// Create a parallel matrix
		MatCreate(PETSC_COMM_WORLD, &matrix);
		MatSetSizes(matrix, PETSC_DECIDE, PETSC_DECIDE, num_rows, num_cols);
		MatSetFromOptions(matrix);
		MatSetUp(matrix);

		// Set matrix values
		for (int i = 0; i < nnz; i++) { MatSetValue(matrix, rows[i] - 1, cols[i] - 1, values[i], INSERT_VALUES); }
		MatAssemblyBegin(matrix, MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(matrix, MAT_FINAL_ASSEMBLY);
	}
}