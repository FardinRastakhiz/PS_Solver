
// cxdll.h

#define CXDLL_API __declspec(dllexport)

//#include <cstdlib>
//#include <iostream>
//#include <memory>
//#include "Solver.h"



// Solvers for Use-Cases
extern "C" {
	class SolverContainer;

	CXDLL_API void* ses_solve_pressure_cpu(int num_rows, int num_non_zero, int* row_indices, int* col_indices, double* values, double* b, double* x, int iteration_count, double precision, int use_open_mp, int num_threads,int preconditioner);
	CXDLL_API void* ses_solve_pressure_gpu(int num_rows, int num_non_zero, int* row_indices, int* col_indices, double* values, double* b, double* x, int target_lib, int iteration_count, double precision, int platform, int device, int preconditioner);
	CXDLL_API void* ses_solve_begin_density_cpu(int num_rows, int num_non_zero, int* row_indices, int* col_indices, double* values, double* b, double* x, int iteration_count, double precision, int use_open_mp, int num_threads, int preconditioner);
	CXDLL_API void* ses_solve_begin_density_gpu(int num_rows, int num_non_zero, int* row_indices, int* col_indices, double* values, double* b, double* x, int target_lib, int iteration_count, int precision, int platform, int device, int preconditioner);
	CXDLL_API int ses_solve_next(SolverContainer * solver_container, double* rhs, double* x, int iteration_count, int precision);
	CXDLL_API int ses_write_devices_to_file();
	CXDLL_API int ses_build_initial_guess(int numRows, int numRowsAct, double* locX, double* locY, double* locZ, double* locActX, double* locActY, double* locActZ, int* bnd, double* x);


}

// Solver for Technical Aspect
// Not required for now


// This class is exported from the dll

//extern "C" CXDLL_API double* ses_get_next_result(int numRows, int numNonzero, int* rowIndices, int* colIndices, double* values, double* rhs);

//extern "C" CXDLL_API void solve_matrix2(int numRows, int numCols, int numNonzero, int* rowIndices, int* colIndices, double* values, double* rhs, double* result);
//
//extern "C" CXDLL_API double* ses_symmetric_cpu_cg(int numRows, int numNonzero, int* rowIndices, int* colIndices, double* values, double* rhs, double* result);