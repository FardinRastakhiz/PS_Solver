
// cxdll.h

#define CXDLL_API __declspec(dllexport)

// Solvers for Use-Cases
extern "C" CXDLL_API int ses_solve_pressure_cpu(int num_rows, int num_non_zero, int* row_indices, int* col_indices, double* values, double* b, double* x, int iteration_count, double precision, int use_open_mp, int num_threads);
extern "C" CXDLL_API int ses_solve_pressure_gpu(int num_rows, int num_non_zero, int* row_indices, int* col_indices, double* values, double* b, double* x,int target_lib, int iteration_count, double precision, int platform, int device);
extern "C" CXDLL_API int ses_solve_begin_density_cpu(int num_rows, int num_non_zero, int* row_indices, int* col_indices, double* values, double* b, double* x, int iteration_count, double precision, int use_open_mp, int num_threads);
extern "C" CXDLL_API int ses_solve_begin_density_gpu(int num_rows, int num_non_zero, int* row_indices, int* col_indices, double* values, double* b, double* x, int target_lib, int iteration_count, int precision, int platform, int device);
extern "C" CXDLL_API int ses_solve_next(double* rhs, double* x, int iteration_count, int precision);



// Solver for Technical Aspect
// Not required for now


// This class is exported from the dll

//extern "C" CXDLL_API double* ses_get_next_result(int numRows, int numNonzero, int* rowIndices, int* colIndices, double* values, double* rhs);

//extern "C" CXDLL_API void solve_matrix2(int numRows, int numCols, int numNonzero, int* rowIndices, int* colIndices, double* values, double* rhs, double* result);
//
//extern "C" CXDLL_API double* ses_symmetric_cpu_cg(int numRows, int numNonzero, int* rowIndices, int* colIndices, double* values, double* rhs, double* result);