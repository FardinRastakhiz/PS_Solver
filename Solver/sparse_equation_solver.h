
// cxdll.h

#define CXDLL_API __declspec(dllexport)

// This class is exported from the dll

extern "C" CXDLL_API void solve_matrix(int numRows, int numNonzero, int* rowIndices, int* colIndices, double* values, double* rhs, double* result);

extern "C" CXDLL_API double* ses_symmetric_cpu_cg(int numRows, int numNonzero, int* rowIndices, int* colIndices, double* values, double* rhs, double* result);