// cxdll.cpp : Defines the exported functions for the DLL.
//

#include "pch.h"
#include "framework.h"
#include "cxdll.h"
#include <iostream>



    // This is an example of an exported function.
    CXDLL_API void solve_matrix(int numRows, int numNonzero, int* rowIndices, int* colIndices, double* values, double* rhs, double* result) {
        std::cout << numNonzero;
        std::cout << "this is from new compiled";

    }
