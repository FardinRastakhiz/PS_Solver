// cxdll.cpp : Defines the exported functions for the DLL.
//

#include "pch.h"
#include "sparse_equation_solver.h"



//#include "framework.h"
//#include <fstream>
//#include <vector>
//#include <stdlib.h>
//#include <stdio.h>
//#include <chrono>
//#include <future>
//#include <ctype.h>
//#include <iostream>
//#include <string>
//#include <time.h>
//#include <memory>

#include "Solver.h"
#include "GPUSolver.h"
#include "SimpleGPUSolver.h"
#include "SimpleCPUSolver.h"
#include "SequentialGPUSolver.h"
#include "SequentialCPUSolver.h"
#include "utilities.h"
#include "vector_factory.h"

//#include "SimpleCPUSolver.h"
//#include "SequentialGPUSolver.h"
//#include "SequentialCPUSolver.h"
//#include "matrix_factory.h"
//#include "vector_factory.h"
//#include "Algorithm.h"
//#include "converters.h"
//#include "IPreconditioner.h"
//#include "utilities.h"


using namespace std;
//using namespace std::chrono;
using namespace ses;

std::unique_ptr<ISolver> solver;


CXDLL_API void ses_solve_pressure_gpu(int num_rows, int num_cols, int nnz, int* row_indices, int* col_indices, double * values, double* b, double* x) {
	
	SolverArgs args(num_rows, num_cols, nnz, row_indices, col_indices, values, b, GMRES);

	// create solvers and solve the matrix
	solver = std::make_unique<SimpleGPUSolver<VI_SELL_MAT, VI_VEC>> (args);
	solver->Solve(1000, 1.0e-10);

	// Get the result
	x = ses::cast_to<double>(solver->GetResult(), num_cols);

	auto new_b = ses::cast_to<double>(solver->CalculateB(), num_rows);
	//std::cout << x1[0] << std::endl;
	//std::cout << x1[1] << std::endl;
	save_vector_pointer(x, num_cols, "output_x.txt");
	save_vector_pointer(new_b, num_rows, "output_b.txt");
}


CXDLL_API void ses_solve_pressure_cpu(int num_rows, int num_cols, int nnz, int* row_indices, int* col_indices, double* values, double* b, double* x) {

	SolverArgs args(num_rows, num_cols, nnz, row_indices, col_indices, values, b, GMRES);

	// create solvers and solve the matrix
	solver = std::make_unique<SimpleCPUSolver<PETSC_MAT,PETSC_VEC >>(args);
	solver->Solve(1000, -1.0);
	// Get the result
	x = ses::cast_to<double>(solver->GetResult(), num_cols);

	save_vector_pointer(x, num_cols, "output_x.txt");
	// save x and b
	if (SimpleCPUSolver<PETSC_MAT, PETSC_VEC >* c = dynamic_cast<SimpleCPUSolver<PETSC_MAT, PETSC_VEC >*>(solver.get()))
	{
		//c->PrintX();
		//c->PrintResultB();
		c->Finalize();
	}
	
}
extern "C" CXDLL_API void ses_solve_begin_density_cpu(int num_rows, int num_cols, int num_non_zero, int* row_indices, int* col_indices, double* values, double* b, double* x) {
	SolverArgs args(num_rows, num_cols, num_non_zero, row_indices, col_indices, values, b, GMRES);

	// create solvers and solve the matrix
	solver = std::make_unique<SequentialCPUSolver<PETSC_MAT, PETSC_VEC >>(args);
	solver->Solve(1000, -1.0);
	// Get the result
	x = ses::cast_to<double>(solver->GetResult(), num_cols);

	//save_vector_pointer(x, num_cols, "output_x0.txt");
}

extern "C" CXDLL_API void ses_solve_next(int size, double* b, double* x) {
	if (SequentialCPUSolver<PETSC_MAT, PETSC_VEC >* c = dynamic_cast<SequentialCPUSolver<PETSC_MAT, PETSC_VEC >*>(solver.get()))
	{
		c->SetNewB(b);
		c->Solve(c->b,1000 , -1.0);
		// Get the result
		x = ses::cast_to<double>(solver->GetResult(), size);

		//save_vector_pointer(x, size, "output_x.txt");
	}
	
}
//
CXDLL_API void ses_solve_begin_density_gpu(
	int num_rows, int num_cols, int num_non_zero,
	int* row_indices, int* col_indices, double* values, double* b, double* x)
{
	SolverArgs args(num_rows, num_cols, num_non_zero, row_indices, col_indices, values, b, GMRES);

	solver = std::make_unique<SequentialGPUSolver<VI_SELL_MAT, VI_VEC>>(args);
	solver->Solve(1000, 0.1);
	x = solver->GetResult();
}

CXDLL_API void ses_solve_next(double* rhs, double* x) {
	/*SequentialGPUSolver<VI_SELL_MAT, VI_VEC>* gpu_seq_solver =
		dynamic_cast<SequentialGPUSolver<VI_SELL_MAT, VI_VEC>*>(solver.get());

	SequentialCPUSolver<VI_SELL_MAT, VI_VEC>* cpu_seq_solver =
		dynamic_cast<SequentialCPUSolver<VI_SELL_MAT, VI_VEC>*>(solver.get());

	assert((gpu_seq_solver  || cpu_seq_solver), "It is not a Sequential Solver");*/

	////SolverArgs args(num_rows, num_cols, num_non_zero, row_indices, col_indices, values, b, GMRES);

	//if (gpu_seq_solver) {
	//	VI_VEC vec;
	//	create_vector(solver->args.num_rows, rhs, vec);
	//	gpu_seq_solver->Solve(vec, 1000, 0.1);
	//}
	//else if (cpu_seq_solver) {
	//	cpu_seq_solver->Solve(vec, 1000, 0.1);
	// }

	//x = solver->GetResult();
}
