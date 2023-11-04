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
#include "ViennaSolver.h"
#include "SimpleViennaSolver.h"
#include "SimplePetscSolver.h"
#include "SequentialViennaSolver.h"
#include "SequentialPetscSolver.h"
#include "utilities.h"
#include "vector_factory.h"

//#include "SimplePetscSolver.h"
//#include "SequentialViennaSolver.h"
//#include "SequentialPetscSolver.h"
//#include "matrix_factory.h"
//#include "vector_factory.h"
//#include "algorithms.h"
//#include "converters.h"
//#include "IPreconditioner.h"
//#include "utilities.h"


using namespace std;
//using namespace std::chrono;
using namespace ses;

std::unique_ptr<ISolver> solver;
enum TargetLibrary {
	PETSC_CPU = 1,
	PETSC_GPU = 2,
	VIENNA_CL_GPU = 3
} target_library;
int default_iteration_count = 1000;
int default_precision = 0.1;

CXDLL_API void ses_solve_pressure_gpu(int num_rows, int num_cols, int nnz, int* row_indices, int* col_indices, double * values, double* b, double* x, int target_lib, int iteration_count, double precision, int platform , int device) {
	target_library = (TargetLibrary)target_lib;
	if (target_library == TargetLibrary::VIENNA_CL_GPU) {
		SolverArgs args(num_rows, num_cols, nnz, row_indices, col_indices, values, b, GMRES);

		// create solvers and solve the matrix
		solver = std::make_unique<SimpleViennaSolver<VI_SELL_MAT, VI_VEC>>(args);
		solver->Solve(iteration_count == -1 ? default_iteration_count : iteration_count, precision == -1 ? default_precision : precision);

		// Get the result
		x = ses::cast_to<double>(solver->GetResult(), num_cols);

		auto new_b = ses::cast_to<double>(solver->CalculateB(), num_rows);
		//std::cout << x1[0] << std::endl;
		//std::cout << x1[1] << std::endl;
		save_vector_pointer(x, num_cols, "output_x.txt");
		save_vector_pointer(new_b, num_rows, "output_b.txt");
	}
	if (target_library == TargetLibrary::PETSC_GPU) {
		SolverArgs args(num_rows, num_cols, nnz, row_indices, col_indices, values, b, GMRES);

		// create solvers and solve the matrix
		solver = std::make_unique<SimplePetscSolver<PETSC_MAT, PETSC_VEC >>(args);
		

		save_vector_pointer(x, num_cols, "output_x.txt");
		// save x and b
		if (SimplePetscSolver<PETSC_MAT, PETSC_VEC >* c = dynamic_cast<SimplePetscSolver<PETSC_MAT, PETSC_VEC >*>(solver.get()))
		{
			// TODO : should get platform and device as parameter
			c->SetOptions(PetscBackend::OPENCL, platform, device);
			c->Solve(iteration_count, precision);
			// Get the result
			x = ses::cast_to<double>(c->GetResult(), num_cols);
			c->Finalize();
		}
	}
	

}


CXDLL_API void ses_solve_pressure_cpu(int num_rows, int num_cols, int nnz, int* row_indices, int* col_indices, double* values, double* b, double* x,int iteration_count, double precision,int use_open_mp ,int num_threads) {

	SolverArgs args(num_rows, num_cols, nnz, row_indices, col_indices, values, b, GMRES);

	// create solvers and solve the matrix
	solver = std::make_unique<SimplePetscSolver<PETSC_MAT,PETSC_VEC >>(args);
	

	save_vector_pointer(x, num_cols, "output_x.txt");
	// save x and b
	if (SimplePetscSolver<PETSC_MAT, PETSC_VEC >* c = dynamic_cast<SimplePetscSolver<PETSC_MAT, PETSC_VEC >*>(solver.get()))
	{
		c->SetOptions(use_open_mp == 0 ? PetscBackend::NORMAL : PetscBackend::OPENMP, 0 , 0 , num_threads );
		c->Solve(iteration_count, precision);
		// Get the result
		x = ses::cast_to<double>(c->GetResult(), num_cols);
		c->Finalize();
	}
	
}
CXDLL_API void ses_solve_begin_density_cpu(int num_rows, int num_cols, int num_non_zero, int* row_indices, int* col_indices, double* values, double* b, double* x , int iteration_count, double precision, int use_open_mp, int num_threads) {
	SolverArgs args(num_rows, num_cols, num_non_zero, row_indices, col_indices, values, b, GMRES);

	// create solvers and solve the matrix
	solver = std::make_unique<SequentialPetscSolver<PETSC_MAT, PETSC_VEC >>(args);

	if (SequentialPetscSolver<PETSC_MAT, PETSC_VEC >* c = dynamic_cast<SequentialPetscSolver<PETSC_MAT, PETSC_VEC >*>(solver.get()))
	{
		
		c->SetOptions(use_open_mp == 0 ? PetscBackend::NORMAL : PetscBackend::OPENMP, 0, 0, num_threads);
		c->Solve(iteration_count, precision);
		// Get the result
		x = ses::cast_to<double>(c->GetResult(), num_cols);
		target_library = PETSC_CPU;
	}
}


//
CXDLL_API void ses_solve_begin_density_gpu(
	int num_rows, int num_cols, int num_non_zero,
	int* row_indices, int* col_indices, double* values, double* b, double* x,int target_lib, int iteration_count , int precision , int platform , int device){
	target_library = (TargetLibrary)target_lib;
	if (target_library == VIENNA_CL_GPU) {
		SolverArgs args(num_rows, num_cols, num_non_zero, row_indices, col_indices, values, b, GMRES);

		solver = std::make_unique<SequentialViennaSolver<VI_SELL_MAT, VI_VEC>>(args);
		solver->Solve(iteration_count == -1 ? default_iteration_count : iteration_count, precision == -1 ? default_precision : precision);
		x = solver->GetResult();
	}
	if (target_library == PETSC_GPU) {
		SolverArgs args(num_rows, num_cols, num_non_zero, row_indices, col_indices, values, b, GMRES);

		// create solvers and solve the matrix
		solver = std::make_unique<SimplePetscSolver<PETSC_MAT, PETSC_VEC >>(args);

		if (SequentialPetscSolver<PETSC_MAT, PETSC_VEC >* c = dynamic_cast<SequentialPetscSolver<PETSC_MAT, PETSC_VEC >*>(solver.get()))
		{
			c->SetOptions(PetscBackend::OPENCL , platform , device);
			c->Solve(iteration_count, precision);
			// Get the result
			x = ses::cast_to<double>(c->GetResult(), num_cols);
		}
	}
	
}

CXDLL_API void ses_solve_next(double* rhs, double* x, int iteration_count, int precision) {
	switch (target_library)
	{
	case PETSC_GPU:
	case PETSC_CPU:
		if (SequentialPetscSolver<PETSC_MAT, PETSC_VEC >* c = dynamic_cast<SequentialPetscSolver<PETSC_MAT, PETSC_VEC >*>(solver.get()))
		{
			c->SetNewB(rhs);
			c->Solve(c->b, iteration_count, precision);

		}
	case VIENNA_CL_GPU:
		SequentialViennaSolver<VI_SELL_MAT, VI_VEC>* gpu_seq_solver =
			dynamic_cast<SequentialViennaSolver<VI_SELL_MAT, VI_VEC>*>(solver.get());
		assert((gpu_seq_solver), "It is not a Sequential Solver");
		gpu_seq_solver->Solve(rhs, iteration_count == -1 ? default_iteration_count : iteration_count, precision == -1 ? default_precision : precision);
		break;
	}
	x = solver->GetResult();
}
