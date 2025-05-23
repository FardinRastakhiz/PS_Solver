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
#include "InitialGuessBuilder.cpp"

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

float sample_f = -1.0f;
enum TargetLibrary {
	PETSC_CPU = 1,
	PETSC_GPU = 2,
	VIENNA_CL_GPU = 3
};

class SolverContainer {
public:
	TargetLibrary target_library;

	int default_iteration_count = 1000;
	int default_precision = 0.1;
	std::unique_ptr<ISolver> solver;
};

const std::string devicesFilePath = "devices.json";

CXDLL_API int ses_write_devices_to_file() {
	std::ofstream file(devicesFilePath);
	if (file.is_open())
	{
		file << "{" << endl;
		file << "	\"platforms\" : {" << endl;
		int max_len = viennacl::ocl::get_platforms().size();
		for (int i = 0; i <= max_len; i++) {
			std::vector<viennacl::ocl::device> devices;
			try {
				if (i == max_len) {
					file << endl;
					break;
				 }
				 devices = viennacl::ocl::get_platforms()[i].devices();
				 if (i != 0 && i < max_len) {
					 file << ',';
				 }
				}
				catch (exception) {
					continue;
				}
				file << "		\"" << viennacl::ocl::get_platforms()[i].info() << "\" : [" << endl;
				for (int j = 0; j < devices.size(); j++) {
					file << "			\"" << devices[j].name() << "\"";
					if (j < devices.size() - 1) {
						file << ',';
					}
						file << endl;
				}
				file << "		]";
		}
		file << "	}" << endl;
		file << "}";
		file.close();
	}
	return 0;
}

CXDLL_API int ses_build_initial_guess(int numRows, int numRowsAct, double* locX, double* locY, double* locZ, double* locActX, double* locActY, double* locActZ, int* bnd, double* x, double** x_out) {
	InitialGuessBuilder* builder = new InitialGuessBuilder();
	// set any options here
	//builder->COORDINATES_MIN_NODES = 1000;
	std::cout << "Using Intital Guess ..." << std::endl;
	builder->build_with_coordinates(numRows, numRowsAct, locX, locY, locZ, locActX, locActY, locActZ, bnd, x);
	*x_out = x;
	return 0;
}

CXDLL_API void* ses_solve_pressure_gpu(int num_rows, int nnz, int* row_indices, int* col_indices, double * values, double* b, double* x, double** x_out, int target_lib, int iteration_count, double precision, int platform , int device , int preconditioner) {
	SolverContainer* solver_container = new SolverContainer();
	solver_container->target_library = (TargetLibrary)target_lib;
	if (solver_container->target_library == TargetLibrary::VIENNA_CL_GPU) {
		SolverArgs args(num_rows, num_rows, nnz, row_indices, col_indices, values, b , x, GMRES);

		// create solvers and solve the matrix
		solver_container->solver = std::make_unique<SimpleViennaSolver<VI_SELL_MAT, VI_VEC>>(args);
		solver_container->solver->Solve(iteration_count == -1 ? solver_container->default_iteration_count : iteration_count, precision == -1 ? solver_container->default_precision : precision);

		// Get the result
		x = ses::cast_to<double>(solver_container->solver->GetResult(), num_rows);
		auto new_b = ses::cast_to<double>(solver_container->solver->CalculateB(), num_rows);
		//std::cout << x1[0] << std::endl;
		//std::cout << x1[1] << std::endl;
		save_vector_pointer(x, num_rows, "output_x.txt");
		
	}
	if (solver_container->target_library == TargetLibrary::PETSC_GPU) {
		SolverArgs args(num_rows, num_rows, nnz, row_indices, col_indices, values, b , x, GMRES);

		// create solvers and solve the matrix
		solver_container->solver = std::make_unique<SimplePetscSolver<PETSC_MAT, PETSC_VEC >>(args);
		

		save_vector_pointer(x, num_rows, "output_x.txt");
		// save x and b
		if (SimplePetscSolver<PETSC_MAT, PETSC_VEC >* c = dynamic_cast<SimplePetscSolver<PETSC_MAT, PETSC_VEC >*>(solver_container->solver.get()))
		{
			
			c->SetOptions(PetscBackend::OPENCL, platform, device,8, iteration_count, precision,preconditioner);
			c->Initialize();
			c->Solve(iteration_count, precision);
			// Get the result
			x = ses::cast_to<double>(c->GetResult(), num_rows);
			c->Finalize();
		}
	}
	*x_out = x;
	
	return (void*)solver_container;
}


CXDLL_API void* ses_solve_pressure_cpu(int num_rows, int nnz, int* row_indices, int* col_indices, double* values, double* b, double* x, double** x_out,int iteration_count, double precision,int use_open_mp ,int num_threads, int preconditioner) {
	//cout << "here is the initialized x last 100 rows" << endl;

	SolverArgs args(num_rows, num_rows, nnz, row_indices, col_indices, values, b , x, CG);

	SolverContainer* solver_container = new SolverContainer();
	// create solvers and solve the matrix
	solver_container->solver = std::make_unique<SimplePetscSolver<PETSC_MAT,PETSC_VEC >>(args);

	if (SimplePetscSolver<PETSC_MAT, PETSC_VEC >* c = dynamic_cast<SimplePetscSolver<PETSC_MAT, PETSC_VEC >*>(solver_container->solver.get()))
	{
		c->SetOptions(use_open_mp == 0 ? PetscBackend::NORMAL : PetscBackend::OPENMP, 0 , 0 , num_threads, iteration_count, precision, preconditioner);
		c->Initialize();
		c->Solve(iteration_count, precision);
		// Get the result
		x = ses::cast_to<double>(c->GetResult(), num_rows);
		save_vector_pointer(c->GetResult(), num_rows, "output_x.txt");
		c->PrintX();
		c->PrintResultB();
		c->Finalize();
	}
	*x_out = x;
	//for (int i = num_rows - 2000; i < num_rows; i++) {
	//	cout << x[i] << endl;
	//}
	return (void*)solver_container;
}
CXDLL_API void* ses_solve_begin_density_cpu(int num_rows, int num_non_zero, int* row_indices, int* col_indices, double* values, double* b, double* x, double** x_out , int iteration_count, double precision, int use_open_mp, int num_threads, int preconditioner) {
	SolverArgs args(num_rows, num_rows, num_non_zero, row_indices, col_indices, values, b, x, GMRES);

	SolverContainer* solver_container = new SolverContainer();
	// create solvers and solve the matrix
	solver_container->solver = std::make_unique<SequentialPetscSolver<PETSC_MAT, PETSC_VEC >>(args);

	if (SequentialPetscSolver<PETSC_MAT, PETSC_VEC >* c = dynamic_cast<SequentialPetscSolver<PETSC_MAT, PETSC_VEC >*>(solver_container->solver.get()))
	{
		
		c->SetOptions(use_open_mp == 0 ? PetscBackend::NORMAL : PetscBackend::OPENMP, 0, 0, num_threads, iteration_count, precision , preconditioner);
		c->Initialize();
		c->Solve(iteration_count, precision);
		// Get the result
		x = ses::cast_to<double>(c->GetResult(), num_rows);
		solver_container->target_library = PETSC_CPU;
	}
	*x_out = x;

	return (void*)solver_container;
}


//
CXDLL_API void* ses_solve_begin_density_gpu(int num_rows, int num_non_zero,	int* row_indices, int* col_indices, double* values, double* b, double* x, double** x_out,int target_lib, int iteration_count , int precision , int platform , int device, int preconditioner){

	SolverContainer* solver_container = new SolverContainer();
	solver_container->target_library = (TargetLibrary)target_lib;
	sample_f = 15.0f;
	if (solver_container->target_library == VIENNA_CL_GPU) {
		SolverArgs args(num_rows, num_rows, num_non_zero, row_indices, col_indices, values, b , x, GMRES);

		solver_container->solver = std::make_unique<SequentialViennaSolver<VI_SELL_MAT, VI_VEC>>(args);
		solver_container->solver->Solve(iteration_count == -1 ? solver_container->default_iteration_count : iteration_count, precision == -1 ? solver_container->default_precision : precision);
		x = solver_container->solver->GetResult();
	}
	if (solver_container->target_library == PETSC_GPU) {
		SolverArgs args(num_rows, num_rows, num_non_zero, row_indices, col_indices, values, b , x, GMRES);

		// create solvers and solve the matrix
		solver_container->solver = std::make_unique<SequentialPetscSolver<PETSC_MAT, PETSC_VEC >>(args);

		if (SequentialPetscSolver<PETSC_MAT, PETSC_VEC >* c = dynamic_cast<SequentialPetscSolver<PETSC_MAT, PETSC_VEC >*>(solver_container->solver.get()))
		{
			c->SetOptions(PetscBackend::OPENCL , platform , device, 8, iteration_count, precision , preconditioner);
			c->Initialize();
			c->Solve(iteration_count, precision);
			// Get the result
			x = ses::cast_to<double>(c->GetResult(), num_rows);
		}
	}
	*x_out = x;

	return (void*)solver_container;
	
}

CXDLL_API int ses_solve_next(SolverContainer* solver_container, double* rhs, double* x, double** x_out, int iteration_count, int precision) {
	switch ((TargetLibrary)(solver_container->target_library))
	{
	//case PETSC_GPU:

	case PETSC_CPU:
		if (SequentialPetscSolver<PETSC_MAT, PETSC_VEC >* c = dynamic_cast<SequentialPetscSolver<PETSC_MAT, PETSC_VEC >*>(solver_container->solver.get()))
		{
			c->SetNewB(rhs);
			c->Solve(iteration_count, precision);
		}
		break;
	case PETSC_GPU:
		if (SequentialPetscSolver<PETSC_MAT, PETSC_VEC >* c = dynamic_cast<SequentialPetscSolver<PETSC_MAT, PETSC_VEC >*>(solver_container->solver.get()))
		{
			c->SetNewB(rhs);
			c->Solve(iteration_count, precision);
		}
		break;
	case VIENNA_CL_GPU:
		SequentialViennaSolver<VI_SELL_MAT, VI_VEC>* gpu_seq_solver =
			dynamic_cast<SequentialViennaSolver<VI_SELL_MAT, VI_VEC>*>(solver_container->solver.get());
		assert((gpu_seq_solver), "It is not a Sequential Solver");
		gpu_seq_solver->Solve(rhs, iteration_count == -1 ? solver_container->default_iteration_count : iteration_count, precision == -1 ? solver_container->default_precision : precision);
		break;
	}
	x = solver_container->solver->GetResult();
	*x_out = x;
	return 0;
}
