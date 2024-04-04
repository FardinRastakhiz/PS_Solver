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

const std::string devicesFilePath = "devices.txt";

CXDLL_API int ses_write_devices_to_file() {
	std::ofstream file(devicesFilePath);
	if (file.is_open())
	{
		int max_len = viennacl::ocl::get_platforms().size();
		for (int i = 0; i < max_len; i++) {
			try {
				std::vector<viennacl::ocl::device> devices = viennacl::ocl::get_platforms()[i].devices();
				file << viennacl::ocl::get_platforms()[i].info() << "-----" << devices.size() << endl;
				for (int j = 0; j < devices.size(); j++) {
					file << devices[j].name() << endl;
				}
			}
			catch(exception){
				//do nothing. must continue
			}
		}
		file.close();
	}
	return 0;
}

CXDLL_API int ses_build_initial_guess(int numRows, int numRowsAct, double* locX, double* locY, double* locZ, double* locActX, double* locActY, double* locActZ, int* bnd, double* x) {
	std::cout << "hello world";
	InitialGuessBuilder* builder = new InitialGuessBuilder();
	// set any options here
	//builder->COORDINATES_MIN_NODES = 1000;
	builder->build_with_coordinates(numRows, numRowsAct, locX, locY, locZ, locActX, locActY, locActZ, bnd, x);
	return 0;
}

CXDLL_API void* ses_solve_pressure_gpu(int num_rows, int nnz, int* row_indices, int* col_indices, double * values, double* b, double* x, int target_lib, int iteration_count, double precision, int platform , int device) {
	std::cout << '1' << std::endl;
	std::cout << '2' << std::endl;
	SolverContainer* solver_container = new SolverContainer();
	solver_container->target_library = (TargetLibrary)target_lib;
	std::cout << '3' << std::endl;
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
	std::cout << '4' << std::endl;
	if (solver_container->target_library == TargetLibrary::PETSC_GPU) {
		std::cout << '5' << std::endl;
		SolverArgs args(num_rows, num_rows, nnz, row_indices, col_indices, values, b , x, GMRES);

		std::cout << '6' << std::endl;
		// create solvers and solve the matrix
		solver_container->solver = std::make_unique<SimplePetscSolver<PETSC_MAT, PETSC_VEC >>(args);
		

		std::cout << '7' << std::endl;
		save_vector_pointer(x, num_rows, "output_x.txt");
		std::cout << '8' << std::endl;
		// save x and b
		if (SimplePetscSolver<PETSC_MAT, PETSC_VEC >* c = dynamic_cast<SimplePetscSolver<PETSC_MAT, PETSC_VEC >*>(solver_container->solver.get()))
		{
			std::cout << '9' << std::endl;
			
			c->SetOptions(PetscBackend::OPENCL, platform, device, iteration_count, precision);
			std::cout << '10' << std::endl;
			c->Initialize();
			std::cout << '11' << std::endl;
			c->Solve(iteration_count, precision);
			std::cout << '12' << std::endl;
			// Get the result
			x = ses::cast_to<double>(c->GetResult(), num_rows);
			std::cout << '13' << std::endl;
			c->Finalize();
			std::cout << '14' << std::endl;
		}
	}
	std::cout << '15' << std::endl;
	
	return (void*)solver_container;
}


CXDLL_API void* ses_solve_pressure_cpu(int num_rows, int nnz, int* row_indices, int* col_indices, double* values, double* b, double* x,int iteration_count, double precision,int use_open_mp ,int num_threads) {
	//cout << "here is the initialized x last 100 rows" << endl;

	SolverArgs args(num_rows, num_rows, nnz, row_indices, col_indices, values, b , x, CG);

	SolverContainer* solver_container = new SolverContainer();
	// create solvers and solve the matrix
	solver_container->solver = std::make_unique<SimplePetscSolver<PETSC_MAT,PETSC_VEC >>(args);

	std::cout << "3333333333333333: " << typeid(solver_container->solver.get()).name() << std::endl;
	if (SimplePetscSolver<PETSC_MAT, PETSC_VEC >* c = dynamic_cast<SimplePetscSolver<PETSC_MAT, PETSC_VEC >*>(solver_container->solver.get()))
	{
		std::cout << "3333333333333333: " << typeid(solver_container->solver.get()).name() << std::endl;
		c->SetOptions(use_open_mp == 0 ? PetscBackend::NORMAL : PetscBackend::OPENMP, 0 , 0 , num_threads, iteration_count, precision);
		c->Initialize();
		c->Solve(iteration_count, precision);
		// Get the result
		//x = ses::cast_to<double>(c->GetResult(), num_rows);
		save_vector_pointer(c->GetResult(), num_rows, "output_x.txt");
		c->PrintX();
		c->PrintResultB();
		c->Finalize();
	}
	
	//for (int i = num_rows - 2000; i < num_rows; i++) {
	//	cout << x[i] << endl;
	//}
	return (void*)solver_container;
}
CXDLL_API void* ses_solve_begin_density_cpu(int num_rows, int num_non_zero, int* row_indices, int* col_indices, double* values, double* b, double* x , int iteration_count, double precision, int use_open_mp, int num_threads) {
	SolverArgs args(num_rows, num_rows, num_non_zero, row_indices, col_indices, values, b, x, GMRES);

	SolverContainer* solver_container = new SolverContainer();
	// create solvers and solve the matrix
	solver_container->solver = std::make_unique<SequentialPetscSolver<PETSC_MAT, PETSC_VEC >>(args);

	if (SequentialPetscSolver<PETSC_MAT, PETSC_VEC >* c = dynamic_cast<SequentialPetscSolver<PETSC_MAT, PETSC_VEC >*>(solver_container->solver.get()))
	{
		
		c->SetOptions(use_open_mp == 0 ? PetscBackend::NORMAL : PetscBackend::OPENMP, 0, 0, num_threads);
		c->Solve(iteration_count, precision);
		// Get the result
		x = ses::cast_to<double>(c->GetResult(), num_rows);
		solver_container->target_library = PETSC_CPU;
	}

	return (void*)solver_container;
}


//
CXDLL_API void* ses_solve_begin_density_gpu(int num_rows, int num_non_zero,	int* row_indices, int* col_indices, double* values, double* b, double* x,int target_lib, int iteration_count , int precision , int platform , int device){

	SolverContainer* solver_container = new SolverContainer();
	std::cout << "1" << std::endl;
	solver_container->target_library = (TargetLibrary)target_lib;
	std::cout << "2" << std::endl;
	sample_f = 15.0f;
	std::cout << "3: " << sample_f << std::endl;
	if (solver_container->target_library == VIENNA_CL_GPU) {
		SolverArgs args(num_rows, num_rows, num_non_zero, row_indices, col_indices, values, b , x, GMRES);

		solver_container->solver = std::make_unique<SequentialViennaSolver<VI_SELL_MAT, VI_VEC>>(args);
		solver_container->solver->Solve(iteration_count == -1 ? solver_container->default_iteration_count : iteration_count, precision == -1 ? solver_container->default_precision : precision);
		x = solver_container->solver->GetResult();
	}
	if (solver_container->target_library == PETSC_GPU) {
		SolverArgs args(num_rows, num_rows, num_non_zero, row_indices, col_indices, values, b , x, GMRES);

		// create solvers and solve the matrix
		solver_container->solver = std::make_unique<SimplePetscSolver<PETSC_MAT, PETSC_VEC >>(args);

		std::cout << "44444444444: " << typeid(solver_container->solver.get()).name() << std::endl;
		if (SequentialPetscSolver<PETSC_MAT, PETSC_VEC >* c = dynamic_cast<SequentialPetscSolver<PETSC_MAT, PETSC_VEC >*>(solver_container->solver.get()))
		{
			std::cout << "44444444444: " << typeid(solver_container->solver.get()).name() << std::endl;
			c->SetOptions(PetscBackend::OPENCL , platform , device);
			c->Solve(iteration_count, precision);
			// Get the result
			x = ses::cast_to<double>(c->GetResult(), num_rows);
		}
	}
	std::cout << "4: " << sample_f << std::endl;
	return (void*)solver_container;
	
}

CXDLL_API int ses_solve_next(SolverContainer* solver_container, double* rhs, double* x, int iteration_count, int precision) {
	std::cout << "a" << std::endl;
	std::cout << "b: " << typeid(solver_container->solver.get()).name() << std::endl;
	switch ((TargetLibrary)(solver_container->target_library))
	{
	//case PETSC_GPU:

	case PETSC_CPU:
		std::cout << "aaaa: " << std::endl;
		if (SequentialPetscSolver<PETSC_MAT, PETSC_VEC >* c = dynamic_cast<SequentialPetscSolver<PETSC_MAT, PETSC_VEC >*>(solver_container->solver.get()))
		{
			std::cout << "bbbb: " << std::endl;
			c->SetNewB(rhs);
			c->Solve(c->b, iteration_count, precision);

		}
		break;
	case PETSC_GPU:
		std::cout << "eeee: " << std::endl;
		if (SequentialPetscSolver<PETSC_MAT, PETSC_VEC >* c = dynamic_cast<SequentialPetscSolver<PETSC_MAT, PETSC_VEC >*>(solver_container->solver.get()))
		{
			std::cout << "ffff: " << std::endl;
			c->SetNewB(rhs);
			c->Solve(c->b, iteration_count, precision);

		}
		break;
	case VIENNA_CL_GPU:
		std::cout << "ccccc: " << std::endl;
		SequentialViennaSolver<VI_SELL_MAT, VI_VEC>* gpu_seq_solver =
			dynamic_cast<SequentialViennaSolver<VI_SELL_MAT, VI_VEC>*>(solver_container->solver.get());
		assert((gpu_seq_solver), "It is not a Sequential Solver");
		gpu_seq_solver->Solve(rhs, iteration_count == -1 ? solver_container->default_iteration_count : iteration_count, precision == -1 ? solver_container->default_precision : precision);
		break;
	}
	std::cout << "ddddd: " << std::endl;
	x = solver_container->solver->GetResult();
	return 0;
}
