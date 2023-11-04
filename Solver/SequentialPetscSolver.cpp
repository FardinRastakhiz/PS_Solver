#include "pch.h"
#include "SequentialPetscSolver.h"
#include "algorithms.h"
#include "IPreconditioner.h"
#include "matrix_factory.h"
#include "vector_factory.h"

#include "petscksp.h"
#include "chrono"
using namespace std::chrono;
namespace ses {

	template<typename mat_T, typename vec_T>
	SequentialPetscSolver<mat_T, vec_T>::SequentialPetscSolver(SolverArgs args) :
		PetscSolver<mat_T, vec_T>(args) {
		PetscErrorCode ierr;
		ierr = PetscInitialize(PETSC_NULL, PETSC_NULL, PETSC_NULL, PETSC_NULL);
		create_matrix(args.num_rows, args.num_cols, args.nnz, args.row_indices, args.col_indices, this->s_values, this->A);
		create_petsc_vector(args.num_rows, this->b);
		create_petsc_vector(args.num_rows, this->x);
		fill_petsc_vector(args.num_rows, this->s_b, this->b);
		this->algorithm = args.algorithm;
	
	}


	template<typename mat_T, typename vec_T>
	void SequentialPetscSolver<mat_T, vec_T>::Solve(int iteration_count, LocalType precision) {
		PetscErrorCode ierr;
		PetscMPIInt    size, rank;
		if (iteration_count != -1)
			PetscOptionsSetValue(NULL, "-ksp_max_it", std::to_string(iteration_count).c_str());
		if (precision != -1.0)
			PetscOptionsSetValue(NULL, "-ksp_atol", std::to_string(precision).c_str());
		ierr = PetscInitialize(PETSC_NULL, PETSC_NULL, PETSC_NULL, PETSC_NULL);
		PetscPrintf(PETSC_COMM_WORLD, "PETSC Initialized \n");

		KSP            ksp;           /* linear solver context */
		PetscReal      norm;          /* norm of solution error */
		PC             pc;



		ierr = KSPCreate(PETSC_COMM_WORLD, &ksp);
		ierr = KSPSetType(ksp, GetKSPType(this->algorithm));
		ierr = KSPSetOperators(ksp, this->A, this->A);
		ierr = KSPSetFromOptions(ksp);
		auto start = high_resolution_clock::now();
		ierr = KSPSolve(ksp, this->b, this->x);
		auto stop = high_resolution_clock::now();
		auto duration = duration_cast<microseconds>(stop - start);
		PetscPrintf(PETSC_COMM_WORLD, "Time taken by solver: %d microseconds \n", duration.count());
		PetscInt its;
		KSPGetIterationNumber(ksp, &its);
		PetscPrintf(PETSC_COMM_WORLD, "Solved With %d Iterations \n", its);
		this->ksp = ksp;
	}
	// call below function always before solve function
	template<class mat_T, class vec_T>
	void SequentialPetscSolver<mat_T, vec_T>::SetOptions(PetscBackend backend, int platform, int device, int num_threads) {
		if (backend == PetscBackend::OPENMP) {
			PetscOptionsSetValue(NULL, "-mat_type", "aijviennacl");
			PetscOptionsSetValue(NULL, "-vec_type", "viennacl");
			PetscOptionsSetValue(NULL, "-viennacl_backend", "openmp");
			PetscOptionsSetValue(NULL, "-omp_num_threads", std::to_string(num_threads).c_str());
		}
		if (backend == PetscBackend::OPENCL) {
			PetscOptionsSetValue(NULL, "-mat_type", "aijviennacl");
			PetscOptionsSetValue(NULL, "-vec_type", "viennacl");
			PetscOptionsSetValue(NULL, "-viennacl_backend", "opencl");
			PetscOptionsSetValue(NULL, "-viennacl_opencl_device", std::to_string(device).c_str());
			PetscOptionsSetValue(NULL, "-viennacl_opencl_platform", std::to_string(platform).c_str());
		}
	}
	template<typename mat_T, typename vec_T>
	void SequentialPetscSolver<mat_T, vec_T>::Solve(vec_T b, int iteration_count, LocalType precision) {
		PetscErrorCode ierr;
		if (iteration_count != -1)
			PetscOptionsSetValue(NULL, "-ksp_max_it", std::to_string(iteration_count).c_str());
		if (precision != -1.0)
			PetscOptionsSetValue(NULL, "-ksp_atol", std::to_string(precision).c_str());
		this->b = b;
		auto start = high_resolution_clock::now();
		ierr = KSPSolve(ksp, this->b, this->x);
		auto stop = high_resolution_clock::now();
		auto duration = duration_cast<microseconds>(stop - start);
		PetscPrintf(PETSC_COMM_WORLD, "Time taken by solver: %d microseconds \n", duration.count());
		PetscInt its;
		KSPGetIterationNumber(ksp, &its);
		PetscPrintf(PETSC_COMM_WORLD, "Solved With %d Iterations \n", its);
		
	}
	template<typename mat_T, typename vec_T>
	void SequentialPetscSolver<mat_T, vec_T>::SetNewB(LocalType* b) {
		create_petsc_vector(this->args.num_rows, this->b);
		fill_petsc_vector(this->args.num_rows, b, this->b);
	}
	template<class mat_T, class vec_T>
	KSPType SequentialPetscSolver<mat_T, vec_T>::GetKSPType(Algorithm alg) {
		KSPType type;
		switch (this->algorithm) {
		case GMRES:
			type = KSPGMRES;
		default:
			type = KSPCG;
		}
		return type;
	}
	/*			constructor			*/
	template SequentialPetscSolver<PETSC_MAT, PETSC_VEC>;
}