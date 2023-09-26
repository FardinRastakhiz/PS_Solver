#include "pch.h"
#include "SequentialCPUSolver.h"
#include "petscksp.h"
namespace ses {

	template<typename mat_T, typename vec_T>
	SequentialCPUSolver<mat_T, vec_T>::SequentialCPUSolver(SolverArgs args) :
		CPUSolver<mat_T, vec_T>(args) {
		PetscErrorCode ierr;
		ierr = PetscInitialize(PETSC_NULL, PETSC_NULL, PETSC_NULL, PETSC_NULL);
		create_matrix(args.num_rows, args.num_cols, args.nnz, args.row_indices, args.col_indices, this->s_values, this->A);
		create_petsc_vectors(args.num_rows, this->bs);
		create_petsc_vectors(args.num_rows, this->xs);
		fill_petsc_vector(args.num_rows, this->s_bs, this->bs);
		this->algorithm = args.algorithm;
	
	}


	template<typename mat_T, typename vec_T>
	void SequentialCPUSolver<mat_T, vec_T>::Solve(int iteration_count, LocalType precision) {
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

		for (int i = 0; i < iteration_count; i++) {
			PetscPrintf(PETSC_COMM_WORLD, "Iteration : %d \n", i);
			auto start = high_resolution_clock::now();
			ierr = KSPSolve(ksp, this->bs[i], this->xs[i]);
			auto stop = high_resolution_clock::now();
			auto duration = duration_cast<microseconds>(stop - start);
			PetscPrintf(PETSC_COMM_WORLD, "Time taken by solver: %d microseconds \n", duration.count());
			PetscInt its;
			KSPGetIterationNumber(ksp, &its);
			PetscPrintf(PETSC_COMM_WORLD, "Solved With %d Iterations \n", its);
		}

		KSPDestroy(&ksp);
	}

	template<typename mat_T, typename vec_T>
	void SequentialCPUSolver<mat_T, vec_T>::Solve(vec_T b, int iteration_count, LocalType precision) {
		throw "Not Implemented";
	}


	/*			constructor			*/
	template SequentialCPUSolver<PETSC_MAT, PETSC_VEC>;
}