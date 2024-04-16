#pragma once
#include "pch.h"
#include "SimplePetscSolver.h"



// IMPORTANT: Must be set prior to any ViennaCL includes if you want to use ViennaCL algorithms on Eigen objects
#define VIENNACL_WITH_UBLAS 1

// System headers
#include <iostream>
#include <vector>
#include <string>



#include "algorithms.h"
#include "IPreconditioner.h"
#include "matrix_factory.h"
#include "vector_factory.h"

#include "petscksp.h"
#include "chrono"
using namespace std::chrono;
namespace ses {
	template<class mat_T, class vec_T>
	SimplePetscSolver<mat_T, vec_T>::SimplePetscSolver(SolverArgs args) :
	PetscSolver<mat_T, vec_T>(args) {
		SetLocalTypes(args);
	}
	template<class mat_T, class vec_T>
	void SimplePetscSolver<mat_T, vec_T>::SetLocalTypes(SolverArgs args) {
		this->args = args;
		
	}
	template<class mat_T, class vec_T>
	void SimplePetscSolver<mat_T, vec_T>::Initialize() {
		PetscErrorCode ierr;
		ierr = PetscInitialize(PETSC_NULL, PETSC_NULL, PETSC_NULL, PETSC_NULL);
		PetscPrintf(PETSC_COMM_WORLD, "PETSC Initialized \n");
		//PetscPrintf(PETSC_COMM_WORLD, "error code \n", ierr);
		create_matrix(this->args.num_rows, this->args.num_cols, this->args.nnz, this->args.row_indices, this->args.col_indices, this->s_values, this->A);
		create_petsc_vector(this->args.num_rows, this->b);
		create_petsc_vector(this->args.num_rows, this->y);
		create_petsc_vector(this->args.num_rows, this->u);
		create_petsc_vector(this->args.num_rows, this->x);
		fill_petsc_vector(this->args.num_rows, this->s_b, this->b);
		fill_petsc_vector(this->args.num_rows, this->s_x, this->x);
		this->algorithm = this->args.algorithm;
	}
	template<class mat_T, class vec_T>
	void SimplePetscSolver<mat_T, vec_T>::Solve(int iteration_count, LocalType precision) {

		PetscErrorCode ierr;
		PetscMPIInt    size, rank;
		KSP            ksp;           /* linear solver context */
		PetscReal      norm;          /* norm of solution error */
		PC             pc;
		
		auto start = high_resolution_clock::now();
		ierr = KSPCreate(PETSC_COMM_WORLD, &ksp);
		//KSPGetPC(ksp, &pc);
		//PCSetType(pc, PCBJACOBI);
		ierr = KSPSetType(ksp, GetKSPType(this->algorithm));
		//for inital guess
		//KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
		ierr = KSPSetOperators(ksp, this->A, this->A);
		ierr = KSPSetFromOptions(ksp);
		ierr = KSPSolve(ksp, this->b, this->x);
		auto stop = high_resolution_clock::now();
		auto duration = duration_cast<microseconds>(stop - start);
		PetscPrintf(PETSC_COMM_WORLD, "Time taken by solver: %d microseconds \n", duration.count());
		PetscInt its;
		KSPGetIterationNumber(ksp, &its);
		PetscPrintf(PETSC_COMM_WORLD, "Solved With %d Iterations \n", its);
		ierr = MatMult(this->A, this->x, this->y); // y = A*x
		ierr = VecAXPY(this->y, -1.0, this->b); // y = y - b
		ierr = VecNorm(this->y, NORM_2, &norm); // Calculate the 2-norm of y

		PetscPrintf(PETSC_COMM_WORLD, "norm of error %g\n", norm);
		KSPDestroy(&ksp);

	}
	// call below function always before solve function
	template<class mat_T, class vec_T>
	void SimplePetscSolver<mat_T, vec_T>::SetOptions(PetscBackend backend, int platform, int device , int num_threads,int iteration_count, LocalType precision, int preconditioner) {
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
		//setting the preconditioner
		// acceptable values : 1 = jacobi, 2 = ILU, 3 = gasm, 4 = icc, 5 = ksp, 6 = bjacobi, 7 = sor, 8 = asm, 9 = cholesky
		switch (preconditioner) {
			case 1:
				PetscOptionsSetValue(NULL, "-pc_type", "jacobi");
				break;
			case 2:
				PetscOptionsSetValue(NULL, "-pc_type", "ilu");
				break;
			case 3:
				PetscOptionsSetValue(NULL, "-pc_type", "gasm");
				break;
			case 4:
				PetscOptionsSetValue(NULL, "-pc_type", "icc");
				break;
			case 5:
				PetscOptionsSetValue(NULL, "-pc_type", "ksp");
				break;
			case 6:
				PetscOptionsSetValue(NULL, "-pc_type", "bjacobi");
				break;
			case 7:
				PetscOptionsSetValue(NULL, "-pc_type", "sor");
				break;
			case 8:
				PetscOptionsSetValue(NULL, "-pc_type", "asm");
				break;
			case 9:
				PetscOptionsSetValue(NULL, "-pc_type", "cholesky");
				break;
		}
		
		if (iteration_count != -1)
			PetscOptionsSetValue(NULL, "-ksp_max_it", std::to_string(iteration_count).c_str());
		if (precision != -1.0)
			PetscOptionsSetValue(NULL, "-ksp_atol", std::to_string(precision).c_str());
	}
	template<class mat_T, class vec_T>
	LocalType* SimplePetscSolver<mat_T, vec_T>::GetResult() {
		PetscScalar* a;
		VecGetArray(this->x, &a);
		LocalType* res = (LocalType*)malloc(sizeof(LocalType) * this->args.num_rows);
		for (int i = 0; i < this->args.num_rows; i++) {
			res[i] = (LocalType)a[i];
		}
		return res;
	}

	template<class mat_T, class vec_T>
	void SimplePetscSolver<mat_T, vec_T>::PrintX() {
		// Optional: Print the solution
		PetscViewer viewer;
		PetscViewerASCIIOpen(PETSC_COMM_SELF, "petsc_x.txt", &viewer); // open an ASCII file for writing
		VecView(this->x, viewer); // save the vector to the file
		PetscViewerDestroy(&viewer); // close the file
	}

	template<class mat_T, class vec_T>
	void SimplePetscSolver<mat_T, vec_T>::PrintResultB() {

		/* Now, multiply A by x and check if it equals b */
		MatMult(this->A, this->x, this->y); // y = A*x
		// Optional: Print the rhs
		PetscViewer vrhs;
		PetscViewerASCIIOpen(PETSC_COMM_SELF, "petsc_rhs.txt", &vrhs); // open an ASCII file for writing
		VecView(this->y, vrhs); // save the vector to the file
		PetscViewerDestroy(&vrhs); // close the file
	}

	template<class mat_T, class vec_T>
	void SimplePetscSolver<mat_T, vec_T>::Finalize() {
		VecDestroy(&this->x);
		VecDestroy(&this->b);
		MatDestroy(&this->A);
		PetscFinalize();
	}
	template<class mat_T, class vec_T>
	KSPType SimplePetscSolver<mat_T, vec_T>::GetKSPType(Algorithm alg) {
		KSPType type;
		switch (this->algorithm) {
		case GMRES:
			type = KSPGMRES;
		default:
			type = KSPCG;
		}
		return type;
	}

	/*				implementation options				*/
	/*--------------------------------------------------*/


	/*			constructor			*/
	template SimplePetscSolver<PETSC_MAT, PETSC_VEC>;

}
