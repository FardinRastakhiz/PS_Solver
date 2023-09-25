#pragma once
#include "pch.h"
#include "SimpleCPUSolver.h"



// IMPORTANT: Must be set prior to any ViennaCL includes if you want to use ViennaCL algorithms on Eigen objects
#define VIENNACL_WITH_UBLAS 1

// System headers
#include <iostream>
#include <vector>
#include <string>



#include "Algorithm.h"
#include "IPreconditioner.h"
#include "matrix_factory.h"
#include "vector_factory.h"

#include "petscksp.h"
#include "chrono"
using namespace std::chrono;
namespace ses {
	template<class mat_T, class vec_T>
	SimpleCPUSolver<mat_T, vec_T>::SimpleCPUSolver(SolverArgs args) :
		CPUSolver<mat_T, vec_T>(args) {
		SetLocalTypes(args);
	}
	template<class mat_T, class vec_T>
	void SimpleCPUSolver<mat_T, vec_T>::SetLocalTypes(SolverArgs args) {
		PetscErrorCode ierr;
		ierr = PetscInitialize(PETSC_NULL, PETSC_NULL, PETSC_NULL, PETSC_NULL);
		create_matrix(args.num_rows, args.num_cols, args.nnz, args.row_indices, args.col_indices, this->s_values, this->A);
		create_petsc_vector(args.num_rows, this->b);
		create_petsc_vector(args.num_rows, this->y);
		create_petsc_vector(args.num_rows, this->u);
		create_petsc_vector(args.num_rows, this->x);
		fill_petsc_vector(args.num_rows, this->s_b, this->b);
		this->algorithm = args.algorithm;
	}
	template<class mat_T, class vec_T>
	void SimpleCPUSolver<mat_T, vec_T>::Solve(int iteration_count, LocalType precision) {
		PetscErrorCode ierr;
		PetscMPIInt    size, rank;
		if(iteration_count != -1)
		PetscOptionsSetValue(NULL, "-ksp_max_it", std::to_string(iteration_count).c_str());
		if(precision != -1.0)
		PetscOptionsSetValue(NULL, "-ksp_atol", std::to_string(precision).c_str());
		ierr = PetscInitialize(PETSC_NULL, PETSC_NULL, PETSC_NULL, PETSC_NULL);
		PetscPrintf(PETSC_COMM_WORLD, "PETSC Initialized \n");

		KSP            ksp;           /* linear solver context */
		PetscReal      norm;          /* norm of solution error */
		PC             pc;

		
		auto start = high_resolution_clock::now();
		ierr = KSPCreate(PETSC_COMM_WORLD, &ksp);
		ierr = KSPSetType(ksp, GetKSPType(this->algorithm));
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
		double* newX = (double*)this->x;
		std::cout << newX[1000000];

	}

	template<class mat_T, class vec_T>
	LocalType* SimpleCPUSolver<mat_T, vec_T>::GetResult() {
		PetscScalar* a;
		VecGetArray(this->x, &a);
		return (LocalType*)a;
	}

	template<class mat_T, class vec_T>
	void SimpleCPUSolver<mat_T, vec_T>::PrintX() {
		// Optional: Print the solution
		PetscViewer viewer;
		PetscViewerASCIIOpen(PETSC_COMM_SELF, "petsc_x.txt", &viewer); // open an ASCII file for writing
		VecView(this->x, viewer); // save the vector to the file
		PetscViewerDestroy(&viewer); // close the file
	}

	template<class mat_T, class vec_T>
	void SimpleCPUSolver<mat_T, vec_T>::PrintResultB() {

		/* Now, multiply A by x and check if it equals b */
		MatMult(this->A, this->x, this->y); // y = A*x
		// Optional: Print the rhs
		PetscViewer vrhs;
		PetscViewerASCIIOpen(PETSC_COMM_SELF, "petsc_rhs.txt", &vrhs); // open an ASCII file for writing
		VecView(this->y, vrhs); // save the vector to the file
		PetscViewerDestroy(&vrhs); // close the file
	}

	template<class mat_T, class vec_T>
	void SimpleCPUSolver<mat_T, vec_T>::Finalize() {
		VecDestroy(&this->x);
		VecDestroy(&this->b);
		MatDestroy(&this->A);
		PetscFinalize();
	}
	template<class mat_T, class vec_T>
	KSPType SimpleCPUSolver<mat_T, vec_T>::GetKSPType(Algorithm alg) {
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
	template SimpleCPUSolver<PETSC_MAT, PETSC_VEC>;

}
