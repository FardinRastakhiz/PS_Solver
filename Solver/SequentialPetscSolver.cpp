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
		SetLocalTypes(args);
	}
	template<class mat_T, class vec_T>
	void SequentialPetscSolver<mat_T, vec_T>::SetLocalTypes(SolverArgs args) {
		this->args = args;
		this->initialized = 0;
	}
	template<class mat_T, class vec_T>
	void SequentialPetscSolver<mat_T, vec_T>::Initialize() {
		this->iter_count = 0;
		PetscErrorCode ierr;
		ierr = PetscInitialize(PETSC_NULL, PETSC_NULL, PETSC_NULL, PETSC_NULL);
		PetscPrintf(PETSC_COMM_WORLD, "PETSC Initialized \n");
		create_matrix(this->args.num_rows, this->args.num_cols, this->args.nnz, this->args.row_indices, this->args.col_indices, this->s_values, this->A);
		create_petsc_vector(this->args.num_rows, this->b);
		create_petsc_vector(this->args.num_rows, this->x);
		fill_petsc_vector(this->args.num_rows, this->s_b, this->b);
		this->algorithm = this->args.algorithm;
	}


	template<typename mat_T, typename vec_T>
	void SequentialPetscSolver<mat_T, vec_T>::Solve(int iteration_count, LocalType precision) {
		PetscErrorCode ierr;
		PetscMPIInt    size, rank;
		KSP            ksp = NULL;           /* linear solver context */
		PetscReal      norm;          /* norm of solution error */
		PC             pc;


		if (this->initialized == 0) {
			ierr = KSPCreate(PETSC_COMM_WORLD, &ksp);
			ierr = KSPSetType(ksp, GetKSPType(this->algorithm));
			KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
			ierr = KSPSetOperators(ksp, this->A, this->A);
			ierr = KSPSetFromOptions(ksp);
		}
		else {
			ksp = this->ksp;
		}
		auto start = high_resolution_clock::now();
		ierr = KSPSolve(ksp, this->b, this->x);
		auto stop = high_resolution_clock::now();
		auto duration = duration_cast<microseconds>(stop - start);
		PetscPrintf(PETSC_COMM_WORLD, "Time taken by solver: %d microseconds \n", duration.count());
		PetscInt its;
		KSPGetIterationNumber(ksp, &its);
		PetscPrintf(PETSC_COMM_WORLD, "Solved With %d Iterations \n", its);
		this->ksp = ksp;
		if(this->initialized == 0) PetscPrintf(PETSC_COMM_WORLD, "Solved For The First Time");
		this->initialized = 1;
		this->iter_count += 1;
	/*	this->PrintX();
		this->PrintA();
		this->PrintB();*/
		PetscPrintf(PETSC_COMM_WORLD, "Solved %d times \n", iter_count);
		
	}
	// call below function always before solve function
	template<class mat_T, class vec_T>
	void SequentialPetscSolver<mat_T, vec_T>::SetOptions(PetscBackend backend, int platform, int device, int num_threads, int iteration_count, LocalType precision , int preconditioner) {
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
	template<typename mat_T, typename vec_T>
	void SequentialPetscSolver<mat_T, vec_T>::SetNewB(LocalType* b) {
		create_petsc_vector(this->args.num_rows, this->b);
		fill_petsc_vector(this->args.num_rows, b, this->b);
	}
	template<class mat_T, class vec_T>
	LocalType* SequentialPetscSolver<mat_T, vec_T>::GetResult() {
		PetscScalar* a;
		VecGetArray(this->x, &a);
		LocalType* res = (LocalType*)malloc(sizeof(LocalType) * this->args.num_rows);
		for (int i = 0; i < this->args.num_rows; i++) {
			res[i] = (LocalType)a[i];
		}
		return res;
	}

	template<class mat_T, class vec_T>
	void SequentialPetscSolver<mat_T, vec_T>::PrintX() {
		// Optional: Print the solution
		PetscViewer viewer;
		PetscViewerASCIIOpen(PETSC_COMM_SELF, "petsc_x.txt", &viewer); // open an ASCII file for writing
		if(this->iter_count == 1)
			PetscViewerASCIIOpen(PETSC_COMM_SELF, "petsc_x_1.txt", &viewer); // open an ASCII file for writing
		if (this->iter_count == 2)
			PetscViewerASCIIOpen(PETSC_COMM_SELF, "petsc_x_2.txt", &viewer); // open an ASCII file for writing
		if (this->iter_count == 3)
			PetscViewerASCIIOpen(PETSC_COMM_SELF, "petsc_x_3.txt", &viewer); // open an ASCII file for writing
		if (this->iter_count == 4)
			PetscViewerASCIIOpen(PETSC_COMM_SELF, "petsc_x_4.txt", &viewer); // open an ASCII file for writing
		if (this->iter_count == 5)
			PetscViewerASCIIOpen(PETSC_COMM_SELF, "petsc_x_5.txt", &viewer); // open an ASCII file for writing
		if (this->iter_count == 6)
			PetscViewerASCIIOpen(PETSC_COMM_SELF, "petsc_x_6.txt", &viewer); // open an ASCII file for writing
		if (this->iter_count == 7)
			PetscViewerASCIIOpen(PETSC_COMM_SELF, "petsc_x_7.txt", &viewer); // open an ASCII file for writing
		VecView(this->x, viewer); // save the vector to the file
		PetscViewerDestroy(&viewer); // close the file
	}

	template<class mat_T, class vec_T>
	void SequentialPetscSolver<mat_T, vec_T>::PrintA() {
		// Optional: Print the solution
		PetscViewer viewer;
		PetscViewerASCIIOpen(PETSC_COMM_SELF, "petsc_a.txt", &viewer); // open an ASCII file for writing
		if (this->iter_count == 1)
			PetscViewerASCIIOpen(PETSC_COMM_SELF, "petsc_a_1.txt", &viewer); // open an ASCII file for writing
		if (this->iter_count == 2)
			PetscViewerASCIIOpen(PETSC_COMM_SELF, "petsc_a_2.txt", &viewer); // open an ASCII file for writing
		if (this->iter_count == 3)
			PetscViewerASCIIOpen(PETSC_COMM_SELF, "petsc_a_3.txt", &viewer); // open an ASCII file for writing
		if (this->iter_count == 4)
			PetscViewerASCIIOpen(PETSC_COMM_SELF, "petsc_a_4.txt", &viewer); // open an ASCII file for writing
		if (this->iter_count == 5)
			PetscViewerASCIIOpen(PETSC_COMM_SELF, "petsc_a_5.txt", &viewer); // open an ASCII file for writing
		if (this->iter_count == 6)
			PetscViewerASCIIOpen(PETSC_COMM_SELF, "petsc_a_6.txt", &viewer); // open an ASCII file for writing
		if (this->iter_count == 7)
			PetscViewerASCIIOpen(PETSC_COMM_SELF, "petsc_a_7.txt", &viewer); // open an ASCII file for writing
		MatView(this->A, viewer); // save the vector to the file
		PetscViewerDestroy(&viewer); // close the file
	}

	template<class mat_T, class vec_T>
	void SequentialPetscSolver<mat_T, vec_T>::PrintB() {
		// Optional: Print the solution
		PetscViewer viewer;
		PetscViewerASCIIOpen(PETSC_COMM_SELF, "petsc_b.txt", &viewer); // open an ASCII file for writing
		if (this->iter_count == 1)
			PetscViewerASCIIOpen(PETSC_COMM_SELF, "petsc_b_1.txt", &viewer); // open an ASCII file for writing
		if (this->iter_count == 2)
			PetscViewerASCIIOpen(PETSC_COMM_SELF, "petsc_b_2.txt", &viewer); // open an ASCII file for writing
		if (this->iter_count == 3)
			PetscViewerASCIIOpen(PETSC_COMM_SELF, "petsc_b_3.txt", &viewer); // open an ASCII file for writing
		if (this->iter_count == 4)
			PetscViewerASCIIOpen(PETSC_COMM_SELF, "petsc_b_4.txt", &viewer); // open an ASCII file for writing
		if (this->iter_count == 5)
			PetscViewerASCIIOpen(PETSC_COMM_SELF, "petsc_b_5.txt", &viewer); // open an ASCII file for writing
		if (this->iter_count == 6)
			PetscViewerASCIIOpen(PETSC_COMM_SELF, "petsc_b_6.txt", &viewer); // open an ASCII file for writing
		if (this->iter_count == 7)
			PetscViewerASCIIOpen(PETSC_COMM_SELF, "petsc_b_7.txt", &viewer); // open an ASCII file for writing
		VecView(this->b, viewer); // save the vector to the file
		PetscViewerDestroy(&viewer); // close the file
	}
	template<class mat_T, class vec_T>
	void SequentialPetscSolver<mat_T, vec_T>::Finalize() {
		VecDestroy(&this->x);
		VecDestroy(&this->b);
		MatDestroy(&this->A);
		PetscFinalize();
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