#pragma once

#include "algorithms.h"
#include "IPreconditioner.h"

#include "PetscSolver.h"

namespace ses {

	template<typename mat_T, typename vec_T>
	class SequentialPetscSolver : public PetscSolver<mat_T, vec_T>
	{
	public:
		SequentialPetscSolver(SolverArgs args);
		void Initialize();
		void SetLocalTypes(SolverArgs args) override;
		void Solve(int iteration_count = 100, LocalType precision = 1e-4) override;
		void SetOptions(PetscBackend backend, int platform = 0, int device = 0, int num_thread = 4, int iteration_count = 100, LocalType precision = 1e-4, int preconditioner = 1);
		void SetNewB(LocalType* b);
		void PrintX();
		void PrintA();
		void PrintB();
		void Finalize();
		LocalType* GetResult() override;
		int initialized;
		int iter_count;
		KSP ksp;
	private:
		KSPType GetKSPType(Algorithm alg);
	};
}