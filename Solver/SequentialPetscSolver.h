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
		void Solve(int iteration_count = 100, LocalType precision = 1e-4) override;
		void SetOptions(PetscBackend backend, int platform = 0, int device = 0, int num_thread = 4);
		void Solve(vec_T b, int iteration_count = 100, LocalType precision = 1e-4);
		void SetNewB(LocalType* b);
		vec_T b;
		vec_T x;
		LocalType* s_b;
		KSP ksp;
	private:
		KSPType GetKSPType(Algorithm alg);
	};
}