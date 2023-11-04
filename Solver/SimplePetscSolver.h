#pragma once

#include "algorithms.h"
#include "IPreconditioner.h"

#include "PetscSolver.h"

namespace ses {
	template<typename mat_T, typename vec_T>
	class SimplePetscSolver : public PetscSolver<mat_T, vec_T>
	{
	public:
		SimplePetscSolver(SolverArgs args);
		void Solve(int iteration_count = 100, LocalType precision = 1e-4) override;
		void SetLocalTypes(SolverArgs args) override;
		void PrintX();
		void PrintResultB();
		void Finalize();
		LocalType* GetResult() override;
		vec_T u;
		vec_T y;
	private:
		KSPType GetKSPType(Algorithm alg);
	};
}