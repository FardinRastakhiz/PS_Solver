#pragma once

#include "CPUSolver.h"


namespace ses {


	template<typename mat_T, typename vec_T>
	class SequentialCPUSolver : public CPUSolver<mat_T, vec_T>
	{
	public:
		SequentialCPUSolver(SolverArgs args);
		void Solve(int iteration_count = 100, LocalType precision = 1e-4) override;
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