#pragma once

#include "GPUSolver.h"


namespace ses {


	template<typename mat_T, typename vec_T>
	class SequentialGPUSolver : public GPUSolver<mat_T, vec_T>
	{
	public:
		SequentialGPUSolver(SolverArgs args);
		void Solve(int iteration_count = 100, LocalType precision = 1e-4) override;
		void Solve(vec_T b, int iteration_count = 100, LocalType precision = 1e-4);
	};
}