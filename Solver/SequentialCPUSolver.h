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
		std::vector<vec_t> bs;
		std::vector<vec_t> xs;
		LocalType** s_bs;

	};
}