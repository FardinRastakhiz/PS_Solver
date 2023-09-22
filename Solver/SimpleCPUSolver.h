#pragma once

#include "Algorithm.h"
#include "IPreconditioner.h"

#include "CPUSolver.h"

namespace ses {
	template<typename mat_T, typename vec_T>
	class SimpleCPUSolver : public CPUSolver<mat_T, vec_T>
	{
	public:
		SimpleCPUSolver(SolverArgs args);
		void Solve(int iteration_count = 100, LocalType precision = 1e-4) override;
		void SetLocalTypes(SolverArgs args) override;
		void PrintX();
		void PrintResultB();
		void Finalize();
		LocalType* GetResult() override;
		vec_T u;
		vec_T y;
	};
}