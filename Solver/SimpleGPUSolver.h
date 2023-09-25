#pragma once
#include "GPUSolver.h"

namespace ses {
	template<class mat_T, class vec_T>
	class SimpleGPUSolver : public GPUSolver<mat_T, vec_T>
	{
	public:
		SimpleGPUSolver(SolverArgs args);
		void SetPlatform() override;
		void SetLocalTypes(SolverArgs args) override;
		void Solve(int iteration_count = 100, LocalType precision = 1e-4) override;
		LocalType* GetResult() override;
		LocalType* CalculateB() override;
		void PrintActiveDevice();
	};
}