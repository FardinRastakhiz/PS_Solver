#pragma once

#include "GPUSolver.h"

#include "solver_factory.h"

namespace ses {


	template<typename mat_T, typename vec_T>
	class SequentialGPUSolver : public GPUSolver<mat_T, vec_T>
	{
	private:
		SolverFactory<vec_T> solverFactory;
		void helper_solver();
	public:
		SequentialGPUSolver(SolverArgs args);
		void SetPlatform() override;
		void SetLocalTypes(SolverArgs args) override;
		void Solve(int iteration_count = 100, LocalType precision = 1e-4) override;
		void Solve(LocalType* b, int iteration_count = 100, LocalType precision = 1e-4);
		LocalType* GetResult() override;
		LocalType* CalculateB() override;
		void PrintActiveDevice();

	};
}