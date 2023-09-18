#pragma once

#include "Algorithm.h"
#include "IPreconditioner.h"

#include "CPUSolver.h"

namespace ses {
	template<typename mat_T, typename vec_T>
	class SimpleCPUSolver : public CPUSolver<mat_T, vec_T>
	{
	public:
		SimpleCPUSolver(mat_T A, vec_T b, Algorithm algorithm, IPreconditioner preconditioner = DummyPreconditioner())
			: CPUSolver<mat_T, vec_T>(A, b, algorithm, preconditioner)
		{
		}

		vec_T Solve(int iteration_count = 100, ScalarType precision = 1e-4) override {
			this->x = CPUSolver<mat_T, vec_T>.Solve(this->A, this->b, this->algorithm, this->preconditioner);
			return this->x;
		}

	};
}