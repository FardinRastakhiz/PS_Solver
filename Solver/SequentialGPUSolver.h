#pragma once

#include "Algorithm.h"
#include "IPreconditioner.h"

#include "GPUSolver.h"


namespace ses {


	template<typename mat_T, typename vec_T>
	class SequentialGPUSolver : public GPUSolver<mat_T, vec_T>
	{

	public:
		SequentialGPUSolver(mat_T A, vec_T b, Algorithm algorithm, IPreconditioner preconditioner = DummyPreconditioner())
			: GPUSolver<mat_T, vec_T>(A, b, algorithm, preconditioner)
		{

		}

		void Solve(int iteration_count = 100, ScalarType precision = 1e-4) override {
			this->x = GPUSolver<mat_T, vec_T>::Solve(this->A, this->b, this->algorithm, this->preconditioner);
		}

		void Solve(vec_T b, int iteration_count = 100, ScalarType precision = 1e-4) {
			throw "Not Implemented";
		}
	};
}