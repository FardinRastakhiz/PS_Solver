#pragma once

#include "Algorithm.h"
#include "IPreconditioner.h"

#include "Solver.h"

namespace ses 
{
	template<typename mat_T, typename vec_T>
	class GPUSolver : public Solver<mat_T, vec_T>
	{
	public:

		GPUSolver(Algorithm algorithm, IPreconditioner preconditioner = DummyPreconditioner())
			: Solver<mat_T, vec_T>(algorithm, preconditioner) {}

		GPUSolver(mat_T A, vec_T b, Algorithm algorithm, IPreconditioner preconditioner = DummyPreconditioner())
			: Solver<mat_T, vec_T>(A, b, algorithm, preconditioner) {}

		static vec_T Solve(mat_T A, vec_T b, Algorithm algorithm, IPreconditioner preconditioner = DummyPreconditioner()) {
			throw "Not Implemented Exception";
		}
	};
}