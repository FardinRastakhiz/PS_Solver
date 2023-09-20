#pragma once

#include "Algorithm.h"
#include "IPreconditioner.h"

#include "Solver.h"

namespace ses
{
	template<typename mat_T, typename vec_T>
	class CPUSolver : public Solver<mat_T, vec_T>
	{
	public:
		CPUSolver(SolverArgs args) :
			Solver<mat_T, vec_T>(args) {}

		//CPUSolver(Algorithm algorithm, IPreconditioner preconditioner = DummyPreconditioner())
		//	: Solver<mat_T, vec_T>(algorithm, preconditioner) {}

		//CPUSolver(mat_T A, vec_T b, Algorithm algorithm, IPreconditioner preconditioner = DummyPreconditioner())
		//	: Solver<mat_T, vec_T>(A, b, algorithm, preconditioner) {}

		static double* Solve(mat_T A, vec_T b, Algorithm algorithm, IPreconditioner preconditioner = DummyPreconditioner()) {
			throw "Not Implemented Exception";
		}
	};
}