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
		GPUSolver(SolverArgs args) :
			Solver<mat_T, vec_T>(args) {	}
	};
}