#pragma once

#include "algorithms.h"
#include "IPreconditioner.h"

#include "Solver.h"

namespace ses
{
	template<typename mat_T, typename vec_T>
	class ViennaSolver : public Solver<mat_T, vec_T>
	{
	public:
		ViennaSolver(SolverArgs args) :
			Solver<mat_T, vec_T>(args) {	}


	};
}