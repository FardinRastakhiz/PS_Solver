#pragma once

#include "algorithms.h"
#include "IPreconditioner.h"

#include "Solver.h"
enum PetscBackend {
	OPENCL = 1,
	OPENMP = 2,
	NORMAL = 3
};
namespace ses
{
	template<typename mat_T, typename vec_T>
	class PetscSolver : public Solver<mat_T, vec_T>
	{
	public:
		PetscSolver(SolverArgs args) :
			Solver<mat_T, vec_T>(args) {}
	};
}