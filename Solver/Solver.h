#pragma once

#include "Algorithm.h"
#include "IPreconditioner.h"
#include "SharedTypes.h"

namespace ses {

	class ISolver {
	public:
		size_t nnz = 0;
		size_t num_rows = 0;
		size_t num_cols = 0;
		ISolver(){}
		virtual void Solve(int iteration_count = 100, ScalarType precision = 1e-4) { throw "Not Implemented Exception"; }
		virtual ScalarType* GetResult() { throw "Not Implemented Exception"; }
	};



	template<typename mat_T, typename vec_T>
	class Solver: public ISolver
	{
	protected:
		mat_T A;
		vec_T b;
		vec_T x;
		Algorithm algorithm;
		IPreconditioner preconditioner;
		bool Solved = false;
	public:
		Solver(Algorithm algorithm, IPreconditioner preconditioner = DummyPreconditioner()) {
			algorithm = algorithm;
			preconditioner = preconditioner;
		}

		Solver(mat_T A_arg, vec_T b_arg, Algorithm algorithm, IPreconditioner preconditioner = DummyPreconditioner()) {
			A = A_arg;
			b = b_arg;
			algorithm = algorithm;
			preconditioner = preconditioner;
		}


		/*virtual void Set(mat_T A);
		virtual void Set(vec_T b);
		virtual void Set(Algorithm algorithm);
		virtual void Set(IPreconditioner preconditioner);

		~Solver();*/
	};



}
