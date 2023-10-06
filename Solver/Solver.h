#pragma once
#include "Algorithm.h"
#include "IPreconditioner.h"
#include "utilities.h"
namespace ses {
	// Argument Types
	struct SolverArgs {
	public:
		int num_rows, num_cols, nnz;
		int* row_indices, *col_indices;
		double* values, * b, * x;
		Algorithm algorithm;
		IPreconditioner preconditioner;
		SolverArgs();
		SolverArgs(int num_rows, int num_cols, int nnz, int* row_indices,
			int* col_indices, double* values, double* b, Algorithm algorithm,
			IPreconditioner preconditioner = DummyPreconditioner());
	};

	class ISolver {
	protected:
		SolverArgs args;
		AlgorithmClass algorithmClass;
	public:
		ISolver(SolverArgs args);
		SolverArgs GetArgs();
		AlgorithmClass GetAlgorithmClass();
		virtual void Solve(int iteration_count = 100, LocalType precision = 1e-4);
		virtual LocalType* GetResult();
		virtual LocalType* CalculateB();
		virtual void SetPlatform();
		virtual void SetLocalTypes(SolverArgs args);
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
		LocalType* s_values;
		LocalType* s_b;
		LocalType* s_x;
	public:
		Solver(SolverArgs args);
		void SetLocalTypes(SolverArgs args) override;
	};
}
