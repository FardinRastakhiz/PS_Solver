#include "pch.h"
#include "Solver.h"

#include "algorithms.h"
#include "IPreconditioner.h"
#include "utilities.h"
namespace ses {
	SolverArgs::SolverArgs(){}
	SolverArgs::SolverArgs(int num_rows, int num_cols, int nnz, int* row_indices,
		int* col_indices, double* values, double* b, Algorithm algorithm,
		IPreconditioner preconditioner) {
		this->num_rows = num_rows;
		this->num_cols = num_cols;
		this->nnz = nnz;
		this->row_indices = row_indices;
		this->col_indices = col_indices;
		this->values = values;
		this->b = b;
		this->x = x;
		this->algorithm = algorithm;
		this->preconditioner = preconditioner;
	}

	ISolver::ISolver(SolverArgs args) { 
		this->args = args; 
	}

	SolverArgs ISolver::GetArgs() {
		return this->args;
	}

	void ISolver::Solve(int iteration_count, LocalType precision) { throw std::exception("Not Implemented Exception"); }
	LocalType* ISolver::GetResult() { throw std::exception("Not Implemented Exception"); }

	LocalType* ISolver::CalculateB() { throw std::exception("Not Implemented Exception"); }
	void ISolver::SetPlatform() { }
	void ISolver::SetLocalTypes(SolverArgs args) { this->args = args; }

	template<typename mat_T, typename vec_T>
	Solver<mat_T, vec_T>::Solver(SolverArgs args) : ISolver(args) {
		SetLocalTypes(args);
	}
	template<typename mat_T, typename vec_T>
	void Solver<mat_T, vec_T>::SetLocalTypes(SolverArgs args) {
		s_values = cast_to_local(args.values, args.nnz);
		s_b = cast_to_local(args.b, args.num_rows);
		s_x = cast_to_local(args.x, args.num_cols);

		std::cout << s_values[0] << std::endl;
		std::cout << s_b[1] << std::endl;
		//std::cout << s_values[0] << std::endl;
		//std::cout << s_b[10099] << std::endl;
	}


	/*			constructor			*/
	template Solver<VI_COMP_Mat, VI_VEC>;
	template Solver<VI_COO_Mat, VI_VEC>;
	template Solver<VI_SELL_MAT, VI_VEC>;
	template Solver<PETSC_MAT, PETSC_VEC>;
}
