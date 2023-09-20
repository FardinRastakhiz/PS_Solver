#pragma once

#include "Algorithm.h"
#include "IPreconditioner.h"

#include "CPUSolver.h"


namespace ses {


	template<typename mat_T, typename vec_T>
	class SequentialCPUSolver : public CPUSolver<mat_T, vec_T>
	{
	private:
		//template<typename MatrixT, typename VectorT>
		//struct monitor_user_data
		//{
		//	monitor_user_data(MatrixT const& A, VectorT const& b, VectorT const& guess) : A_ptr(&A), b_ptr(&b), guess_ptr(&guess) {}
		//	MatrixT const* A_ptr;
		//	VectorT const* b_ptr;
		//	VectorT const* guess_ptr;
		//};

		//template<typename VectorT, typename NumericT, typename MatrixT>
		//bool my_custom_monitor(VectorT const& current_approx, NumericT residual_estimate, void* user_data)
		//{
		//	// Extract residual:
		//	monitor_user_data<MatrixT, VectorT> const* data = reinterpret_cast<monitor_user_data<MatrixT, VectorT> const*>(user_data);
		//	// Form residual r = b - A*x, taking an initial guess into account: r = b - A * (current_approx + x_initial)
		//	VectorT x = current_approx + *data->guess_ptr;
		//	VectorT residual = *data->b_ptr - viennacl::linalg::prod(*data->A_ptr, x);
		//	VectorT initial_residual = *data->b_ptr - viennacl::linalg::prod(*data->A_ptr, *data->guess_ptr);
		//	//std::cout << "Residual estimate vs. true residual: " << residual_estimate << " vs. " << viennacl::linalg::norm_2(residual) / viennacl::linalg::norm_2(initial_residual) << std::endl;
		//	return false; // no termination of iteration
		//}

	public:
		SequentialCPUSolver(SolverArgs args) :
			CPUSolver<mat_T, vec_T>(args) {}

		//SequentialCPUSolver(mat_T A, vec_T b, Algorithm algorithm, IPreconditioner preconditioner = DummyPreconditioner())
		//	: CPUSolver<mat_T, vec_T>(algorithm, preconditioner)
		//{
		//	this.A = A;
		//}

		vec_T Solve(vec_T b);

		void Solve(int iteration_count = 100, LocalType precision = 1e-4) override {
			this->x = CPUSolver<mat_T, vec_T>::Solve(this->A, this->b, this->algorithm, this->preconditioner);
		}

		void Solve(vec_T b, int iteration_count = 100, LocalType precision = 1e-4) {
			throw "Not Implemented";
		}
	};
}