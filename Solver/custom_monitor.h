#pragma once


template<typename MatrixT, typename VectorT>
struct monitor_user_data
{
	monitor_user_data(MatrixT const& A, VectorT const& b, VectorT const& guess) : A_ptr(&A), b_ptr(&b), guess_ptr(&guess) {}
	MatrixT const* A_ptr;
	VectorT const* b_ptr;
	VectorT const* guess_ptr;
};

template<typename VectorT, typename NumericT, typename MatrixT>
bool my_custom_monitor(VectorT const& current_approx, NumericT residual_estimate, void* user_data)
{
	// Extract residual:
	monitor_user_data<MatrixT, VectorT> const* data = reinterpret_cast<monitor_user_data<MatrixT, VectorT> const*>(user_data);
	// Form residual r = b - A*x, taking an initial guess into account: r = b - A * (current_approx + x_initial)
	VectorT x = current_approx + *data->guess_ptr;
	VectorT residual = *data->b_ptr - viennacl::linalg::prod(*data->A_ptr, x);
	VectorT initial_residual = *data->b_ptr - viennacl::linalg::prod(*data->A_ptr, *data->guess_ptr);
	//std::cout << "Residual estimate vs. true residual: " << residual_estimate << " vs. " << viennacl::linalg::norm_2(residual) / viennacl::linalg::norm_2(initial_residual) << std::endl;
	return false; // no termination of iteration
}