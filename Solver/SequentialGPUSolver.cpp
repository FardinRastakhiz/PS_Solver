#pragma once

#include "pch.h"
#include "SequentialGPUSolver.h"

//#include "Solver.h"
//#include "GPUSolver.h"




//#ifndef NDEBUG
//#define NDEBUG
//#endif

// Necessary to obtain a suitable performance in ublas
#ifndef NDEBUG
#define BOOST_UBLAS_NDEBUG
#endif

//
//#ifndef NDEBUG
//#define BOOST_UBLAS_NDEBUG
//#endif

#ifndef VIENNACL_WITH_OPENCL
#define VIENNACL_WITH_OPENCL
#endif


//
// ublas includes
//
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/operation_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/lu.hpp>


// IMPORTANT: Must be set prior to any ViennaCL includes if you want to use ViennaCL algorithms on Eigen objects
#define VIENNACL_WITH_UBLAS 1

// System headers
#include <iostream>
#include <vector>



#include "viennacl/forwards.h"
#include "viennacl/scalar.hpp"
#include "viennacl/vector.hpp"
#include "viennacl/matrix.hpp"
#include "viennacl/coordinate_matrix.hpp"
#include "viennacl/compressed_matrix.hpp"
#include "viennacl/ell_matrix.hpp"
#include "viennacl/hyb_matrix.hpp"
#include "viennacl/sliced_ell_matrix.hpp"

#include "viennacl/meta/result_of.hpp"
#include "viennacl/traits/size.hpp"

#include "viennacl/tools/timer.hpp"
#include "viennacl/tools/tools.hpp"

#include "viennacl/io/matrix_market.hpp"
#include "viennacl/ocl/backend.hpp"

#include "viennacl/linalg/cg.hpp"
#include "viennacl/linalg/bicgstab.hpp"
#include "viennacl/linalg/gmres.hpp"
#include "viennacl/linalg/norm_2.hpp"
#include "viennacl/linalg/power_iter.hpp"
#include "viennacl/linalg/prod.hpp"
#include "viennacl/linalg/ilu.hpp"
#include "viennacl/linalg/jacobi_precond.hpp"
#include "viennacl/linalg/matrix_operations.hpp"
#include "viennacl/linalg/sparse_matrix_operations.hpp"
#include <viennacl/linalg/opencl/kernels/compressed_matrix.hpp>

#ifndef VIENNACL_WITH_UBLAS
#define VIENNACL_WITH_UBLAS
#endif


#define BENCHMARK_RUNS          10


#include "Algorithm.h"
#include "IPreconditioner.h"
#include "matrix_factory.h"
#include "vector_factory.h"


namespace ses {


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


	template<typename mat_T, typename vec_T>
	SequentialGPUSolver<mat_T, vec_T>::SequentialGPUSolver(SolverArgs args) :
		GPUSolver<mat_T, vec_T>(args) {
		SetPlatform();
		SetLocalTypes(args);


	}

	template<class mat_T, class vec_T>
	void SequentialGPUSolver<mat_T, vec_T>::SetPlatform() {
		std::cout << "SimpleGPUSolver.cpp - SetPlatform()" << std::endl;
		std::vector<viennacl::ocl::platform> platforms = viennacl::ocl::get_platforms();
		std::vector<viennacl::ocl::device> devices = platforms[2].devices(CL_DEVICE_TYPE_GPU);
		viennacl::ocl::setup_context(0, devices[0]);
	}

	template<class mat_T, class vec_T>
	void SequentialGPUSolver<mat_T, vec_T>::SetLocalTypes(SolverArgs args) {
		std::cout << "SimpleGPUSolver.cpp - SetLocalTypes(SolverArgs args)" << std::endl;
		// create target library vectors and matrices
		Solver<mat_T, vec_T>::SetLocalTypes(args);
		/*VI_SELL_MAT mat; VI_VEC vec;*/
		create_matrix(args.num_rows, args.num_cols, args.nnz, args.row_indices, args.col_indices, this->s_values, this->A);
		create_vector(args.num_rows, this->s_b, this->b);
	}

	template<typename mat_T, typename vec_T>
	void SequentialGPUSolver<mat_T, vec_T>::helper_solver() {
		throw std::exception("Not Implemented Exception");
	}


	template<typename mat_T, typename vec_T>
	void SequentialGPUSolver<mat_T, vec_T>::Solve(int iteration_count, LocalType precision) {
		TagFactory tagFactory = ses::TagFactory(ISolver::GetAlgorithmClass(), 1.0E-10, 1000, 50U);
		this->solverFactory = ses::SolverFactory<vec_T>(ISolver::GetAlgorithmClass(), tagFactory);
		this->x = this->solverFactory.GetSolver(ISolver::GetAlgorithmClass())(this->A, this->b);
		//viennacl::linalg::gmres_tag my_tag(1.0E-10, 1000, 50U);
		//viennacl::linalg::gmres_solver<viennacl::vector<LocalType>> _my_solver(my_tag);

		//this->x = viennacl::zero_vector<LocalType>(this->b.size(), viennacl::traits::context(this->b));

		//this->x = _my_solver(this->A, this->b);
	}

	template<typename mat_T, typename vec_T>
	void SequentialGPUSolver<mat_T, vec_T>::Solve(LocalType* b, int iteration_count, LocalType precision) {
		//viennacl::linalg::gmres_solver<viennacl::vector<LocalType>> _my_solver(my_tag);
		monitor_user_data<mat_T, vec_T>	my_monitor_data(this->A, this->b, this->x);
		this->solverFactory.GetSolver(ISolver::algorithmClass)
			.set_monitor(my_custom_monitor<vec_T, LocalType, mat_T >, &my_monitor_data);
		this->solverFactory.GetSolver(ISolver::algorithmClass).set_initial_guess(this->x);
		this->x = this->solverFactory.GetSolver(ISolver::algorithmClass)(this->A, this->b);
		//this->x = this->solverFactory.GetGMRESSolver()(this->A, this->b);
	}

	std::vector<LocalType> vec_result;
	template<class mat_T, class vec_T>
	LocalType* SequentialGPUSolver<mat_T, vec_T>::GetResult() {
		vec_result = std::vector<LocalType>(this->x.size());
		viennacl::copy(this->x.begin(), this->x.end(), vec_result.begin());
		//LocalType* vec2 = &vec[0];
		return &vec_result[0];
	}

	std::vector<LocalType> b_result;
	template<class mat_T, class vec_T>
	LocalType* SequentialGPUSolver<mat_T, vec_T>::CalculateB() {
		VI_VEC new_b = viennacl::linalg::prod(this->A, this->x);
		b_result = std::vector<LocalType>(new_b.size());
		viennacl::copy(new_b.begin(), new_b.end(), b_result.begin());
		//LocalType* vec2 = &vec[0];
		return &b_result[0];
	}

	template<class mat_T, class vec_T>
	void SequentialGPUSolver<mat_T, vec_T>::PrintActiveDevice() {
		viennacl::ocl::device current_device = viennacl::ocl::current_device();
		std::cout << "current device is:" << current_device.name() << std::endl;
	}



	/*			constructor			*/
	template SequentialGPUSolver<VI_COMP_Mat, VI_VEC>;
	template SequentialGPUSolver<VI_COO_Mat, VI_VEC>;
	template SequentialGPUSolver<VI_SELL_MAT, VI_VEC>;
}










