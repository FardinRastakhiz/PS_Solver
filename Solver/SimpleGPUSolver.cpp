#pragma once
#include "pch.h"
#include "SimpleGPUSolver.h"
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
	template<class mat_T, class vec_T>
	SimpleGPUSolver<mat_T, vec_T>::SimpleGPUSolver(SolverArgs args) :
		GPUSolver<mat_T, vec_T>(args) {
		SetPlatform();
		SetLocalTypes(args);
	}

	template<class mat_T, class vec_T>
	void SimpleGPUSolver<mat_T, vec_T>::SetPlatform() {
		std::cout << "SimpleGPUSolver.cpp - SetPlatform()" << std::endl;
		std::vector<viennacl::ocl::platform> platforms = viennacl::ocl::get_platforms();
		std::vector<viennacl::ocl::device> devices = platforms[2].devices(CL_DEVICE_TYPE_GPU);
		viennacl::ocl::setup_context(0, devices[0]);
	}

	template<class mat_T, class vec_T>
	void SimpleGPUSolver<mat_T, vec_T>::SetLocalTypes(SolverArgs args) {
		std::cout << "SimpleGPUSolver.cpp - SetLocalTypes(SolverArgs args)" << std::endl;
		// create target library vectors and matrices
		Solver<mat_T, vec_T>::SetLocalTypes(args);
		/*VI_SELL_MAT mat; VI_VEC vec;*/
		create_matrix(args.num_rows, args.num_cols, args.nnz, args.row_indices, args.col_indices, this->s_values, this->A);
		create_vector(args.num_rows, this->s_b, this->b);
		/*
		auto b_val = reinterpret_cast<VI_VEC>(this->b)(0);
		auto a_val = reinterpret_cast<VI_SELL_MAT>(this->A)(0, 0);
		std::cout << b_val << std::endl;
		std::cout << a_val << std::endl;*/
	}

	template<class mat_T, class vec_T>
	void SimpleGPUSolver<mat_T, vec_T>::Solve(int iteration_count, LocalType precision) {
		viennacl::linalg::gmres_tag my_tag(1.0E-10, 1000, 50U);
		viennacl::linalg::gmres_solver<viennacl::vector<LocalType>> _my_solver(my_tag);

		this->x = viennacl::zero_vector<LocalType>(this->b.size(), viennacl::traits::context(this->b));
		this->x = _my_solver(this->A, this->b);
	}
	std::vector<LocalType> vec_result;
	template<class mat_T, class vec_T>
	LocalType* SimpleGPUSolver<mat_T, vec_T>::GetResult() {
		vec_result = std::vector<LocalType>(this->x.size());
		viennacl::copy(this->x.begin(), this->x.end(), vec_result.begin());
		//LocalType* vec2 = &vec[0];
		return &vec_result[0];
	}

	std::vector<LocalType> b_result;
	template<class mat_T, class vec_T>
	LocalType* SimpleGPUSolver<mat_T, vec_T>::CalculateB() {
		VI_VEC new_b = viennacl::linalg::prod(this->A, this->x);
		b_result = std::vector<LocalType>(new_b.size());
		viennacl::copy(new_b.begin(), new_b.end(), b_result.begin());
		//LocalType* vec2 = &vec[0];
		return &b_result[0];
	}

	template<class mat_T, class vec_T>
	void SimpleGPUSolver<mat_T, vec_T>::PrintActiveDevice() {
		viennacl::ocl::device current_device = viennacl::ocl::current_device();
		std::cout << "current device is:" << current_device.name() << std::endl;
	}


	/*				implementation options				*/
	/*--------------------------------------------------*/


	/*			constructor			*/
	template SimpleGPUSolver<VI_COMP_Mat, VI_VEC>;
	template SimpleGPUSolver<VI_COO_Mat, VI_VEC>;
	template SimpleGPUSolver<VI_SELL_MAT, VI_VEC>;

}
