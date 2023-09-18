#pragma once

#include "Algorithm.h"
#include "IPreconditioner.h"

#include "GPUSolver.h"

namespace ses {
	template<class mat_T, class vec_T>
	class SimpleGPUSolver : public GPUSolver<mat_T, vec_T>
	{
	public:
		SimpleGPUSolver(mat_T A, vec_T b, Algorithm algorithm, IPreconditioner preconditioner = DummyPreconditioner())
			: GPUSolver<mat_T,vec_T>(A, b, algorithm, preconditioner) {	}

		void Solve(int iteration_count = 100, ScalarType precision = 1e-4) override {
			this->x = GPUSolver<mat_T, vec_T>::Solve(this->A, this->b, this->algorithm, this->preconditioner);
		}

		//ScalarType* GetResult() override {
		//	std::vector<ScalarType> vec(this->x.size());
		//	//viennacl::copy(this->x.begin(), this->x.end(), vec.begin());
		//	ScalarType* vec2 = &vec[0];
		//	return vec2;
		//}

		/*~SimpleGPUSolver() {

		}*/

	};
}