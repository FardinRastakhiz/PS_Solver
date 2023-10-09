#pragma once
#ifndef ALGORITHMS_CLASSES

#include "IPreconditioner.h"
#include <exception>

namespace ses {
	enum Algorithm
	{
		CG,
		PCG,
		BIPCG,
		GMRES
	};

	struct AlgorithmClass {};

	struct CGAlgorithm : AlgorithmClass {};
	struct BIPCGAlgorithm : AlgorithmClass {};
	struct GMRESAlgorithm : AlgorithmClass {};
	struct PCGAlgorithm : AlgorithmClass {};


	AlgorithmClass get_algorithm_class(Algorithm algorithm);
}

#define ALGORITHMS_CLASSES
#endif