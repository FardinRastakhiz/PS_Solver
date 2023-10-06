#pragma once

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


	AlgorithmClass get_algorithm_class(Algorithm algorithm) {

		switch (algorithm)
		{
		case Algorithm::CG:
			return CGAlgorithm();
		case Algorithm::GMRES:
			return GMRESAlgorithm();
		case Algorithm::BIPCG:
			return BIPCGAlgorithm();
		case Algorithm::PCG:
			return PCGAlgorithm();
		}
		throw std::exception("Not Implemented Exception");
	}
}

