#include "pch.h"
#include "algorithms.h"

namespace ses {


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
