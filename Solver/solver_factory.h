#pragma once

#include "Algorithm.h"
#include "utilities.h"

#include "viennacl/linalg/cg.hpp"
#include "viennacl/linalg/bicgstab.hpp"
#include "viennacl/linalg/gmres.hpp"

namespace ses {
	class TagFactory {
	private:
		viennacl::linalg::gmres_tag gmres_tag;
		viennacl::linalg::cg_tag cg_tag;
		viennacl::linalg::bicgstab_tag bicgstab_tag;
		Algorithm algorithm;
	public:
		TagFactory();
		TagFactory(GMRESAlgorithm algorithm, LocalType tolerance,
			unsigned int iteration, unsigned int krilov_dim = 50U);

		void CreateTag(GMRESAlgorithm algorithmIdentifier, LocalType tolerance,
			unsigned int iteration, unsigned int krilov_dim);
		void CreateTag(CGAlgorithm algorithmIdentifier, LocalType tolerance,
			unsigned int iteration, unsigned int krilov_dim);
		void CreateTag(BIPCGAlgorithm algorithmIdentifier, LocalType tolerance,
			unsigned int iteration, unsigned int krilov_dim);
		void CreateTag(PCGAlgorithm algorithmIdentifier, LocalType tolerance,
			unsigned int iteration, unsigned int krilov_dim);

		viennacl::linalg::gmres_tag GetTag(GMRESAlgorithm algorithmIdentifier);
		viennacl::linalg::cg_tag GetTag(CGAlgorithm algorithmIdentifier);
		viennacl::linalg::bicgstab_tag GetTag(BIPCGAlgorithm algorithmIdentifier);
		viennacl::linalg::cg_tag GetTag(PCGAlgorithm algorithmIdentifier);
	};

	template<typename vec_T>
	class SolverFactory {
	private:
		viennacl::linalg::gmres_solver<vec_T> gmres_solver;
		viennacl::linalg::cg_solver<vec_T> cg_solver;
		viennacl::linalg::bicgstab_solver<vec_T> bicgstab_solver;
		Algorithm algorithm;
		TagFactory tag_factory;
	public:
		SolverFactory();
		SolverFactory(AlgorithmClass algorithm, TagFactory tagFactory);


		void CreateSolver(GMRESAlgorithm algorithm, TagFactory tagFactory);
		void CreateSolver(CGAlgorithm algorithm, TagFactory tagFactory);
		void CreateSolver(BIPCGAlgorithm algorithm, TagFactory tagFactory);
		void CreateSolver(PCGAlgorithm algorithm, TagFactory tagFactory);

		viennacl::linalg::gmres_solver<vec_T> GetSolver(GMRESAlgorithm algorithmIdentifier);
		viennacl::linalg::cg_solver<vec_T> GetSolver(CGAlgorithm algorithmIdentifier);
		viennacl::linalg::bicgstab_solver<vec_T> GetSolver(BIPCGAlgorithm algorithmIdentifier);
		viennacl::linalg::cg_solver<vec_T> GetSolver(PCGAlgorithm algorithmIdentifier);
	};
}