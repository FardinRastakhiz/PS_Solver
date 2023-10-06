#pragma once
#include "pch.h"
#include "solver_factory.h"

namespace ses {

	/*			Tag Factory			*/
	TagFactory::TagFactory() {}

	TagFactory::TagFactory(GMRESAlgorithm algorithm, LocalType tolerance,
		unsigned int iteration, unsigned int krilov_dim) {
		this->CreateTag(algorithm, tolerance, iteration, krilov_dim);
	}

	void TagFactory::CreateTag(GMRESAlgorithm algorithmIdentifier, LocalType tolerance,
		unsigned int iteration, unsigned int krilov_dim) {
		this->gmres_tag = viennacl::linalg::gmres_tag(tolerance, iteration, krilov_dim);
	}
	void TagFactory::CreateTag(CGAlgorithm algorithmIdentifier, LocalType tolerance,
		unsigned int iteration, unsigned int krilov_dim) {
		this->cg_tag = viennacl::linalg::cg_tag(tolerance, iteration);
	}
	void TagFactory::CreateTag(BIPCGAlgorithm algorithmIdentifier, LocalType tolerance,
		unsigned int iteration, unsigned int krilov_dim) {
		this->bicgstab_tag = viennacl::linalg::bicgstab_tag(tolerance, iteration);
	}
	void TagFactory::CreateTag(PCGAlgorithm algorithmIdentifier, LocalType tolerance,
		unsigned int iteration, unsigned int krilov_dim) {
		this->cg_tag = viennacl::linalg::cg_tag(tolerance, iteration);
	}

	viennacl::linalg::gmres_tag TagFactory::GetTag(GMRESAlgorithm algorithmIdentifier) {
		return this->gmres_tag;
	}
	viennacl::linalg::cg_tag TagFactory::GetTag(CGAlgorithm algorithmIdentifier) {
		return this->cg_tag;
	}
	viennacl::linalg::bicgstab_tag TagFactory::GetTag(BIPCGAlgorithm algorithmIdentifier) {
		return this->bicgstab_tag;
	}
	viennacl::linalg::cg_tag TagFactory::GetTag(PCGAlgorithm algorithmIdentifier) {
		return this->cg_tag;
	}

	/*			Solver Factory			*/
	template<typename vec_T>
	SolverFactory<vec_T>::SolverFactory() {}

	template<typename vec_T>
	SolverFactory<vec_T>::SolverFactory(AlgorithmClass algorithm, TagFactory tagFactory) {
		this->tag_factory = TagFactory(algorithm, tolerance, iteration, krilov_dim);
		this->CreateSolver(algorithm, tagFactory);
	}

	
	template<typename vec_T>
	void SolverFactory<vec_T>::CreateSolver(GMRESAlgorithm algorithm, TagFactory tagFactory) {
		this->gmres_solver = viennacl::linalg::gmres_solver(this->tag_factory.GetTag(algorithm));
	}
	template<typename vec_T>
	void SolverFactory<vec_T>::CreateSolver(CGAlgorithm algorithm, TagFactory tagFactory) {
		this->cg_solver = viennacl::linalg::cg_solver(this->tag_factory.GetTag(algorithm));
	}
	template<typename vec_T>
	void SolverFactory<vec_T>::CreateSolver(BIPCGAlgorithm algorithm, TagFactory tagFactory) {
		this->bicgstab_solver = viennacl::linalg::bicgstab_solver(this->tag_factory.GetTag(algorithm));
	}
	template<typename vec_T>
	void SolverFactory<vec_T>::CreateSolver(PCGAlgorithm algorithm, TagFactory tagFactory) {
		this->cg_solver = viennacl::linalg::cg_solver(this->tag_factory.GetTag(algorithm));
	}



	template<typename vec_T>
	viennacl::linalg::gmres_solver<vec_T> SolverFactory<vec_T>::GetSolver(GMRESAlgorithm algorithmIdentifier) {
		return this->gmres_solver;
	}
	template<typename vec_T>
	viennacl::linalg::cg_solver<vec_T> SolverFactory<vec_T>::GetSolver(CGAlgorithm algorithmIdentifier) {
		return this->cg_solver;
	}
	template<typename vec_T>
	viennacl::linalg::bicgstab_solver<vec_T> SolverFactory<vec_T>::GetSolver(BIPCGAlgorithm algorithmIdentifier) {
		return this->bicgstab_solver;
	}

	template<typename vec_T>
	viennacl::linalg::cg_solver<vec_T> SolverFactory<vec_T>::GetSolver(PCGAlgorithm algorithmIdentifier) {
		return this->cg_solver;
	}


	/*			constructor			*/
	template SolverFactory<VI_VEC>;
}