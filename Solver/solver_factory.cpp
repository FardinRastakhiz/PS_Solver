#include "pch.h"
#include "solver_factory.h"

namespace ses {

	/*			Tag Factory			*/
	TagFactory::TagFactory() {}

	template<typename algo_T>
	TagFactory::TagFactory(algo_T algorithm, LocalType tolerance,
		unsigned int iteration, unsigned int krilov_dim) {
		this->CreateTag(algorithm, tolerance, iteration, krilov_dim);
	}

	void TagFactory::CreateTag(GMRESAlgorithm algorithmIdentifier, LocalType tolerance,
		unsigned int iteration, unsigned int krilov_dim) {
		this->gmres_tag = new viennacl::linalg::gmres_tag(tolerance, iteration, krilov_dim);
	}
	void TagFactory::CreateTag(CGAlgorithm algorithmIdentifier, LocalType tolerance,
		unsigned int iteration, unsigned int krilov_dim) {
		this->cg_tag = new viennacl::linalg::cg_tag(tolerance, iteration);
	}
	void TagFactory::CreateTag(BIPCGAlgorithm algorithmIdentifier, LocalType tolerance,
		unsigned int iteration, unsigned int krilov_dim) {
		this->bicgstab_tag = new viennacl::linalg::bicgstab_tag(tolerance, iteration);
	}
	void TagFactory::CreateTag(PCGAlgorithm algorithmIdentifier, LocalType tolerance,
		unsigned int iteration, unsigned int krilov_dim) {
		this->cg_tag = new viennacl::linalg::cg_tag(tolerance, iteration);
	}

	viennacl::linalg::gmres_tag TagFactory::GetTag(GMRESAlgorithm algorithmIdentifier) {
		return *this->gmres_tag;
	}
	viennacl::linalg::cg_tag TagFactory::GetTag(CGAlgorithm algorithmIdentifier) {
		return *this->cg_tag;
	}
	viennacl::linalg::bicgstab_tag TagFactory::GetTag(BIPCGAlgorithm algorithmIdentifier) {
		return *this->bicgstab_tag;
	}
	viennacl::linalg::cg_tag TagFactory::GetTag(PCGAlgorithm algorithmIdentifier) {
		return *this->cg_tag;
	}

	/*			Solver Factory			*/
	template<typename vec_T>
	SolverFactory<vec_T>::SolverFactory() {}

	template<typename vec_T>
	template<typename algo_T>
	SolverFactory<vec_T>::SolverFactory(algo_T algorithm, TagFactory tag_factory): tag_factory(tag_factory){
		this->tag_factory = tag_factory;

		this->CreateSolver(algorithm);
	}


	template<typename vec_T>
	void SolverFactory<vec_T>::CreateSolver(GMRESAlgorithm algorithm) {
		this->gmres_solver = new viennacl::linalg::gmres_solver<vec_T>(this->tag_factory.GetTag(algorithm));
	}
	template<typename vec_T>
	void SolverFactory<vec_T>::CreateSolver(CGAlgorithm algorithm) {
		this->cg_solver = new viennacl::linalg::cg_solver<vec_T>(this->tag_factory.GetTag(algorithm));
	}
	template<typename vec_T>
	void SolverFactory<vec_T>::CreateSolver(BIPCGAlgorithm algorithm) {
		this->bicgstab_solver = new viennacl::linalg::bicgstab_solver<vec_T>(this->tag_factory.GetTag(algorithm));
	}
	template<typename vec_T>
	void SolverFactory<vec_T>::CreateSolver(PCGAlgorithm algorithm) {
		this->cg_solver = new viennacl::linalg::cg_solver<vec_T>(this->tag_factory.GetTag(algorithm));
	}



	template<typename vec_T>
	viennacl::linalg::gmres_solver<vec_T> SolverFactory<vec_T>::GetSolver(GMRESAlgorithm algorithmIdentifier) {
		return *this->gmres_solver;
	}
	template<typename vec_T>
	viennacl::linalg::cg_solver<vec_T> SolverFactory<vec_T>::GetSolver(CGAlgorithm algorithmIdentifier) {
		return *this->cg_solver;
	}
	template<typename vec_T>
	viennacl::linalg::bicgstab_solver<vec_T> SolverFactory<vec_T>::GetSolver(BIPCGAlgorithm algorithmIdentifier) {
		return *this->bicgstab_solver;
	}

	template<typename vec_T>
	viennacl::linalg::cg_solver<vec_T> SolverFactory<vec_T>::GetSolver(PCGAlgorithm algorithmIdentifier) {
		return *this->cg_solver;
	}


	/*			constructor			*/


	template SolverFactory<VI_VEC>;

	template SolverFactory<VI_VEC>::SolverFactory(GMRESAlgorithm algorithm, TagFactory tag_factory);
	template SolverFactory<VI_VEC>::SolverFactory(CGAlgorithm algorithm, TagFactory tag_factory);
	template SolverFactory<VI_VEC>::SolverFactory(BIPCGAlgorithm algorithm, TagFactory tag_factory);
	template SolverFactory<VI_VEC>::SolverFactory(PCGAlgorithm algorithm, TagFactory tag_factory);

	template TagFactory::TagFactory(GMRESAlgorithm algorithm, LocalType tolerance, unsigned int iteration, unsigned int krilov_dim);
	template TagFactory::TagFactory(CGAlgorithm algorithm, LocalType tolerance, unsigned int iteration, unsigned int krilov_dim);
	template TagFactory::TagFactory(BIPCGAlgorithm algorithm, LocalType tolerance, unsigned int iteration, unsigned int krilov_dim);
	template TagFactory::TagFactory(PCGAlgorithm algorithm, LocalType tolerance, unsigned int iteration, unsigned int krilov_dim);
	//template SolverFactory<GMRESAlgorithm>;
	//template SolverFactory<CGAlgorithm>;
	//template SolverFactory<BIPCGAlgorithm>;
	//template SolverFactory<PCGAlgorithm>;


}