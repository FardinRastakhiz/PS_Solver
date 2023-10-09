#include "pch.h"
#include "vector_factory.h"

//#include "viennacl/forwards.h"
//#include "viennacl/vector.hpp"
//#include "viennacl/tools/tools.hpp"
//#include "viennacl/linalg/sparse_matrix_operations.hpp"

#include "petscksp.h"
#include "viennacl/scalar.hpp"
#include "viennacl/vector.hpp"
#include <vector>

namespace ses {
	void create_vector(int size, LocalType* values, VI_VEC& out_vec) {
		out_vec = VI_VEC(size);
		//std::vector<int, PreAllocator<int>> cpu_vec(0, PreAllocator<int>(&values2[0], size));
		std::vector<LocalType> cpu_vec(size);
		for (size_t i = 0; i < size; i++)
		{
			cpu_vec[i] = values[i];
		}
		viennacl::copy(cpu_vec.begin(), cpu_vec.end(), out_vec.begin());

		/*std::vector<LocalType> vec(out_vec.size());
		viennacl::copy(out_vec.begin(), out_vec.end(), vec.begin());
		LocalType* vec2 = &vec[0];
		std::cout << "vec2: " << vec2[10099] << std::endl;*/
	}
	void fill_petsc_vector(int size, LocalType* values, PETSC_VEC& out_vec) {
		// Set rhs values
		for (int i = 0; i < size; i++) { VecSetValue(out_vec, i, values[i], INSERT_VALUES); }
		VecAssemblyBegin(out_vec);
		VecAssemblyEnd(out_vec);
	}
	void create_petsc_vector(int size, PETSC_VEC& out_vec) {
		// Create vectors
		VecCreate(PETSC_COMM_WORLD, &out_vec);
		VecSetSizes(out_vec, PETSC_DECIDE, size);
		VecSetFromOptions(out_vec);
	}
	//void create_petsc_vectors(int iterations ,int size, std::vector<PETSC_VEC> out_vecs) {
	//	for (int i = 0; i < iterations; i++) {
	//		VecCreate(PETSC_COMM_WORLD, &(out_vecs[i]));
	//		VecSetSizes(out_vecs[i], PETSC_DECIDE, size);
	//		VecSetFromOptions(out_vecs[i]);
	//	}
	//}
	//void fill_petsc_vectors(int iterations, int size, LocalType** values, std::vector<PETSC_VEC> out_vecs) {
	//	for (int i = 0; i < iterations; i++) {
	//		for (int j = 0; j < size; j++) { VecSetValue(out_vecs[i], j, values[i][j], INSERT_VALUES); }
	//		VecAssemblyBegin(out_vecs[i]);
	//		VecAssemblyEnd(out_vecs[i]);
	//	}
}

////
////#ifndef PRE_ALLOCATOR
//
//	//
//	//template <typename T>
//	//class PreAllocator
//	//{
//	//private:
//	//	T*const memory_ptr;
//	//	std::size_t memory_size;
//	//
//	//public:
//	//	typedef std::size_t     size_type;
//	//	typedef T* pointer;
//	//	typedef T               value_type;
//	//
//	//	PreAllocator(T* memory_ptr, std::size_t memory_size) : memory_ptr(memory_ptr), memory_size(memory_size) {}
//	//
//	//	PreAllocator(const PreAllocator& other) throw() : memory_ptr(other.memory_ptr), memory_size(other.memory_size) {};
//	//
//	//	template<typename U>
//	//	friend class PreAllocator;
//	//
//	//	template<typename U>
//	//	PreAllocator(const PreAllocator<U>& other) throw() : memory_ptr(other.memory_ptr), memory_size(other.memory_size) {};
//	//
//	//	template<typename U>
//	//	PreAllocator& operator = (const PreAllocator<U>& other) { return *this; }
//	//	PreAllocator<T>& operator = (const PreAllocator& other) { return *this; }
//	//	~PreAllocator() {}
//	//
//	//
//	//	pointer allocate(size_type n, const void* hint = 0) { return memory_ptr; }
//	//	void deallocate(T* ptr, size_type n) {}
//	//
//	//	size_type max_size() const { return memory_size; }
//	//};
//
////
////#define PRE_ALLOCATOR
////#endif
