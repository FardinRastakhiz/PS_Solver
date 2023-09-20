#include "pch.h"
//
//#include "vector_factory.h"
//
//namespace ses {
//	void create_vector(int size, LocalType* values, VI_VEC& out_vec) {
//		//std::vector<int, PreAllocator<int>> cpu_vec(0, PreAllocator<int>(&values2[0], size));
//		std::vector<int> cpu_vec(size);
//		for (size_t i = 0; i < size; i++)
//		{
//			cpu_vec[i] = values[i];
//		}
//		viennacl::copy(cpu_vec.begin(), cpu_vec.end(), out_vec.begin());
//	}
//}
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
