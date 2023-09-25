#include "pch.h"
#include "utilities.h"
#include <stdio.h>

namespace ses {


	/*				Save To Files			*/
	// Saving Vector!!
	template<typename T>
	void save_vector(viennacl::vector<T> x, const char* file_name) {
		size_t size = x.size();
		FILE* fp;
		fp = fopen(file_name, "w");
		for (unsigned i = 0; i < size; i++) {
			float output = x[i];
			fprintf(fp, "%d, %f\n", i, output);
		}
		fclose(fp);
	}
	
	void save_vector_pointer(double* x, size_t size, const char* file_name) {
		FILE* fp;
		fp = fopen(file_name, "w");
		for (unsigned i = 0; i < size; i++) {
			float output = (float)x[i];
			fprintf(fp, "%d, %f\n", i, output);
		}
		fclose(fp);
	}

	void save_vector_pointer(float* x, size_t size, const char* file_name) {
		FILE* fp;
		fp = fopen(file_name, "w");
		for (unsigned i = 0; i < size; i++) {
			float output = x[i];
			fprintf(fp, "%d, %f\n", i, output);
		}
		fclose(fp);
	}

	// implementing equivalent of "is" is c# for type checking
	template < class T, class U >
	bool is_inst_of(U u) {
		return dynamic_cast<T>(u) != nullptr;
	}

	/*				Converters				*/
	LocalType* cast_to_local(double* arr, size_t size) {

		if (typeid(double) == typeid(LocalType)) {
			return arr;
		}

		LocalType* new_arr = (LocalType*)malloc(sizeof(LocalType) * size);
		for (size_t i = 0; i < size; i++)
		{
			new_arr[i] = (LocalType)arr[i];
		}
		return new_arr;
	}

	template <typename out_T>
	out_T* cast_to(LocalType* arr, size_t size) {
		if (typeid(out_T) == typeid(LocalType)) {
			return arr;
		}

		out_T* new_arr = (out_T*)malloc(sizeof(out_T) * size);
		for (size_t i = 0; i < size; i++)
		{
			new_arr[i] = (out_T)arr[i];
		}
		return new_arr;
	}


#if USE_FLOAT
	template float* cast_to<float>(LocalType* arr, size_t size);
	template void save_vector<float>(viennacl::vector<float> x, const char* file_name);
	//template float save_vector_pointer<float>(float* x, size_t size, const char* file_name);
#else
	template double* cast_to<double>(LocalType* arr, size_t size);
	template void save_vector<double>(viennacl::vector<double> x, const char* file_name);
	//template double save_vector_pointer<double>(double* x, size_t size, const char* file_name);
#endif

	//void fill_triplet(vector<T>& tripletList, int numNonzero, int* rowIndices, int* colIndices, double* values) {
	//	auto start = high_resolution_clock::now();
	//	tripletList.reserve(numNonzero);
	//	for (int i = 0; i < numNonzero; i++)
	//	{
	//		tripletList.push_back(T(rowIndices[i] - 1, colIndices[i] - 1, values[i]));
	//	}
	//	auto stop = high_resolution_clock::now();
	//	auto duration = duration_cast<microseconds>(stop - start);
	//	cout << "Time taken by Reader: "
	//		<< duration.count() << " microseconds" << endl;
	//}

	//bool is_symmetric(int numNonzero, int* rowIndices, int* colIndices, double* values) {
	//	bool symmetric = true;
	//	for (int i = 0; i < numNonzero; i++)
	//	{
	//		for (int j = 0; j < numNonzero; j++) {
	//			if (rowIndices[i] == colIndices[j] && rowIndices[j] == colIndices[i]) {
	//				if (values[i] != values[j]) {
	//					return false;
	//				}
	//			}
	//		}
	//	}
	//	return true;
	//}

}