#pragma once



#include "viennacl/scalar.hpp"
#include "viennacl/vector.hpp"

#include <stdio.h>
//#include <iostream>
//#include <fstream>


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


// implementing equivalent of "is" is c# for type checking
template < class T, class U >
bool is_inst_of(U u) {
    return dynamic_cast<T>(u) != nullptr;
}