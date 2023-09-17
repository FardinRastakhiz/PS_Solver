#pragma once



#include "viennacl/scalar.hpp"
#include "viennacl/vector.hpp"

#include <stdio.h>
//#include <iostream>
//#include <fstream>
template<typename T>
void save_vector(viennacl::vector<T> x, const char* file_name) {
    int size = x.size();
    FILE* fp;
    fp = fopen(file_name, "w");
    for (unsigned i = 0; i < size; i++) {
        float output = x[i];
        fprintf(fp, "%d, %f\n", i, output);
    }
    fclose(fp);
}