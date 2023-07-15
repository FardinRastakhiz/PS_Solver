// cxdll.cpp : Defines the exported functions for the DLL.
//
#define _CRT_SECURE_NO_WARNINGS
#include "pch.h"
#include "framework.h"
#include "sparse_equation_solver.h"
#include <iostream>
#include <omp.h>
#include "./Eigen/Eigen"
#include <fstream>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <chrono>
#include <future>
#include <ctype.h>

#include "IAlgorithm.h"
#include "IPreconditioner.h"


#include "CgAlgorithm.h"



CXDLL_API double* ses_symmetric_cpu_cg(int numRows, int numNonzero, int* rowIndices,
	int* colIndices, double* values, double* rhs, double* result) {
	ses::IAlgorithm algorithm = ses::CgAlgorithm();

	return nan;
}




using namespace std;
using namespace std::chrono;
typedef Eigen::Triplet<double> T;






const string X_FILE_NAME = "X";
const string SOLVEDB_FILE_NAME = "SOLVEDB";
const string READING_X_FILENAME = "Slatec_X";
const string SLATECB_FILENAME = "SLATECB";

void fill_triplet(vector<T>& tripletList, int numNonzero, int* rowIndices, int* colIndices, double* values) {
	auto start = high_resolution_clock::now();
	tripletList.reserve(numNonzero);
	for (int i = 0; i < numNonzero; i++)
	{
		tripletList.push_back(T(rowIndices[i] - 1, colIndices[i] - 1, values[i]));
	}
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);
	cout << "Time taken by Reader: "
		<< duration.count() << " microseconds" << endl;
}
bool is_symmetric(int numNonzero, int* rowIndices, int* colIndices, double* values) {
	bool symmetric = true;
	for (int i = 0; i < numNonzero; i++)
	{
		for (int j = 0; j < numNonzero; j++) {
			if (rowIndices[i] == colIndices[j] && rowIndices[j] == colIndices[i]) {
				if (values[i] != values[j]) {
					return false;
				}
			}
		}
	}
	return true;
}

bool is_positive_difinite(Eigen::SparseMatrix<double>& A) {
	try {
		// Perform Cholesky decomposition
		Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > cholesky(A);
		if (cholesky.info() == Eigen::Success) {
			return true;
		}
		else {
			return false;
		}
	}
	catch (int _) {
		return false;
	}
}

void fill_vector(Eigen::VectorXd& B, int numRows, double* rhs) { //here rhs means right hand side
	for (int i = 0; i < numRows; i++)
	{
		B(i) = rhs[i];
	}
}
Eigen::VectorXd buildGuessedX(int numRows) {
	Eigen::VectorXd result(numRows);
	double pack = (double)1 / (double)numRows;
	for (int i = 0; i < numRows; i++) {
		result(i) = i * pack;
	}
	return result;
}


Eigen::VectorXd solve(Eigen::SparseMatrix<double>& A, Eigen::VectorXd& B, int numRows) {
	auto start = high_resolution_clock::now();
	Eigen::VectorXd x(numRows);
	Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> solver;
	// uncomment the code bellow to use all cpu cores
	//Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper > solver;
	solver.compute(A);
	x = solver.solve(B);
	//x = solver.solveWithGuess(B, buildGuessedX(numRows));
	cout << "number of threads:";
	cout << Eigen::nbThreads() << endl;
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);
	cout << "Time taken by solver: "
		<< duration.count() << " microseconds" << endl;
	return x;

}
void write_to_file(Eigen::VectorXd X, string filename) {
	std::ofstream file(filename);
	if (file.is_open())
	{
		for (int i = 0; i < X.size(); i++) {
			file << X[i] << endl;
		}
	}
}
Eigen::VectorXd spmv(int numRows, int numNonzero, int* rowIndices, int* colIndices, double* values, Eigen::VectorXd x) {
	Eigen::VectorXd b(numRows);
	// Initialize b to 0
	for (int i = 0; i < numRows; i++) {
		b[i] = 0.0;
	}
	// Perform sparse matrix-vector multiplication
	for (int i = 0; i < numRows; i++) {
		for (int j = 0; j < numNonzero; j++) {
			if (rowIndices[j] - 1 == i) {
				b[i] += x(colIndices[j] - 1) * values[j];
			}
		}
	}
	return b;
}

Eigen::VectorXd spmv2(int numRows, int numNonzero, int* rowIndices, int* colIndices, double* values, Eigen::VectorXd x) {
	Eigen::VectorXd b(numRows);
	// Initialize b to 0
	for (int i = 0; i < numRows; i++) {
		b[i] = 0.0;
	}
	// Perform sparse matrix-vector multiplication
	for (int j = 0; j < numNonzero; j++) {
		b[rowIndices[j] - 1] += values[j] * x(colIndices[j] - 1);
	}
	return b;
}

void fill_vector(Eigen::VectorXd& X) {
	ifstream infile(READING_X_FILENAME);
	int i = 0;
	for (string line; getline(infile, line); i++)
	{
		double value = stod(line);
		X(i) = value;
	}
}

// This is an example of an exported function.
CXDLL_API void solve_matrix(int numRows, int numNonzero, int* rowIndices, int* colIndices, double* values, double* rhs, double* result) {
	if (is_symmetric) {
		cout << "the matrix is symmetric" << endl;
	}
	else {
		cout << "the matrix is not symmetric" << endl;
	}
	Eigen::initParallel();
	int m = 8;
	omp_set_num_threads(m);
	Eigen::setNbThreads(m);
	int n = numRows;
	cout << "Reading A ..." << endl;
	Eigen::VectorXd B(n);
	Eigen::VectorXd solvedBVector(n);
	Eigen::SparseMatrix<double> A(n, n);
	std::vector<T> tripletList;
	fill_triplet(tripletList, numNonzero, rowIndices, colIndices, values);
	A.setFromTriplets(tripletList.begin(), tripletList.end());
	//if (is_positive_difinite(A))
	//{
	//    cout << "Matrix is positive difinite" << endl;
	//}
	//else {
	//    cout << "Matrix is not positive definite" << endl;
	//}
	cout << "Reading A Finished" << endl;
	cout << "Reading B ..." << endl;
	fill_vector(B, numRows, rhs);
	cout << "Reading B Finished" << endl;
	cout << "Solving ..." << endl;
	Eigen::VectorXd X = solve(A, B, numRows);
	Eigen::VectorXd solvedB = spmv2(numRows, numNonzero, rowIndices, colIndices, values, X);
	write_to_file(X, X_FILE_NAME);
	// write solved b to file
	write_to_file(solvedB, SOLVEDB_FILE_NAME);

	//uncomment below lines for multiplying slatec solver
	//int n = numRows;
	//Eigen::VectorXd X(numRows);
	//fill_vector(X);
	//Eigen::VectorXd solvedB = spmv(numRows, numNonzero, rowIndices, colIndices, values, X);
	//// write solved b to file
	//write_to_file(solvedB , SLATECB_FILENAME);

}
