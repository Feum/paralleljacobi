#pragma once

#include <vector>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <functional>

using namespace std;

class Matrix
{
	vector<vector<double>> _mat;
	size_t N;

public:

	Matrix(int32_t n)
	{
		N = n;
		srand(123);
		vector<vector<double>> m (n, vector<double>(n)); 
		for(size_t i = 0; i<n; i++)
			for (size_t j = 0; j<n; j++)
				m[i][j] = rand() % 100 + 1; 
		_mat = ref (m);
	}
	
	//generate a Diagonally Dominant random matrix
	Matrix(int32_t n, bool dominant)
	{	
		N = n;
		srand(123);
		vector<vector<double>> m (n, vector<double>(n)); 
		for(size_t i = 0; i<n; i++)
			for (size_t j = 0; j<n; j++)
				m[i][j] = (i == j) ? (rand() % 100 + 10*n + 1) :  (rand() % 10 + 1); 
		_mat = ref (m);
	}

	inline bool isDiagonallyDominant() const
	{
		double nonDiagSum;
		double diagElem;
		bool DD;
		for (size_t i = 0; i<N; i++)
		{
			nonDiagSum = 0;
			diagElem = abs(_mat[i][i]);
			for (size_t j = 0; j<N; j++)
			{
				nonDiagSum += abs(_mat[i][j]);
			}
			nonDiagSum -= diagElem;

			DD = (diagElem > nonDiagSum);
		}

		return DD;

	}

	vector<double> operator[] (size_t i) const
	{
		return _mat[i];
	}

	~Matrix(void)
	{
		vector<vector<double>>().swap(_mat); //free mat
	}
};

