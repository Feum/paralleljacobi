#pragma once

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <string>
#include <ff/utils.hpp>
#include "Matrix.h"

using namespace std;
using namespace ff;

class Jacobi_Seq
{
	
	__declspec(align(8)) vector<double> B;
	__declspec(align(8)) vector<double> oldX;
	__declspec(align(8)) vector<double> newX; 
	size_t N;
	size_t nIt;
	bool showConvergence;
	
	__inline void fillRandVector(vector<double> & B, const size_t N) {
		srand(345);
		for (size_t i = 0; i < N; i++)
			B[i] = rand() % 100 + 1;
}

	
inline void jacobiIteration(const Matrix & 			__restrict__ A, 
							const vector<double> & 	__restrict__ B,
							const vector<double> & 	__restrict__ X,
							vector <double> & 		__restrict__ newX,
														  size_t N) 
{
	double sumAx;

  for (size_t i = 0; i < N; i++) {
	sumAx = 0;
	

	for (size_t j = 0; j < N; j++) 
	{
		sumAx += ((i == j) ? 0 : (A[i][j] * X[j]));
	}
	
	newX[i] = (B[i] - sumAx) / A[i][i];
  }

}
	
	public:
		
		//Prepare the equation system
		Jacobi_Seq(size_t N, size_t nIt, bool showConvergence = false):N(N), nIt(nIt), showConvergence(showConvergence)
		
		{
			B = vector<double>(N);
			oldX = vector<double>(N);
			newX = vector<double>(N);
			fillRandVector(B, N);
			fillRandVector(oldX, N); //Starting vector (random) 
			
		}
		
		//Run Jacobi in sequential with nWor == 0
		//Precondition: N > 0
		__inline void startJacobi()
		{
				Matrix A (N, true);

				if (showConvergence)
				cout << (A.isDiagonallyDominant()? "The Coefficient Matrix is Diagonally Dominant. Granted Convergence." :
									   "The Coefficient Matrix is NOT Diagonally Dominant. Method MAY not Converge") << endl;
			
				ffTime(START_TIME);
				#pragma loop_count min(2), max (100), avg (10)
				for (size_t t = 0; t < nIt; t++)
				{	
					
						jacobiIteration(A, B, oldX, newX, N);
						vector<double>().swap(oldX); //free old x
						oldX = ref (newX);
				
						vector<double> temp (N);
						newX = ref (temp); 
				}
			
				ffTime(STOP_TIME);	
				
			
		}
		
		//take the time and return the results of the equation system
		__inline double measureTime(bool showResults) 
		
		{		
		
				double time = ffTime(GET_TIME);
				if (showResults) 
				{
				cout << "Aproximated solution with " << nIt << " iterations:" << endl;
				cout << "X = [ ";
				for (size_t i = 0; i<N; i++)
					cout << oldX[i] << ", ";
				cout << "] " << endl; 
				
				cout << "Time: " << time << " (ms)" << endl;
				
				}


				
				return time;
			
			
			
		}
		
		//Free the memory
		~Jacobi_Seq(void)
			{
				vector<double>().swap(B); //free B
				vector<double>().swap(oldX); //free old X 
				vector<double>().swap(newX); //free new X
			}

	
	};
