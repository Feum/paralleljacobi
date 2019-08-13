#pragma once

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <string>
#include <ff/pipeline.hpp>
#include <ff/map.hpp>
#include <ff/parallel_for.hpp>
#include "Matrix.h"

using namespace std;
using namespace ff;

class Jacobi_Par
{
	__declspec(align(8)) vector<double> B;
	__declspec(align(8)) vector<double> oldX;
	__declspec(align(8)) vector<double> newX; 
	__declspec(align(16)) vector<int32_t> M;
	size_t N;
	size_t nWor;
	size_t nIt;
	bool showConvergence;
	
	__inline void fillRandVector(vector<double> & B, const size_t N) {
		srand(345);
		
		# pragma simd
		for (size_t i = 0; i < N; i++)
			B[i] = rand() % 100 + 1;
	}

	__inline void preprocessModules(vector<int32_t> & __restrict__ M,
									ParallelForReduce<double> & pfr,
									int32_t N)
	{
			size_t incr = ((N & (N - 1)) == 0) ? N : pow(2,floor(log(N)/log(2)));
			M = vector<int32_t> (N);
			
			pfr.parallel_for(0, N, 1, 1, [&](const int32_t i){ //counting the spurious elements
				if (((i*incr) % N)  >= incr) 
					M[i] = 1;
					
			}); 
	

			for (size_t j = 1; j<M.size(); j++) //prefix sums
			{
				M[j]=M[j]+M[j-1];
			}
				
		
	}
	
	__inline void jacobiIteration(  const Matrix & __restrict__ A, 
									const vector<double> &  __restrict__ B,
									const vector<double> &  __restrict__ X,
										  vector <double> & __restrict__ newX,
														ParallelForReduce<double> & pfr,
														size_t N,
														size_t nWor,
									const vector<int32_t> & __restrict__ M) 
{					
		size_t k = 0, sqN = N*N;
			__declspec(align(8)) vector<double> T (sqN);
								
			pfr.parallel_for(0, sqN, N, [&](const size_t j){
					size_t i = j/N, k = j%N; 			
					T[j] = i == k ? 0 : A[i][k] * X[k]; 
								
			});
				
							
			size_t toBeProcessed = sqN; //Number of AijXj yet to be summed in TOTAL
			size_t step = 2; //Step of the inner par_for
			int32_t incr = 1; //where to find next/previous AijXj
			size_t perRawProcessable; //Number of AijXj yet to be summed PER RAW
			
			
			while (toBeProcessed > N) //While there's something to sum...
			{
				perRawProcessable = toBeProcessed/N;
						
				srand(123); 	
				pfr.parallel_for(0, sqN, step, N-(N/10), [&](const size_t i){ //Sum in place couples in parallel

					size_t position = i % N;
					size_t  next = i+incr, nextNext = i+2*incr, prec=i-incr;
					double succ;

					if ((perRawProcessable % 2) == 1) //If i have to sum an ODD number of things per RAW...
						
						{	if (position/incr == 1) //If i'm here i've skipped an element. Sum also the previous element.
									succ = T[next] + T[prec]; 
									
					 		else if ((perRawProcessable - position/incr) == 3) //Sum 3 things instead of 2 if only 3 things are left to sum in the RAW
									succ = T[next] + T[nextNext]; 
									
							else if ((perRawProcessable - position/incr) <= 1) //Do nothing if i end up in the last element. Already considered before
								succ = 0; 									   //Will be passed over as an actual result, but i'll skip it later
							else	
								succ = T[next];	//In any other case, sum two non-zero adiacent things.
						}
					else //If i have to sum an EVEN number of things per RAW...
						{
						
							if ((perRawProcessable - position/incr)%2 == 1) //If there are an odd number of element before me, i've skipped one
								succ = T[prec]; 
								
							else if ((perRawProcessable - position/incr) <= 1) //Do nothing if i end up in the last element. Already considered before
								succ = 0;									   //Will be passed over as an actual result, but i'll skip it later
							else	
								succ = T[next];	//In any other case, sum two non-zero adiacent things.
						}
				
					
					T[i] = T[i] + succ;
				});	
				
				toBeProcessed = toBeProcessed - N*((perRawProcessable)/2 + (perRawProcessable % 2)); //Update toBeProcessed counter
				step *=2; //double the step to skip already summed values.
				incr *=2 ; //double the incr according to the step
														
			}
							

			pfr.parallel_for(0, sqN, incr, [&](const int32_t i){				
						
								
						
						if (((int32_t)(i % N) - incr) < 0) { //Skip spurious results where i had to "do nothing"
							
							//formula to fit in the [0, N) vectors. Uses the preprocessed M vector to count the n of spurious before every position
							
						 	
							int32_t closedCycles = (i/(N*incr))*(N - incr); //every N increments there are a fixed number of spurious elements
							int32_t midCycle = M[(i/incr) % N]; //the number of spurious in the middle of a cycle of N increments
							
							int32_t nOfSpuriousBeforeMe = closedCycles + midCycle;  //the number of elements before me i mustn't consider
							    
							    
						
							int32_t p = (i/incr) - nOfSpuriousBeforeMe; //my position in the vector excluding the elements that will not go in the vector
						
							//Building new X vector
							newX[p] = (B[p]-T[i]) / A[p][p];
							
					}
							
							
			});					
					     
}
	
	public:
		
		//Prepare the equation system
		Jacobi_Par(size_t N, size_t nIt, size_t nWor, bool showConvergence=true):N(N), nWor(nWor), nIt(nIt), showConvergence(showConvergence)
		
		{
			B = vector<double>(N);
			oldX = vector<double>(N);
			newX = vector<double>(N);
			fillRandVector(B, N);
			ParallelForReduce<double> pfr(nWor);
			preprocessModules(M,pfr,N);
			fillRandVector(oldX, N); //Starting vector (random) 
			

			
		}
		
		//Run Jacobi in parallel with nWor > 0; 
		//Precondition: nWor > 0, N > 0.
		__inline void startJacobi()
		{
				ParallelForReduce<double> pfr(nWor,true);
				Matrix A (N, true);
			
				if (showConvergence) cout << (A.isDiagonallyDominant()? "The Coefficient Matrix is Diagonally Dominant. Granted Convergence." :
									   "The Coefficient Matrix is NOT Diagonally Dominant. Method MAY not Converge") << endl;
			

				ffTime(START_TIME);	
				#pragma loop_count min(2), max (100), avg (10)
				for (size_t t = 0; t < nIt; t++)
				{	
							
						jacobiIteration(A, B, oldX, newX, pfr, N, nWor, M);
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
		~Jacobi_Par(void)
			{
				vector<double>().swap(B); //free B
				vector<double>().swap(oldX); //free old X 
				vector<double>().swap(newX); //free new X
				vector<int32_t>().swap(M); //free M
			}

	
	};
