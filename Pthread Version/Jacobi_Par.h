#pragma once

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <string>
#include <pthread.h>
#include <deque>
#include <atomic>
#include <chrono>
#include <ff/utils.hpp>
#include "Matrix.h"

using namespace std;
using namespace ff;

enum ops {mult,sum,divis,nop};


struct task { //operation task: granularity = raw
	ops operation = nop; //type of operation (sum, multiplication or division)
	size_t raw; //raw involved in the operation
	__declspec(align(8)) vector<double> Ai; //A vector initially storing the concerned raw of the matrix, then storing the A[i][j]*X[j] products
	__declspec(align(8)) vector<double> X; //(pointer to) a vector storing the iteration vector
	double partial; //Initialized with the value of B[raw], then contains the partial results.
	
	
};

__declspec(align(64)) vector<deque<task>> taskRepos; //Operations ready to be done
vector<pthread_t> workers;
uint32_t ID;
vector<pthread_mutex_t> mutex_front;
vector<pthread_mutex_t> mutex_back;
int32_t nOfRepos;
__declspec(align(8)) vector<double> newX; 
atomic<int32_t> nOfTaskLeft (0);
atomic<bool> endOfAlgorithm (false);
pthread_mutex_t endOfIt_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t endOfIt_cond = PTHREAD_COND_INITIALIZER;


void *worker(void*){

	srand(pthread_self());
	while (!endOfAlgorithm) 
	{
		
		task r; 
		bool locked = false;
		size_t max = 0;		
		int32_t idx =  rand() % nOfRepos;

		//Thread safe task pop begin
		if (!(taskRepos[idx]).empty()) {

			if (pthread_mutex_trylock(&(mutex_back[idx])) == 0)
			{
				pthread_mutex_lock(&(mutex_front[idx]));
				if ((taskRepos[idx]).size() > 1) pthread_mutex_unlock(&(mutex_front[idx]));
				else locked = true;
					
				if (!(taskRepos[idx]).empty()) {r = (taskRepos[idx]).back(); (taskRepos[idx]).pop_back();}
				pthread_mutex_unlock(&(mutex_back[idx]));
				
				if (locked){ pthread_mutex_unlock(&(mutex_front[idx])); locked = false; }
			}
			else
			{		
				pthread_mutex_lock(&(mutex_front[idx])); 
				if ((taskRepos[idx]).size() > 1) {if (!taskRepos[idx].empty()) r = (taskRepos[idx]).front(); if (!taskRepos[idx].empty()) (taskRepos[idx]).pop_front(); else r.operation = nop;}
				pthread_mutex_unlock(&(mutex_front[idx]));
			}
		//thread safe task pop end

			switch (r.operation) 
			{
				
				case mult:
					# pragma ivdep
					for(size_t j = 0; j<(r.Ai).size(); j++)
						(r.Ai)[j] *= (r.raw == j)? 1 : (r.X)[j];
					r.operation = sum;	
					break;
				case sum:
					for (size_t j = 0; j<(r.Ai).size(); j++)
						r.partial += (r.raw == j)? 0 : -(r.Ai)[j];
					r.operation = divis; 
					break;
				case divis:
					newX[r.raw] = r.partial/((r.Ai)[r.raw]);  
					vector<double>().swap(r.Ai); //free Ai copy
					r.operation = nop;
					nOfTaskLeft--; 
					break;
				default:
					break;
					
			}
			//thread safe task push begin
			if (r.operation != nop)
				if (pthread_mutex_trylock(&(mutex_back[idx])) == 0)
					{
						pthread_mutex_lock(&(mutex_front[idx]));
							if ((taskRepos[idx]).size() > 1) pthread_mutex_unlock(&(mutex_front[idx]));
							else locked = true;
							
						if (r.operation != nop)(taskRepos[idx]).push_back(r);
						pthread_mutex_unlock(&(mutex_back[idx]));
						
						if (locked) {pthread_mutex_unlock(&(mutex_front[idx])); locked = false;}
					}
					else
					{		
						pthread_mutex_lock(&(mutex_front[idx]));
						if (r.operation != nop)(taskRepos[idx]).push_front(r);
						pthread_mutex_unlock(&(mutex_front[idx]));
					}
			 //thread safe task push end

		}
						
			 if (nOfTaskLeft == 0) 
			 {
				 pthread_mutex_lock(&endOfIt_mutex);
				 pthread_cond_signal(&endOfIt_cond);
				 pthread_mutex_unlock(&endOfIt_mutex);
		     }

	}
   			 
	pthread_exit(nullptr);

}
class Jacobi_Par
{
	__declspec(align(8)) vector<double> B;
	__declspec(align(8)) vector<double> oldX;
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

	
	
	public:
		
		//Prepare the equation system
		Jacobi_Par(size_t N, size_t nIt, size_t nWor, bool showConvergence=true):N(N), nWor(nWor), nIt(nIt), showConvergence(showConvergence)
		
		{
			B = vector<double>(N);
			oldX = vector<double>(N);
			newX = vector<double>(N);
			fillRandVector(B, N);
			fillRandVector(oldX, N); //Starting vector (random) 
			workers = vector<pthread_t>(nWor);
			ID = 0;
			nOfRepos = floor((double)sqrt(N));
			taskRepos = vector<deque<task>>(nOfRepos, deque<task>()); 
			mutex_front = vector<pthread_mutex_t>(nOfRepos);
			mutex_back = vector<pthread_mutex_t>(nOfRepos);
			nOfTaskLeft = N;
			for (auto t : mutex_front)
				pthread_mutex_init(&t, NULL);
			
			for (auto t : mutex_back)
				pthread_mutex_init(&t, NULL);			
			
		}
		
		//Run Jacobi in parallel with nWor > 0; 
		//Precondition: nWor > 0, N > 0.
		__inline void startJacobi()
		{
				Matrix A (N, true);
		
			
				if (showConvergence) cout << (A.isDiagonallyDominant()? "The Coefficient Matrix is Diagonally Dominant. Granted Convergence." :
									   "The Coefficient Matrix is NOT Diagonally Dominant. Method MAY not Converge") << endl;


				ffTime(START_TIME);


				for (size_t t=0; t < nWor; t++)
					pthread_create(&(workers[t]),nullptr,&worker,nullptr);
			
				for (size_t i = 0; i < nIt; i++)
				{
								
					for(size_t t = 0; t < N; t++)
						{
							task r;
					
							r.operation = mult;
							r.raw = t;
						
							r.partial = B[t];
							r.Ai = A[t];
					
							r.X = ref (oldX);
							size_t idx = t % nOfRepos;
										
							pthread_mutex_lock(&(mutex_front[idx]));
							(taskRepos[idx]).push_front(r);
							pthread_mutex_unlock(&(mutex_front[idx]));
							
						}
						
					pthread_mutex_lock(&endOfIt_mutex);
					while (nOfTaskLeft > 0)
						pthread_cond_wait(&endOfIt_cond, &endOfIt_mutex); 
					pthread_mutex_unlock(&endOfIt_mutex);		
					
					vector<double>().swap(oldX); //free oldX
					oldX = ref (newX);
				
					newX = vector<double>(N);
					
					nOfTaskLeft = N;				
		
				}
				
			ffTime(STOP_TIME);

	

	
				
			endOfAlgorithm = true;
				
			for (size_t t = 0; t < nWor; t++)
				pthread_join(workers[t], nullptr);	
				
							
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
				vector<double>().swap(newX); //free new X
			}

	
	};
