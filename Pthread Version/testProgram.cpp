#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <string>
#include "Jacobi_Par.h"
#include "Jacobi_Seq.h"

using namespace std;

int main(int argc, char* argv[]) {

	bool benchmarkMode = false;
	//./Jacobi [coeff square matrix dim] [numb of iterat] [nworkers]
	if ((argc == 2) && ((string)argv[1])=="-benchmark")
		benchmarkMode = true;
	else if (argc != 4) {cout << "Use: " <<argv[0] << "[coeff square matrix dim] [num of iterat] [nworkers]"; return -1;}
		
		
	//Ad Hoc Mode
	if (!benchmarkMode)
	
	{
		size_t N = stoi(argv[1]);
		size_t nIt = stoi(argv[2]);
		size_t nWor = stoi(argv[3]);
		double time;
	
		//PARALLEL
		if (nWor > 0)	
		
		{
		
		Jacobi_Par jacobi = Jacobi_Par(N,nIt,nWor,false);
		
		jacobi.startJacobi();
		time = jacobi.measureTime(false);
		cout << "N=" << N <<";nWorkers=" <<nWor<<";time(ms): " << time <<endl;	
		}
	
		//SEQUENTIAL
		else
		
		{
		Jacobi_Seq jacobi = Jacobi_Seq(N,nIt);
		
		jacobi.startJacobi();
		time = jacobi.measureTime(true);
		cout << "N=" << N <<";Seq;time(ms): " << time <<endl;		
			
		}
	}
	//Benchmark Mode
	
	else
		{
			//Theory max par degree n^2 - n

			int i = 7;
			int N = 400;
			double time;
			
			Jacobi_Par jacobi = Jacobi_Par(N, 10, i,false);
			jacobi.startJacobi();
			time = jacobi.measureTime(false);


			
		}
		
}
