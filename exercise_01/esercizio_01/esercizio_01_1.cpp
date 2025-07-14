#include <iostream>
#include <fstream>
#include <cmath>

#include "../ParallelRandomNumberGenerator/random.h"

using namespace std;

int main() {

	int M = 10000;			//total number of points
	int N = 100;			//Number of blocks
	int L = M / N;			//point increment per block

	Random rand;
	rand.Initialize();

	//----------------file output:
	string filename1 = "mean.dat";
	ofstream out(filename1);
	
	if (!out.is_open()) {
		cerr << "Errore nell'apertura del file " << filename1 << endl;
		return -1;
	}

	//---------------------------------compute averages:
	
	double* r = new double[M];
	double* mysum = new double[N];						//accumulators for the total block average of the first integral
	double* sigmasquared = new double[N];				//accumulators for the total block average of the second integral
	double* mysum_squares = new double[N];				//accumulators for the sum of the squared averages of the first integral
														//(to get another average!)
	double* mysum_squares_sigma = new double[N];		//accumulators for the sum of the squared averages of the second integral
	double* myerr_mean = new double[N];					//error after each N block for the first integral
	double* myerr_sigma = new double[N];				//error after each N block for the second integral

	//first for i=0:
	mysum[0] = 0;
	sigmasquared[0] = 0;
	for (int j = 0; j < L; j++) {
		r[j] = rand.Rannyu();						// generate a random number array
		mysum[0] += r[j];							// accumulate the sum for mean
		sigmasquared[0] += pow(r[j] - 0.5, 2);		// accumulate the sum for variance
	}
	mysum[0] /= L;
	sigmasquared[0] /= L;
	mysum_squares[0] = mysum[0] * mysum[0];
	mysum_squares_sigma[0] = sigmasquared[0] * sigmasquared[0];
	//output results:
	out << mysum[0] << " " << 0 << " " << sigmasquared[0] << " " << 0 << "\n";	// output file mean and error 0 for first block

	//then the same for all the other blocks:
	for (int i = 1; i < N; i++) {			//iteration on the i-th block
		mysum[i] = 0;
		sigmasquared[i] = 0;

		for (int j = 0; j < L; j++) {
			r[i * L + j] = rand.Rannyu();
			mysum[i] += r[i * L + j] / L;
			sigmasquared[i] += pow(r[i * L + j] - 0.5, 2) / L;
		}

		//average and accumulate (by de-averaging the sum of the previous step, adding it to the current one, and computing
		//average on the current nuber of blocks (step)
		mysum_squares[i] = (mysum[i] * mysum[i] + mysum_squares[i - 1] * i) / (i + 1);
		mysum_squares_sigma[i] = (sigmasquared[i] * sigmasquared[i] + mysum_squares_sigma[i - 1] * i) / (i + 1);
		mysum[i] = (mysum[i] + mysum[i - 1] * i)/(i + 1);
		sigmasquared[i] = (sigmasquared[i] + sigmasquared[i - 1] * i) / (i + 1);

		//compute statistical uncertainty:
		myerr_mean[i] = sqrt((mysum_squares[i] - mysum[i] * mysum[i]) / i);
		myerr_sigma[i] = sqrt((mysum_squares_sigma[i] - sigmasquared[i] * sigmasquared[i]) / i);

		//output results:
		out << mysum[i] << " " << myerr_mean[i] << " " << sigmasquared[i] << " " << myerr_sigma[i] << "\n"; // output file mean and its error
	}

	out.close();

	//-----------------------------Chi2 Test:
	
	//----------------file output:
	string filename2 = "chi2.dat";
	out.open(filename2);

	if (!out.is_open()) {
		cerr << "Errore nell'apertura del file " << filename2 << endl;
		return -1;
	}

	double n[N];			//counts the observed values in the N-th interval
	double chi2;			//chi2 of the n-th iteration (n<100)

	for (int tries = 0; tries < 100; tries ++){
		chi2 = 0;

		for (int i = 0; i < N; i++) {
			n[i] = 0;		

			for (int j = 0; j < M; j++) {
				r[j] = rand.Rannyu();
				if (r[j] >= double(i)/N and r[j] < double(i + 1)/N){	//if r falls into the bin
					n[i]++;
				}
			}
			chi2 += pow(n[i] - double(M / N), 2) / double(M / N);		//compute chi2 with expected value for uniform distribution
		}

		out << chi2 << "\n";
	}

	out.close();

	delete[] r;
	delete[] mysum;
	delete[] sigmasquared;
	delete[] mysum_squares;
	delete[] mysum_squares_sigma;
	delete[] myerr_mean;
	delete[] myerr_sigma;

	return 0;
}