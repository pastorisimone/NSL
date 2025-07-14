#include <iostream>
#include <fstream>
#include <cmath>

#include "../ParallelRandomNumberGenerator/random.h"

using namespace std;

double my_cos(double x) {
	
	return M_PI * 0.5 * cos(M_PI * x * 0.5);

}

void MC_integrate(double (*my_funct)(double), Random& rand, int M=10000, int N=100, string filename = "mean.dat") {

	int L = M / N;			//point increment per block

	//----------------file output:
	ofstream out(filename);

	if (!out.is_open()) {
		cerr << "Errore nell'apertura del file " << filename << endl;
	}
	else{

		//----------------compute averages:

		double r;
		double* mysum = new double[N];
		double* mysum_squares = new double[N];
		double* myerr_mean = new double[N];

		//first for i=0:
		mysum[0] = 0;
		for (int j = 0; j < L; j++) {
			r = rand.Rannyu();							// generate a random number array
			mysum[0] += my_funct(r);					// accumulate the sum for mean
		}
		mysum[0] /= L;
		mysum_squares[0] = mysum[0] * mysum[0];
		//output results:
		out << mysum[0] << " " << 0 << "\n";	// output file mean and error 0 for first block

		//then the same for all the other blocks:
		for (int i = 1; i < N; i++) {			//iteration on the i-th block
			mysum[i] = 0;

			for (int j = 0; j < L; j++) {
				r = rand.Rannyu();
				mysum[i] += my_funct(r) / L;
			}

			mysum_squares[i] = (mysum[i] * mysum[i] + mysum_squares[i - 1] * i) / (i + 1);//cfr
			mysum[i] = (mysum[i] + mysum[i - 1] * i) / (i + 1);

			//compute statistical uncertainty:
			myerr_mean[i] = sqrt((mysum_squares[i] - mysum[i] * mysum[i]) / i);

			//output results:
			out << mysum[i] << " " << myerr_mean[i] << "\n"; // output file mean and its error

		}

			delete[] mysum;
			delete[] mysum_squares;
			delete[] myerr_mean;
	}

		out.close();

}

double TAY(double x) {

	return 1.5 - 1.5 * x * x;//2/3
}

double integranda(double x) {

	return my_cos(x) / TAY(x);
}

double funct_sample(double (*my_funct)(double), Random& rand, double funct_max) {

	while (true) {
		double x = rand.Rannyu();		//point on the described interval [0,1)
		double u = rand.Rannyu();

		if (u < my_funct(x) / funct_max) {
			return x;
		}
	}
}


void MC_integrate_2(double (*my_funct)(double), Random& rand, int M = 10000, int N = 100, string filename = "mean.dat") {

	int L = M / N;			//point increment per block

	//----------------file output:
	ofstream out(filename);

	if (!out.is_open()) {
		cerr << "Errore nell'apertura del file " << filename << endl;
	}
	else {

		//----------------compute averages:

		double r;
		double* mysum = new double[N];
		double* mysum_squares = new double[N];
		double* myerr_mean = new double[N];
		double this_max = M_PI * 0.5;

		//first for i=0:
		mysum[0] = 0;
		for (int j = 0; j < L; j++) {
			r = funct_sample(TAY, rand, this_max);		// generate a random number array
			mysum[0] += my_funct(r);					// accumulate the sum for mean
		}
		mysum[0] /= L;
		mysum_squares[0] = mysum[0] * mysum[0];
		//output results:
		out << mysum[0] << " " << 0 << "\n";	// output file mean and error 0 for first block

		//then the same for all the other blocks:
		for (int i = 1; i < N; i++) {			//iteration on the i-th block
			mysum[i] = 0;

			for (int j = 0; j < L; j++) {
				r = funct_sample(TAY, rand, this_max);
				mysum[i] += my_funct(r) / L;
			}

			mysum_squares[i] = (mysum[i] * mysum[i] + mysum_squares[i - 1] * i) / (i + 1);//cfr
			mysum[i] = (mysum[i] + mysum[i - 1] * i) / (i + 1);

			//compute statistical uncertainty:
			myerr_mean[i] = sqrt((mysum_squares[i] - mysum[i] * mysum[i]) / i);

			//output results:
			out << mysum[i] << " " << myerr_mean[i] << "\n"; // output file mean and its error

		}

		delete[] mysum;
		delete[] mysum_squares;
		delete[] myerr_mean;
	}

	out.close();

}

int main() {

	int M = 10000;			//total number of points
	int N = 100;			//Number of blocks

	Random rand;
	rand.Initialize();

	MC_integrate(my_cos, rand, M, N);
	MC_integrate_2(integranda, rand, M, N, "importance.dat");


return 0;
}