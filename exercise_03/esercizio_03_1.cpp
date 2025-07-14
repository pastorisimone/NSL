#include <iostream>
#include <fstream>
#include <cmath>

#include "../ParallelRandomNumberGenerator/random.h"

using namespace std;

//statistical error
double error(double avg, double avg2, int n) {
	if (n == 0) return 0;
	return sqrt((avg2 - avg * avg) / n);
}

//Brownian motion:
double BM_step(double mu, double sigma, double step, double z) {

	return mu * step + sigma * z * sqrt(step);
}

double* build_BM(Random& rand, int N, double* t, double mu = 0, double sigma = 1, double IC = 0) {

	double* result = new double[N];
	result[0] = IC;
	double z;
	double step;

	for (int i = 1; i < N; i++) {
		z = rand.Gauss(0, 1);				//													<|--------------------<<<<<<
		step = t[i] - t[i - 1];
		//result[i] = result[i - 1] + mu * (t[i] - t[i - 1]) + sigma * z * sqrt(t[i] - t[i - 1]);
		result[i] = result[i - 1] + BM_step(mu, sigma, step, z);
	}

	return result;
}

//Geometric Brownian motion:
double GBM_step(double mu, double sigma, double step, double z) {

	return exp((mu - 0.5 * sigma * sigma) * step + sigma * z * sqrt(step));
}

//Sample GMB at t with initial condition IC:
double direct_GBM(Random& rand, double t, double mu, double sigma, double IC) {

	double BM = rand.Gauss(0, t);
	
	return IC * exp((mu - 0.5 * sigma * sigma) * t + sigma * BM);
}

double* build_GBM(Random& rand, int N, double* t, double mu, double sigma, double IC) {

	double* result = new double[N];
	result[0] = IC;
	double z;
	double step;

	for (int i = 1; i < N; i++) {
		z = rand.Gauss(0, 1);				
		step = t[i] - t[i - 1];
		result[i] = result[i - 1] * GBM_step(mu, sigma, step, z);
	}

	return result;
}

int main() {

	Random rand;
	rand.Initialize();

	int n_steps = 100;
	double T = 1;
	double* t = new double[n_steps];
	for (int i = 0; i < n_steps; i++) t[i] = i*T/ n_steps;
	double r = 0.1;
	double sigma = 0.25;
	double S0 = 100;

	//------------------------test the generators:
	string filename1 = "BM_simulation.dat";
	ofstream out(filename1);

	if (!out.is_open()) {
		cerr << "Errore nell'apertura del file " << filename1 << endl;
		return -1;
	}

	double* my_BM = build_BM(rand, n_steps, t);
	double* my_GBM = build_GBM(rand, n_steps, t, r, sigma, S0);

	for (int i = 0; i < n_steps; i++) {
		out << t[i] << " " << my_BM[i] << " " << my_GBM[i] << endl;
	}

	out.close();

	//------------------------compute call and put option price with Monte Carlo method:
	int M = 10000;			//total number of points
	int N = 100;			//Number of blocks
	int L = M / N;			//point increment per block

	//params:
	double K = 100;

	double avg[4] = { 0 };  //[CALL_direct, PUT_direct, CALL_discrete, PUT_discrete]
	double av2[4] = { 0 };
	double temp = -777.;

	//file output:
	string filename2 = "EuropeanOptions.dat";
	out.open(filename2);

	if (!out.is_open()) {
		cerr << "Errore nell'apertura del file " << filename2 << endl;
		return -1;
	}


	for (int i = 0; i < N; i++) {	

		double accu[4] = { 0 };

		for (int j = 0; j < L; j++) {

			//direct GBM:
			temp = direct_GBM(rand, T, r, sigma, S0);
			accu[0] += exp(-r * T) * max(temp - K, 0.) / L;
			accu[1] += exp(-r * T) * max(K - temp, 0.) / L;

			//discretized GBM:
			my_GBM = build_GBM(rand, n_steps, t, r, sigma, S0);
			temp = my_GBM[n_steps - 1];
			accu[2] += exp(-r * T) * max(temp - K, 0.) / L;
			accu[3] += exp(-r * T) * max(K - temp, 0.) / L;
		}

		av2[0] = (av2[0] * i + accu[0] * accu[0]) / (i + 1);
		av2[1] = (av2[1] * i + accu[1] * accu[1]) / (i + 1);
		av2[2] = (av2[2] * i + accu[2] * accu[2]) / (i + 1);
		av2[3] = (av2[3] * i + accu[3] * accu[3]) / (i + 1);

		avg[0] = (avg[0] * i + accu[0]) / (i + 1);
		avg[1] = (avg[1] * i + accu[1]) / (i + 1);
		avg[2] = (avg[2] * i + accu[2]) / (i + 1);
		avg[3] = (avg[3] * i + accu[3]) / (i + 1);

		out << avg[0] << " " << error(avg[0], av2[0], i) << " ";
		out << avg[1] << " " << error(avg[1], av2[1], i) << " ";
		out << avg[2] << " " << error(avg[2], av2[2], i) << " ";
		out << avg[3] << " " << error(avg[3], av2[3], i) << endl;
	}



	out.close();

	delete[] t;
	delete[] my_BM;
	delete[] my_GBM;

	return 0;
}