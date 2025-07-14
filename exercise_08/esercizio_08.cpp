#include <iostream>
#include <fstream>
#include <cmath>
#include <functional>
#include <iomanip>

#include "../ParallelRandomNumberGenerator/random.h"

using namespace std;


struct Vect {		//3d vector

	double x, y, z;

	Vect() : x(0), y(0), z(0) {}
	Vect(double X, double Y, double Z) : x(X), y(Y), z(Z) {}

	Vect operator+(Vect& v2) const {
		return Vect(x + v2.x, y + v2.y, z + v2.z);
	}

	Vect operator-(Vect& v2) const {
		return Vect(x - v2.x, y - v2.y, z - v2.z);
	}

	Vect operator*(double scalar) const {
		return Vect(x * scalar, y * scalar, z * scalar);
	}

	double norm2() const {
		return x * x + y * y + z * z;
	}

};

Vect Step(Random& rand, Vect X0, double delta_x = 2, double delta_y = 0.1, double delta_z = 0.1) {

	Vect x_new;
	x_new.x = (rand.Rannyu() * delta_x - delta_x / 2.);
	x_new.y = X0.y + (rand.Rannyu() * delta_y - delta_y / 2.);
	x_new.z = X0.z + (rand.Rannyu() * delta_z - delta_z / 2.);
	return x_new;
}
/*Vect* Metropolis_3D(Random& rand, int n_step, double (*p)(double, double, double), Vect X0, double equilibrate = 0.5) {

	Vect* X = new Vect[n_step];
	X[0] = X0;
	Vect x_new;
	double a;
	int n_step_eq = int(double(n_step) * equilibrate);

	double delta = 0.0529;										//------------------------------scala ----- acceptance????

	//equilibrate:
	for (int n = 0; n < n_step_eq; n++) {

		Vect step(rand.Rannyu() * 2 * delta - delta, rand.Rannyu() * 2 * delta - delta, rand.Rannyu() * 2 * delta - delta);
		x_new = X[n] + step;

		a = min(1., p(x_new.x, x_new.y, x_new.z) / p(X[n].x, X[n].y, X[n].z));

		double r = rand.Rannyu();
		if (r <= a) X[n + 1] = x_new;					//accept
		else X[n + 1] = X[n];							//reject
	}
	X[0] = X[n_step_eq];

	//actually start generating points after reaching equilibrium
	for (int n = 0; n < n_step - 1; n++) {

		Vect step(rand.Rannyu() * 2 * delta - delta, rand.Rannyu() * 2 * delta - delta, rand.Rannyu() * 2 * delta - delta);
		x_new = X[n] + step;

		a = min(1., p(x_new.x, x_new.y, x_new.z) / p(X[n].x, X[n].y, X[n].z));

		double r = rand.Rannyu();
		if (r <= a) X[n + 1] = x_new;					//accept
		else X[n + 1] = X[n];							//reject
	}

	return X;

}*/


double* Metropolis_1D(Random& rand, int n_step, double (*p)(double, double, double), Vect x0, double equilibrate = 0.5) {


	double* x = new double[n_step];
	x[0] = x0.x;
	double sigma = x0.y;
	double mu = x0.z;
	double x_new;
	double a;
	int n_step_eq = int(double(n_step) * equilibrate);

	double delta = 0.0529;										//------------------------------scala ----- acceptance????

	//equilibrate:
	for (int n = 0; n < n_step_eq; n++) {

		x_new = x[n] + (rand.Rannyu() * 2 * delta - delta);

		a = min(1., p(x_new, sigma, mu) / p(x[n], sigma, mu));

		double r = rand.Rannyu();
		if (r <= a) x[n + 1] = x_new;					//accept
		else x[n + 1] = x[n];							//reject
	}
	x[0] = x[n_step_eq];

	//actually start generating points after reaching equilibrium
	for (int n = 0; n < n_step - 1; n++) {

		x_new = x[n] + (rand.Rannyu() * 2 * delta - delta);

		a = min(1., p(x_new, sigma, mu) / p(x[n], sigma, mu));

		double r = rand.Rannyu();
		if (r <= a) x[n + 1] = x_new;					//accept
		else x[n + 1] = x[n];							//reject
	}


	return x;

}


double MC_integrate(function<double(double, double, double)> my_funct, Random& rand, double* funct_sample, double sigma, double mu, 
					double& error_out, int M = 10000, int N = 100, string filename = "no_file") {

	int L = M / N;			//point increment per block				<|-------------<<<<< cambiare in base all'array?

	double r;
	double* mysum = new double[N];
	double* mysum_squares = new double[N];
	double* myerr_mean = new double[N];
	double result;													//
	myerr_mean[0] = 0;

	ofstream result_out;
	if (filename != "no_file") {
		result_out.open(filename, ios::out);
		if (!result_out.is_open()) {
			cerr << "Errore nell'apertura del file " << "output_dump_SA.dat" << endl;
		}
	}


	//first for i=0:
	mysum[0] = 0;
	for (int j = 0; j < L; j++) {
		r = funct_sample[j];						// generated sample array		// <|------<<<<<< oppure estrarre a caso da f_sample?
		mysum[0] += my_funct(r, sigma, mu);					// accumulate the sum for mean
	}
	mysum[0] /= L;
	mysum_squares[0] = mysum[0] * mysum[0];
	if (filename != "no_file") result_out << mysum[0] << setw(12) << 0 << endl;

	//then the same for all the other blocks:
	for (int i = 1; i < N; i++) {			//iteration on the i-th block
		mysum[i] = 0;

		for (int j = 0; j < L; j++) {
			r = funct_sample[j + L * i];
			mysum[i] += my_funct(r, sigma, mu) / L;
		}
		mysum_squares[i] = (mysum[i] * mysum[i] + mysum_squares[i - 1] * i) / (i + 1);
		mysum[i] = (mysum[i] + mysum[i - 1] * i) / (i + 1);
		myerr_mean[i] = sqrt((mysum_squares[i] - mysum[i] * mysum[i]) / i);
		if (filename != "no_file") result_out << mysum[i] << setw(12) << myerr_mean[i] << endl;
	}

	error_out = myerr_mean[N - 1];
	result = mysum[N - 1];

	result_out.close();
	delete[] mysum_squares;
	delete[] myerr_mean;
	delete[] mysum;

	return result;

}

double V_pot(double x) {

	return pow(x, 4) - 2.5 * pow(x, 2);
}

double E_kin(double x, double sigma, double mu) {

	double arg_forw = pow(x + mu, 2) / pow(sigma, 2);
	double arg_back = pow(x - mu, 2) / pow(sigma, 2);

	return -0.5 * (exp(-0.5 * arg_back) * (arg_back - 1) + exp(-0.5 * arg_forw) * (arg_forw - 1)) / pow(sigma, 2);
}

double trial_wavefunct(double x, double sigma, double mu) {

	return exp(-0.5 * pow(x - mu, 2) / (sigma * sigma)) + exp(-0.5 * pow(x + mu, 2) / (sigma * sigma));
}

double p_density(double x, double sigma, double mu) {

	double psi = trial_wavefunct(x, sigma, mu);
	return psi * psi;
}

double target(double x, double sigma, double mu) {

	return (E_kin(x, sigma, mu) / trial_wavefunct(x, sigma, mu)) + V_pot(x);
}

double boltzmann_delta(double L0, double L1, double T) {

	return exp((L0 - L1)/T);
}

int main() {

	Random rand;
	rand.Initialize();
	
	//----------------file output:
	ofstream out;
	out.open("output_dump_SA.dat", ios::out);
	if (!out.is_open()) {
		cerr << "Errore nell'apertura del file " << "output_dump_SA.dat" << endl;
	}
	out << "T" << setw(12) << "sigma" << setw(12) << "mu" << setw(12) << "H_value" << setw(12) << "error" << setw(12) << "acceptance" << endl;

	//variables:
	int n_schedule = 100;				//dimension of the annealing schedule 
	int n_steps = 100;					//metropolis steps at a single temperature
	int n_samples_int = 10000;			//# samples for integration
	double T[n_schedule];				//array with temperature time series for simulated annealing
	Vect X0(0.0529, 2, 3);			//initial position

	double error;						//integration error to print
	//double cost;
	//Vect state_samples[n_samples];		//samples explored in SA
	//int accept_reject[n_samples];		//step acceptance
	//state_samples[0] = X0;

	Vect s_old = X0;		//old state
	Vect s_new;							//new state
	double L0 = 0;						//old state cost function
	double L1 = 1;						//new state cost function
	double* x_samples;					//sample storage
	

	//begin simulated annealing:
	x_samples = Metropolis_1D(rand, n_samples_int, p_density, s_old);			//sample from psi**2(x)|(sigma(L0), mu(L0)) for integration
	L0 = MC_integrate(target, rand, x_samples, s_old.y, s_old.z, error);
	//save data:
	out << 1 << setw(12) << s_old.y << setw(12) << s_old.z << setw(12) << L0 << setw(12) << error << setw(12) << 1 << endl;
	for (int i = 0; i < n_schedule; i++) {	//iterate on the annealing schedule
		T[i] = 1. - (1. / n_schedule) * i;

		for (int j = 0; j < n_steps; j++) {	
			s_new = Step(rand, s_old);

			x_samples = Metropolis_1D(rand, n_samples_int, p_density, s_new);	//sample from psi**2(x)|(sigma(L0), mu(L0)) for integration
			L1 = MC_integrate(target, rand, x_samples, s_new.y, s_new.z, error);
			//save data:
			out << T[i] << setw(12) << s_new.y << setw(12) << s_new.z
				<< setw(12) << L0 << setw(12) << error;

			if (L1 > L0) {
				if (rand.Rannyu() < boltzmann_delta(L0, L1, T[i])) {
					s_old = s_new;
					L0 = L1;
					//accept_reject[i * n_steps + j] = 1;
					out << setw(12) << 1 << endl;			//accept
				}
				else {
					//accept_reject[i * n_steps + j] = 0;
					out << setw(12) << 0 << endl;			//reject
				}
			}
			else {
				s_old = s_new;
				L0 = L1;
				//accept_reject[i * n_steps + j] = 1;
				out << setw(12) << 1 << endl;			//accept
			}
			
		}
		

	}

	out.close();

	//compute integral at minimum values:
	n_samples_int = 1000000;
	x_samples = Metropolis_1D(rand, n_samples_int, p_density, s_old);	//sample from psi**2(x)|(sigma(L0), mu(L0)) for integration
	L1 = MC_integrate(target, rand, x_samples, s_old.y, s_old.z, error, 10000, 100, "output_integral.dat");

	//save sampling:
	out.open("samples.dat", ios::out);
	if (!out.is_open()) {
		cerr << "Errore nell'apertura del file " << "samples.dat" << endl;
	}
	for (int i = 0; i < n_samples_int; i++) out << x_samples[i] << endl;

	delete[] x_samples;

	return 0;
}