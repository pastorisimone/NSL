#include <iostream>
#include <fstream>
#include <cmath>
#include <functional>

#include "../ParallelRandomNumberGenerator/random.h"

using namespace std;

struct Vect {		//3d vector

	double x, y, z;

	Vect(): x(0), y(0), z(0) {}
	Vect(double X, double Y, double Z): x(X), y(Y), z(Z) {}

	Vect operator+(Vect& v2) const {
		return Vect(x + v2.x, y + v2.y, z + v2.z);
	}

	Vect operator-(Vect& v2) const {
		return Vect(x - v2.x, y - v2.y, z - v2.z);
	}

	Vect operator*(double scalar) const {
		return Vect(x*scalar, y*scalar, z*scalar);
	}

	double norm2() const {
		return x * x + y * y + z * z;
	}

};

double my_Gauss3D(Vect X, Vect mean, double stdev) {

	double coeff = 1.0 / (pow(stdev * sqrt(2.0 * M_PI), 3));
	double exponent = -(X - mean).norm2() / (2.0 * stdev * stdev);

	return coeff * exp(exponent);

}

double H_ground_state(Vect x) {

	double a0 = 0.0529;
	double coeff = 1. / (pow(a0, 1.5) * sqrt(M_PI));
	double wave = coeff * exp(-sqrt(x.norm2()) / a0);

	return abs(wave) * abs(wave);

}

double H_2p_state(Vect x) {

	double a0 = 0.0529;
	double coeff = sqrt(2/M_PI)/ (pow(a0, 2.5) * 8);
	double r = sqrt(x.norm2());
	double costheta = x.z/r;			
	double wave = coeff * costheta * r * exp(-r / (2 * a0));

	return abs(wave) * abs(wave);

}

Vect* Metropolis(Random& rand, int n_step, double (*p)(Vect), Vect X0, double delta = 0.0529, string T = "Uniform", double equilibrate = 0.5) {

	Vect* X = new Vect[n_step];
	X[0] = X0;
	Vect x_new;
	double a;
	int n_step_eq = int(double(n_step) * equilibrate);
	int accepted = 0;

	if (T == "Uniform") {								//if outside for loop for efficency

		//equilibrate:
		for (int n = 0; n < n_step_eq; n++) {

			Vect step(rand.Rannyu() * 2 * delta - delta, rand.Rannyu() * 2 * delta - delta, rand.Rannyu() * 2 * delta - delta);
			x_new = X[n] + step;

			a = min(1., p(x_new) / p(X[n]));

			double r = rand.Rannyu();
			if (r <= a) X[n + 1] = x_new;					//accept
			else X[n + 1] = X[n];							//reject
		}
		X[0] = X[n_step_eq];

		//actually start generating points after reaching equilibrium
		for (int n = 0; n < n_step - 1; n++) {

			Vect step(rand.Rannyu() * 2 * delta - delta, rand.Rannyu() * 2 * delta - delta, rand.Rannyu() * 2 * delta - delta);
			x_new = X[n] + step;

			a = min(1., p(x_new) / p(X[n]));

			double r = rand.Rannyu();
			if (r <= a) {
				X[n + 1] = x_new;							//accept
				accepted++;
			}
			else X[n + 1] = X[n];							//reject
		}

		cout << endl << T << ": acceptance = " << double(accepted)/double(n_step) << endl;
	}
	else if (T == "Gauss") {

		double stdev = delta;

		//equilibrate:
		for (int n = 0; n < n_step_eq; n++) {

			x_new = Vect(rand.Gauss(X[n].x, stdev), rand.Gauss(X[n].y, stdev), rand.Gauss(X[n].z, stdev));

			a = min(1., my_Gauss3D(X[n], x_new, stdev) * p(x_new) / (my_Gauss3D(x_new, X[n], stdev) * p(X[n])));

			double r = rand.Rannyu();
			if (r <= a) X[n + 1] = x_new;					//accept
			else X[n + 1] = X[n];							//reject
		}
		X[0] = X[n_step_eq];

		//actually start generating points after reaching equilibrium
		for (int n = 0; n < n_step - 1; n++) {

			x_new = Vect(rand.Gauss(X[n].x, stdev), rand.Gauss(X[n].y, stdev), rand.Gauss(X[n].z, stdev));

			a = min(1., my_Gauss3D(X[n], x_new, stdev) * p(x_new) / (my_Gauss3D(x_new, X[n], stdev) * p(X[n])));

			double r = rand.Rannyu();
			if (r <= a) {
				X[n + 1] = x_new;							//accept
				accepted++;
			}
			else X[n + 1] = X[n];							//reject
		}

		cout << endl << T << ": acceptance = " << double(accepted) / double(n_step) << endl;
	}
	else {
		cout << endl << "Transition Probability Error: input types 'Uniform', 'Gauss'" << endl;
		exit(1);
	}

	return X;

}

void MC_integrate_3D(function<double(Vect)> my_funct, Random& rand, Vect* funct_sample, int M = 10000, int N = 100, string filename = "mean.dat") {

	int L = M / N;			//point increment per block			

	//----------------file output:
	ofstream out(filename);

	if (!out.is_open()) {
		cerr << "Errore nell'apertura del file " << filename << endl;
	}
	else {

		//----------------compute averages:

		Vect r;
		double* mysum = new double[N];
		double* mysum_squares = new double[N];
		double* myerr_mean = new double[N];

		//first for i=0:
		mysum[0] = 0;
		for (int j = 0; j < L; j++) {
			r = funct_sample[j];						// generated sample array
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
				r = funct_sample[j+L*i];
				mysum[i] += my_funct(r) / L;
			}

			mysum_squares[i] = (mysum[i] * mysum[i] + mysum_squares[i - 1] * i) / (i + 1);
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

	Random rand;
	rand.Initialize();


	int n_samples = 1000000;
	int n_blk = 100;
	int n_samples_to_print = 0.03*n_samples;
	Vect X0(0.0529, 0.0, 0.0);

	double delta = 0.0529;
	//sampling the ground state:
	Vect* ground_samples = Metropolis(rand, n_samples, H_ground_state, X0, delta);
	ofstream ground_out("H_ground_samples.dat");
	for (int i = 0; i < n_samples_to_print; ++i) {
		ground_out << ground_samples[i].x << " "
			<< ground_samples[i].y << " "
			<< ground_samples[i].z << "\n";
	}
	ground_out.close();

	MC_integrate_3D([](Vect r) { return sqrt(r.norm2()); }, rand, ground_samples, n_samples, n_blk, "r_ground_state_uniform.dat");


	delta = delta * 2.5;
	//sampling the 2p state:
	Vect* p2_samples = Metropolis(rand, n_samples, H_2p_state, X0, delta);
	ofstream p2_out("H_2p_samples.dat");
	for (int i = 0; i < n_samples_to_print; ++i) {
		p2_out << p2_samples[i].x << " "
			<< p2_samples[i].y << " "
			<< p2_samples[i].z << "\n";
	}
	p2_out.close();

	MC_integrate_3D([](Vect r) { return sqrt(r.norm2()); }, rand, p2_samples, n_samples, n_blk, "r_2p_state_uniform.dat");

	//Gaussian T:
	delta = 0.0529 * 0.7;
	Vect* ground_samples_gauss = Metropolis(rand, n_samples, H_ground_state, X0, delta, "Gauss");
	MC_integrate_3D([](Vect r) { return sqrt(r.norm2()); }, rand, ground_samples_gauss, n_samples, n_blk, "r_ground_state_gauss.dat");

	delta = delta * 2;
	Vect* p2_samples_gauss = Metropolis(rand, n_samples, H_2p_state, X0, delta, "Gauss");
	MC_integrate_3D([](Vect r) { return sqrt(r.norm2()); }, rand, p2_samples_gauss, n_samples, n_blk, "r_2p_state_gauss.dat");


	return 0;
}
