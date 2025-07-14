#include <iostream>
#include <fstream>
#include <cmath>

#include "../ParallelRandomNumberGenerator/random.h"

using namespace std;

//function that generates the sine of a uniformly distributed angle theta
//as the inclination of the segment between a random point inside the unitary circle and the origin
double rand_sintheta(Random& rand) {
	
	//generate the point in a square around the unitary circle:
	double X1 = rand.Rannyu() * 2 - 1;		//X coordinate of the point, uniformly generated in [-1;1)
	double Y1 = rand.Rannyu() * 2 - 1;		//Y coordinate of the point, uniformly generated in [-1;1)
	double l_squared = X1 * X1 + Y1 * Y1;

	while (l_squared >= 1) {				//reject points outside the circle
		X1 = rand.Rannyu() * 2 - 1;
		Y1 = rand.Rannyu() * 2 - 1;	
		l_squared = X1 * X1 + Y1 * Y1;
	}

	return Y1 / sqrt(l_squared);
}

int main() {

	int M = 10000;			//total number of throws
	int N = 100;			//Number of blocks
	int n = M / N;			//increment per block

	double L = 0.5;			//length of the needle
	double d = 0.8;			//distance between the lines

	double YL;				//stores the projection of the length L of the needle along the Y axis
	double Y1;				//Y coordinate of first generated point
	int n_hit[N] = { 0 };	//stores number of needles that are hitting the lines
	double my_PI[N];
	double squares[N];
	double err[N];

	Random rand;
	rand.Initialize();


	//----------------file output:
	string filename1 = "theta.dat";
	string filename2 = "PI.dat";
	ofstream out(filename1);

	if (!out.is_open()) {
		cerr << "Errore nell'apertura del file " << filename1 << endl;
		return -1;
	}

	//check function generating needles:
	for (int j = 0; j < M; j++) {
		out << rand_sintheta(rand) << "\n";
	}

	out.close();
	

	//----------------compute PI:
	out.open(filename2);

	if (!out.is_open()) {
		cerr << "Errore nell'apertura del file " << filename2 << endl;
		return -1;
	}

	//first for i = 0:
	for (int j = 0; j < n; j++) {
		Y1 = rand.Rannyu() * 10 - 5;									//throw a needle on a 10*10 square centered on the origin
		YL = L * rand_sintheta(rand);									//random needle orientation, compute y coordinate
		if (int(10 + (YL + Y1) / d) != int(10 + Y1 / d)) n_hit[0]++;	//check for hits, the 10 value added is needed to check negative numbers
	}

	if (n_hit[0] == 0) n_hit[0] = 1;	//avoid infinities... not the best but necessary!

	//compute PI from the probability of a needle hitting a line
	my_PI[0] = (n * 2 * L) / (n_hit[0] * d);
	squares[0] = my_PI[0] * my_PI[0];
	out << my_PI[0] << " " << 0 << "\n"; // output file mean and its error (0 for first block)

	//then for each i up to N:
	for (int i = 1; i < N; i++) {
		n_hit[i] = n_hit[i - 1];

		for (int j = 0; j < n; j++) {
			Y1 = rand.Rannyu() * 10 - 5;
			YL = L * rand_sintheta(rand);
			if (int(10 + (YL + Y1) / d) != int(10 + Y1 / d)) n_hit[i]++;	//check for hits
		}

		//compute PI from the probability of a needle hitting a line
		my_PI[i] = (2 * L * n*(i+1)) / (n_hit[i] * d);										
		squares[i] = (my_PI[i] * my_PI[i] + squares[i - 1] * i) / (i + 1);
		my_PI[i] = (my_PI[i] + my_PI[i - 1] * i) / (i + 1);		//accumulate and compute average

		//compute statistical uncertainty:
		err[i] = sqrt((squares[i] - my_PI[i] * my_PI[i]) / i);

		out << my_PI[i] << " " << err[i] << "\n"; // output file mean and its error
	}

	out.close();

	return 0;
}