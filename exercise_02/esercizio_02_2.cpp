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

//distance from origin
double squared_distance(double* pos) {
	return pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2];
}

//discrete random walk:
void DRW_update(double* pos, Random& rand, double RW_step = 1.) {	//pos = [x, y, z]

	double r = rand.Rannyu();

	if (r < 1. / 6) pos[0] += RW_step;
	else if (r < 2. / 6) pos[1] += RW_step;
	else if (r < 3. / 6) pos[2] += RW_step;
	else if (r < 4. / 6) pos[0] -= RW_step;
	else if (r < 5. / 6) pos[1] -= RW_step;
	else pos[2] -= RW_step;

}

//continuous random walk:
void CRW_update(double* pos, Random& rand, double RW_step = 1.) {	//pos = [x, y]

	double theta = rand.Rannyu() * 2 * M_PI;						//angle with x axis and y axis
	double phi = rand.Rannyu() * 2 * M_PI;							//angle with z axis and x axis

	pos[0] += RW_step * cos(phi) * cos(theta);
	pos[1] += RW_step * cos(phi) * sin(theta);
	pos[2] += RW_step * sin(phi);

}


int main() {

	int M = 10000;			//total number of points
	int N = 100;			//Number of blocks
	int L = M / N;			//point increment per block

	Random rand;
	rand.Initialize();

	double* pos_DRW = new double[3]{0, 0, 0};
	double* pos_CRW = new double[3]{0, 0, 0};

	//----------------file output:
	string filename1 = "RW_simulation.dat";
	ofstream out(filename1);

	if (!out.is_open()) {
		cerr << "Errore nell'apertura del file " << filename1 << endl;
		return -1;
	}


	//simulate random walks
	for (int i = 0; i < M; i++) {
		DRW_update(pos_DRW, rand);
		CRW_update(pos_CRW, rand);

		out << pos_DRW[0] << " " << pos_DRW[1] << " " << pos_DRW[2] << " " << pos_CRW[0] << " " << pos_CRW[1] << " " << pos_CRW[2] << "\n";
	}

	out.close();
	
	//------------------------------------------COMPUTE RW PROPERTIES

	int n_steps = 100;
	double accu_DRW[n_steps] = { 0 };
	double accu_CRW[n_steps] = { 0 };
	double accu_squares_DRW[n_steps] = { 0 };
	double accu_squares_CRW[n_steps] = { 0 };


	for (int i = 0; i < N; i++) {

		double accu_DRW_block[n_steps] = { 0 };
		double accu_CRW_block[n_steps] = { 0 };

		for (int j = 0; j < L; j++) {

			double* pos_DRW = new double[3] {0, 0, 0};
			double* pos_CRW = new double[3] {0, 0, 0};

			if (j != L-1) {
				for (int step = 0; step < n_steps; step++) {
					DRW_update(pos_DRW, rand);
					CRW_update(pos_CRW, rand);

					accu_DRW_block[step] += squared_distance(pos_DRW) / L;
					accu_CRW_block[step] += squared_distance(pos_CRW) / L;
				}
			}
			else {	//if at the end of a block
				for (int step = 0; step < n_steps; step++) {
					DRW_update(pos_DRW, rand);
					CRW_update(pos_CRW, rand);

					accu_DRW_block[step] += squared_distance(pos_DRW) / L;
					accu_CRW_block[step] += squared_distance(pos_CRW) / L;

					accu_squares_DRW[step] += accu_DRW_block[step] * accu_DRW_block[step] / N;	//only intersted in the last block, so /N
					accu_squares_CRW[step] += accu_CRW_block[step] * accu_CRW_block[step] / N;

					accu_DRW[step] += accu_DRW_block[step] / N;
					accu_CRW[step] += accu_CRW_block[step] / N;

					//if (i % (N / 10) == 0 && step == n_steps - 1) cout << i << "%" << endl;
				}
			}
		}	
	}

	//file output:
	string filename2 = "RW_props.dat";
	out.open(filename2);

	if (!out.is_open()) {
		cerr << "Errore nell'apertura del file " << filename2 << endl;
		return -1;
	}

	//output data:
	for (int step = 0; step < n_steps; step++) {
		double err_DRW = error(accu_DRW[step], accu_squares_DRW[step], N - 1);
		double err_CRW = error(accu_CRW[step], accu_squares_CRW[step], N - 1);

		out << sqrt(accu_DRW[step]) << " " << err_DRW << " " << sqrt(accu_CRW[step]) << " " << err_CRW << endl;
	}

	delete[] pos_DRW;
	delete[] pos_CRW;

	out.close();

	return 0;
}