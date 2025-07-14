#include <iostream>
#include <fstream>
#include <cmath>
#include <functional>
#include <iomanip>
#include <utility>
#include "mpi.h"

#include "esercizio_09_lib.h"
#include "../ParallelRandomNumberGenerator/random.h"

using namespace std;

int main(int argc, char* argv[]) {

	Random rand;
	rand.Initialize();

	ofstream out_results;
	out_results.open("notable_results.dat");
	if (!out_results) {
		cerr << "ERROR: can't open output results file" << endl;
		return -1;
	}
	out_results << "NODE" << setw(30) << "MODE" << setw(30) << "COST" << endl;
	out_results.close();

	int size, rank;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	
	int n_ind = 300;				//individuals
	int n_send = 7;					//individuals to send on average to each node
	int gens = 8000;				//generations
	int N_migr = 50;				//time step before migration
	int n_epoch = gens / N_migr;	//how many migrations
	int n_best = (size - 1)*n_send;	//best individuals of a "continent"
	Path best;

	int* targets = new int[size];	//list of nodes to which the i-node should broadcast to
	
	if (rank == 0) {
		for (int i = 0; i < size; i++) targets[i] = i;
	}

	//define population:
	string file_input = "cap_prov_ita.dat";
	Population pop_0(file_input, n_ind, rand);

	pop_0.check_bonds();

	int to, from;				//nodes in the exchange: who to send to and the sender

	for (int epoch = 0; epoch < n_epoch; epoch++) {
		for (int i = 0; i < rank * gens; i++) {
			int r = rand.Rannyu();
		}
		pop_0.Evolve(N_migr);

		//select best individuals to migrate:
		pop_0.sort_by_cost();

		for (int b = 0; b < n_best; b++) {
			best = pop_0.get_individual(b);
			
			//convert into an array of doubles to send:
			int length = best.get_length();
			double* send_buffer = new double[3 * length];
			for (int i = 0; i < length; ++i) {
				send_buffer[3 * i] = best.get_city(i).x;
				send_buffer[3 * i + 1] = best.get_city(i).y;
				send_buffer[3 * i + 2] = best.get_city(i).id;
			}

			double* recv_buffer = new double[3 * length];

			if (rank == 0) {
				//find casual target to send to:
				for (int i = 0; i < size; i++) {
					int r = int(rand.Rannyu() * (i + 1));
					int temp = targets[i];
					targets[i] = targets[r];
					targets[r] = temp;
				}
			}

			MPI_Bcast(targets, size, MPI_INT, 0, MPI_COMM_WORLD);

			to = targets[rank];
			from = -1;	//to make the compiler happy
			for (int i = 0; i < size; i++) {
				if (targets[i] == rank) {
					from = i;
					break;
				}
			}

			if (to == rank) {
				delete[] send_buffer;
				delete[] recv_buffer;
				continue;
			}

			//send and wait for receive:
			MPI_Sendrecv(send_buffer, 3 * length, MPI_DOUBLE, to, 0,
				recv_buffer, 3 * length, MPI_DOUBLE, from, 0,
				MPI_COMM_WORLD, MPI_STATUS_IGNORE);

			//rebuild Path object from the array received:
			City* recv_cities = new City[length];
			for (int i = 0; i < length; ++i) {
				recv_cities[i] = City(recv_buffer[3 * i], recv_buffer[3 * i + 1], int(recv_buffer[3 * i + 2]));
			}
			Path migrant(recv_cities, length);

			pop_0.set_individual(migrant, b);

			delete[] send_buffer;
			delete[] recv_buffer;
			delete[] recv_cities;


		}

	}

	double local_min = pop_0.get_minimum();
	cout << "Node " << rank << "/" << size << " estimation of TSP shortest path with migration: " << local_min << endl;
	for (int i = 0; i < size; i++) {
		if (rank == i) {
			out_results.open("notable_results.dat", ios::app);
			if (!out_results) {
				cerr << "ERROR: can't open output results file" << endl;
				return -1;
			}
			out_results << rank << setw(30) << "MIGR" << setw(30) << local_min << endl;
			out_results.close();
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}

	//find the node with global minimum:
	double global_min;
	int best_node;

	MPI_Allreduce(&local_min, &global_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);	//send local min to all nodes,
																					//apply MPI_MIN to it and saves result in
																					//global_min, who is sent back to all nodes

	if (local_min == global_min) best_node = rank;
	else best_node = -1;

	int global_best_node;
	MPI_Allreduce(&best_node, &global_best_node, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

	if (rank == global_best_node) {
		pop_0.printf_ind("TSP_result_MIGR.dat", true);  
	}

	//------------------------------------start with independent GA searches:
	for (int i = 0; i < rank * 20; i++) {
		int r = rand.Rannyu();
	}

	Population pop_1(file_input, n_ind, rand);

	pop_1.check_bonds();
	pop_1.Evolve(gens);
	local_min = pop_1.get_minimum();
	cout << "Node " << rank << "/" << size << " estimation of TSP shortest path without migration: " << local_min << endl;
	for (int i = 0; i < size; i++) {
		if (rank == i) {
			out_results.open("notable_results.dat", ios::app);
			if (!out_results) {
				cerr << "ERROR: can't open output results file" << endl;
				return -1;
			}
			out_results << rank << setw(30) << "NO_MIGR" << setw(30) << local_min << endl;
			out_results.close();
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}

	//find the node with global minimum:
	MPI_Allreduce(&local_min, &global_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);	

	if (local_min == global_min) best_node = rank;
	else best_node = -1;

	MPI_Allreduce(&best_node, &global_best_node, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

	if (rank == global_best_node) {
		pop_1.printf_ind("TSP_result_NO_MIGR.dat", true);
	}

	MPI_Finalize();
	return 0;
}