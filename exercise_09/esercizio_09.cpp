#include <iostream>
#include <fstream>
#include <cmath>
#include <functional>
#include <iomanip>
#include <utility>

#include "../ParallelRandomNumberGenerator/random.h"

using namespace std;

//******************************************************************************************************************************
// ----------------------------- Travelling Salesman Problem: ------------------------------------------------------------------
struct City {
	double x;	//x coord.
	double y;	//y coord.
	int id;		//city identifier

	City() : x(0.), y(0.), id(0.) {}

	City(double _x, double _y, int _id) : x(_x), y(_y), id(_id) {}

	double distance(City& B) {	//compute distance from another city "B"
		if (id == B.id) {
			cerr << "error: distance from city to itself" << endl;
			exit(-1);
		}
		double dx = x - B.x;
		double dy = y - B.y;

		return sqrt(dx * dx + dy * dy);
	}

};

City* write_map(Random& rand, int length = 34, double delim_x = 100., double delim_y = 100., bool circle = false) { 

	City* path = new City[length];
	bool unique;
	double _x;
	double _y;


	if (!circle) {						//generate cities in a square
										//delim_x*delim_y centered in (0;0)
		for (int i = 0; i < length; i++) {
			unique = false;

			while (!unique) {		//check that all cities are different (very bad luck!)
				_x = rand.Rannyu() * delim_x - delim_x/2;
				_y = rand.Rannyu() * delim_y - delim_y/2;

				unique = true;
				for (int j = 0; j < i; j++) {
					if ((_x == path[j].x) && (_y == path[j].y)) {
						unique = false;
						break;
					}
				}
			}

			path[i] = City(_x, _y, i);
		}
	}
	else {
		double R = delim_x / 2.; //sqrt(pow(delim_x, 2) + pow(delim_y, 2)) / 2; circoscritto
		double theta;

		for (int i = 0; i < length; i++) {
			unique = false;

			while (!unique) {		//check that all cities are different (very bad luck!)
				theta = rand.Rannyu() * 2. * M_PI;
				_x = R * cos(theta);
				_y = R * sin(theta);

				unique = true;
				for (int j = 0; j < i; j++) {
					if ((_x == path[j].x) && (_y == path[j].y)) {
						unique = false;
						break;
					}
				}
			}

			path[i] = City(_x, _y, i);
		}
	}

	return path;
}

class Path {
private:
	City* _itinerary;
	int _n;

public:
	//*********builders:
	Path() : _itinerary(nullptr), _n(34) {}														//empty constructor

	Path(City* cities, int n_cities) {															//constructor
		_n = n_cities;
		_itinerary = new City[_n];
		for (int i = 0; i < n_cities; i++) _itinerary[i] = cities[i];
	}
	Path(const Path& copy) {																	//copy constructor
		_n = copy._n;
		_itinerary = new City[_n];
		for (int i = 0; i < _n; i++) _itinerary[i] = copy._itinerary[i];
	}
	Path(Random& rand, int n_cities = 34, bool circle = false, double delim_x = 100., double delim_y = 100.) {	//fast constructor
		_n = n_cities;
		_itinerary = write_map(rand, n_cities, delim_x, delim_y, circle);
	}
	~Path() {
		delete[] _itinerary;
	}


	//*********gets:
	double get_cost(string type = "L1") {						//computes cost function of a path
																//if L1 = sum(abs(X_a - X_b))
																//else L2 = sum((X_a - X_b)**2)
		double L = 0;

		if (type == "L1") {
			for (int i = 0; i < _n - 1; i++) {
				L += _itinerary[i].distance(_itinerary[i + 1]);
			}
			L += _itinerary[_n - 1].distance(_itinerary[0]);
		}
		else if (type=="L2") {
			for (int i = 0; i < _n - 1; i++) {
				L += pow(_itinerary[i].distance(_itinerary[i + 1]), 2);
			}
			L += pow(_itinerary[_n - 1].distance(_itinerary[0]), 2);
		}
		else {
			cerr << "wrong cost function type in Path.get_cost(): supported types are 'L1' (distance) and 'L2' (norm)" << endl;
		}

		return L;
	}

	int get_length() {
		return _n;
	}

	City get_city(int i) {
		if (i < _n) return _itinerary[i];
		else {
			cerr << "index " << i << " out of bounds for Path of length " << _n << endl;
			exit(1);
		}
	}

	int get_cityID(int i) {
		if (i < _n) return _itinerary[i].id;
		else {
			cerr << "index " << i << " out of bounds for Path of length " << _n << endl;
			exit(1);
		}
	}

	City* get_itinerary() {
		City* itinerary = new City[_n];
		for (int i = 0; i < _n; i++) itinerary[i] = _itinerary[i];
		return itinerary;
	}

	//*********operator =:
	Path& operator=(const Path& other) {
		if (this == &other)
			return *this;

		delete[] _itinerary;

		_n = other._n;
		_itinerary = new City[_n];
		for (int i = 0; i < _n; i++)
			_itinerary[i] = other._itinerary[i];

		return *this;
	}

	//********Swap
	void Swap(int a, int b) {
		City temp = _itinerary[a];
		_itinerary[a] = _itinerary[b];
		_itinerary[b] = temp;
	}

	//*********Utility
	void set_city(int i, const City& city) {
		if (i < _n) _itinerary[i] = city;
		else {
			cerr << "index " << i << " out of bounds for Path of length " << _n << endl;
			exit(1);
		}
	}

	void printIDs() {
		cout << endl << "[ ";
		for (int i = 0; i < _n - 1; i++)
			cout << _itinerary[i].id << ", ";
		cout << _itinerary[_n-1].id << "]" << endl;
	}
};

// *****************************************************************************************************************************
// -----------------------------      Genetic Algorithm:      ------------------------------------------------------------------
class Population {
private:

	Random& _rand;
	Path* _individuals;						
	int  _n;								//total population
	double _p_pair, _p_shift, _p_inversion;	//probability of mutation
	double _p_cross;						//probability of crossover
	bool _ordered;	


	void merge(int left, int mid, int right, string cost) {			//merge for mergeSort ordering of every individual cost
		int n1 = mid - left + 1;
		int n2 = right - mid;

		Path* L = new Path[n1];
		Path* R = new Path[n2];

		for (int i = 0; i < n1; i++) L[i] = _individuals[left + i];
		for (int i = 0; i < n2; i++) R[i] = _individuals[mid + i + 1];

		int i = 0, j = 0, k = left;
		while ((i < n1) && (j < n2)) {
			if (L[i].get_cost(cost) <= R[j].get_cost(cost)) {		//più efficace se cost già salvato..    <|----------<<<<<
				_individuals[k] = L[i];								//AGGIUNGI COST TYPE
				k++;
				i++;
			}
			else {
				_individuals[k] = R[j];
				k++;
				j++;
			}
		}

		while (i < n1) {
			_individuals[k] = L[i];
			k++;
			i++;
		}
		while (j < n2) {
			_individuals[k] = R[j];
			k++;
			j++;
		}

		delete[] L;
		delete[] R;
	}

	void merge_sort(int left, int right, string cost) {				//mergeSort algorithm based on cost function
		if (left >= right) {
			return;
		}

		int mid = (left + right) / 2;
		merge_sort(left, mid, cost);
		merge_sort(mid + 1, right, cost);
		merge(left, mid, right, cost);

	}

public:
	//*********builders:
	Population(int n, Random& rand, int n_cities = 34, 
		bool circle = false, double delim_x = 100., double delim_y = 100.): _rand(rand) {		//builder
		_ordered = false;
		_n = n;
		_p_pair = 0.1;
		_p_shift = 0.1;
		_p_inversion = 0.1;
		_p_cross = 0.5;
		City* original = write_map(_rand, n_cities, delim_x, delim_y, circle);	

		_individuals = new Path[_n];
		for (int i = 0; i < _n; i++) {
			City* permutation = new City[n_cities];
			for (int j = 0; j < n_cities; j++) permutation[j] = original[j];
			for (int j = 2; j < n_cities; j++) {
				int index = 1 + int(rand.Rannyu() * (j));					//Fisher-Yates leaving the first element unchanged
				swap(permutation[j], permutation[index]);
			}

			_individuals[i] = Path(permutation, n_cities);
			delete[] permutation;
		}

		delete[] original;
	}
	~Population() {
		delete[] _individuals;
	}

	void set_p(double p_pair = 0.1, double p_shift = 0.1, double p_inversion = 0.1, double p_cross = 0.5) {
		if (p_pair + p_shift + p_inversion + p_cross > 1) {
			cerr << "ERROR: sum of mutation probabilities is grater than 1" << endl;
			exit(-1);
		}

		_p_pair = p_pair;
		_p_shift = p_shift;
		_p_inversion = p_inversion;
		_p_cross = p_cross;
	}

	//*********gets:
	int get_size() {
		return _n;
	}

	double sum_cost(string cost = "L1", double size = 1., bool squared = false) {		//computes sum of cost of first size*n elements
		if (size > 1.) {
			cerr << "ERROR: requested size > pop size" << endl;
			return -1;
		}
		double sum = 0;
		for (int i = 0; i < int(size * _n); i++) {
			if (!squared) sum += _individuals[i].get_cost(cost);
			else sum += _individuals[i].get_cost(cost) * _individuals[i].get_cost(cost);;
		}	
		return sum;
	}

	double get_minimum(string cost = "L1") {
		if (!_ordered) this->sort_by_cost(cost);
		return _individuals[0].get_cost(cost);
	}

	Path get_individual(int i) {
		if (i >= _n) {
			cerr << "index out of bounds for Population of size " << _n << endl;
		}
		else return _individuals[i];
	}

	//*********utility:
	void print(string cost = "false") {
		cout << endl;
		for (int i = 0; i < _n; i++) {
			cout << "[";
			for (int j = 0; j < _individuals[i].get_length()-1; j++) {
				cout << " " << _individuals[i].get_cityID(j) << ",";
			}
			cout << " " << _individuals[i].get_cityID(_individuals[i].get_length() - 1) << "]";
			if (cost != "false"){
				cout << " cost: " << cost << " " << _individuals[i].get_cost(cost);
			}
			cout << endl;
		}
		cout << endl;
	}

	void printf_ind(string output, bool order = false, int index = 0, string cost = "L1") {	//saves individual to file
		if (order == true && _ordered == false) this->sort_by_cost(cost);
		Path ind = _individuals[index];

		ofstream out(output);
		if (!out) {
			cerr << "ERROR: can't open " << output << endl;
			exit(1);
		}

		for (int i = 0; i < ind.get_length(); i++) {
			City cit = ind.get_city(i);
			out << cit.id << setw(12) << cit.x << setw(12) << cit.y << endl;
		}
		
		out.close();
	}

	void head(int h, string cost = "false") {	//prints the first h elements
		if (h > _n) h = _n;
		cout << endl;
		for (int i = 0; i < h; i++) {
			cout << "[";
			for (int j = 0; j < _individuals[i].get_length() - 1; j++) {
				cout << " " << _individuals[i].get_cityID(j) << ",";
			}
			cout << " " << _individuals[i].get_cityID(_individuals[i].get_length() - 1) << "]";
			if (cost != "false") {
				cout << " cost: " << cost << " " << _individuals[i].get_cost(cost);
			}
			cout << endl;
		}
		cout << endl;
	}

	void check_bonds() {		//checks every individual to have initial city 0 and to visit every city once
		for (int ind = 0; ind < _n; ind++) {
			//City start = _individuals[ind].get_city(0);
			if (_individuals[ind].get_cityID(0) != 0) {
				cerr << "wrong start!" << endl;
				exit(-1);
			}
			for (int i = 0; i < _individuals[ind].get_length(); i++) {
				//City A = _individuals[ind].get_city(i);
				for (int j = i + 1; j < _individuals[ind].get_length(); j++) {
					//City B = _individuals[ind].get_city(j);
					if (_individuals[ind].get_cityID(i) == _individuals[ind].get_cityID(j)) {
						cerr << "not unique!" << endl;
						exit(-1);
					}
				}
			}
		}
		cout << endl << "Bonds correctly checked!" << endl;
	}

	void check_neighbor_order(int ind) {									//to check if mutation changes order:
		if (_ordered) {
			double this_cost = _individuals[ind].get_cost();
			if (ind == 0) {
				if (this_cost > _individuals[ind + 1].get_cost()) _ordered = false;
			}
			else if (ind == _n - 1) {
				if (this_cost < _individuals[ind - 1].get_cost()) _ordered = false;
			}
			else {
				if ((this_cost > _individuals[ind + 1].get_cost()) or (this_cost < _individuals[ind - 1].get_cost())) _ordered = false;
			}
		}
	}

	void sort_by_cost(string cost = "L1") {
		merge_sort(0, _n - 1, cost);
		_ordered = true;
	}

	//*********genetic operators:
	int Select(double p = 2, string cost = "L1") {
		if (!_ordered) this->sort_by_cost(cost);
		double r = _rand.Rannyu();
		return int(_n * pow(r, p));											
	}

	pair<Path, Path> Crossover(int mom, int dad) {						//crossover operator: returns progeny with mixed genetic material
		int length = _individuals[mom].get_length();
		City* son = new City[length];
		City* daughter = new City[length];

		int split = _rand.Rannyu() * (length - 1);
		//if (split == 0 or split == length - 1) return make_pair(_individuals[mom], _individuals[dad]);

		for (int i = 0; i < split; i++) {
			daughter[i] = _individuals[mom].get_city(i);
			son[i] = _individuals[dad].get_city(i);
		}

		City c_mom;
		City c_dad;
		int d_counter = split;
		int s_counter = split;
		for (int i = 0; i < length; i++) {
			bool c_mom_found = false;
			bool c_dad_found = false;
			c_mom = _individuals[mom].get_city(i);
			c_dad = _individuals[dad].get_city(i);
			for (int j = 0; j < split; j++) {
				if (c_mom.id == son[j].id) {
					c_mom_found = true;
				}
				if(c_dad.id == daughter[j].id) {
					c_dad_found = true;
				}
			}
			if (!c_mom_found) {
				son[s_counter] = c_mom;
				s_counter++;
			}
			if (!c_dad_found) {
				daughter[d_counter] = c_dad;
				d_counter++;
			}
		}

		Path S(son, length);
		Path D(daughter, length);

		delete[] son;
		delete[] daughter;

		return make_pair(S, D);
	}

	void Mutation_pair(int ind) {							//pair permutation
		int length = _individuals[ind].get_length();				
		int a = 1 + int(_rand.Rannyu() * (length - 1));
		int b = a;
		while(a==b) b = 1 + int(_rand.Rannyu() * (length - 1));
		_individuals[ind].Swap(a, b);

		//cout << endl << " (a, b, length): " << a << " " << b << " " << length << endl;			<|-----------<<<<<
		this->check_neighbor_order(ind);
	}
	
	void Mutation_shift(int ind) {							//shifts m cities of +n positions
		int l = _individuals[ind].get_length();
		int j = 1 + int(_rand.Rannyu() * (l - 3 ));						//where to mutate (excluded first because fixed, and last because can't be shifted)
		int m = 1 + int(_rand.Rannyu() * ((l - j)/2 - 1));				//how many to shift
		int n = 1 + int(_rand.Rannyu() * ((l - j)/2 - m - 1));			//how much to shift

		for (int k = j; k < j + m; k++) _individuals[ind].Swap(k, k + m + n);

		this->check_neighbor_order(ind);
	}


	void Mutation_inversion(int ind) {						//invertion of the order of m cities
		int l = _individuals[ind].get_length();
		int j = 1 + int(_rand.Rannyu() * (l - 3));			//where to mutate								<|---------<<<< l-2?
		int m = 1 + int(_rand.Rannyu() * (l - j - 1));		//how many to invert

		for (int k = 0; k < m / 2; k++) _individuals[ind].Swap(k + j, j + m - k - 1);

		this->check_neighbor_order(ind);
	}

	//*********genetic operations:

	//random search for optimum using only casual mutations to explore solution space
	double random_search(int generations, string cost = "L1", double p_pair = 0.3, double p_shift = 0.3, double p_inversion = 0.3) {
		this->set_p(p_pair, p_shift, p_inversion, 0);
		double r;
		this->sort_by_cost();

		for (int gen = 0; gen < generations; gen++) {
			for (int i = 1; i < _n; i++) {														//leave the first element, (gen n-1's optimal 
				r = _rand.Rannyu();																//solution), unchanged
				if (r < _p_pair) this->Mutation_pair(i);										//randomly select a mutation, if any
				else if (r < _p_pair + _p_shift) this->Mutation_shift(i);
				else if (r < _p_pair + _p_shift + _p_inversion) this->Mutation_inversion(i);
			}

			if (!_ordered) this->sort_by_cost();
		}

		this->set_p();
		return _individuals[0].get_cost();

	}

	//complete evolution of the population
	void Evolve(int generations, string cost = "L1", string output = "none") {
		int a, b, r;

		ofstream out(output);
		if (!out) {
			cerr << "ERROR: can't open output file" << endl;
			exit(1);
		}

		for (int gen = 0; gen < generations; gen++){
			Path* new_gen = new Path[_n];
	
			//crossover:
			if (!_ordered) this->sort_by_cost("L1");

			for (int i = 0; i < _n/2; i++) {		//each steps creates two scions
				a = this->Select();					//select random pair with probability depending on ordering (=cost function)
				b = this->Select();
				if (_rand.Rannyu() < _p_cross){		//crossover with probability p_cross
					auto scions = Crossover(a, b);
					new_gen[i] = scions.first;
					new_gen[i + _n/2] = scions.second;
				}
				else {
					new_gen[i] = _individuals[a];
					new_gen[i + _n / 2] = _individuals[b];
				}
			}
			if (_n % 2 != 0) {						//odd population case
				a = this->Select();
				b = this->Select();
				if (_rand.Rannyu() < _p_cross) {
					auto scions = Crossover(a, b);
					new_gen[_n - 1] = scions.first;
				}
				else {
					new_gen[_n - 1] = _individuals[a];
				}
			}

			_ordered = false;
			delete[] _individuals;
			_individuals = new_gen;

			//mutation:
			for (int i = 0; i < _n; i++) {					
				r = _rand.Rannyu();																
				if (r < _p_pair) this->Mutation_pair(i);										
				else if (r < _p_pair + _p_shift) this->Mutation_shift(i);
				else if (r < _p_pair + _p_shift + _p_inversion) this->Mutation_inversion(i);
			}

			if (output != "none") {
				if (!_ordered) this->sort_by_cost(cost);
				out << this->sum_cost(cost, 0.5) / (0.5*_n) << setw(20) << (sum_cost(cost, 0.5, true) - sum_cost(cost, 0.5) * sum_cost(cost, 0.5)) / (0.5 * _n) << setw(20) << this->get_minimum() << endl;
			}
		}

		out.close();
	}

};

// *****************************************************************************************************************************
// -----------------------------             Main:            ------------------------------------------------------------------

int main() {

	Random rand;
	rand.Initialize();

	ofstream out_results("notable_results.dat");	
	if (!out_results) {
		cerr << "ERROR: can't open output results file" << endl;
		return -1;
	}
	out_results << "ALG_TYPE" << setw(30) << "CONFIG" << setw(30) << "GENS" << setw(30) << "COST" << endl;

	Population pop_0(10, rand);
	pop_0.print("L1");
	pop_0.check_bonds();

	pop_0.sort_by_cost();
	pop_0.print("L1");

	pop_0.Mutation_pair(5);
	pop_0.Mutation_shift(0);
	pop_0.Mutation_inversion(8);
	pop_0.Mutation_inversion(9);
	pop_0.print("L1");
	pop_0.check_bonds();

	cout << endl << "********************************** RANDOM SEARCH **********************************" << endl;

	int gens = 100;
	double result = pop_0.random_search(gens);
	pop_0.print("L1");
	pop_0.check_bonds();
	cout << "result of random search after " << gens << " generations: cost = " << result << endl;
	out_results << "RANDOM_SEARCH" << setw(30) << "SQUARE" << setw(30) << gens << setw(30) << result << endl;

	cout << endl << "********************************** EVOLUTION **********************************" << endl;

	int n_paths = 200;
	gens = 10000;
	Population pop_1(n_paths, rand);
	pop_1.check_bonds();
	//save initial solution
	pop_1.printf_ind("square_TSP_init.dat");

	pop_1.sort_by_cost();
	if (n_paths > 20) pop_1.head(10);
	else pop_1.print("L1");

	pop_1.Evolve(gens, "L1", "TSP_GA_cost_square.dat");

	pop_1.sort_by_cost();
	cout << "after " << gens << " generations of evolution: " << endl;
	if (n_paths > 20) pop_1.head(10);
	else pop_1.print("L1");
	cout << "result of genetic algorithm after " << gens << " generations: cost = " << pop_1.get_minimum() << endl;
	out_results << "GENETIC_ALG" << setw(30) << "SQUARE " << setw(30) << gens << setw(30) << pop_1.get_minimum() << endl;
	
	//save last result:
	pop_1.printf_ind("square_TSP_best.dat", true);

	cout << endl << "**************************** EVOLUTION ON A CIRCLE ***************************" << endl;

	gens = 10000;

	Population pop_c(n_paths, rand, 34, true);
	//save initial solution
	pop_c.printf_ind("circle_TSP_init.dat");
	pop_c.check_bonds();

	pop_c.sort_by_cost();
	if (n_paths > 20) pop_c.head(10);
	else pop_c.print("L1");

	//save cost evolution:
	pop_c.Evolve(gens, "L1", "TSP_GA_cost_circle.dat");

	pop_c.sort_by_cost();
	cout << "after " << gens << " generations of evolution: " << endl;
	if (n_paths > 20) pop_c.head(10);
	else pop_c.print("L1");
	cout << "result of genetic algorithm after " << gens << " generations: cost = " << pop_c.get_minimum() << endl;
	out_results << "GENETIC_ALG" << setw(30) << "CIRCLE " << setw(30) << gens << setw(30) << pop_c.get_minimum() << endl;
	cout << "circonference = " << 100 * M_PI << endl;

	Population pop_f(n_paths, rand, 34, true);
	result = pop_f.random_search(gens);
	cout << "result of random search after " << gens << " generations: cost = " << result << endl;
	out_results << "RANDOM_SEARCH" << setw(30) << "CIRCLE " << setw(30) << gens << setw(30) << result << endl;

	//save best result:
	pop_c.printf_ind("circle_TSP_best.dat", true);

	return 0;
}
