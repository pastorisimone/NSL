#ifndef ESERCIZIO_09_H
#define ESERCIZIO_09_H

#include <vector>
#include <cmath>
#include <iostream>
#include <utility>
#include <string>

#include "../ParallelRandomNumberGenerator/random.h"

using namespace std;

struct City {
    double x, y;    //city coord.
    int id;         //city identifier

    City();                                 //void constructor
    City(double _x, double _y, int _id);    //constructor
    double distance(City& B);               //computes distance from another city
};

City* write_map(Random& rand, int length = 34, double delim_x = 100., double delim_y = 100., bool circle = false); 

class Path {
public:
    Path();                                 //void constructor
    Path(City* cities, int n_cities);       //constructor
    Path(const Path& copy);                 //copy constructor
    Path(Random& rand, int n_cities = 34, bool circle = false, double delim_x = 100., double delim_y = 100.);   //random constructor
    Path(string filename);                  //construct from file
    ~Path();

    double get_cost(std::string type = "L1");   //get cost of path, supports L1 and L2
    int get_length();                           //get path length
    City get_city(int i);                       //get i-th city
    int get_cityID(int i);                      //get i-th city ID
    City* get_itinerary();                      //get path itinerary as a pointer to an array of cities
    Path& operator=(const Path& other);         //defines assignment operator
    void Swap(int a, int b);                    //swaps city in position a with one in position b
    void set_city(int i, const City& city);     //modify the i-th city
    void printIDs();                            //prints to screen the path as a collection of city IDs

private:
    City* _itinerary;   //collection of cities visited in path order
    int _n;             //dimension
};

class Population {
public:
    Population(int n, Random& rand, int n_cities = 34, bool circle = false, double delim_x = 100., double delim_y = 100.);  //random constructor
    Population(string filename, int n, Random& rand);   //build population from file as permutation of file itinerary
    ~Population();

    void set_p(double p_pair = 0.1, double p_shift = 0.1, double p_inversion = 0.1, double p_cross = 0.5);  //set probability for mutations (seap, shift, inversion) and crossbreeding
    int get_size();                                     //get _n size of population
    double sum_cost(std::string cost = "L1");           //computes sum of cost of all individuals
    double get_minimum(std::string cost = "L1");        //returns the minimum cost of all population
    Path get_individual(int i);                         //get the i-th individual
    void set_individual(Path migrant, int i);           //modify the i-th individual
    void print(std::string cost = "false");             //print to screen the population
    void printf_ind(string output, bool order = false, int index = 0, string cost = "L1");  //print to "output" file the individual in the 
                                                                                            //index position, option to orfer pop. based on cost
    void head(int h, std::string cost = "false");       //print to screen the first h elements, can print the chosen cost if not "false"
    void check_bonds();                                 //checks that population respects the conditions of the problem (same first spot, no duplicates, for each ind.)
    void check_neighbor_order(int ind);                 //checks that, for the ind. position, [ind-1:ind+1] ar in order. if not updates
                                                        //_ordered. Used in mutation operators
    void sort_by_cost(std::string cost = "L1");         //sorts population by given cost

    int Select(double p = 2);                           //selection operator
    pair<Path, Path> Crossover(int mom, int dad);       //crossover operator
    void Mutation_pair(int ind);                        //mutation operator: pair swap
    void Mutation_shift(int ind);                       //mutation operator: sequence shift
    void Mutation_inversion(int ind);                   //mutation operator: block inversion

    //genetic operations:
    //random search with no crossbreeding:
    double random_search(int generations, std::string cost = "L1", double p_pair = 0.3, double p_shift = 0.3, double p_inversion = 0.3);
    //genetic algorithm search of optimum:
    void Evolve(int generations, std::string cost = "L1", string output = "none");

private:
    Random& _rand;
    Path* _individuals;                                             
    int _n;                                                         //total population
    double _p_pair, _p_shift, _p_inversion, _p_cross;               //probability of mutation and of crossover
    bool _ordered;                                                  //if the paths are ordered by cost

    void merge(int left, int mid, int right, std::string cost);     //methods to order the population by cost
    void merge_sort(int left, int right, std::string cost);
};

#endif // ESERCIZIO_09_H
