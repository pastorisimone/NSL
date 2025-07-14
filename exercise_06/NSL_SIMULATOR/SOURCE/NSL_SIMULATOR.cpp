/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "system.h"

using namespace std;

void write_input(double sim_type, double J, double H, double T,
        int n_spins = 50, int N_blk = 20, int N_stp = 10000) {              //writes file input.dat for the simulation

    /*example:
    SIMULATION_TYPE        2/3    1.0    0.0
    RESTART                0
    TEMP                   1.0
    NPART                  50
    RHO                    1.0
    R_CUT                  0.0
    DELTA                  0.0
    NBLOCKS                20
    NSTEPS                 20000
    */

    ofstream input("../INPUT/input.dat");       //write input file for this simulation
    input << "SIMULATION_TYPE " << sim_type << " " << J << " " << H << endl;
    input << "RESTART " << 0 << endl;
    input << "TEMP " << T << endl;
    input << "NPART " << n_spins << endl;
    input << "RHO " << 1.0 << endl;
    input << "R_CUT " << 0. << endl;
    input << "DELTA " << 0. << endl;
    input << "NBLOCKS " << N_blk << endl;
    input << "NSTEPS " << N_stp << endl;
    input << endl << "ENDINPUT" << endl;
    input.close();
}

string* get_last_elements(string readName) {    //reads the last 2 elements from an output file (global average and its error)
    string* result = new string[2];
    ifstream readFile(readName);
    if (!readFile.is_open()) {
        cerr << "ERROR: can't open " << readName << endl;
        return result;
    }
    string lastLine;        //stores the last line (last block) of output.dat
    string temp;

    while (getline(readFile, temp)) {
        if (!temp.empty()) {
            lastLine = temp;
        }
    }
    readFile.close();
    if (lastLine.empty()) {
        cerr << "ERROR: " << readName << " is empty" << endl;
        return result;
    }

    istringstream iss(lastLine);
    int count = 0;

    while (iss >> temp) {
        result[0] = result[1];       //I only need global average and error, the line includes also blk and blk average
        result[1] = temp;
        count++;
    }

    if (count < 2) {
        std::cerr << "Errore: l'ultima riga contiene meno di 2 elementi" << std::endl;
        return result;
    }

    return result;
}

void start_simulation(double input_sim_type, double input_J, double input_H, double input_T,
                    string properties_template, ofstream& writeFile, int eq_steps = 0, bool only_Energy = false) {

    string readName;        //name of the file to read from
    string* lastElements;   //temporary stores last elements of an output

    write_input(input_sim_type, input_J, input_H, input_T);       //write input file for this simulation
    
    //set input properties
    ifstream file1;
    ofstream file2;
    file1.open(properties_template, ios::in | ios::binary);
    file2.open("../INPUT/properties.dat", ios::out | ios::binary | ios::trunc);
    if (!file1.is_open() || !file2.is_open()) {
        cerr << "Error copying " << properties_template << endl;
        exit(1);
    }
    file2 << file1.rdbuf();
    file1.close();
    file2.close();


    //start simulation:
    System SYS;
    SYS.initialize();
    SYS.initialize_properties();
    SYS.block_reset(0);
    //equilibrate:
    for (int i = 0; i < eq_steps; i++) SYS.step();
    //simulate:
    for (int i = 0; i < SYS.get_nbl(); i++) {
        for (int j = 0; j < SYS.get_nsteps(); j++) {
            SYS.step();
            SYS.measure();
        }
        SYS.averages(i + 1);
        SYS.block_reset(i + 1);
    }
    SYS.finalize();

    //save results
    if (only_Energy) {
        string new_filename = "../OUTPUT/total_energy_equilibrium" + to_string(input_sim_type) + ".dat";
        rename("../OUTPUT/total_energy.dat", new_filename.c_str());
    }
    else {  //appends last results of outout file to a target file to get temperature dependence
        if (input_H != 0.) {
            //magnetization:
            lastElements = get_last_elements("../OUTPUT/magnetization.dat");
            writeFile << lastElements[0] << setw(20) << lastElements[1] << endl;
        }
        else {
            //energy:
            lastElements = get_last_elements("../OUTPUT/total_energy.dat");
            writeFile << lastElements[0] << setw(20) << lastElements[1] << setw(20);
            //specific heat:
            lastElements = get_last_elements("../OUTPUT/specific_heat.dat");
            writeFile << lastElements[0] << setw(20) << lastElements[1] << setw(20);
            //susceptibility:
            lastElements = get_last_elements("../OUTPUT/susceptibility.dat");
            writeFile << lastElements[0] << setw(20) << lastElements[1] << setw(20);
        }
    }

}

int main(int argc, char* argv[]) {

    int N_blk = 10;
    int N_stp = 10000;
    int n_spins = 50;
    int n_temps = 16;
    double t_min = 0.5;
    double t_max = 2.;

    double temperatures[n_temps];
    for (int i = 0; i < n_temps; ++i) temperatures[i] = t_min + i * (t_max - t_min) / (n_temps - 1);

    //open file stream to save results to:
    //gibbs:
    ofstream writeFile_gibbs("../OUTPUT/data_vs_temperature_gibbs.dat");
    if (!writeFile_gibbs.is_open()) {
        cerr << "ERROR: can't open ../OUTPUT/data_vs_temperature_gibbs.dat.dat" << endl;
        return 1;
    }
    writeFile_gibbs << "T" << setw(20) << "energy" << setw(20) << "energy_error"
        << setw(20) << "cv" << setw(20) << "cv_error"
        << setw(20) << "suscept" << setw(20) << "suscept_error"
        << setw(20) << "magn" << setw(20) << "magn_error" << endl;
    //metropolis:
    ofstream writeFile_MRT2("../OUTPUT/data_vs_temperature_MRT2.dat");
    if (!writeFile_MRT2.is_open()) {
        cerr << "ERROR: can't open ../OUTPUT/data_vs_temperature.dat_MRT2.dat" << endl;
        return 1;
    }
    writeFile_MRT2 << "T" << setw(20) << "energy" << setw(20) << "energy_error"
        << setw(20) << "cv" << setw(20) << "cv_error"
        << setw(20) << "suscept" << setw(20) << "suscept_error"
        << setw(20) << "magn" << setw(20) << "magn_error" << endl;

    //-------------- first, plot energy to check when the simulation has reached equilibrium ----------------

    start_simulation(3., 1., 0., 1., "../INPUT/properties_6_01.dat", writeFile_gibbs, 0, true);
    start_simulation(2., 1., 0., 1., "../INPUT/properties_6_01.dat", writeFile_MRT2, 0, true);

    // ------------------------------- now get temperature dependence -------------------------------------

    for (int t = 0; t < n_temps; ++t) {
        double T = temperatures[t];
        cout << "temperature: " << T << endl;
        double beta = 1.0 / T;

        //write temperature:
        writeFile_gibbs << T << setw(20);
        writeFile_MRT2 << T << setw(20);

        //**************************** simulation h = 0: U, C, chi (gibbs) ***********************************
        start_simulation(3., 1., 0., T, "../INPUT/properties_6_01.dat", writeFile_gibbs, 3000);

        //******************************* simulation h = 0.02: M (gibbs) *************************************
        start_simulation(3., 1., 0.02, T, "../INPUT/properties_6_02.dat", writeFile_gibbs, 3000);

        //simulations using metropolis samplings
        //**************************** simulation h = 0: U, C, chi (MRT^2) ***********************************
        start_simulation(2., 1., 0., T, "../INPUT/properties_6_01.dat", writeFile_MRT2, 3000);

        //******************************* simulation h = 0.02: M (MRT^2) *************************************
        start_simulation(2., 1., 0.02, T, "../INPUT/properties_6_02.dat", writeFile_MRT2, 3000);

    }

    writeFile_gibbs.close();
    writeFile_MRT2.close();

    return 0;


}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
