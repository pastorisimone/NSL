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
#include "system.h"

using namespace std;

int main(int argc, char* argv[]) {


    ifstream file1;
    ofstream file2;
    int n_data = 10000;   //10^5 troppo alto                                                       <|-------------<<<<<<

    //set initial config
    file1.open("../INPUT/CONFIG/config.fcc", ios::in | ios::binary);
    file2.open("../INPUT/CONFIG/config.xyz", ios::out | ios::binary | ios::trunc);
    if (!file1.is_open() || !file2.is_open()) {
        cerr << "Error copying config.fcc" << endl;
        exit(1);
    }
    file2 << file1.rdbuf();
    file1.close();
    file2.close();

    //-----------------------------   es   7.02   --------------------------------
    //----------------------------- MD simulation --------------------------------
    //set input parameters
    file1.open("../INPUT/input_7_02_MD.dat", ios::in | ios::binary);
    file2.open("../INPUT/input.dat", ios::out | ios::binary | ios::trunc);
    if (!file1.is_open() || !file2.is_open()) {
        cerr << "Error copying input_7_02_MD.dat" << endl;
        exit(1);
    }
    file2 << file1.rdbuf();
    file1.close();
    file2.close();


    cout << endl << "Es 7.02 (MD)" << endl;
    //Es 7.02 MD
    System SYS;
    SYS.initialize();
    SYS.initialize_properties();
    SYS.block_reset(0);

    SYS.equilibration(0.01);        //consider the system as "at equilibrium" when T fluctuations are under 1% (arbitrary)
    SYS.inst_reset();

    for (int j = 0; j < n_data; j++) {
        SYS.step();
        SYS.measure();
        SYS.output_U_inst();
        SYS.output_T_inst();
        if (j % (n_data / 10) == 0) cout << endl << j / (n_data / 100) << "%";
    }
    SYS.block_reset(0);
    
    //save results
    rename("../OUTPUT/U_instant.dat", "../OUTPUT/U_instant_7_02_MD.dat");
    rename("../OUTPUT/T_instant.dat", "../OUTPUT/T_instant_7_02_MD.dat");

    //****************  print U, P, GOFR ********************
    for (int i = 0; i < SYS.get_nbl(); i++) { //loop over blocks
        for (int j = 0; j < SYS.get_nsteps(); j++) { //loop over steps in a block
            SYS.step();
            SYS.measure();
        }
        SYS.averages(i + 1);
        SYS.block_reset(i + 1);
    }
    SYS.finalize();

    //save results
    rename("../OUTPUT/potential_energy.dat", "../OUTPUT/potential_energy_MD.dat");
    rename("../OUTPUT/pressure.dat", "../OUTPUT/pressure_MD.dat");
    rename("../OUTPUT/gofr.dat", "../OUTPUT/gofr_MD.dat");



    //----------------------------- MC simulation --------------------------------
    //set input parameters
    file1.open("../INPUT/input_7_02_MC.dat", ios::in | ios::binary);
    file2.open("../INPUT/input.dat", ios::out | ios::binary | ios::trunc);
    if (!file1.is_open() || !file2.is_open()) {
        cerr << "Error copying input_7_02_MC.dat" << endl;
        exit(1);
    }
    file2 << file1.rdbuf();
    file1.close();
    file2.close();

    cout << endl << "Es 7.02 (MC)" << endl;
    //Es 7.02 MC
    System SYS2;
    SYS2.initialize();
    SYS2.initialize_properties();
    SYS2.block_reset(0);

    //SYS2.equilibration(500.);        //consider the system as "at equilibrium" after 500 steps
    for (int j = 0; j < 500; j++) {
        SYS2.step();
        SYS2.measure();
        SYS2.output_U_inst();
        if (j % (n_data / 10) == 0) cout << endl << j / (n_data / 100) << "%";
    }
    SYS2.block_reset(0);

    //save results
    rename("../OUTPUT/U_instant.dat", "../OUTPUT/U_equilibration_MC.dat");
    SYS2.inst_reset();

    for (int j = 0; j < n_data; j++) {
        SYS2.step();
        SYS2.measure();
        SYS2.output_U_inst();
        if (j % (n_data / 10) == 0) cout << endl << j / (n_data / 100) << "%";
    }
    SYS2.block_reset(0);

    //save results
    rename("../OUTPUT/U_instant.dat", "../OUTPUT/U_instant_7_02_MC.dat");

    //****************  print U, P, GOFR ********************
    for (int i = 0; i < SYS2.get_nbl(); i++) { //loop over blocks
        for (int j = 0; j < SYS2.get_nsteps(); j++) { //loop over steps in a block
            SYS2.step();
            SYS2.measure();
        }
        SYS2.averages(i + 1);
        SYS2.block_reset(i + 1);
    }
    SYS2.finalize();

    //save results
    rename("../OUTPUT/potential_energy.dat", "../OUTPUT/potential_energy_MC.dat");
    rename("../OUTPUT/pressure.dat", "../OUTPUT/pressure_MC.dat");
    rename("../OUTPUT/gofr.dat", "../OUTPUT/gofr_MC.dat");


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
