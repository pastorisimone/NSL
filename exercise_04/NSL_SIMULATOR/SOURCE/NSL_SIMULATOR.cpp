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
#include <fstream>
#include <sstream>
#include <iomanip>
#include "system.h"

using namespace std;

int main(int argc, char* argv[]) {

    ifstream file1;
    ofstream file2;

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

    //set input parameters
    file1.open("../INPUT/input_4_01.dat", ios::in | ios::binary);
    file2.open("../INPUT/input.dat", ios::out | ios::binary | ios::trunc);
    if (!file1.is_open() || !file2.is_open()) {
        cerr << "Error copying input_4_01.dat" << endl;
        exit(1);
    }
    file2 << file1.rdbuf();
    file1.close();
    file2.close();

    int nconf = 1;
    //---------------------------------------------------------Es 4.01----------------------------------------------------------------
    cout << endl << "Es 4.01" << endl;
    //Es 4.01
    System SYS;
    SYS.initialize();
    SYS.initialize_properties();
    SYS.block_reset(0);

    cout << "blocks: " << SYS.get_nbl() << " steps: " << SYS.get_nsteps() << endl;
    for (int i = 0; i < SYS.get_nbl(); i++) {
        for (int j = 0; j < SYS.get_nsteps(); j++) {
            SYS.step();
            SYS.measure();
            if (j % 50 == 0) {
                // SYS.write_XYZ(nconf); // opzionale
                nconf++;
            }
        }
        SYS.averages(i + 1);
        SYS.block_reset(i + 1);
    }
    SYS.finalize();

    //save results
    rename("../OUTPUT/temperature.dat", "../OUTPUT/temperature_4_01.dat");
    rename("../OUTPUT/pofv.dat", "../OUTPUT/pofv_4_01.dat");

    cout << endl << "Es 4.02" << endl;
    //---------------------------------------------------------Es 4.02----------------------------------------------------------------
    // modify config.xyz to use only half box for every coord.
    double scale_factor = 0.5;

    ifstream infile("../INPUT/CONFIG/config.xyz");
    ofstream outfile("temp.xyz");
    if (!infile.is_open() || !outfile.is_open()) {
        cerr << "Error modifying config.xyz" << endl;
        exit(1);
    }

    string line;
    getline(infile, line); // part. number
    outfile << line << endl;

    getline(infile, line); // comment
    outfile << line << endl;

    string temp;
    double x, y, z;
    while (infile >> temp >> x >> y >> z) {
        outfile << temp << "  "
            << fixed << setprecision(12)
            << x * scale_factor << "  "
            << y * scale_factor << "  "
            << z * scale_factor << endl;
    }

    infile.close();
    outfile.close();

    // overwrite file name
    remove("../INPUT/CONFIG/config.xyz");
    rename("temp.xyz", "../INPUT/CONFIG/config.xyz");

    // set input parameters
    file1.open("../INPUT/input_4_02.dat", ios::in | ios::binary);
    file2.open("../INPUT/input.dat", ios::out | ios::binary | ios::trunc);
    if (!file1.is_open() || !file2.is_open()) {
        cerr << "Error copying input_4_02.dat" << endl;
        exit(1);
    }
    file2 << file1.rdbuf();
    file1.close();
    file2.close();

    // start simulation
    cout << "blocks: " << SYS.get_nbl() << " steps: " << SYS.get_nsteps() << endl;

    System SYS2;
    SYS2.initialize();
    SYS2.initialize_properties();
    SYS2.initialize_velocities_delta();
    SYS2.block_reset(0);
    nconf = 1;

    SYS2.measure();
    SYS2.averages_delta(1);//AGGIUNTO PER PROVARE A VEDERE DISTR. ORIGINALE........
    SYS2.block_reset(1);
    for (int i = 0; i < SYS2.get_nbl(); i++) {
        for (int j = 0; j < SYS2.get_nsteps(); j++) {
            SYS2.step();
            SYS2.measure();
            if (j % 50 == 0) {
                // SYS2.write_XYZ(nconf);
                nconf++;
            }
        }
        SYS2.averages(i + 1);
        SYS2.block_reset(i + 1);
    }
    SYS2.finalize();

    //save results
    rename("../OUTPUT/temperature.dat", "../OUTPUT/temperature_4_02.dat");
    rename("../OUTPUT/potential_energy.dat", "../OUTPUT/potential_energy_4_02.dat");
    rename("../OUTPUT/pofv.dat", "../OUTPUT/pofv_4_02.dat");

    cout << endl << "Es 4.03" << endl;
    //---------------------------------------------------------Es 4.03----------------------------------------------------------------

    //set input parameters
    file1.open("../INPUT/input_4_03.dat", ios::in | ios::binary);
    file2.open("../INPUT/input.dat", ios::out | ios::binary | ios::trunc);
    if (!file1.is_open() || !file2.is_open()) {
        cerr << "Error copying input_4_03.dat" << endl;
        exit(1);
    }
    file2 << file1.rdbuf();
    file1.close();
    file2.close();

    //set initial config
    file1.open("../INPUT/CONFIG/config_4_02.xyz", ios::in | ios::binary);
    file2.open("../INPUT/CONFIG/config.xyz", ios::out | ios::binary | ios::trunc);
    if (!file1.is_open() || !file2.is_open()) {
        cerr << "Error copying config.fcc" << endl;
        exit(1);
    }
    file2 << file1.rdbuf();
    file1.close();
    file2.close();

    file1.open("../INPUT/CONFIG/conf_4_02-1.xyz", ios::in | ios::binary);
    file2.open("../INPUT/CONFIG/conf-1.xyz", ios::out | ios::binary | ios::trunc);
    if (!file1.is_open() || !file2.is_open()) {
        cerr << "Error copying conf-1.xyz" << endl;
        exit(1);
    }
    file2 << file1.rdbuf();
    file1.close();
    file2.close();


    // start simulation
    System SYS3;
    SYS3.initialize();
    SYS3.initialize_properties();
    SYS3.block_reset(0);
    nconf = 1;

    cout << "blocks: " << SYS2.get_nbl() << " steps: " << SYS2.get_nsteps() << endl;

    SYS3.measure();
    SYS3.block_reset(1);
    for (int i = 0; i < SYS3.get_nbl(); i++) {
        for (int j = 0; j < SYS3.get_nsteps(); j++) {
            SYS3.step_back();
            SYS3.measure();
            if (j % 50 == 0) {
                // SYS3.write_XYZ(nconf);
                nconf++;
            }
        }
        SYS3.averages(i + 1);
        SYS3.block_reset(i + 1);
    }
    SYS3.finalize();

    //save results
    rename("../OUTPUT/temperature.dat", "../OUTPUT/temperature_4_03.dat");
    rename("../OUTPUT/potential_energy.dat", "../OUTPUT/potential_energy_4_03.dat");
    rename("../OUTPUT/pofv.dat", "../OUTPUT/pofv_4_03.dat");

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
