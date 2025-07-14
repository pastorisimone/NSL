#include <iostream>
#include <fstream>
#include <cmath>

#include "../ParallelRandomNumberGenerator/random.h"

using namespace std;

double Sample_Exp(Random& rand, double lambda) {
    double u = rand.Rannyu();
    return -log(1 - u) / lambda;                        //inverse cdf, u is F(x)
}

double Sample_CauchyLorentz(Random& rand, double mu, double gamma) {
    double u = rand.Rannyu();
    return mu + gamma * tan(M_PI * (u - 0.5));          //inverse cdf, u is F(x)
}

int main() {
    Random rand;
    rand.Initialize();

    int Nsamples = 10000;  // #datapoints
    
    double lambda = 1.0;  // exponential
    double mu = 0.0;      // average for Lorentz distribution
    double gamma = 1.0;   // width of Lorentz distribution

    ofstream outExp("exp_distribution.dat");
    ofstream outCauchy("cauchy_distribution.dat");
    ofstream outCLT("CLT.dat");

    if (!outExp.is_open()) {
        cerr << "Errore nell'apertura del file exp_distribution.dat" << endl;
        return -1;
    }
    if (!outCauchy.is_open()) {
        cerr << "Errore nell'apertura del file cauchy_distribution.dat" << endl;
        return -1;
    }
    if (!outCLT.is_open()) {
        cerr << "Errore nell'apertura del file CLT.dat" << endl;
        return -1;
    }

    //sample probability distributions:
    for (int i = 0; i < Nsamples; i++) {
        outExp << Sample_Exp(rand, lambda) << endl;
        outCauchy << Sample_CauchyLorentz(rand, mu, gamma) << endl;
    }

    //-----------------------------test the CLT:
    int N[] = { 1, 2, 10, 100 }; //

    Nsamples = 10000;  // #datapoints

    double SN_unif[4];  //accumulators
    double SN_exp[4];
    double SN_CL[4];

    //title line output file
    outCLT << "uniform_1" << " " << "uniform_2" << " " << "uniform_10" << " " << "uniform_100" << " "
        << "exp_1" << " " << "exp_2" << " " << "exp_10" << " " << "exp_100" << " "
        << "CL_1" << " " << "CL_2" << " " << "CL_10" << " " << "CL_100"
        << "\n";

    for (int i = 0; i < Nsamples; i++) {
        for (int steps = 0; steps < 4; steps++) {
            SN_unif[steps] = 0;
            SN_exp[steps] = 0;
            SN_CL[steps] = 0;
            for (int draw = 0; draw < N[steps]; draw++) {
                SN_unif[steps] += rand.Rannyu();
                SN_exp[steps] += Sample_Exp(rand, lambda);
                SN_CL[steps] += Sample_CauchyLorentz(rand, mu, gamma);
            }
        }

        //print output
        outCLT << SN_unif[0] / N[0] << " " << SN_unif[1] / N[1] << " " << SN_unif[2] / N[2] << " " << SN_unif[3] / N[3] << " "
            << SN_exp[0] / N[0] << " " << SN_exp[1] / N[1] << " " << SN_exp[2] / N[2] << " " << SN_exp[3] / N[3] << " "
            << SN_CL[0] / N[0] << " " << SN_CL[1] / N[1] << " " << SN_CL[2] / N[2] << " " << SN_CL[3] / N[3] << "\n";
            
    }

    outExp.close();
    outCauchy.close();

    return 0;
}