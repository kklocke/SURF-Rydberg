#include "rydberg.hpp"
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_histogram.h>

using namespace std;

double Gamma = -0.9; // -0.84701;
double T = 20000.;
double dt = .25;
double b = -5e-3;
double po = .0;
double dx = sqrt(1.5);
double kappa = 5e-4;
double exRate = 5e-4;
int Nsites = 1e5;

static int seed1= time(NULL);// 437;
static int seed2= time(NULL); //4357;
static boost::mt19937 rnd_gen1(seed1);
static boost::mt19937 rnd_gen2(seed2);


vector<int> excitation_sites(int Tmax, double dT, double r, int Nsites, double offset) {
    random_device rd;
    mt19937 gen(rd());
    vector<int> sites(int(double(Tmax) / dT),-1);
    uniform_int_distribution<> uniDist(0,Nsites);
    geometric_distribution<> geoDist(r*dT);
    // Generate times
    cout << "Excitation Sites\n";
    double myT = offset;
    while (myT < Tmax) {
        int nextT = geoDist(gen);
        myT += float(nextT)*dT;
        if (myT > float(Tmax)) {
            break;
        }
        sites[int(myT/dt)] = uniDist(gen);
        cout << "Current Time: " << myT << "\tSite: " << sites[int(myT/dt)] << "\n";
    }
    cout << "---------------------------------------\n";
    return sites;
}

int main() {
    /* Part 1:
     * Simulate the system on a 1-D lattice with dx --> infinity so that the sites are decoupled
     * Every few time steps, compute the histogram and plot in gnuplot
     */
    FILE * datFile;
    datFile = fopen("hist_data.txt", "w");
    gsl_histogram *h = gsl_histogram_alloc (1000);
    // gsl_histogram_set_ranges_uniform(h, -10., 3.);
    gsl_histogram_set_ranges_uniform(h, 1, 10);
    Lattice1D L(T, dt, Gamma, b, po, Nsites, dx, kappa);
    // Lattice1D L;
    vector<vector<double> > trackP;
    vector<vector<double> > trackH;
    int nsteps = int(T / dt);
    vector<double> myP = L.getP();
    vector<double> myH = L.getH();
    vector<double> myM;
    vector<double> myT;

    vector<int> toExcite = excitation_sites(T, dt, exRate, Nsites, 0);

    for (int i = 0; i < nsteps; i++) {
        if (i % 1 == 0) {
            // gsl_histogram_reset(h);
            double M = 0;
            for (int j = 0; j < Nsites; j++) {
                // gsl_histogram_increment(h, log10(myP[j] + 1e-9));
                if (i > 40000) {
                    gsl_histogram_increment(h, -Gamma*b*myH[j] + 1 - Gamma + kappa*dt*float(i));
                }
                M += myH[j];
            }
            M /= float(Nsites);
            M *= -Gamma*b;
            M += 1 - Gamma + kappa*dt*float(i);
            myM.push_back(M);
            myT.push_back(float(i)*dt);
            // gsl_histogram_fprintf(datFile, h, "%f", "%f");
            // fprintf(datFile, "TIME: %f\n", float(i)*dt);
        }
        if (i % 4000 == 0) {
            cout << "T: " << float(i)*dt << "\n";
        }
        // if ((i % 1 == 0) && (i < 32000)) {
        //     vector<double> tmpP;
        //     vector<double> tmpH;
        //     for (int j = 5000; j < 6000; j++) {
        //         tmpP.push_back(myP[j]);
        //         tmpH.push_back(myH[j]);
        //     }
        //     trackP.push_back(tmpP);
        //     trackH.push_back(tmpH);
        //     // myP = L.update();
        //     myP = L.getP();
        //     // myP = L.euler_update();
        //     myH = L.getH();
        // }
        myP = L.update();
        myH = L.getH();
        // if (fabs(i * dt - 1600) < .01) {
        //     L.exciteSite(5500, 1.0);
        //     L.exciteSite(5501, 1.0);
        //     L.exciteSite(5502, 1.0);
        //     L.exciteSite(5503, 1.0);
        //     L.exciteSite(5504, 1.0);
        // }
        if (toExcite[i] != -1) {
            int s = toExcite[i];
            for (int j = -5; j < 6; j++) {
                if ((s + j < 0) || (s + j > Nsites)) {
                    continue;
                }
                L.exciteSite(s+j, 1.);
            }
        }
    }
    gsl_histogram_fprintf(datFile, h, "%f", "%f");
    // fprintf(datFile, "TIME: %f\n", float(i)*dt);
    gsl_histogram_free(h);
    fclose(datFile);
    ofstream mFile;
    mFile.open("mean_data.txt");
    for (int i = 0; i < int(myM.size()); i++) {
        mFile << myT[i] << "\t" << myM[i] << "\n";
    }
    mFile.close();
    // ofstream pFile;
    // pFile.open("rho_data.txt");
    // for (int i = 0; i < int(trackP.size()); i++) {
    //     for (int j = 0; j < 1000; j++) {
    //         pFile << trackP[i][j] << " ";
    //     }
    //     pFile << "\n";
    // }
    // pFile.close();
    // ofstream hFile;
    // hFile.open("h_data.txt");
    // for (int i = 0; i < int(trackH.size()); i++) {
    //     for (int j = 0; j < 1000; j++) {
    //         hFile << trackH[i][j] << " ";
    //     }
    //     hFile << "\n";
    // }
    // hFile.close();

    return 0;
}
