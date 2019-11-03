#include "rydberg.hpp"
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_histogram.h>
#include <string>

using namespace std;

double Gamma = -.9; // -0.84701;
double T = 300000.;
double dt = .25;
double b = 0.; // -5e-3;
double po = .0;
double dx = 1.; // sqrt(1.5)
double kappa = 0.; // 2.82e-3;
double exRate = 5e-4; // 5e-4
double useTau = 0.;
int myL = 400;
int Nsites = myL*myL; // 5e3;
// int myL = 500;
// int Nsites = myL;
double trapDepth = 1.;

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

int main(int argc, char* argv[]) {
    /* Part 1:
     * Simulate the system on a 1-D lattice with dx --> infinity so that the sites are decoupled
     * Every few time steps, compute the histogram and plot in gnuplot
     */

    if (argc >= 2) {
        kappa = atof(argv[1]);
    }
    string hFname = "Results/hist_data_";
    string mFname = "Results/mean_data_";
    string rhoHM_fname = "Results/rho_data_";
    string hHM_fname = "Results/h_data_";
    string projx_fname = "Results/h_x_data_";
    string projy_fname = "Results/h_y_data_";
    string std_fname = "Results/std_data_";

    if (argc >= 3) {
        hFname += string(argv[2]) + "_2D.txt";
        mFname += string(argv[2]) + "_2D.txt";
        rhoHM_fname += string(argv[2]) + "_2D.txt";
        hHM_fname += string(argv[2]) + "_2D.txt";
        projx_fname += string(argv[2]) + "_2D.txt";
        projy_fname += string(argv[2]) + "_2D.txt";
        std_fname += string(argv[2]) + "_2D.txt";
    }
    int cutOffTime = 6000;
    if (argc >= 5) {
        trapDepth = atof(argv[3]);
        exRate = atof(argv[4]);
        // exRate = .0000001;
        // cutOffTime = atoi(argv[4]);
    }

    FILE * datFile;
    datFile = fopen(hFname.c_str(), "w");
    gsl_histogram *h = gsl_histogram_alloc (4000);
    // gsl_histogram_set_ranges_uniform(h, -10., 3.);
    gsl_histogram_set_ranges_uniform(h, 1., 3.5);
    if (kappa < .0009) {
        gsl_histogram_set_ranges_uniform(h, 1., 10.);
    }
    // Lattice1D L(T, dt, Gamma, b, po, Nsites, dx, kappa);
    Lattice2D L(T, dt, Gamma, b, po, myL, dx, kappa);
    vector<vector<double> > trackP;
    vector<vector<double> > trackH;
    vector<vector<double> > trackHx;
    vector<vector<double> > trackHy;
    int nsteps = int(T / dt);
    vector<double> myP = L.getP();
    vector<double> myH = L.getH();
    vector<double> myM;
    vector<double> myMH;
    vector<double> myT;
    vector<double> myS;

    // vector<int> toExcite = excitation_sites(T, dt, exRate, Nsites, 0);

    for (int i = 0; i < nsteps; i++) {
        if (i == (int)(4000. / dt)) {
            cout << "setting b\n";
            L.set_b(.005);
            useTau = exRate;
        }
        /* if (i == (int)(float(cutOffTime) / dt)) {
            L.set_b(.0);
            useTau = 0.;
        } */
        if (i % 1 == 0) {
            // gsl_histogram_reset(h);
            double M = 0;
            double Mh = 0;
            for (int j = 0; j < Nsites; j++) {
                // gsl_histogram_increment(h, log10(myP[j] + 1e-9));
                if (i > 4*30000) {
                    gsl_histogram_increment(h, myH[j]);
                }
                M += myP[j];
                Mh += myH[j];
            }
            M /= float(Nsites);
            Mh /= float(Nsites);
            myM.push_back(M);
            myMH.push_back(Mh);
            myT.push_back(float(i)*dt);
            // gsl_histogram_fprintf(datFile, h, "%f", "%f");
            // fprintf(datFile, "TIME: %f\n", float(i)*dt);
        }
        if (i % 4000 == 0) {
            cout << "T: " << float(i)*dt << "\n";
        }
        if ((i % 8 == 0) && (i >= 4*4000) & (i < 4*1000000)) {
           vector<double> tmpH;
           vector<double> tmpHx(myL,0.);
           vector<double> tmpHy(myL,0.);
           if (L.getPmean() > 1e-4) {
               vector<double> tmpP;
               tmpP.push_back(i*dt);
               for (int i1 = 100; i1 < 300; i1++) {
                   for (int i2 = 100; i2 < 300; i2++) {
                       tmpP.push_back(myP[myL*i1+i2]);
                   }
               }
               trackP.push_back(tmpP);
           }
           double tmpSS;
           double tmpS;
           // for (int j = 0; j < Nsites; j++) {
           //     // tmpP.push_back(myP[j]);
           //     // tmpH.push_back(myH[j]);
           //     tmpS += myH[j];
           //     tmpSS += myH[j]*pow(float(j) - float(Nsites / 2.), 2);
           // }
           for (int j1 = 0; j1 < myL; j1++) {
               for (int j2 = 0; j2 < myL; j2++) {
                   tmpS += myH[j1*myL + j2];
                   tmpSS += myH[j1*myL + j2]*(pow(float(j1) - float(myL)/2., 2) + pow(float(j2) - float(myL)/2.,2));
               }
           }
           myS.push_back(sqrt(tmpSS / tmpS));

           /* for (int i1 = 0; i1 < myL; i1++) {
               for (int i2 = 0; i2 < myL; i2++) {
                   int tmpInd = myL*i1 + i2;
                   tmpHx[i1] += myH[tmpInd] / float(myL);
                   tmpHy[i2] += myH[tmpInd] / float(myL);
               }
           } */

           // trackHx.push_back(tmpHx);
           // trackHy.push_back(tmpHy);

           // trackH.push_back(tmpH);
           // myP = L.update();
           // myP = L.getP();
           // myP = L.euler_update();
           // myH = L.getH();
        }
        myP = L.trap_update(useTau, trapDepth);
        myH = L.getH();
        // if (fabs(i * dt - 1600) < .01) {
        //     L.exciteSite(5500, 1.0);
        //     L.exciteSite(5501, 1.0);
        //     L.exciteSite(5502, 1.0);
        //     L.exciteSite(5503, 1.0);
        //     L.exciteSite(5504, 1.0);
        // }
        // if (toExcite[i] != -1) {
        //     int s = toExcite[i];
        //     for (int j = -5; j < 6; j++) {
        //         if ((s + j < 0) || (s + j > Nsites)) {
        //             continue;
        //         }
        //         L.exciteSite(s+j, 1.);
        //     }
        // }
    }
    /* gsl_histogram_fprintf(datFile, h, "%f", "%f");
    // fprintf(datFile, "TIME: %f\n", float(i)*dt);
    // gsl_histogram_free(h);
    fclose(datFile); */
    ofstream mFile;
    mFile.open(mFname);
    for (int i = 0; i < int(myM.size()); i++) {
        mFile << myT[i] << "\t" << myM[i] << "\t" << myMH[i] << "\n";
    }
    mFile.close();
    ofstream pFile;
    pFile.open(rhoHM_fname.c_str());
    for (int i = 0; i < int(trackP.size()); i++) {
        // for (int j = 100; j < 400; j++) {
        for (int j = 0; j < int(trackP[i].size()); j++) {
            pFile << trackP[i][j] << " ";
        }
        pFile << "\n";
    }
    pFile.close();
    /* ofstream hFile;
    hFile.open(hHM_fname.c_str());
    for (int i = 0; i < int(trackH.size()); i++) {
       for (int j = 0; j < int(trackH[i].size()); j++) {
            hFile << trackH[i][j] << " ";
        }
        hFile << "\n";
    }
    hFile.close(); */
    /* vector<double> sx(trackHx.size(),0.);
    vector<double> sy(trackHy.size(),0.);
    ofstream hxFile;
    hxFile.open(projx_fname.c_str());
    for (int i = 0; i < int(trackHx.size()); i++) {
        double mySum = 0.;
        double myStdSum = 0.;
        for (int j = 0; j < int(trackHx[i].size()); j++) {
            hxFile << trackHx[i][j] << " ";
            double locX = float(j) - float(myL)/2.;
            mySum += trackHx[i][j];
            myStdSum += trackHx[i][j] * locX * locX;
        }
        sx[i] = sqrt(myStdSum / mySum);
        hxFile << "\n";
    }
    hxFile.close();
    ofstream hyFile;
    hyFile.open(projy_fname.c_str());
    for (int i = 0; i < int(trackHy.size()); i++) {
        double mySum = 0.;
        double myStdSum = 0.;
        for (int j = 0; j < int(trackHy[i].size()); j++) {
            hyFile << trackHy[i][j] << " ";
            double locX = float(j) - float(myL)/2.;
            mySum += trackHx[i][j];
            myStdSum += trackHx[i][j] * locX * locX;
        }
        sy[i] = sqrt(myStdSum / mySum);
        hyFile << "\n";
    }
    hyFile.close(); */
    ofstream stdFile;
    stdFile.open(std_fname.c_str());
    for (int i = 0; i < int(myS.size()); i++) {
        stdFile << myS[i] << "\n";
        // stdFile << sx[i] << " " << sy[i] << "\n";
    }
    stdFile.close();
    return 0;
}
