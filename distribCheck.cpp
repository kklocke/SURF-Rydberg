#include "rydberg.hpp"
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_histogram.h>
// #include <mgl2/mgl.h>

using namespace std;

double Gamma = .93;
double T = 700.1;
double dt = .7;
double b = 1e-2;
double po = 1e-2;
double dx = 1e2;
int Nsites = 1e5;

int main() {
    /* Part 1:
     * Simulate the system on a 1-D lattice with dx --> infinity so that the sites are decoupled
     * Every few time steps, compute the histogram and plot in gnuplot
     */

    FILE * datFile;
    datFile = fopen("hist_data.txt", "w");
    gsl_histogram *h = gsl_histogram_alloc (100);
    gsl_histogram_set_ranges_uniform(h, -7., 0.);
    Lattice1D L(T, dt, Gamma, b, po, Nsites, dx);
    int nsteps = int(T / dt);
    vector<double> myP;
    for (int i = 0; i < nsteps; i++) {
        myP = L.update();
        if (i % 20 == 0) {
            gsl_histogram_reset(h);
            for (int j = 0; j < Nsites; j++) {
                gsl_histogram_increment(h, log10(myP[j] + 1e-6));
            }
            gsl_histogram_fprintf(datFile, h, "%f", "%f");
            fprintf(datFile, "TIME: %f\n", float(i)*dt);
        }
    }
    gsl_histogram_free(h);
    fclose(datFile);

    /* Part 2:
     * Simulate the system by propagating the distribution at each time step
     * Start with working on a single site
     */
    OneSite S(T, dt, Gamma, b, po);
    return 0;
}
