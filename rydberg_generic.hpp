#include <vector>
#include <time.h>
#include <algorithm>
#include <math.h>
#include <cmath>
#include <random>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <cstdio>

using namespace std;

class Rydberg
{
public:
    Rydberg();
    Rydberg(vector<double> config);
    ~Rydberg();
    vector<double> getConfig();
    void reset();
    void update();
    vector<vector<double> > simulation();
    vector<vector<vector<double> > > simulation(int N);
    vector<vector<double> > sim_avg(int N);
private:
    double T, Tmax, Gamma, b, dt, dx, po;
    int dims, len, Nsites;
    vector<double> p, h;
};
