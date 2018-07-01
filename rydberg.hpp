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
#include <boost/random.hpp>
#include <boost/random/gamma_distribution.hpp>
#include <boost/random/poisson_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/mersenne_twister.hpp>

using namespace std;

class OneSite
{
public:
    OneSite();
    OneSite(double Tmax, double dt, double Gamma, double b, double po);
    ~OneSite();
    vector<double> getConfig();
    void reset();
    double update();
    vector<vector<double> > simulation();
    vector<vector<double> > simulation(int N);
private:
    double T, Tmax, Gamma, b, dt, po, p, h;
};

class Lattice1D
{
public:
    Lattice1D();
    Lattice1D(double myT, double myDT, double myGamma, double myB, double myPo, int myNsites, double myDx, double myKappa);
    ~Lattice1D();
    vector <double> getConfig();
    void reset();
    vector<double> update();
    bool is_zero();
    vector<vector<double> > simulation();
    vector<vector<vector<double> > > simulation(int N);
    vector<vector<double> > sim_avg(int N);
    void exciteSite(int site, double seed);
    vector<double> getP();
    vector<double> getH();
    vector<double> euler_update();
private:
    double T, Tmax, Gamma, dt, po, dx, b, kappa;
    int Nsites;
    vector<double> p, h;
};

class Lattice2D
{
public:
    Lattice2D();
    Lattice2D(double myT, double myDT, double myGamma, double myB, double myPo, int L, double myDx);
    ~Lattice2D();
    vector <double> getConfig();
    void reset();
    vector<double> update();
    bool is_zero();
    vector<vector<double> > simulation();
    vector<vector<vector<double> > > simulation(int N);
    vector<vector<double> > sim_avg(int N);
private:
    double T, Tmax, Gamma, dt, po, dx, b;
    int Nsites, L;
    vector<double> p, h;
};

class Lattice3D
{
public:
    Lattice3D();
    Lattice3D(double myT, double myDT, double myGamma, double myB, double myPo, int L, double myDx);
    ~Lattice3D();
    vector <double> getConfig();
    void reset();
    vector<double> update();
    bool is_zero();
    vector<vector<double> > simulation();
    vector<vector<vector<double> > > simulation(int N);
    vector<vector<double> > sim_avg(int N);
private:
    double T, Tmax, Gamma, dt, po, dx, b;
    int Nsites, L;
    vector<double> p, h;
};

double gamma_rand(double shape, double scale, boost::mt19937& rng);
double poiss_rand(double shape, boost::mt19937& rng);
