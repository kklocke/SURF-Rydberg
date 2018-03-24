#include <iostream>
#include <fstream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <algorithm>
#include <math.h>
#include <random>

using namespace std;

default_random_engine generator;

class OneSite
{
public:
    OneSite() {
        T = 0.;
        Tmax = 500.;
        dt = .1;
        delta  = -.2;
        mu = 1.0;
        l = .03;
        gamma = .03;
        po = .01;
        p = po;
        h = 0.;
    }

    OneSite(double Tmax, double dt, double delta, double mu, double l, double gamma, double po) {
        T = 0.;
        Tmax = Tmax;
        dt = dt;
        delta  = delta;
        mu = mu;
        l = l;
        gamma = gamma;
        po = po;
        p = po;
        h = 0.;
    }

    ~OneSite() {}

    vector<double> getConfig() {
        vector<double> my_config;
        my_config.push_back(dt);
        my_config.push_back(delta);
        my_config.push_back(mu);
        my_config.push_back(l);
        my_config.push_back(gamma);
        my_config.push_back(po);
        return my_config;
    }

    void reset() {
        T = 0.;
        p = po;
        h = 0.;
    }

    double update() {
        double beta = -delta - l * h;
        double lambda = 2.*beta / (pow(gamma, 2) * (exp(beta * dt) - 1.));
        poisson_distribution<int> PoissonDist(lambda * p * exp(beta*dt));
        gamma_distribution<double> GammaDist(PoissonDist(generator), 1.0);
        double pStar = GammaDist(generator) / lambda;
        pStar /= (1 + pStar * dt * mu);
        p = pStar;
        T += dt;
        h += pStar * dt;
        return pStar;
    }

    vector< vector<double> > simulation() {
        int N = int(Tmax / dt);
        vector<double> Tlist(N, 0.);
        vector<double> Plist(N, 0.);
        reset();

        for (int i = 0; i <= N; i++) {
            Tlist[i]  = T;
            Plist[i] = p;
            update();
        }
        vector< vector<double> > res;
        res.push_back(Tlist);
        res.push_back(Plist);
        return res;
    }

    vector<vector < double>> sim_w_scale() {
        vector<double> Tlist, Plist;
        reset();
        int i = 0;
        while (T < Tmax) {
            Tlist[i] = T;
            Plist[i] = p;
            update();
            dt *= sqrt(p/ Plist[i]);
            i++;
        }
        vector<vector<double> > res;
        res.push_back(Tlist);
        res.push_back(Plist);
        return res;
    }

    vector< vector<double> > simulation(int N) {
        vector< vector<double> > res;
        if (N <= 0) {
            return res;
        }
        res = simulation();
        for (int i = 1; i < N; i++) {
            vector< vector<double> > tmp = simulation();
            res.push_back(tmp[1]);
        }
        // return res;
        int Nsteps = int(Tmax / dt);
        vector<double> means(Nsteps, 0.);
        vector<double> stdevs(Nsteps, 0.);
        // Work on speeding this part up later
        for (int i = 0; i < Nsteps; i++) {
            double m = 0.;
            for (int j = 0; j < N; j++) {
                m += res[j+1][i];
            }
            m /= double(N);
            double stdev = 0.;
            for (int j = 0; j < N; j++) {
                stdev += pow(res[j+1][i] - m, 2);
            }
            stdev = sqrt(stdev / (N-1.));
            means[i] = m;
            stdevs[i] = stdev;
        }
        vector<vector<double> > toRet;
        toRet.push_back(res[0]);
        toRet.push_back(means);
        toRet.push_back(stdevs);
        return toRet;
    }


private:
    double T, Tmax, gamma, l, mu, delta, dt, h, p, po;
};


vector<double> mean_stdev(vector<double> dat) {
    vector<double> res;
    double sum = accumulate(dat.begin(), dat.end(), 0.0);
    double m = sum / dat.size();
    double accum = 0.0;
    for_each (dat.begin(), dat.end(), [&](const double d) {accum += (d-m)*(d-m);});
    double stdev = sqrt(accum / (dat.size()-1));
    res.push_back(m);
    res.push_back(stdev);
    return res;
}

class Lattice
{
public:
    Lattice() {
        T = 0.;
        Tmax = 200;
        dt = .2;
        delta  = -.1;
        mu = 1.0;
        l = .1;
        gamma = .03;
        po = .01;
        Nsites = 1e5;
        dx = 1.;
        D = .5;
        p = vector<double>(Nsites, po);
        h = vector<double>(Nsites, 0.);
    }

    Lattice(double Tmax, double dt, double delta, double mu, double l, double gamma, double po, int Nsites, double dx, double D) {
        T = 0.;
        Tmax = Tmax;
        dt = dt;
        delta  = delta;
        mu = mu;
        l = l;
        gamma = gamma;
        po = po;
        Nsites = Nsites;
        dx = dx;
        D = D;
        p = vector<double>(Nsites, po);
        h = vector<double>(Nsites, 0.);
    }

    ~Lattice() {}

    vector<double> getConfig() {
        vector<double> my_config;
        my_config.push_back(dt);
        my_config.push_back(delta);
        my_config.push_back(mu);
        my_config.push_back(l);
        my_config.push_back(gamma);
        my_config.push_back(po);
        my_config.push_back(Nsites);
        my_config.push_back(D);
        my_config.push_back(dx);
        return my_config;
    }

    void reset() {
        T = 0.;
        p = vector<double>(Nsites,po);
        h = vector<double>(Nsites,0.);
    }

    vector<double> update() {
        vector<double> new_p (Nsites, 0.);
        vector<double> new_h (Nsites, 0.);

        for (int i = 0; i < Nsites; i++) {
            double alpha = 0.;
            if (i != 0) {
                alpha += p[i-1];
            }
            if (i != Nsites - 1) {
                alpha += p[i+1];
            }
            alpha *= D / (pow(dx, 2));
            double beta = -delta - l*h[i] - 2*D/pow(dx, 2);
            double lambda = 2*beta / (pow(gamma, 2) * (exp(beta*dt) - 1));
            poisson_distribution<int> PoissonDist(lambda * p[i] * exp(beta*dt));
            gamma_distribution<double> GammaDist(alpha + PoissonDist(generator), 1.0);
            double pStar = GammaDist(generator) / lambda;
            pStar /= (1 + pStar * dt * mu);
            new_p[i] = pStar;
            new_h[i] = h[i] + dt*pStar;
        }
        p = new_p;
        h = new_h;
        T += dt;
        return new_p;
    }

    vector< vector<double> > simulation() {
        int N = int(Tmax / dt);
        vector<double> Tlist(N, 0.);
        vector<vector<double> > Plist;
        reset();

        for (int i = 0; i <= N; i++) {
            Tlist[i] = T;
            Plist.push_back(p);
            update();
        }
        Plist.insert(Plist.begin(), Tlist);
        // [[T's], [p's at T=0], [p's at T=1] ...]
        return Plist;
    }

    vector<vector<vector<double> > > simulation(int N) {
        vector<vector<vector<double> > > res;
        if (N <= 0) {
            return res;
        }
        vector<vector<double> > tmp = simulation();
        vector<vector<double>> Tvec;
        Tvec.push_back(tmp[0]);
        tmp.erase(tmp.begin());
        res.push_back(Tvec);
        res.push_back(tmp);
        for (int i = 1; i < N; i++) {
            tmp = simulation();
            tmp.erase(tmp.begin());
            res.push_back(tmp);
        }
        int NTime = Tvec[0].size();
        vector<vector<double>> allM;
        vector<vector<double>> allS;
        for (int i = 0; i < NTime; i++) {
            vector<double> mVec (Nsites, 0.);
            vector<double> sVec (Nsites, 0.);
            for (int j = 0; j < Nsites; j++) {
                double m = 0.;
                for (int k = 0; k < N; k++) {
                    m += res[k+1][i][j];
                }
                m /= double(N);
                double s = 0.;
                for (int k = 0; k < N; k++) {
                    s += pow(res[k+1][i][j] - m, 2);
                }
                s = sqrt(s / (N-1));
                mVec[j] = m;
                sVec[j] = s;
            }
            allM.push_back(mVec);
            allS.push_back(sVec);
        }
        vector<vector<vector<double> > > toRet;
        toRet.push_back(Tvec);
        toRet.push_back(allM);
        toRet.push_back(allS);
        return toRet;
    }


private:
    double T, Tmax, gamma, l, mu, delta, dt, po, dx, D;
    int Nsites;
    vector<double> p, h;
};







int main() {
    srand(time(NULL));
    time_t t1, t2, t3;
    time(&t1);
    cout << "One site with 5000 simulations\n";
    ofstream datFile;
    datFile.open("rydberg_data_lpt03_deltapt2.txt");
    OneSite R;
    int numSims = 5000;
    vector< vector<double> > dat = R.simulation(numSims);
    for (int i = 0; i < int(dat[0].size()); i++) {
        datFile << dat[0][i] << " " << dat[1][i] << " " << dat[2][i] << "\n";
    }
    datFile.close();
    ofstream configFile;
    configFile.open("rydberg_config_lpt03_deltapt2.txt");
    vector<double> mc = R.getConfig();
    for (int i = 0; i < int(mc.size()); i++) {
        configFile << mc[i] << "\n";
    }
    configFile.close();
    time(&t2);
    cout << "Time elapsed: " << difftime(t2,t1) << "\n";
    // cout << "Lattice with 1e5 sites and 1 simulation\n";
    // ofstream datFile2;
    // datFile2.open("rydberg_lattice_data_1e5.txt");
    // Lattice L;
    // numSims = 1;
    // vector<vector<vector<double> > > dat2 = L.simulation(numSims);
    // int NTimes = dat2[0][0].size();
    // int NSites = dat2[1][0].size();
    // for (int i = 0; i < NSites; i++) {
    //     datFile2 << "# SITE NUMBER: " << i << "\n";
    //     for (int j = 0; j < NTimes; j++) {
    //         datFile2 << dat2[0][0][j] << " " << dat2[1][j][i] << " " << dat2[2][j][i] << "\n";
    //     }
    // }
    // datFile2.close();
    // vector<double> mc2 = L.getConfig();
    // ofstream configFile2;
    // configFile2.open("rydberg_lattice_config_1e5.txt");
    // for (int i = 0; i < int(mc2.size()); i++) {
    //     configFile2 << mc2[i] << "\n";
    // }
    // configFile2.close();
    // time(&t3);
    // cout << "Time elapsed: " << difftime(t3, t2) << "\n";
    return 0;
}
