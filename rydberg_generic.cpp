#include "rydberg_generic.hpp"

default_random_engine generator;

#define ndx(obj,i,j,k,n)    (obj[(i) + (j)*n + (k)*n*n])
#define inds2n(i,j,k,l)     (i + j*l + k*l*l)

vector<int> map_to_inds(int N, int dim) {
    int k = int(N / (dim * dim));
    int i_plus_j = N % (dim * dim);
    int j = int(i_plus_j / dim);
    int i = i_plus_j % dim;
    vector<int> ret(3, 0.);
    ret[0] = i;
    ret[1] = j;
    ret[2] = k;
    return ret;
}

Rydberg::Rydberg() {
    Gamma = .4;
    dt = .02;
    dx = 1.;
    po = .01;
    b = .001;
    T = 0.;
    Tmax = 200.;
    dims = 1;
    len = 1;
    Nsites = int(pow(len, dims));
    p = vector<double> (Nsites, po);
    h = vector<double> (Nsites, 0.);
}

Rydberg::Rydberg(vector<double> config) {
    Tmax = config[0];
    dt = config[1];
    dx = config[2];
    po = config[3];
    dims = config[4];
    len = config[5];
    Gamma = config[6];
    b = config[7];
    Nsites = int(pow(len, dims));
    T = 0.;
    p = vector<double> (Nsites, po);
    h = vector<double> (Nsites, 0.);
}

Rydberg::~Rydberg() {}

vector<double> Rydberg::getConfig() {
    vector<double> ret(8, 0.);
    ret[0] = Tmax;
    ret[1] = dt;
    ret[2] = dx;
    ret[3] = double(dims);
    ret[4] = double(len);
    ret[5] = Gamma;
    ret[6] = b;
    ret[7] = po;
    return ret;
}

void Rydberg::reset() {
    T = 0.;
    for (int i = 0; i < Nsites; i++) {
        p[i] = po;
        h[i] = 0.;
    }
}

void Rydberg::update() {
    vector<double> new_p (Nsites, 0.);
    vector<double> new_h (Nsites, 0.);
    int i = 0;
    int j = 0;
    int k = 0;
    for (int ind = 0; ind < Nsites; ind++) {
        vector<int> my_inds = map_to_inds(ind, len);
        i = my_inds[0];
        j = my_inds[1];
        k = my_inds[2];
        double pSum = 0.;
        vector<int> adjs;
        if (dims >= 1) {
            adjs.push_back(inds2n(i-1, j, k, len));
            adjs.push_back(inds2n(i+1, j, k, len));
        }
        else if (dims >= 2) {
            adjs.push_back(inds2n(i, j-1, k, len));
            adjs.push_back(inds2n(i, j+1, k, len));
        }
        else if (dims == 3) {
            adjs.push_back(inds2n(i, j, k-1, len));
            adjs.push_back(inds2n(i, j, k+1, len));
        }
        for (int adj_ind = 0; adj_ind < int(adjs.size()); adj_ind++) {
            if ((adjs[adj_ind] > 0) && (adjs[adj_ind] < Nsites)) {
                pSum += p[adj_ind];
            }
        }
        double alpha = pSum / pow(dx, 2);
        double beta = 1 - Gamma - Gamma*b*h[ind] - 2*dims / pow(dx, 2);
        double mu = -1 + 2*alpha;
        double lambda = 2*beta / (pow(Gamma, 2) * (exp(beta*dt) - 1));
        poisson_distribution<int> PoissonDist(lambda * p[ind] * exp(beta*dt));
        gamma_distribution<double> GammaDist(mu + 1 + PoissonDist(generator), 1.0);
        double pStar = GammaDist(generator) / lambda;
        pStar /= (1 + pStar * dt * 2);
        new_p[ind] = pStar;
        new_h[ind] = h[ind] + dt*pStar;
    }
    p = new_p;
    h = new_h;
    T += dt;
}


vector<vector<double> > Rydberg::simulation() {
    int N = int(Tmax / dt);
    vector<double> Tlist(N, 0.);
    vector<vector<double> > Plist;
    reset();
    for (int i = 0; i < N; i++) {
        Tlist[i] = T;
        Plist.push_back(p);
        if ((Nsites == 1) && (p[0] = 0.)) {
            T += dt;
        }
        else {
            update();
        }
    }
    Plist.insert(Plist.begin(), Tlist);
    return Plist;
}

vector<vector<vector<double> > > Rydberg::simulation(int N) {
    vector<vector<vector<double> > > res;
    if (N <= 0) {
        return res;
    }
    for (int i = 0; i < N; i++) {
        res.push_back(simulation());
    }
    // take average and stdev
    int Nsteps = int(Tmax / dt);
    vector<vector<double>> all_means;
    vector<vector<double>> all_stds;
    vector<vector<double>> all_T;
    all_T.push_back(res[0][0]);

    for (int j = 0; j < Nsteps; j++) {
        vector<double> means(Nsites, 0.);
        vector<double> stds(Nsites, 0.);
        for (int i = 0; i < Nsites; i++) {
            double myMean = 0.;
            for (int k = 0; k < N; k++) {
                myMean += res[k][j+1][i];
            }
            myMean /= double(N);
            double myStd = 0.;
            for (int k = 0; k < N; k++) {
                myStd += pow(res[k][j+1][i], 2);
            }
            myStd = sqrt(myStd / (N-1));
            means[i] = myMean;
            stds[i] = myStd;
        }
        all_means.push_back(means);
        all_stds.push_back(stds);
    }
    vector<vector<vector<double> > > ret;
    ret.push_back(all_T);
    ret.push_back(all_means);
    ret.push_back(all_stds);
    return ret;
}

vector<vector<double> > Rydberg::sim_avg(int N) {
    vector<vector<vector<double> > > sim_res = simulation(N);
    vector<double> myT = sim_res[0][0];
    int Nsteps = int(Tmax / dt);
    vector<double> myMeans(Nsteps, 0.);
    vector<double> myStds(Nsteps, 0.);
    for (int i = 0; i < Nsteps; i++) {
        double m = 0.;
        double s = 0.;
        for (int j = 0; j < Nsites; j++) {
            m += sim_res[1][i][j];
            s += pow(sim_res[2][i][j],2);
        }
        m /= double(Nsites);
        s = sqrt(s) / double(Nsites);
        myMeans[i] = m;
        myStds[i] = s;
    }
    vector<vector<double> > ret;
    ret.push_back(myT);
    ret.push_back(myMeans);
    ret.push_back(myStds);
    return ret;
}
