#include "rydberg.hpp"

default_random_engine generator (time(NULL));

OneSite::OneSite() {
    T = 0.;
    Tmax = 500.;
    dt = .1;
    Gamma = .3;
    b = .005;
    po = .01;
    p = po;
    h = 0.;
}


OneSite::OneSite(double myT, double myDT, double myGamma, double myB, double myPo) {
    T = 0.;
    Tmax = myT;
    dt = myDT;
    Gamma = myGamma;
    b = myB;
    po = myPo;
    p = po;
    h = 0.;
}

OneSite::~OneSite() {}

vector<double> OneSite::getConfig() {
    vector<double> my_config;
    my_config.push_back(dt);
    my_config.push_back(Gamma);
    my_config.push_back(b);
    my_config.push_back(po);
    return my_config;
}

void OneSite::reset() {
    T = 0.;
    p = po;
    h = 0.;
}

double OneSite::update() {
    double beta = 1. - Gamma - Gamma * b * h;
    double lambda = 2.*beta / (exp(beta*dt) - 1.);
    poisson_distribution<int> PoissonDist(lambda * p * exp(beta*dt));
    gamma_distribution<double> GammaDist(PoissonDist(generator), 1.0);
    double pStar = GammaDist(generator) / lambda;
    pStar /= (1+2*dt*pStar);
    h += .5*dt*(pStar + p);
    p = pStar;
    T += dt;
    return pStar;
}

vector< vector<double> > OneSite::simulation() {
    int N = int(Tmax / dt);
    vector<double> Tlist(N, 0.);
    vector<double> Plist(N, 0.);
    vector<double> Hlist(N, 0.);
    reset();

    for (int i = 0; i <= N; i++) {
        Tlist[i]  = T;
        Plist[i] = p;
        Hlist[i] = h;
        if (p == 0.) {
            T += dt;
        }
        else {
            update();
        }
    }
    vector< vector<double> > res;
    res.push_back(Tlist);
    res.push_back(Plist);
    res.push_back(Hlist);
    return res;
}

vector< vector<double> > OneSite::simulation(int N) {
    vector< vector<double> > res;
    if (N <= 0) {
        return res;
    }
    res = simulation();
    for (int i = 1; i < N; i++) {
        vector< vector<double> > tmp = simulation();
        res.push_back(tmp[1]);
        res.push_back(tmp[2]);
    }
    // return res;
    int Nsteps = int(Tmax / dt);
    vector<double> means(Nsteps, 0.);
    vector<double> stdevs(Nsteps, 0.);
    vector<double> Hmeans(Nsteps, 0.);
    for (int i = 0; i < Nsteps; i++) {
        double m = 0.;
        double h = 0.;
        for (int j = 0; j < N; j++) {
            m += res[2*j+1][i];
            h += res[2*j + 2][i];
        }
        m /= double(N);
        h /= double(N);
        double stdev = 0.;
        for (int j = 0; j < N; j++) {
            stdev += pow(res[2*j+1][i] - m, 2);
        }
        stdev = sqrt(stdev / (N-1.));
        means[i] = m;
        stdevs[i] = stdev;
        Hmeans[i] = h;
    }
    vector<vector<double> > toRet;
    toRet.push_back(res[0]);
    toRet.push_back(means);
    toRet.push_back(stdevs);
    toRet.push_back(Hmeans);
    return toRet;
}

//
// vector<double> mean_stdev(vector<double> dat) {
//     vector<double> res;
//     double sum = accumulate(dat.begin(), dat.end(), 0.0);
//     double m = sum / dat.size();
//     double accum = 0.0;
//     for_each (dat.begin(), dat.end(), [&](const double d) {accum += (d-m)*(d-m);});
//     double stdev = sqrt(accum / (dat.size()-1));
//     res.push_back(m);
//     res.push_back(stdev);
//     return res;
// }

Lattice1D::Lattice1D() {
    T = 0.;
    Tmax = 200;
    dt = .2;
    dx = 1.;
    Gamma = .7;
    b = .01;
    po = .01;
    Nsites = 5e1;
    p = vector<double>(Nsites, po);
    h = vector<double>(Nsites, 0.);
}

Lattice1D::Lattice1D(double my_Tmax, double my_dt, double my_Gamma, double my_b, double my_po, int my_Nsites, double my_dx) {
    T = 0.;
    Tmax = my_Tmax;
    dt = my_dt;
    Gamma = my_Gamma;
    b = my_b;
    po = my_po;
    Nsites = my_Nsites;
    dx = my_dx;
    p = vector<double>(Nsites, po);
    h = vector<double>(Nsites, 0.);
}

Lattice1D::~Lattice1D() {}

vector<double> Lattice1D::getConfig() {
    vector<double> my_config;
    my_config.push_back(dt);
    my_config.push_back(Gamma);
    my_config.push_back(b);
    my_config.push_back(po);
    my_config.push_back(dx);
    my_config.push_back(Nsites);
    return my_config;
}

void Lattice1D::reset() {
    T = 0.;
    for (int i = 0; i < Nsites; i ++) {
        p[i] = po;
        h[i] = 0.;
    }
}

vector<double> Lattice1D::update() {
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
        alpha /= (dx*dx);
        double beta = 1.-Gamma-Gamma*b*h[i] - 2./(dx*dx);
        double scale = 1 / ((1.-Gamma) * (1.-Gamma));
        scale = 10000.;
        double lambda = 2.*beta * scale / (exp(beta*dt) - 1.);
        poisson_distribution<int> PoissonDist(lambda * p[i] * exp(beta*dt));
        gamma_distribution<double> GammaDist(2.*alpha + PoissonDist(generator), 1.0);
        double pStar = GammaDist(generator) / lambda;
        pStar /= (1. + pStar * dt * 2.);
        new_p[i] = pStar;
        new_h[i] = h[i] + .5 * dt * (pStar + p[i]);
    }
    p = new_p;
    h = new_h;
    T += dt;
    return new_p;
}

bool Lattice1D::is_zero() {
    for (int i = 0; i < Nsites; i++) {
        if (p[i] != 0.) {
            return false;
        }
    }
    return true;
}

vector< vector<double> > Lattice1D::simulation() {
    int N = int(Tmax / dt);
    vector<double> Tlist(N, 0.);
    vector<vector<double> > dataList;
    reset();
    for (int i = 0; i < N; i++) {
        Tlist[i] = T;
        dataList.push_back(p);
        dataList.push_back(h);
        bool checkZero = is_zero();
        if (checkZero) {
            T += dt;
        }
        else {
                update();
        }
    }
    dataList.insert(dataList.begin(), Tlist);
    return dataList;
}

vector<vector<vector<double> > > Lattice1D::simulation(int N) {
    vector<vector<vector<double> > > res;
    if (N <= 0) {
        return res;
    }
    vector<vector<double> > tmp = simulation();
    vector<vector<double> > Tvec;
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
    vector<vector<double> > allM;
    vector<vector<double> > allS;
    vector<vector<double> > allH;
    for (int i = 0; i < NTime; i++) {
        vector<double> mVec (Nsites, 0.);
        vector<double> sVec (Nsites, 0.);
        vector<double> hVec (Nsites, 0.);
        for (int j = 0; j < Nsites; j++) {
            double m = 0.;
            double mh = 0.;
            for (int k = 0; k < N; k++) {
                m += res[k+1][2*i][j];
                mh += res[k+1][2*i+1][j];
            }
            m /= double(N);
            mh /= double(N);
            double s = 0.;
            for (int k = 0; k < N; k++) {
                s += pow(res[k+1][2*i][j] - m, 2);
            }
            if (N == 1) {
                s = 0;
            }
            else {
                s = sqrt(s / (N-1));
            }
            mVec[j] = m;
            sVec[j] = s;
            hVec[j] =  mh; //Gamma*(1.+b*mh) - 1.;
        }
        allM.push_back(mVec);
        allS.push_back(sVec);
        allH.push_back(hVec);
    }
    vector<vector<vector<double> > > toRet;
    toRet.push_back(Tvec);
    toRet.push_back(allM);
    toRet.push_back(allS);
    toRet.push_back(allH);
    return toRet;
}

vector<vector<double> > Lattice1D::sim_avg(int N) {
    vector<vector<vector<double> > > res = simulation(N);
    int Nsteps = int(Tmax / dt);
    vector<double> myT(Nsteps, 0.);
    vector<double> myM(Nsteps, 0.);
    vector<double> myS(Nsteps, 0.);
    vector<double> myGap(Nsteps, 0.);
    for (int i = 0; i < Nsteps; i++) {
        double m = 0.;
        double s = 0.;
        double gap = 0.;
        int n = 0.;
        for (int j = 0; j < Nsites; j++) {
            if (res[1][i][j] > 1e-20) { //
                m += res[1][i][j];
                s += pow(res[2][i][j], 2);
                gap += res[3][i][j];
                n++; //
            }
        }
        m /= n; // Nsites
        s = sqrt(s) / n; // Nsites
        gap /= n; // Nsites
        if (N == 1) {
            s = 0.;
            for (int j = 0; j < Nsites; j++) {
                s += pow(res[1][i][j] - m, 2);
            }
            s = sqrt(s / (Nsites - 1));
        }
        myM[i] = m;
        myS[i] = s;
        myT[i] = res[0][0][i];
        myGap[i] = gap;
    }
    vector<vector<double> > ret;
    ret.push_back(myT);
    ret.push_back(myM);
    ret.push_back(myS);
    ret.push_back(myGap);
    return ret;
}

Lattice2D::Lattice2D() {
    T = 0.;
    Tmax = 200;
    dt = .1;
    dx = 1.;
    Gamma = .4;
    b = .1;
    po = .01;
    L = 100;
    Nsites = L * L;
    p = vector<double>(Nsites, po);
    h = vector<double>(Nsites, 0.);
}

Lattice2D::Lattice2D(double myT, double myDT, double myGamma, double myB, double myPo, int myL, double myDx) {
    T = 0.;
    Tmax = myT;
    dt = myDT;
    Gamma = myGamma;
    b = myB;
    po = myPo;
    L = myL;
    Nsites = L * L * L;
    dx = myDx;
    p = vector<double>(Nsites, po);
    h = vector<double>(Nsites, 0.);
}

Lattice2D::~Lattice2D() {}

vector<double> Lattice2D::getConfig() {
    vector<double> my_config;
    my_config.push_back(dt);
    my_config.push_back(Gamma);
    my_config.push_back(b);
    my_config.push_back(po);
    my_config.push_back(dx);
    my_config.push_back(L);
    return my_config;
}

void Lattice2D::reset() {
    T = 0.;
    p = vector<double>(Nsites,po);
    h = vector<double>(Nsites,0.);
}

vector<double> Lattice2D::update() {
    vector<double> new_p (Nsites, 0.);
    vector<double> new_h (Nsites, 0.);

    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            int ind = i + L*j;
            double alpha = 0.;
            if (i != 0) {
                alpha += p[ind -1];
            }
            if (i != L-1) {
                alpha += p[ind + 1];
            }
            if (j != 0) {
                alpha += p[ind - L];
            }
            if (j != L-1) {
                alpha += p[ind + L];
            }
            alpha /= (dx*dx);
            double beta = 1.-Gamma-Gamma*b*h[ind] - 4./(dx*dx);
            double lambda = 2*beta / (Gamma*Gamma * (exp(beta*dt) - 1));
            poisson_distribution<int> PoissonDist(lambda * p[ind] * exp(beta*dt));
            gamma_distribution<double> GammaDist(2*alpha + PoissonDist(generator), 1.0);
            double pStar = GammaDist(generator) / lambda;
            pStar /= (1 + pStar * dt * 2);
            new_p[ind] = pStar;
            new_h[ind] = h[ind] + .5*dt*(pStar+p[ind]);
        }
    }
    p = new_p;
    h = new_h;
    T += dt;
    return new_p;
}

bool Lattice2D::is_zero() {
    for (int i = 0; i < Nsites; i++) {
        if (p[i] != 0.) {
            return false;
        }
    }
    return true;
}

vector< vector<double> > Lattice2D::simulation() {
    int N = int(Tmax / dt);
    vector<double> Tlist(N, 0.);
    vector<vector<double> > Plist;
    reset();

    for (int i = 0; i < N; i++) {
        Tlist[i] = T;
        Plist.push_back(p);
        bool checkZero = is_zero();
        if (checkZero) {
            T += dt;
        }
        else {
            update();
        }
    }
    Plist.insert(Plist.begin(), Tlist);
    return Plist;
}

vector<vector<vector<double> > > Lattice2D::simulation(int N) {
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
            if (N == 1) {
                s = 0;
            }
            else {
                s = sqrt(s / (N-1));
            }
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

vector<vector<double> > Lattice2D::sim_avg(int N) {
    vector<vector<vector<double> > > res = simulation(N);
    int Nsteps = int(Tmax / dt);
    vector<double> myT(Nsteps, 0.);
    vector<double> myM(Nsteps, 0.);
    vector<double> myS(Nsteps, 0.);

    for (int i = 0; i < Nsteps; i++) {
        double m = 0.;
        double s = 0.;
        for (int j = 0; j < Nsites; j++) {
            m += res[1][i][j];
            s += pow(res[2][i][j], 2);
        }
        m /= Nsites;
        s = sqrt(s) / Nsites;
        if (N == 1) {
            s = 0.;
            for (int j = 0; j < Nsites; j++) {
                s += pow(res[1][i][j] - m, 2);
            }
            s = sqrt(s / (Nsites - 1));
        }
        myM[i] = m;
        myS[i] = s;
        myT[i] = res[0][0][i];
    }
    vector<vector<double> > ret;
    ret.push_back(myT);
    ret.push_back(myM);
    ret.push_back(myS);
    return ret;
}


Lattice3D::Lattice3D() {
    T = 0.;
    Tmax = 200;
    dt = .1;
    dx = 1000.;
    Gamma = .3;
    b = .01;
    po = .01;
    L = 30;
    Nsites = L * L * L;
    p = vector<double>(Nsites, po);
    h = vector<double>(Nsites, 0.);
}

Lattice3D::Lattice3D(double myT, double myDT, double myGamma, double myB, double myPo, int myL, double myDx) {
    T = 0.;
    Tmax = myT;
    dt = myDT;
    Gamma = myGamma;
    b = myB;
    po = myPo;
    L = myL;
    Nsites = L * L * L;
    dx = myDx;
    p = vector<double>(Nsites, po);
    h = vector<double>(Nsites, 0.);
}

Lattice3D::~Lattice3D() {}

vector<double> Lattice3D::getConfig() {
    vector<double> my_config;
    my_config.push_back(dt);
    my_config.push_back(Gamma);
    my_config.push_back(b);
    my_config.push_back(po);
    my_config.push_back(dx);
    my_config.push_back(L);
    return my_config;
}

void Lattice3D::reset() {
    T = 0.;
    p = vector<double>(Nsites,po);
    h = vector<double>(Nsites,0.);
}

vector<double> Lattice3D::update() {
    vector<double> new_p (Nsites, 0.);
    vector<double> new_h (Nsites, 0.);

    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            for (int k = 0; k < L; k++) {
                int ind = i + L*j + L*L*k;
                double alpha = 0.;
                if (i != 0) {
                    alpha += p[ind - 1];
                }
                if (i != L - 1) {
                    alpha += p[ind + 1];
                }
                if (j != 0) {
                    alpha += p[ind - L];
                }
                if (j != L - 1) {
                    alpha += p[ind + L];
                }
                if (k != 0) {
                    alpha += p[ind - L*L];
                }
                if (k != L - 1) {
                    alpha += p[ind + L*L];
                }
                alpha /= (dx * dx);
                double beta = 1. - Gamma - Gamma*b*h[ind] - 6./(dx * dx);
                double lambda = 2.*beta / (Gamma * Gamma * (exp(beta*dt) - 1));
                poisson_distribution<int> PoissonDist(lambda * p[ind] * exp(beta*dt));
                gamma_distribution<double> GammaDist(2*alpha + PoissonDist(generator), 1.0);
                double pStar = GammaDist(generator) / lambda;
                pStar /= (1 + pStar * dt * 2);
                new_p[ind] = pStar;
                new_h[ind] = h[ind] + .5*dt*(pStar + p[ind]);
            }
        }
    }
    p = new_p;
    h = new_h;
    T += dt;
    return new_p;
}

bool Lattice3D::is_zero() {
    for (int i = 0; i < Nsites; i++) {
        if (p[i] != 0.) {
            return false;
        }
    }
    return true;
}

vector< vector<double> > Lattice3D::simulation() {
    int N = int(Tmax / dt);
    vector<double> Tlist(N, 0.);
    vector<vector<double> > Plist;
    reset();

    for (int i = 0; i < N; i++) {
        Tlist[i] = T;
        Plist.push_back(p);
        bool checkZero = is_zero();
        if (checkZero) {
            T += dt;
        }
        else {
            update();
        }
    }
    Plist.insert(Plist.begin(), Tlist);
    return Plist;
}

vector<vector<vector<double> > > Lattice3D::simulation(int N) {
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
            if (N == 1) {
                s = 0;
            }
            else {
                s = sqrt(s / (N-1));
            }
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

vector<vector<double> > Lattice3D::sim_avg(int N) {
    vector<vector<vector<double> > > res = simulation(N);
    int Nsteps = int(Tmax / dt);
    vector<double> myT(Nsteps, 0.);
    vector<double> myM(Nsteps, 0.);
    vector<double> myS(Nsteps, 0.);

    for (int i = 0; i < Nsteps; i++) {
        double m = 0.;
        double s = 0.;
        for (int j = 0; j < Nsites; j++) {
            m += res[1][i][j];
            s += pow(res[2][i][j], 2);
        }
        m /= Nsites;
        s = sqrt(s) / Nsites;
        if (N == 1) {
            s = 0.;
            for (int j = 0; j < Nsites; j++) {
                s += pow(res[1][i][j] - m, 2);
            }
            s = sqrt(s / (Nsites - 1));
        }
        myM[i] = m;
        myS[i] = s;
        myT[i] = res[0][0][i];
    }
    vector<vector<double> > ret;
    ret.push_back(myT);
    ret.push_back(myM);
    ret.push_back(myS);
    return ret;
}
