#include "rydberg.hpp"

default_random_engine generator (time(NULL));

static int seed1= time(NULL); //437;
static int seed2= time(NULL); //4357;
static boost::mt19937 rnd_gen1(seed1);
static boost::mt19937 rnd_gen2(seed2);

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
    double pShape = lambda * p * exp(beta * dt);
    double pOut = poiss_rand(pShape, rnd_gen2);
    // poisson_distribution<int> PoissonDist(lambda * p * exp(beta*dt));
    double pStar = gamma_rand(pOut, 1.0, rnd_gen1) / lambda;
    // gamma_distribution<double> GammaDist(PoissonDist(generator), 1.0);
    // double pStar = GammaDist(generator) / lambda;
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

Lattice1D::Lattice1D() {
    T = 0.;
    Tmax = 200;
    dt = .1;
    dx = 2.;
    Gamma = -.63;
    b = 0.;
    po = 1e-1;
    Nsites = 1e5;
    p = vector<double>(Nsites, po);
    kappa = 1e-3;
    // p = vector<double>(Nsites, 0.);
    // int i = int(po * Nsites);
    // for (int j = 0; j < i; j++) {
    //     p[j] = .99;
    // }
    // random_shuffle(p.begin(), p.end());
    h = vector<double>(Nsites, 1.5);
    pMean = po;
}

Lattice1D::Lattice1D(double my_Tmax, double my_dt, double my_Gamma, double my_b, double my_po, int my_Nsites, double my_dx, double myKappa) {
    T = 0.;
    Tmax = my_Tmax;
    dt = my_dt;
    Gamma = my_Gamma;
    b = my_b;
    po = my_po;
    Nsites = my_Nsites;
    dx = my_dx;
    p = vector<double>(Nsites, po);
    kappa = myKappa;
    // p = vector<double>(Nsites, 0.);
    // for (int i = 0; i < int(Nsites / 100); i++) {
    //     p[i] = po * 100.;
    // }
    // random_shuffle(p.begin(), p.end());
    h = vector<double>(Nsites, 0.);
    for (int i = int(Nsites/3); i < int(2*Nsites/3);i++) {
        h[i] = 1.5;
        if (i % 2 == 0) {
            h[i] += pow(-1.,int(i/2))*.1;
        }
    }
    pMean = po;
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
    // int j = int(po * Nsites);
    for (int i = 0; i < Nsites; i ++) {
        p[i] = po;
        h[i] = 0.;
        // if (i < j) {
        //     p[i] = 0.99;
        // }
    }
    random_shuffle(p.begin(), p.end());
    pMean = po;
}

vector<double> Lattice1D::getP() {
    return p;
}

double Lattice1D::getPmean() {
    return pMean;
}

void Lattice1D::exciteSite(int site, double seed) {
    p[site] = seed;
}

vector<double> Lattice1D::getH() {
    return h;
}

vector<double> Lattice1D::update() {
    double sigma = sqrt(2.);
    double Dtilde = 3.;
    vector<double> new_p (Nsites, 0.);
    vector<double> new_h (Nsites, 0.);
    pMean = 0.;
    for (int i = 0; i < Nsites; i++) {
        double alpha = 0.;
        if (i != 0) {
            alpha += p[i-1];
        }
        else {
            alpha += p[Nsites - 1];
        }
        if (i != Nsites - 1) {
            alpha += p[i+1];
        }
        else {
            alpha += p[0];
        }
        alpha /= (dx*dx);
        double beta = 1. + (kappa * T) -Gamma-Gamma*b*h[i] - 2./(dx*dx);
        double lambda = 2.*beta / (exp(beta*dt) - 1.);
        lambda /= (sigma * sigma);
        double pShape = lambda * p[i] * exp(beta * dt);
        double pOut = poiss_rand(pShape, rnd_gen2);
        double gShape = pOut + 2.*alpha / (sigma * sigma);
        double gOut = gamma_rand(gShape, 1.0, rnd_gen1) / lambda;
        double pStar = gOut;
        // pStar /= (1. + pStar * dt * 2.);
        pStar /= (1. + pStar * dt);
        new_p[i] = pStar;
        pMean += pStar;
        double hStar = h[i] + 0.5*dt*(pStar + p[i]);
        if (i == 0) {
            hStar += dt * Dtilde * h[Nsites - 1] / (dx * dx);
            // hStar += 0;
        }
        else {
            hStar += dt * Dtilde * h[i-1] / (dx * dx);
        }
        hStar += dt * Dtilde * h[abs((i+1) % Nsites)] / (dx*dx);
        // if (i == Nsites - 1) {
        //     hStar += 0;
        // }
        // else {
        //     hStar += dt * Dtilde * h[i+1] / (dx*dx);
        // }
        hStar -= 2. * dt * Dtilde*h[i]/(dx*dx);
        new_h[i] = hStar;
    }
    pMean /= double(Nsites);
    p = new_p;
    h = new_h;
    T += dt;
    return new_p;
}

vector<double> Lattice1D::tau_update(double tau) {
    double sigma = sqrt(2.);
    double Dtilde = 3.;
    vector<double> new_p(Nsites, 0.);
    vector<double> new_h(Nsites, 0.);
    pMean = 0.;
    for (int i = 0; i < Nsites; i++) {
        double alpha = 0.;
        if (i != 0) {
            alpha += p[i-1];
        }
        else {
            alpha += p[Nsites - 1];
        }
        if (i != Nsites - 1) {
            alpha += p[i+1];
        }
        else {
            alpha += p[0];
        }
        alpha /= dx*dx;
        double beta = 1. + (kappa * T) - Gamma-Gamma*b*h[i] - 2./(dx*dx);
        alpha -= beta*tau;
        double lambda = 2.*beta / (exp(beta*dt) - 1.);
        lambda /= (sigma*sigma);
        double u = p[i] + tau;
        double uShape = lambda * u * exp(beta * dt);
        double uOut = poiss_rand(uShape, rnd_gen2);
        double gShape = uOut + 2.*alpha / (sigma*sigma);
        double gOut = gamma_rand(gShape, 1.0, rnd_gen1) / lambda;
        double uStar = gOut;
        if (uStar < tau) {
            uStar = tau;
        }
        double pStar = uStar - tau;
        pStar /= (1. + pStar * dt);
        new_p[i] = pStar;
	pMean += pStar;
        double hStar = h[i] + 0.5*dt*(pStar + p[i]);
        if (i == 0) {
            hStar += dt*Dtilde*h[Nsites - 1] / (dx*dx);
        }
        else {
            hStar += dt*Dtilde * h[i-1] / (dx*dx);
        }
        hStar += dt*Dtilde*h[abs((i+1)%Nsites)] / (dx*dx);
        hStar -= 2.*dt*Dtilde*h[i]/(dx*dx);
        new_h[i] = hStar;
    }
    p = new_p;
    h = new_h;
    pMean /= double(Nsites);
    T += dt;
    return new_p;
}

void Lattice1D::set_b(double myB) {
    b = myB;
}

vector<double> potential_grad(vector<double> x) {
    int N = (int)x.size();
    vector<double> res(N,0.);
    vector<double> F(N,0.);
    float mu = (float)(float(N)/2.);
    float sigma = mu/2;
    for (int i = 0; i < N; i++) {
        float pos = (float)i - mu;

        F[i] = (pos / (sigma*sigma)) * exp(-(pos*pos) / (2*sigma*sigma));
    }
    for (int i = 1; i < N-1; i++) {
        res[i] = x[i]*F[i] - x[i+1]*F[i+1];
    }
    res[0] = 2*res[1] - res[2];
    res[N-1] = 2*res[N-2] - res[N-3];
    return res;
}

vector<double> potential_drive(int N) {
    vector<double> res(N,0.);
    float mu = (float)(float(N)/2.);
    float sigma = mu/2;
    for (int i = 0; i < N; i++) {
        float pos = (float)i - mu;
        res[i] = exp(-(pos*pos)/(2*sigma*sigma));
    }
    return res;
}

vector<double> Lattice1D::trap_update(double tau, double depth) {
    double sigma = sqrt(2.);
    double Dtilde = kappa;
    double Dp = 1.;
    vector<double> new_p(Nsites, 0.);
    vector<double> new_h(Nsites, 0.);
    pMean = 0.;
    vector<double> forceTerm = potential_grad(h);
    vector<double> driveTerm = potential_drive(Nsites);
    for (int i = 0; i < Nsites; i++) {
        double alpha = 0.;
        // double beta = 1. - Gamma - h[i];
        double beta = .5*h[i] - 2.*tau - Gamma;
        if (i != 0) {
            alpha += p[i-1];
            beta -= Dp/(dx * dx);
        }
        if (i != Nsites - 1) {
            alpha += p[i+1];
            beta -= Dp/(dx * dx);
        }
        alpha *= Dp/(dx*dx);
        alpha += tau*h[i];
        alpha -= beta * tau;
        double lambda = 2.*beta / (exp(beta * dt) - 1.);
        lambda /= (sigma * sigma);
        double u = p[i] + tau;
        double uShape = lambda * u * exp(beta * dt);
        double uOut = poiss_rand(uShape, rnd_gen2);
        double gShape = uOut + 2.*alpha / (sigma * sigma);
        double gOut = gamma_rand(gShape, 1.0, rnd_gen1) / lambda;
        double uStar = gOut;
        if (uStar < tau) {
            uStar = tau;
        }
        double pStar = uStar - tau;
        pStar /= (1. + pStar * dt);
        new_p[i] = pStar;
        pMean += pStar;
        double hStar = h[i] - b*0.5*dt*(pStar + p[i]);
        if (i != 0) {
            hStar += dt*Dtilde*h[i-1]/(dx*dx);
            hStar -= dt*Dtilde*h[i]/(dx*dx);
        }
        if (i != Nsites - 1) {
            hStar += dt*Dtilde*h[i+1]/(dx*dx);
            hStar -= dt*Dtilde*h[i]/(dx*dx);
        }
        new_h[i] = hStar - depth * forceTerm[i]; //  + dt*kappa*driveTerm[i];
    }
    p = new_p;
    h = new_h;
    pMean /= double(Nsites);
    T += dt;
    return new_p;
}

vector<double> Lattice1D::euler_update() {
    vector<double> new_p(Nsites, 0.);
    vector<double> new_h(Nsites, 0.);
    for (int i = 0; i < Nsites; i++) {
        double myLapl = -2*p[i];
        if (i != 0) {
            myLapl += p[i-1];
        }
        if (i != Nsites - 1) {
            myLapl += p[i+1];
        }
        myLapl /= (dx * dx);
        double dp = myLapl + (1 - Gamma*(1+b*h[i]))*p[i] - 2*(p[i]*p[i]);
        new_p[i] = p[i] + dt*dp;
        new_h[i] = h[i] + 0.5*(p[i] + new_p[i])*dt;
    }
    p = new_p;
    h = new_h;
    T += dt;
    return new_p;
}

bool Lattice1D::is_zero() {
    return (pMean == 0.);
    // for (int i = 0; i < Nsites; i++) {
    //     if (p[i] != 0.) {
    //         return false;
    //     }
    // }
    // return true;
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
            hVec[j] = mh; //Gamma*(1.+b*mh) - 1.;
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
    kappa = 5e-5;
    L = 100;
    Nsites = L * L;
    p = vector<double>(Nsites, po);
    h = vector<double>(Nsites, 0.);
    pMean = po;

}

Lattice2D::Lattice2D(double myT, double myDT, double myGamma, double myB, double myPo, int myL, double myDx, double myKappa) {
    T = 0.;
    Tmax = myT;
    dt = myDT;
    Gamma = myGamma;
    b = myB;
    po = myPo;
    kappa = myKappa;
    L = myL;
    Nsites = L*L;
    dx = myDx;
    p = vector<double>(Nsites, po);
    h = vector<double>(Nsites, 0.);
    for (int x = int(L/3); x < int(2*L/3); x++) {
        for (int y = int(L/3); y < int(2*L/3); y++) {
            h[ind(x,y)] = 1.5;
        }
    }
    pMean = po;
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

void Lattice2D::set_b(double myB) {
    b = myB;
}

void Lattice2D::reset() {
    T = 0.;
    p = vector<double>(Nsites,po);
    h = vector<double>(Nsites,0.);
    pMean = po;
}

vector<double> Lattice2D::getP() {
    return p;
}

double Lattice2D::getPmean() {
    return pMean;
}

double Lattice2D::getHmean() {
    return hMean;
}

void Lattice2D::exciteSite(int site, double seed) {
    p[site] = seed;
}

vector<double> Lattice2D::getH() {
    return h;
}

int Lattice2D::ind(int i, int j) {
    i = (i % L);
    j = (j % L);
    if (i < 0) {
        i += L;
    }
    if (j < 0) {
        j += L;
    }
    int myRes = L*i+j;
    assert(myRes >= 0);
    assert(myRes < L*L);
    return myRes;
}

int ind_2D(int i, int j, int L) {
    i = (i % L);
    j = (j % L);
    if (i < 0) {
        i += L;
    }
    if (j < 0) {
        j += L;
    }
    int myRes = L*i + j;
    assert(myRes >= 0);
    assert(myRes < L * L);
    return myRes;
}

vector<double> Lattice2D::update() {
    vector<double> new_p (Nsites, 0.);
    vector<double> new_h (Nsites, 0.);
    double sigma = sqrt(2.);
    double Dtilde = 3.;
    pMean = 0.;
    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            double alpha = 0.;
            alpha += p[ind(i+1,j)];
	    alpha += p[ind(i-1,j)];
	    alpha += p[ind(i,j+1)];
	    alpha += p[ind(i,j-1)];
            alpha /= (dx*dx);
            double beta = 1.+(kappa*T)-Gamma-Gamma*b*h[ind(i,j)] - 4./(dx*dx);
            double lambda = 2*beta / (sigma*sigma * (exp(beta*dt) - 1));
	    double pShape = lambda * p[ind(i,j)] * exp(beta*dt);
	    double pOut = poiss_rand(pShape, rnd_gen2);
            // poisson_distribution<int> PoissonDist(lambda * p[ind(i,j)] * exp(beta*dt));
            // gamma_distribution<double> GammaDist(2*alpha + PoissonDist(generator), 1.0);
            // double pStar = GammaDist(generator) / lambda;
            double gShape = pOut + 2.*alpha/(sigma*sigma);
	    double pStar = gamma_rand(gShape, 1.0, rnd_gen1) / lambda;
	    pStar /= (1 + pStar * dt);
            // cout << alpha << " " << beta << " " << lambda << " " << pStar << "\n";
            new_p[ind(i,j)] = pStar;
	    pMean += pStar;
	    double hStar = h[ind(i,j)] + 0.5*dt*(pStar + p[ind(i,j)]);
            hStar += dt*Dtilde*h[ind(i+1,j)]/(dx*dx);
            hStar += dt*Dtilde*h[ind(i-1,j)]/(dx*dx);
	    hStar += dt*Dtilde*h[ind(i,j+1)]/(dx*dx);
	    hStar += dt*Dtilde*h[ind(i,j-1)]/(dx*dx);
    	    hStar -= 4.*dt*Dtilde*h[ind(i,j)]/(dx*dx);
            new_h[ind(i,j)] = hStar;
        }
    }
    pMean /= double(Nsites);
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

vector<double> potential_grad_2D(vector<double> myX, int N, int L) {
    vector<double> res(N,0.);
    vector<double> Fx(N,0.);
    vector<double> Fy(N,0.);
    float mu = (float)(float(L) / 2.);
    float sigma = mu / 2.;
    for (int i = 0; i < N; i++) {
        float x = (float)((int) i/L);
        float y = (float)(i % L);
        float posx = x - mu;
        float posy = y - mu;
        float posr = sqrt(posx * posx + posy * posy);
        float expFac = exp(-(posr * posr) / (2*sigma * sigma));
        int ix = (int)x;
        int iy = (int)y;
        float dxh = 0.;
        float dyh = 0.;
        if (ix == 0) {
            dxh = myX[ind_2D(ix+1,iy,L)] - myX[i];
        } else if (ix == (L-1)) {
            dxh = myX[i] - myX[ind_2D(ix-1,iy,L)];
        } else {
            dxh = 0.5*(myX[ind_2D(ix+1,iy,L)] - myX[ind_2D(ix-1,iy,L)]);
        }
        if (iy == 0) {
            dyh = myX[ind_2D(ix,iy+1,L)] - myX[i];
        } else if (iy == (L-1)) {
            dyh = myX[i] - myX[ind_2D(ix,iy-1,L)];
        } else {
            dyh = 0.5*(myX[ind_2D(ix,iy+1,L)] - myX[ind_2D(ix,iy-1,L)]);
        }
        float tmp = dxh*(-posx) + dyh*(-posy) + myX[i]*(-2 + (posr * posr)/(sigma * sigma));
        res[i] = tmp * expFac / (sigma*sigma);
    }
    return res;
}

vector<double> Lattice2D::seeding_update(double tau, double depth) {
    double sigma = sqrt(2.);
    double Dtilde = kappa;
    double Dp = .2;
    vector<double> new_p(Nsites, 0.);
    vector<double> new_h(Nsites, 0.);
    double old_pMean = pMean;
    pMean = 0.;
    vector<double> forceTerm = potential_grad_2D(h,Nsites,L);
    for (int i = 0; i < Nsites; i++) {
        double alpha = 0.;
        double beta = .5*h[i] - Gamma; // - 2.*tau - Gamma;
        int x = int(i / L);
        int y = i % L;
        if (x != 0) {
            alpha += p[ind(x-1, y)];
            beta -= Dp/(dx * dx);
        }
        if (x != L - 1) {
            alpha += p[ind(x+1,y)];
            beta -= Dp/(dx * dx);
        }
        if (y != 0) {
            alpha += p[ind(x, y-1)];
            beta -= Dp/(dx * dx);
        }
        if (y != L - 1) {
            alpha += p[ind(x,y+1)];
            beta -= Dp/(dx * dx);
        }
        alpha *= Dp/(dx * dx);
        // alpha += tau*h[i];
        // alpha -= beta * tau;
        double lambda = 2.*beta / (exp(beta * dt) - 1.);
        lambda /= (sigma * sigma);
        double pShape = lambda * p[i] * exp(beta * dt);
        double pOut = poiss_rand(pShape, rnd_gen2);
        double gShape = pOut + 2.*alpha / (sigma * sigma);
        double gOut = gamma_rand(gShape, 1.0, rnd_gen1) / lambda;
        double pStar = gOut;
        if (pStar < 0) {
            pStar = 0.;
        }
        pStar /= (1. + pStar * dt);
        new_p[i] = pStar;
        if (new_p[i] > 0) {
            avalanche_set.insert(i);
            instantaneous_set.insert(i);
        }
        pMean += pStar;
        double hStar = h[i] - b * 0.5 * dt * (pStar + p[i]);
        if (x != 0) {
            hStar += dt * Dtilde * h[ind(x-1,y)] / (dx * dx);
            hStar -= dt * Dtilde * h[i] / (dx * dx);
        }
        if (x != L-1) {
            hStar += dt * Dtilde * h[ind(x+1,y)] / (dx * dx);
            hStar -= dt * Dtilde * h[i] / (dx * dx);
        }
        if (y != 0) {
            hStar += dt * Dtilde * h[ind(x,y-1)] / (dx * dx);
            hStar -= dt * Dtilde * h[i] / (dx * dx);
        }
        if (y != L-1) {
            hStar += dt * Dtilde * h[ind(x,y+1)] / (dx * dx);
            hStar -= dt * Dtilde * h[i] / (dx * dx);
        }
        new_h[i] = hStar - depth * forceTerm[i];
    }
    // If addSeed flag is true, then randomly seed an excitation
    if (pMean == 0.) {
        double r = double(rand()) / double(RAND_MAX);
        if (avalanche_set.size() > 0) {
            all_avalanche_sets.push_back(avalanche_set.size());
            all_avalanche_sizes.push_back(avalanche_sizes);
            start_times.push_back(avalanche_start_time);
            avalanche_lengths.push_back(T - avalanche_start_time);
            avalanche_set.clear();
            avalanche_sizes.clear();
            instantaneous_set.clear();
        }
        if (r < (tau * float(Nsites))) {
            int addSite = int(r / tau);
            new_p[addSite] = .1;
            cout << "T: " << T << " Adding seed at site " << addSite << " given pMean " << pMean << "\n";
            pMean += .1;
            avalanche_start_time = T;
        }
    }
    else {
        avalanche_sizes.push_back(instantaneous_set.size());
        instantaneous_set.clear();
    }
    p = new_p;
    h = new_h;
    pMean /= double(Nsites);
    T += dt;
    return new_p;
}

vector<double> Lattice2D::trap_update(double tau, double depth) {
    double sigma = sqrt(2.);
    double Dtilde = kappa;
    double Dp = .2;
    vector<double> new_p(Nsites, 0.);
    vector<double> new_h(Nsites, 0.);
    // p[ind(i,j)], use L
    pMean = 0.;
    hMean = 0.;
    vector<double> forceTerm = potential_grad_2D(h, Nsites, L);
    for (int i = 0; i < Nsites; i++) {
        double alpha = 0.;
        double beta = .5*h[i] - 2.*tau - Gamma; // check dimensionality here
        int x = int(i / L);
        int y = i % L;
        if (x != 0) {
            alpha += p[ind(x-1, y)];
            beta -= Dp/(dx * dx);
        }
        if (x != L - 1) {
            alpha += p[ind(x+1,y)];
            beta -= Dp/(dx * dx);
        }
        if (y != 0) {
            alpha += p[ind(x, y-1)];
            beta -= Dp/(dx * dx);
        }
        if (y != L - 1) {
            alpha += p[ind(x,y+1)];
            beta -= Dp/(dx * dx);
        }
        alpha *= Dp/(dx * dx);
        alpha += tau*h[i];
        alpha -= beta * tau;
        double lambda = 2.*beta / (exp(beta * dt) - 1.);
        lambda /= (sigma * sigma);
        double u = p[i] + tau;
        double uShape = lambda * u * exp(beta * dt);
        double uOut = poiss_rand(uShape, rnd_gen2);
        double gShape = uOut + 2.*alpha / (sigma * sigma);
        double gOut = gamma_rand(gShape, 1.0, rnd_gen1) / lambda;
        double uStar = gOut;
        if (uStar < tau) {
            uStar = tau;
        }
        double pStar = uStar - tau;
        pStar /= (1. + pStar * dt);
        new_p[i] = pStar;
        pMean += pStar;
        double hStar = h[i] - b * 0.5 * dt * (pStar + p[i]);
        if (x != 0) {
            hStar += dt * Dtilde * h[ind(x-1,y)] / (dx * dx);
            hStar -= dt * Dtilde * h[i] / (dx * dx);
        }
        if (x != L-1) {
            hStar += dt * Dtilde * h[ind(x+1,y)] / (dx * dx);
            hStar -= dt * Dtilde * h[i] / (dx * dx);
        }
        if (y != 0) {
            hStar += dt * Dtilde * h[ind(x,y-1)] / (dx * dx);
            hStar -= dt * Dtilde * h[i] / (dx * dx);
        }
        if (y != L-1) {
            hStar += dt * Dtilde * h[ind(x,y+1)] / (dx * dx);
            hStar -= dt * Dtilde * h[i] / (dx * dx);
        }
        new_h[i] = hStar - depth * forceTerm[i];
        hMean += new_h[i];
    }
    p = new_p;
    h = new_h;
    pMean /= double(Nsites);
    hMean /= double(Nsites);
    T += dt;
    return new_p;
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
    dx = 1.;
    Gamma = .3;
    b = .01;
    po = .01;
    kappa = 5e-5;
    L = 30;
    Nsites = L * L * L;
    p = vector<double>(Nsites, po);
    h = vector<double>(Nsites, 0.);
    pMean = po;
}

Lattice3D::Lattice3D(double myT, double myDT, double myGamma, double myB, double myPo, int myL, double myDx, double myKappa) {
    T = 0.;
    Tmax = myT;
    dt = myDT;
    Gamma = myGamma;
    b = myB;
    po = myPo;
    kappa = myKappa;
    L = myL;
    Nsites = L * L * L;
    dx = myDx;
    p = vector<double>(Nsites, po);
    h = vector<double>(Nsites, 0.);
    pMean = po;
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
    pMean = po;
}

vector<double> Lattice3D::getP() {
    return p;
}

double Lattice3D::getPmean() {
    return pMean;
}

void Lattice3D::exciteSite(int site, double seed) {
    p[site] = seed;
}

vector<double> Lattice3D::getH() {
    return h;
}

int Lattice3D::ind(int i, int j, int k) {
    i = (i % L);
    j = (j % L);
    k = (k % L);
    if (i < 0) {
        i += L;
    }
    if (j < 0) {
        j += L;
    }
    if (k < 0) {
        k += L;
    }
    int myRes = L*L*i+L*j+k;
    assert(myRes >= 0);
    assert(myRes < L*L*L);
    return myRes;
}

vector<double> Lattice3D::update() {
    double sigma = sqrt(2.);
    double Dtilde = 3.;
    vector<double> new_p (Nsites, 0.);
    vector<double> new_h (Nsites, 0.);
    pMean = 0.;
    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            for (int k = 0; k < L; k++) {
                double alpha = 0.;
                alpha += p[ind(i+1,j,k)];
                alpha += p[ind(i-1,j,k)];
                alpha += p[ind(i,j+1,k)];
                alpha += p[ind(i,j-1,k)];
                alpha += p[ind(i,j,k+1)];
                alpha += p[ind(i,j,k-1)];
                alpha /= (dx*dx);
                double beta = 1.+(kappa*T)-Gamma-Gamma*b*h[ind(i,j,k)]-6./(dx*dx);
                double lambda = 2.*beta / (sigma*sigma *(exp(beta*dt) - 1));
                double pShape = lambda * p[ind(i,j,k)] * exp(beta*dt);
                double pOut = poiss_rand(pShape, rnd_gen2);
                double gShape = pOut + 2.*alpha/(sigma*sigma);
                double pStar = gamma_rand(gShape, 1., rnd_gen1) / lambda;
                pStar /= (1 + pStar * dt);
                new_p[ind(i,j,k)] = pStar;
                pMean += pStar;
                double hStar = h[ind(i,j,k)] + 0.5*dt*(pStar + p[ind(i,j,k)]);
                hStar += dt*Dtilde*h[ind(i+1,j,k)]/(dx*dx);
                hStar += dt*Dtilde*h[ind(i-1,j,k)]/(dx*dx);
                hStar += dt*Dtilde*h[ind(i,j+1,k)]/(dx*dx);
                hStar += dt*Dtilde*h[ind(i,j-1,k)]/(dx*dx);
                hStar += dt*Dtilde*h[ind(i,j,k+1)]/(dx*dx);
                hStar += dt*Dtilde*h[ind(i,j,k-1)]/(dx*dx);
                hStar -= 6.*dt*Dtilde*h[ind(i,j,k)]/(dx*dx);
                new_h[ind(i,j,k)] = hStar;
            }
        }
    }
    pMean /= double(Nsites);
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

/* int ind_3D(int i, int j, int k, int L) {
    i = (i % L);
    j = (j % L);
    k = (k % L);
    if (i < 0) {
        i += L;
    }
    if (j < 0) {
        j += L;
    }
    if (k < 0) {
        k += L;
    }
    int myRes = L*L*i + L*j + k;
    assert (myRes >= 0);
    assert (myRes < L^3);
    return myRes;
}

vector<double> potential_grad_3D(vector<double> myX, int N, int L) {
    vector<double> res(N,0.);
    vector<double> Fx(N,0.);
    vector<double> Fy(N,0.);
    vector<double> Fz(N,0.);
    float mu = float(L) / 2.;
    float sigma = mu / 3.;
    for (int i = 0; i < N; i++) {
        float x = float(i / L^2);
        float y = float((i % L^2)/L);
        float z = float(i % L);
        float posx = x - mu;
        float posy = y - mu;
        float posz = z - mu;
        float posr2 = posx*posx + posy*posy + posz*posz;
        float expFac = exp(-posr2 / (2. * sigma * sigma));
        int ix = (int)x;
        int iy = (int)y;
        int iz = (int)z;

        if (ix == 0) {
            dxh = myX[ind_3D(ix+1,iy,iz,L)] - myX[i];
        } else if (ix == (L-1)) {
            dxh = myX[i] - myX[ind_3D(ix-1,iy,iz,L)];
        } else {
            dxh = 0.5*(myX[ind_3D(ix+1,iy,iz,L)] - myX[ind_3D(ix-1,iy,iz,L)]);
        }

        if (iy == 0) {
            dyh = myX[ind_3D(ix,iy+1,iz,L)] - myX[i];
        } else if (iy == (L-1)) {
            dyh = myX[i] - myX[ind_3D(ix,iy-1,iz,L)];
        } else {
            dyh = 0.5*(myX[ind_3D(ix,iy+1,iz,L)] - myX[ind_3D(ix,iy-1,iz,L)]);
        }

        if (iz == 0) {
            dzh = myX[ind_3D(ix,iy,iz+1,L)] - myX[i];
        } else if (iz == (L-1)) {
            dzh = myX[i] - myX[ind_3D(ix,iy,iz-1,L)];
        } else {
            dzh = 0.5*(myX[ind_3D(ix,iy,iz+1,L)] - myX[ind_3D(ix,iy,iz-1,L)]);
        }

        float tmp = -dxh*posx - dyh*posy - dyz*posz + myX[i]*(-2 + posr2 / (sigma*sigma));
        res[i] = tmp * expFac / (sigma * sigma);
    }
    return res;
}




vector<double> Lattice3D::trap_update(double tau, double depth) {
    double sigma = sqrt(2.);
    double Dtilde = kappa;
    double Dp = .2;
    vector<double> new_p(Nsites, 0.);
    vector<double> new_h(Nsites, 0.);
    pMean = 0.;
    vector<double> forceTerm = potential_grad_2D(h, Nsites, L);
    for (int i = 0; i < Nsites; i++) {
        double alpha = 0.;
        double beta = .5*h[i] - 2.*tau - Gamma; // check dimensionality here
        int x = int(i / L^2);
        int y = int((i % L^2) / L);
        int z = int(i % L);
        if (x != 0) {
            alpha += p[ind_3D(x-1,y,z,L)];
            beta -= Dp/(dx * dx);
        }
        if (x != L - 1) {
            alpha += p[ind_3D(x+1,y,z,L)];
            beta -= Dp/(dx * dx);
        }
        if (y != 0) {
            alpha += p[ind_3D(x,y-1,z,L)];
            beta -= Dp/(dx * dx);
        }
        if (y != L - 1) {
            alpha += p[ind_3D(x,y+1,z,L)];
            beta -= Dp/(dx * dx);
        }
        if (z != 0) {
            alpha += p[ind_3D(x,y,z-1,L)];
            beta -= Dp/(dx * dx);
        }
        if (z != L - 1) {
            alpha += p[ind_3D(x,y,z+1,L)];
            beta -= Dp/(dx * dx);
        }
        alpha *= Dp/(dx * dx);
        alpha += tau*h[i];
        alpha -= beta * tau;
        double lambda = 2.*beta / (exp(beta * dt) - 1.);
        lambda /= (sigma * sigma);
        double u = p[i] + tau;
        double uShape = lambda * u * exp(beta * dt);
        double uOut = poiss_rand(uShape, rnd_gen2);
        double gShape = uOut + 2.*alpha / (sigma * sigma);
        double gOut = gamma_rand(gShape, 1.0, rnd_gen1) / lambda;
        double uStar = gOut;
        if (uStar < tau) {
            uStar = tau;
        }
        double pStar = uStar - tau;
        pStar /= (1. + pStar * dt);
        new_p[i] = pStar;
        pMean += pStar;
        double hStar = h[i] - b * 0.5 * dt * (pStar + p[i]);
        if (x != 0) {
            hStar += dt * Dtilde * h[ind_3D(x-1,y,z,L)] / (dx * dx);
            hStar -= dt * Dtilde * h[i] / (dx * dx);
        }
        if (x != L-1) {
            hStar += dt * Dtilde * h[ind_3D(x+1,y,z,L)] / (dx * dx);
            hStar -= dt * Dtilde * h[i] / (dx * dx);
        }
        if (y != 0) {
            hStar += dt * Dtilde * h[ind_3D(x,y-1,z,L)] / (dx * dx);
            hStar -= dt * Dtilde * h[i] / (dx * dx);
        }
        if (y != L-1) {
            hStar += dt * Dtilde * h[ind_3D(x,y+1,z,L)] / (dx * dx);
            hStar -= dt * Dtilde * h[i] / (dx * dx);
        }
        if (z != 0) {
            hStar += dt * Dtilde * h[ind_3D(x,y,z-1,L)] / (dx * dx);
            hStar -= dt * Dtilde * h[i] / (dx * dx);
        }
        if (z != L-1) {
            hStar += dt * Dtilde * h[ind_3D(x,y,z+1,L)] / (dx * dx);
            hStar -= dt * Dtilde * h[i] / (dx * dx);
        }
        new_h[i] = hStar - depth * forceTerm[i];
    }
    p = new_p;
    h = new_h;
    pMean /= double(Nsites);
    T += dt;
    return new_p;
} */


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


double gamma_rand(double shape, double scale, boost::mt19937& rng) {
    if ((shape <= 0.) || (scale <= 0.) || (shape != shape) || (scale != scale)) {
        // cout << "GAMMA RET 0\n";
        // cout << "Gamma shape, scale: " << shape << "\t" << scale << "\n";
        return 0.;
    }
    // assert(shape > 0.);
    // assert(scale > 0.);
    boost::gamma_distribution<double> gd(shape);
    boost::variate_generator<boost::mt19937&, boost::gamma_distribution<double> > var_gamma(rng, gd);
    double myRet = scale * var_gamma();
    // if (myRet == 0.) {
    //     cout << "GAMMA OUTPUT 0\t" << (myRet * 1.e10) << "\n";
    //     cout << "shape: " << shape << "\n";
    // }
    return myRet;
}

double poiss_rand(double shape, boost::mt19937& rng) {
    if ((shape <= 0) || (shape != shape)) {
        return 0.;
    }
    boost::poisson_distribution<int, double> pd(shape);
    boost::variate_generator<boost::mt19937&, boost::poisson_distribution<int, double> > var_poisson(rng, pd);
    return var_poisson();
}
