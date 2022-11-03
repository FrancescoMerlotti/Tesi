#include "./Thermal.h"

using namespace std;

int main(int argc, char** argv) {

    // Testing if there's the correct number of arguments
    assert(argc == 3);
    // Initializing the truncation
    ntrunc = atoi(argv[1]);
    rphase = atoi(argv[2]);
    cout << endl << "N = " << ntrunc << endl;
    // Reading paramteres from file
    Input();
    // Random states QFI calculation
    ofstream Output;
    Output.open("./graph_thermal/output_thermal_" + to_string(ntrunc) + ".dat");
    for(int igen = 0; igen < ngen; igen++) {
        // Initializing the states
        vector<double> phi = Phase(rphase);
        vector<double> psi = Simplex();
        while(!Check(psi))
            psi = Simplex();
        // Tracing the state with the highest QFI
        if(H(psi, phi) > qfi) {
            qfi = H(psi, phi);
            psi_best = psi;
            phi_best = phi;
        }
        // Output on file
        Output << setw(15) << MeanNumber(psi) << setw(15) << H(psi, phi) << endl;
    }
    Output.close();
    // Best state
    Print("psi", psi_best, phi_best);
    
    return 0;
}

// ---------------------------------------------------

// Input
void Input() {
    ifstream Primes, Seed, Setup;
	Setup.open("Setup.dat");
	Setup >> nseed >> dtime >> ngen >> rphase;
	Setup.close();
	Primes.open("../Random/Primes");
	for(int i = 0; i < nseed; i++)
    	Primes >> p1 >> p2;
	Primes >> p1 >> p2;
    // cout << p1 << " " << p2 << endl;
	Primes.close();
	Seed.open("../Random/Seed");
	Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
	Seed.close();
	// Random setup
	rnd.SetRandom(seed, p1, p2);
}

// ---------------------------------------------------

// Simplex
vector<double> Simplex() {
    vector<double> coeff(ntrunc+1);
    double temp, sum = 0;
    for(int itrunc = 0; itrunc <= ntrunc; itrunc++) {
        Decorrelation();
        temp = rnd.Exponential(1.);
        coeff[itrunc] = temp;
        sum += temp;
    }

    for(int itrunc = 0; itrunc <= ntrunc; itrunc++)
        coeff[itrunc] = sqrt(coeff[itrunc] / sum);
    
    return coeff;
}

// Phases
vector<double> Phase(int r) {
    vector<double> phase(ntrunc+1);
    for(int itrunc = 0; itrunc <= ntrunc; itrunc++) {
        Decorrelation();
        if(r)
            phase[itrunc] = rnd.Rannyu(0., 2 * M_PI);
        else
            phase[itrunc] = 0.;
    }
    return phase;
}

// Best
vector<double> Best() {
    vector<double> best(ntrunc+1, 0.);
    best[0] = sqrt(0.5 * double(ntrunc-1) / ntrunc);
    best[ntrunc] = sqrt(0.5 * double(ntrunc+1) / ntrunc);
    return best;
}

// ---------------------------------------------------

// Norm
double Norm(vector<double> v) {
    double sum = 0.;
    for(int i = 0; i < v.size(); i++)
        sum += v[i] * v[i];
    return sum;
}

// Inner product
double InnerProduct(vector<double> v, vector<double> w) {
    double sum = 0.;
    assert(v.size() == w.size());
    for(int i = 0; i < v.size(); i++)
        sum += v[i] * w[i];
    return sum;
}

// Check real vectors
bool Check(vector<double> v) {
    // Test whether the state is normalized
	return Norm(v) == 1.;
}

// Decorrelation
void Decorrelation() {
    for(int i = 0; i < dtime; i++)
        rnd.Rannyu();
}

// Print vector<double>
void Print(string t, vector<double> v, vector<double> p) {
    ofstream Output;
    Output.open("./best_thermal/output_" + to_string(ntrunc) + ".dat");
    Output << endl << t << ":" << endl;
    for(int i = 0; i < v.size(); i++)
        Output << setprecision(10) << "(" << v[i] <<", " << p[i] <<")" << endl;
    // Number of photons & quantum Fisher information
    Output << endl << "N_tot = " << MeanNumber(psi_best) << ", H = " << H(psi_best, phi_best) << endl;
    // Fidelity
    Output << endl << "Fidelity = " << InnerProduct(v, Best()) << endl;
    Output.close();
}

// ---------------------------------------------------

// Mean number of photos brought by psi
double N1(vector<double> v) {
    double sum = 0.;
    for(int k = 1; k <= ntrunc; k++)
        sum += k * v[k] * v[k];
    return sum;
}

// Mean squared number of photos brought by psi
double N2(vector<double> v) {
    double sum = 0.;
    for(int k = 1; k <= ntrunc; k++)
        sum += k * k * v[k] * v[k];
    return sum;
}

// ---------------------------------------------------

// Sum 1
double S1(vector<double> v, vector<double> p) {
    double sum = 0.;
    for(int k = 0; k <= ntrunc-1; k++)
        sum += sqrt(k+1) * v[k+1] * v[k] * cos(p[k] - p[k+1]);
    return sum;
}

// Sum 2
double S2(vector<double> v, vector<double> p){
    double sum = 0.;
    for(int k = 1; k <= ntrunc-1; k++)
        sum += sqrt(k+1) * v[k+1] * v[k] * sin(p[k] - p[k+1]);
    return sum;
}

// ---------------------------------------------------

// Mean number of photons
double MeanNumber(vector<double> v) {
    return N1(v) + 1. / (exp(1./t) - 1.);
}

// Quatum Fisher Information
double H(vector<double> v, vector<double> p) {
    double n2 = N2(v), n1 = N1(v), s1 = S1(v, p), s2 = S2(v, p), tp = 1. / (exp(1./t) - 1.);
    return n2 - n1*n1 + n1 + tp*(2*n1 + 1.) - 2*(s2 + s1)/sinh(1./t);
}

// Round Half
double RoundHalf(double a) {
    if(a - int(a) < 0.25)
        return int(a);
    if(a - int(a) < 0.75)
        return int(a) + 0.5;
    else
        return int(a) + 1.;
}