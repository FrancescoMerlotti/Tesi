#include "./Vacuum.h"

using namespace std;

int main(int argc, char** argv) {

    // Testing if there's the correct number of arguments
    assert(argc == 2);
    // Initializing the truncation
    ntrunc = atoi(argv[1]);
    // Reading paramteres from file
    Input();
    // Random states QFI calculation
    ofstream Output;
    Output.open("./graph_vacuum/output_vacuum_" + to_string(ntrunc) + ".dat");
    for(int igen = 0; igen < ngen; igen++) {
        vector<double> psi = Simplex();
        if(H(psi) > qfi) {
            qfi = H(psi);
            psi_best = psi;
        }
        Output << setw(15) << N1(psi) << setw(15) << H(psi) << endl;
    }
    Output.close();
    // Best state traced
    Print("psi", psi_best);
    cout << "N_tot = " << N1(psi_best) << ", H = " << H(psi_best) << endl;
    // Analytical best state
    Print("psi_best", Best());
    cout << "N_tot = " << N1(Best()) << ", H = " << H(Best()) << endl;
    // Fidelity
    cout << endl << "Fidelity = " << pow(InnerProduct(psi_best, Best()), 2) << endl << endl;

    return 0;
}

// ---------------------------------------------------

// Input
void Input() {
    ifstream Primes, Seed, Setup;
	Setup.open("Setup.dat");
	Setup >> nseed >> dtime >> ngen;
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
void Print(string t, vector<double> v) {
    cout << endl << t << ":" << endl;
    for(int i = 0; i < v.size(); i++)
        cout << "(" << v[i] << ")" << endl;
    cout << endl;
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

// Quatum Fisher Information
double H(vector<double> v) {
    double n1 = N1(v), n2 = N2(v);
    return n2 - n1*n1 + n1;
}