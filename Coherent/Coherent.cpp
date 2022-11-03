#include "./Coherent.h"

using namespace std;

int main(int argc, char** argv) {
    // Execution time calculation
    time_t time_in = time(NULL);
    // Testing if there's the correct number of arguments
    assert(argc == 3);
    // Initializing the truncation
    ntrunc = atoi(argv[1]);
    rphase = atoi(argv[2]);
    // Reading paramteres from file
    Input();
    // Initializing best and worst state
    psi_best = Vacuum();
    phi_best = Phase(0);
    cout << endl << "N = " << ntrunc << ", P = " << rphase << endl;
    // Random states QFI calculation
    ofstream Output;
    Output.open(rep + "/graph" + subrep + "/output_" + to_string(ntrunc) + ".dat");
    for(int igen = 0; igen < ngen; igen++) {
        // Initializing the states
        vector<double> phi = Phase(rphase);
        vector<double> psi = Simplex();
        while(!Check(psi))
            psi = Simplex();
        // Tracing the state with the highest QFI
        FollowBest(psi, phi);
        // Output on file
        Output << setw(15) << MeanNumber(psi) << setw(15) << H(psi, phi) << endl;
    }
    Output.close();
    // Best state
    Print("best", psi_best, phi_best);
    // Phase study
    vector<double> psi = psi_best;
    vector<double> phi = phi_best;
    Output.open(rep + "/phases" + subrep + "/output_phases_" + to_string(ntrunc) + ".dat");
    for(int igen = 0; igen < 100000; igen++) {
        // Tracing the state with the highest QFI
        FollowBest(psi, phi);
        // Output on file
        for(int i = 1; i <= ntrunc; i++)
            Output << setprecision(10) << setw(20) << phi[i];
        Output << setprecision(10) << setw(20) << H(psi, phi) << endl;
        // Replacing the phases
        phi = Phase(1);
    }
    Output.close();
    // Execution time
    cout << "Execution time: " << difftime(time(NULL), time_in) << " s" << endl << endl;
    
    return 0;
}

// ---------------------------------------------------

// Input
void Input() {
    ifstream Primes, Seed, Setup;
	Setup.open("Setup.dat");
	Setup >> nseed >> dtime >> ngen >> alpha >> theta;
	Setup.close();
    theta = theta * M_PI;
    Subrep();
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
    phase[0] = 0.;
    for(int itrunc = 1; itrunc <= ntrunc; itrunc++) {
        Decorrelation();
        if(r)
            phase[itrunc] = rnd.Rannyu(0., 2 * M_PI);
        else
            phase[itrunc] = 0.;
    }
    return phase;
}

// Vacuum state
vector<double> Vacuum() {
       vector<double> coeff(ntrunc+1, 0.);
       coeff[0] = 1.;
       return coeff;
}

// ---------------------------------------------------

double Norm(vector<double> v) {
    double sum = 0.;
    for(int i = 0; i < v.size(); i++)
        sum += v[i] * v[i];
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

// Follow best
void FollowBest(vector<double> v, vector<double> p) {
    if(H(v, p) > qfi_best) {
        qfi_best = H(v, p);
        psi_best = v;
        phi_best = p;
    }
}

// Subrepository
void Subrep() {
    if(rphase)
        subrep = "/random";
    else
        subrep = "/null";
    
    if(theta == 0.)
        subrep += "/null";
    else
        subrep += "/pi";
    // cout << subrep << endl;
}

// Print vector<double>
void Print(string t, vector<double> v, vector<double> p) {
    ofstream Output;
    int prec = 10, wd = 15;
    Output.open(rep + "/state_" + t + subrep + "/output_" + to_string(ntrunc) + ".dat");
    for(int i = 0; i < v.size(); i++)
        Output << setprecision(prec) << setw(wd) << v[i] << setw(wd) << p[i] << endl;
    Output << endl << "N_tot = " << MeanNumber(v) << ", H = " << H(v, p) << endl;
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
        sum += sqrt(k+1.) * v[k+1] * v[k] * cos(theta + p[k] - p[k+1]);
    return sum;
}

// Sum 2
double S2(vector<double> v, vector<double> p){
    double sum = 0.;
    for(int k = 0; k <= ntrunc-1; k++)
        sum += (k+1.) * sqrt(k+1.) * v[k+1] * v[k] * cos(theta + p[k] - p[k+1]);
    return sum;
}

// Sum 3
double S3(vector<double> v, vector<double> p) {
    double sum = 0.;
    for(int k = 0; k <= ntrunc-2; k++)
        sum += sqrt((k+1.)*(k+2.)) * v[k+2] * v[k] * cos(2*theta + p[k] - p[k+2]);
    return sum;
}

// ---------------------------------------------------

// Mean number of photons
double MeanNumber(vector<double> v) {
    return N1(v) + alpha * alpha;
}

// Quatum Fisher Information
double H(vector<double> v, vector<double> p) {
    double n2 = N2(v), n1 = N1(v), s1 = S1(v, p), s2 = S2(v, p), s3 = S3(v, p);
    return n2 - n1*n1 + n1 + 4*alpha*(n1*s1 - s2) + 2*alpha*alpha*(n1 + s3 + 1- 2*s1*s1);
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
