#include "./Optimization.h"

using namespace std;

int main(int argc, char** argv) {
    // Exectution time calculation
    time_t time_in = time(NULL);
    // Testing if there's the correct number of arguments
    assert(argc == 2);
    // Initializing the truncation
    ntrunc = atoi(argv[1]);
    // Input
    Input();
    // State initialization
    psi_best = Vacuum();
    phi_best = Phase(0);
    phi = phi_best;
    for(int igen = 0; igen < 1e6; igen++) {
        // Generating the state randomly
        psi = Simplex();
        // Checking the generation
        while(!Check(psi))
            psi = Simplex();
        // Following the best state generated
        FollowBest();
    }
    PrintBest(0);
    time_t time_bs = time(NULL);
    cout << endl << "Best state search execution time: " << difftime(time_bs, time_in) << " s" << endl;
    // Simulated annealing setup (starting from the best)
    psi_old = psi_best;
    phi_old = phi_best;
    H_old = H(psi_old, phi_old);
    // Simulated annealing performing
    for(int itemp = 0; itemp < ntemp; itemp++) {
        // Resetting acceptance
        Reset();
        // Setting temperature
        temperature = temp[itemp];
        // Performing Metropolis moves
        for(int istep = 0; istep < nstep; istep++) {
            // Metropolis move
            Move();
            // Tracing the best state
            FollowBest();
        }
        // Print on file
        Print(itemp);
        // Step & acceptance monitoring
        cout << endl << "Step "<< itemp <<":" << endl;
        cout << "Acceptance rate: " << accepted / attempted << endl;
        cout << endl << "-----------------------------" << endl;
    }
    // Print the best state on file
    PrintBest(1);
    time_t time_sa = time(NULL);
    psi = psi_best;
    cout << endl << "Simulated annealing execution time: " << difftime(time_sa, time_bs) << " s" << endl;
    // Best phase research
    for(int igen = 0; igen < 1e6; igen++) {
        // Replacing the phases
        phi = Phase(1);
        // Tracing the state with the highest QFI
        FollowBest();
    }
    // Print the best state on file
    PrintBest(1);
    time_t time_ph = time(NULL);
    cout << endl << "Best phases search execution time: " << difftime(time_ph, time_sa) << " s" << endl;
    // Fidelity
    phi = Phase(0);
    ofstream Fidelity;
    Fidelity.open("./fidelity" + repository + "/output_fidelity_" + to_string(ntrunc) + ".dat");
    for(int igen = 0; igen < 1e6; igen++) {
        // Generating the state randomly
        psi = Simplex();
        // Checking the generation
        while(!Check(psi))
            psi = Simplex();
        // Fidelity
        if(H(psi, phi) >= 8.38) {
            vector<double> ip = InnerProduct(psi_best, phi_best, psi, phi);
            Fidelity << setprecision(10) << setw(15) << H(psi, phi) << setw(15) << pow(ip[0], 2) + pow(ip[1], 2) << endl;
        }
    }
    Fidelity.close();
    cout << endl << "Fidelity calculation execution time: " << difftime(time(NULL), time_ph) << " s" << endl;
    // Execution time estimation
    cout << endl << "Total execution time: " << difftime(time(NULL), time_in) << " seconds" << endl << endl;

    return 0;
}

// --------------------------------------------------------

// Input
void Input() {
    ifstream Primes, Seed, Setup;
	Setup.open("Setup.dat");
    // Input parameters from file
	Setup >> ngen >> nseed >> dtime >> alpha >> theta >> dpsi >> nstep;
	Setup.close();
    if(theta == 1.)
        repository = "/pi";
    if(theta == 0.)
        repository = "/null";
    theta *= M_PI;
    // Setting up random number generator
	Primes.open("../Random/Primes");
	for(int i = 0; i < nseed; i++)
    	Primes >> p1 >> p2;
	Primes >> p1 >> p2;
	Primes.close();
	Seed.open("../Random/Seed");
	Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
	Seed.close();
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

// Norm
double Norm(vector<double> v) {
    double sum = 0.;
    for(int i = 0; i < v.size(); i++)
        sum += v[i] * v[i];
    return sum;
}

// Complex inner product
vector<double> InnerProduct(vector<double> v1, vector<double> p1, vector<double> v2, vector<double> p2) {
    vector<double> sum {0., 0.};
    // assert(v1.size() == p1.size() == v2.size() == p2.size());
    for(int i = 0; i < v1.size(); i++) {
        sum[0] += v1[i] * v2[i] * cos(p1[i] - p2[i]);
        sum[1] += v1[i] * v2[i] * sin(p1[i] - p2[i]);
    }        
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
void FollowBest() {
    if(H(psi, phi) > H_best) {
        H_best = H(psi, phi);
        psi_best = psi;
        phi_best = phi;
    }
}

// Print
void Print(int t) {
    ofstream Output;
    if(t==0)
        Output.open("./output" + repository + "/output_"+to_string(ntrunc)+".dat");
    else
        Output.open("./output" + repository + "/output_"+to_string(ntrunc)+".dat", ios::app);
    for(int i = 0; i <= ntrunc; i++)
        Output << setprecision(10) << setw(25) << psi[i];
    Output << setprecision(10) << setw(25) << MeanNumber(psi) << setw(25) << H(psi, phi) << endl;
}

// Print Best
void PrintBest(int f) {
    ofstream Output;
    int prec = 10, wd = 25;
    if(!f)
        Output.open("./best" + repository + "/output_best_" + to_string(ntrunc) + ".dat");
    else
        Output.open("./best" + repository + "/output_best_" + to_string(ntrunc) + ".dat", ios::app);
    for(int i = 0; i <= ntrunc; i++)
        Output << setprecision(prec) << setw(wd) << psi_best[i] << setw(wd) << phi_best[i] << endl;
    Output << endl << setw(wd) << MeanNumber(psi_best) << setw(wd) << H(psi_best, phi_best) << endl;
    Output << endl << "--------------------------------------------------" << endl << endl;
    Output.close();
}

// --------------------------------------------------------

// Out of Bound
bool OutOfBound(double x) {
    // Test whether the coeffiecient is out of bound
    return x >= 1 or x < 0;
}

// Move
void Move() {
    int o;
    double boltzmann_new, boltzmann_old;
    for(int itrunc = 0; itrunc <= ntrunc; itrunc++) {
        // Choosing a coefficient different from self
        int o = int(rnd.Rannyu(0., ntrunc+1));
        while(o == itrunc)
            o = int(rnd.Rannyu(0., ntrunc+1));
        // Generating an acceptable move
        double move = dpsi * rnd.Rannyu(-0.5, 0.5);
        while(OutOfBound(psi[itrunc] * psi[itrunc] + move) or OutOfBound(psi[o] * psi[o] - move)) {
            move = dpsi * rnd.Rannyu(-0.5, 0.5);
        }
        // Proposed move
        psi[itrunc] = sqrt(psi[itrunc] * psi[itrunc] + move);
        psi[o] = sqrt(psi[o] * psi[o] - move);   
    }
    // Metropolis test for the proposed move
    boltzmann_new = Boltzmann(H(psi, phi));
    boltzmann_old = Boltzmann(H_old);
    if(rnd.Rannyu() <= boltzmann_old / boltzmann_new) {
        psi_old = psi;
        H_old = H(psi, phi);
        accepted += 1.;
    }
    else {
        psi = psi_old;
    }
    attempted += 1.;
}

// Boltzmann weight
double Boltzmann(double energy) {
    // Boltzmann weight (k = 1)
    return exp(- energy / temperature);
}

// Reset
void Reset() {
    // Reset acceptance
    attempted = 0.;
    accepted = 0.;
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
