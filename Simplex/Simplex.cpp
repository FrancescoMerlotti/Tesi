#include "./Simplex.h"

using namespace std;

int main(int argc, char** argv) {
	
	assert(argc == 2);

	ntrunc = atoi(argv[1]);

	Input();

	for(int j = 0; j < ngen; j++) {
		Print("u", Simplex(), j);
		Print("n", NonUniformSimplex(), j);
	}
		

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

// Non-uniform Simplex
vector<double> NonUniformSimplex() {
	vector<double> coeff(ntrunc+1);
    double temp, sum = 0;
    for(int itrunc = 0; itrunc <= ntrunc-1; itrunc++) {
        Decorrelation();
        temp = rnd.Rannyu(0., 1. - sum);
        coeff[itrunc] = temp;
        sum += temp;
    }
	coeff[ntrunc] = 1. - sum;

    for(int itrunc = 0; itrunc <= ntrunc; itrunc++)
        coeff[itrunc] = sqrt(coeff[itrunc]);
    
    return coeff;
}

// Uniform Simplex
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
void Print(string t, vector<double> v, int j) {
	ofstream Output;
	if(j==0)
		Output.open("output_simplex_" + t + ".dat");
	else
		Output.open("output_simplex_" + t + ".dat", ios::app);
    for(int i = 0; i < v.size(); i++)
        Output << setw(15) << v[i];
    Output << endl;
	Output.close();
}