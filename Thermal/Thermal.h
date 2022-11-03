#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cassert>
#include <vector>
#include <string>
#include "../Random/random.h"

using namespace std;

Random rnd;
int seed[4];
int p1, p2;

int nseed, dtime, ngen, rphase;
int ntrunc;

double qfi = 0.;
double t = 0.02;
vector<double> psi_best, phi_best;

void Input();

vector<double> Simplex();
vector<double> Phase(int);
vector<double> Best();

double Norm(vector<double>);
double InnerProduct(vector<double>, vector<double>);
bool Check(vector<double>);
void Decorrelation();
void Print(string, vector<double>, vector<double>);

double N1(vector<double>);
double N2(vector<double>);

double S1(vector<double>, vector<double>);
double S2(vector<double>, vector<double>);

double MeanNumber(vector<double>);
double H(vector<double>, vector<double>);
double RoundHalf(double);
