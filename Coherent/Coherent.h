#ifndef __Coherent_h__
#define __Coherent_h__

#include <iostream>
#include <iomanip>
#include <cmath>
#include <time.h>
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
string rep = "./output_coherent", subrep;

double qfi_best = 0.;
double alpha, theta;
vector<double> psi_best, phi_best;

void Input();

vector<double> Simplex();
vector<double> Phase(int);
vector<double> Vacuum();

double Norm(vector<double>);
bool Check(vector<double>);
void Decorrelation();
void FollowBest(vector<double>, vector<double>);
void Subrep();
void Print(string, vector<double>, vector<double>);

double N1(vector<double>);
double N2(vector<double>);

double S1(vector<double>, vector<double>);
double S2(vector<double>, vector<double>);
double S3(vector<double>, vector<double>);

double MeanNumber(vector<double>);
double H(vector<double>, vector<double>);
double RoundHalf(double);

#endif