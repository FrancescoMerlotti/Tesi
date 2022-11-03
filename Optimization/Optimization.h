#ifndef __Optimization_h__
#define __Optimization_h__

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

// Random number generator setup
Random rnd;
int seed[4];
int p1, p2;

// Input
int nseed, dtime, ngen, nstep;
int ntrunc, ntemp = 5;

// Simulated annealing parameters
double temp[5] {0.1, 0.067, 0.045, 0.03, 0.02};
double dpsi, temperature;
double H_old, H_best;
double attempted, accepted;

// Alpha
double alpha, theta;
vector<double> psi, phi;
vector<double> psi_old, phi_old;
vector<double> psi_best, phi_best;

void Input();

vector<double> Simplex();
vector<double> Phase(int);
vector<double> Vacuum();

double Norm(vector<double>);
bool Check(vector<double>);
void Decorrelation();
void FollowBest();
void Print(int);
void PrintBest(int);

bool OutOfBound(double);
void Move();
double Boltzmann(double);
void Reset();

double N1(vector<double>);
double N2(vector<double>);

double S1(vector<double>, vector<double>);
double S2(vector<double>, vector<double>);
double S3(vector<double>, vector<double>);

double MeanNumber(vector<double>);
double H(vector<double>, vector<double>);
double RoundHalf(double);

#endif