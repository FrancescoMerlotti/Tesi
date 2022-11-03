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

int nseed, dtime, ngen;
int ntrunc;

void Input();

vector<double> NonUniformSimplex();
vector<double> Simplex();

double Norm(vector<double>);
double InnerProduct(vector<double>, vector<double>);
bool Check(vector<double>);
void Decorrelation();
void Print(string, vector<double>, int);