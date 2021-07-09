//
// Created by Daniel Owen on 2/22/21
// Using lambert_multi.m from Dr. Brian Kaplinger
//

#include "helperFuncs.h"
#include <cmath>
#include </usr/local/include/armadillo>

using namespace arma;

#define pi datum::pi

void lambert_multi(vec R1, vec R2, double tof, int m, double mu, vec V[2]);
int mlambert(vec R1, vec R2, double tf, int m, double mu, vec V[2]);
int lambert_LancasterBlanchard(vec r1, vec r2, double tf, int m, double muC, vec V[2]);
void LancasterBlanchard(double x, double q, int m, double Ts[4]);
void sigmax(double y, double sigs[4]);

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}


