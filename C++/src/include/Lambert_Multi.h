//
// Created by Daniel Owen on 2/22/21
// Using lambert_multi.m from Dr. Brian Kaplinger
//

#include "helperFuncs.h"
#include <cmath>

void lambert_multi(vector R1, vector R2, double tof, int m, double mu, vector V[2]);
void mlambert(vector R1, vector R2, double tf, int m, double mu, vector V[2]);
void lambert_LancasterBlanchard(vector r1, vector r2, double tf, int m, double muC, vector V[2]);
void LancasterBlanchard(double x, double q, int m, double Ts[4]);
void sigmax(double y, double sigs[4]);

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}


