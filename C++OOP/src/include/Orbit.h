#include "helperFuncs.h"
#include <math.h>
#include <stdio.h>
#include </usr/local/include/armadillo>

#define pi datum::pi

using namespace arma;

class Orbit {
    public:
        double a;
        double e;
        double inc;
        double RAAN;
        double argPeri;
        double TA;

        Orbit(double);
        Orbit(double, vec, vec);
        Orbit(double, double, double, double, double, double, double);
        void Kepler(double);
        void state2OE(vec, vec);
        vec* getState();

     private:
        double mu;
        vec state[2];
};