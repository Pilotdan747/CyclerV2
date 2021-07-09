#include <stdio.h>
#include "include/cycler.h"

int main() {

    //double phi = 4.2081;
    //double dT1 = 1.1768e+07;
    //double dT2 = 6.9230e+07;
    //double dT3 = 1.8918e+07;

    //Balistic Cycler
    double phi = 1.7671;
    double dT1 = 1.5153e+07;
    double dT2 = 7.4979e+07;
    double dT3 = 1.5153e+07;

    double dV = cycleMulti(dT1, dT2, dT3, phi, false, true);

    return 0;
}