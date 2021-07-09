#include "include/cycler.h"
#include <iostream>

int main() {

    double dT1 = 6*30*24*3600;
    double dT2 = 24*30*24*3600;
    double dT3 = 6*30*24*3600;
    double phi = 30*pi/180;

    printf("Starting\n");

    double dV = cycleMulti(dT1, dT2, dT3, phi);

    printf("It Ran!! dV is: %f\n", dV);

    return 0;
}