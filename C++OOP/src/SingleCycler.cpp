#include <stdio.h>
#include "include/cycler.h"

int main() {

    //double phi = 4.2081;
    //double dT1 = 1.1768e+07;
    //double dT2 = 6.9230e+07;
    //double dT3 = 1.8918e+07;

    //Balistic Cycler
    double phi = 1.759026;
    double dT1 = 15063676.236436;
    double dT2 = 74899160.094651;
    double dT3 = 15341387.382460;

    double dV = cycleMulti(dT1, dT2, dT3, phi, true);
    printf("dV: \n%f\n", dV);

    return 0;
}