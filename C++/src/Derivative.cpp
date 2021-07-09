#include "include/cycler.h"
#include <stdio.h>

int main() {

    FILE *infile = fopen("Output2.csv", "r");

    double test;
    int i, j, k, l;

    FILE *outfile = fopen("Derivate.csv", "w");

    while (true) {
        if (feof(infile)) {
            break;
        }

        fscanf(infile, "%f, %d, %d, %d, %d", &test, &i, &j, &k, &l);
        double phi = 2*pi/100*i;
        double dT1 = (90 + 210.0/100*j) * 24 * 3600;
        double dT2 = (23 + 7.0/100*k) * 30 * 24 * 3600;
        double dT3 = (90 + 210.0/100*l) * 24 * 3600;

        // f(x)
        double dV = cycleMulti(dT1, dT2, dT3, phi);

        // f(x+h)
        double dVdT1 = cycleMulti(dT1 + 1, dT2, dT3, phi);
        double dVdT2 = cycleMulti(dT1, dT2 + 1, dT3, phi);
        double dVdT3 = cycleMulti(dT1, dT2, dT3 + 1, phi);
        double dVPhi = cycleMulti(dT1, dT2, dT3, phi + 1*pi/180.0);

        // f(x+h) - f(x)
        double fpridT1 = (dVdT1 - dV);
        double fpridT2 = (dVdT2 - dV);
        double fpridT3 = (dVdT3 - dV);
        double fpriPhi = (dVPhi - dV)/(1*pi/180.0);

        printf("%f\n", dVdT1 - dV);

        fprintf(outfile, "%f, ", fpridT1);
        fprintf(outfile, "%f, ", fpridT2);
        fprintf(outfile, "%f, ", fpridT3);
        fprintf(outfile, "%f\n", fpriPhi);
    }

    return 0;
}