#include "include/cycler.h"
#include "include/helperFuncs.h"
#include "stdio.h"

int main() {

    FILE *outfile2 = fopen("VinfS1L1.csv", "w");
    fclose(outfile2);
    double phi = 0;
    double dT = 0;
    int dim1 = 1000;
    double dV[dim1];
    int count = 0;
    double min = 1000;
    for (int i = 0; i < dim1; i++) {
        dT = (2.5 + 1.0/dim1*i)*(365.25*24.0*3600.0);

        dV[i] = cycleS1L1(dT, 0);

        if (dV[i] < min) {
            min = dV[i];
        }

        if (fabs(dV[i]) < 1e-3) {
            printf("i is: %d \t T is: %f \t dV is: %f\n", i, dT, dV[i]);
            count++;
        }

        if (i == 1 || i == 1000000 - 1) {
            printf("T2 is: %f\n", dT/365.25/24.0/3600.0);
        }
    }

    printf("Min: %f\n", min);
    double cost = cycleS1L1(2.8276*3.154e7, 0);
    printf("S1L1 Test: %f\n", cost);
    printf("Count: %d\n", count);

    FILE *outfile = fopen("OutputS1L1.csv", "w");

    fprintf(outfile, "%d\n", dim1);

    double num = 1;
    for (int i = 0; i < dim1; i++) {
        num = dV[i];

        fprintf(outfile, "%f, ", num);
    }

    fclose(outfile);

    return 0;
}