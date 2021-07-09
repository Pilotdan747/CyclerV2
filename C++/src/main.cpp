#include "cmath"
#include "include/helperFuncs.h"
#include "include/cycler.h"
#include <time.h>
#include <chrono>
#include <omp.h>
#include <stdio.h>

//Defines pi
#define pi 4*atan(1)

int main() {
    double SynodicT, re, rm, muSun, Te, Tm;

    //Set Up & Initial Conditions
    int dim1, dim2, dim3, dim4; //Number of points in each dimension
    dim1 = 50;
    dim2 = 50;
    dim3 = 50;
    dim4 = 50;

    //Radii of Earth and Mars from sun
    re = 1.495979e8;
    rm = 2.279483e8;

    muSun = 1.32712440e11;

    //Orbital Periods of Earth and Mars
    Te = 2*pi/sqrt(muSun)*pow(re, 3.0/2.0);
    Tm = 2*pi/sqrt(muSun)*pow(rm, 3.0/2.0);

    //Synodic Period
    SynodicT = 1/(fabs(1/Te - 1/Tm));

    //Set up clock for timing
    using namespace std::chrono;
    auto start = high_resolution_clock::now();

    //Main array that stores delta V values

    double ****dV = (double ****) malloc(dim1 * sizeof(double***));

    for (int i = 0; i < dim1; i++) {
        dV[i] = (double ***) malloc(dim2 * sizeof(double**));
        for (int j = 0; j < dim2; j++) {
            dV[i][j] = (double **) malloc(dim3 * sizeof(double*));
            for (int k = 0; k < dim3; k++) {
                dV[i][j][k] = (double *) malloc(dim4 * sizeof(double));
            }
        }
    }


//Main loop region
//Tests all of phi and delta T 1-3 times
//Stores results in the delta V array

int count = 0;
int count2 = 0;

FILE *outfile2 = fopen("Output2.csv", "w");

#pragma omp parallel
    {
        //Testing code
        int threadID = omp_get_thread_num();
        if (threadID == 0) {
            printf("Num threads is: %d\n", omp_get_num_threads());
            printf("Num procs is: %d\n", omp_get_num_procs());
        }
#pragma omp for
        for (int i = 0; i < dim1; i++) {
            for (int j = 0; j < dim2; j++) {
                for (int k = 0; k < dim3; k++) {
                    for (int l = 0; l < dim4; l++) {
                        //Set dT1-3 and phi for each iteration
                        double phi = 2*pi/dim1*i;
                        double dT1 = (90 + 210.0/dim2*j) * 24 * 3600;
                        double dT2 = (23 + 7.0/dim3*k) * 30 * 24 * 3600;
                        double dT3 = (90 + 210.0/dim4*l) * 24 * 3600;


                        double dT4 = SynodicT*2 - (dT1 + dT2 + dT3);

                        double ans;

                        if (dT4 < 0) {
                            ans = 1000;
                        } else {
                            //Calc delta V
                            ans = cycleMulti(dT1, dT2, dT3, phi);
                        }

                        if (isnan(ans)) {
                            count2++;
                            //printf("NAN\n");
                        }

                        if (ans < 1) {
                            count++;
                            fprintf(outfile2, "%f, ", ans);
                            fprintf(outfile2, "%d, ", i);
                            fprintf(outfile2, "%d, ", j);
                            fprintf(outfile2, "%d, ", k);
                            fprintf(outfile2, "%d\n", l);
                        }

                        //Store answer
                        dV[i][j][k][l] = ans;
                    }
                }
            }
        }
    }

    fclose(outfile2);

    //Stops lock and gets run time
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    double time = duration.count();
    printf("It took %f seconds to run\n", time/1e6);
    printf("Count: %d\n", count);
    printf("Count2: %d\n", count2);


    
    //Outputs main array to a CSV file
    //Puts size of each dimentsion in its own line
    //Prints array into 1 huge line
    FILE *outfile = fopen("Output.csv", "w");

    fprintf(outfile, "%d\n", dim1);
    fprintf(outfile, "%d\n", dim2);
    fprintf(outfile, "%d\n", dim3);
    fprintf(outfile, "%d\n", dim4);

    double min = 1000;
    double num = 1;
    for (int i = 0; i < dim1; i++) {
        for (int j = 0; j < dim2; j++) {
            for (int k = 0; k < dim3; k++) {
                for (int l = 0; l < dim4; l++) {
                    num = dV[i][j][k][l];
                    if (num < min) {
                        min = num;
                    }

                    fprintf(outfile, "%f, ", num);
                }
            }
        }
    }

    fclose(outfile);

    printf("Min: %f\n", min);

    return 0;
}
