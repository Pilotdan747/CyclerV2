#include <math.h>
#include "/usr/local/include/nlopt.h"
#include <stdio.h>
#include "include/cycler.h"

//X0 -> phi
//X1 -> dT1
//X2 -> dT2
//X3 -> dT3
int count = 0;
double optFunc(unsigned n, const double *x, double *grad, void *my_func_data) {
    count++;
    double dV = cycleMulti(x[1], x[2], x[3], x[0]);

    if (isnan(dV)) {
        dV = 100;
    }

    return dV;
}

double dT4Constraint(unsigned n, const double *x, double *grad, void *data) {
    double re = 1.495979e8;
    double rm = 2.279483e8;

    double muSun = 1.32712440e11;

    //Orbital Periods of Earth and Mars
    double Te = 2*pi/sqrt(muSun)*pow(re, 3.0/2.0);
    double Tm = 2*pi/sqrt(muSun)*pow(rm, 3.0/2.0);

    //Synodic Period
    double SynodicT = 1/(fabs(1/Te - 1/Tm));
    
    double dT4 = SynodicT*2 - (x[1] + x[2] + x[3]);
    return -dT4; //+ 365.25*24*3600;
 }

int main() {

    double lb[4] = {0, 90*24*3600, 23*30*24*3600, 90*24*3600};
    double ub[4] = {2*pi, 300*24*3600, 30*30*24*3600, 300*24*3600};
    

    nlopt_opt opt;
    //opt = nlopt_create(NLOPT_LN_COBYLA, 4);
    opt = nlopt_create(NLOPT_LN_BOBYQA, 4);
    nlopt_set_lower_bounds(opt, lb);
    nlopt_set_upper_bounds(opt, ub);
    nlopt_set_min_objective(opt, optFunc, NULL);

    nlopt_set_xtol_rel(opt, 1e-10);

    nlopt_add_inequality_constraint(opt, dT4Constraint, NULL, 1e-8);

    FILE *infile = fopen("Output2.csv", "r");
    FILE *outfile = fopen("Optimal.csv", "w");

    FILE *outfile2 = fopen("Vinf.csv", "w");
    fclose(outfile2);

    while (true) {
        double test;
        int i = 0;
        int j = 0;
        int k = 0;
        int l = 0;

        if (feof(infile)) {
            break;
        }

        fscanf(infile, "%f, %d, %d, %d, %d", &test, &i, &j, &k, &l);

        //82, 62, 46, 68

        double phi = 2*pi/50*i;
        double dT1 = (90 + 210.0/50*j) * 24 * 3600;
        double dT2 = (23 + 7.0/50*k) * 30 * 24 * 3600;
        double dT3 = (90 + 210.0/50*l) * 24 * 3600;
        double guess[4] = {phi, dT1, dT2, dT3};
        printf("Guess: %f, %f, %f, %f\n", guess[0], guess[1], guess[2], guess[3]);

        double sol;

        nlopt_optimize(opt, guess, &sol);

        if (!i==0 && !j==0 && !k==0 && !l==0) {
            printf("It took %d iterations\n", count);
            printf("Min is: %f\n", sol);
            printf("Answer: %f, %f, %f, %f\n\n", guess[0], guess[1], guess[2], guess[3]);
            fprintf(outfile, "%f, %f, %f, %f, %f\n", sol, guess[0], guess[1], guess[2], guess[3]);
            cycleMulti(guess[1], guess[2], guess[3], guess[0], true);
        }
        
    }
    fprintf(outfile, "\n");
    fclose(outfile);
    fclose(infile);
    return 0;
}