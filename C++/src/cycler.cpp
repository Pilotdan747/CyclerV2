//
// Created by Daniel Owen on 2019-05-17.
//

#include "include/cycler.h"
#include <chrono>
#include <stdio.h>

// NOT USED ANYMORE
/*double cycle(double dT1, double dT2, double dT3, double phi) {
    double re, rm, muSun, Ve, Vm, Te, Tm, SynodicT, dT4, dThetaM, thetaM, dThetaE, thetaE, dT, dV,
            dV1, dV2, dV3, dV4;
    vector Re1, Rm1, Rm2, Rm3, Re4, Re5, vEarth, vMars;

    re = 1.495979e8;
    rm = 2.279483e8;

    muSun = 1.32712440e11;

    //Assume CRP Model
    Ve = sqrt(muSun/re);
    Vm = sqrt(muSun/rm);

    Te = 2*pi/sqrt(muSun)*pow(re, 3.0/2.0);
    Tm = 2*pi/sqrt(muSun)*pow(rm, 3.0/2.0);

    SynodicT = 1/(fabs(1/Te - 1/Tm));

    //Already check for negative
    dT4 = SynodicT*2 - (dT1 + dT2 + dT3);

    thetaE = 0;
    thetaM = phi;


// Earth to Mars
    //Position Vectors
    Re1.x = re * 1; Re1.y = re * 0; Re1.z = re * 0;
    Rm1.x = rm * cos(phi); Rm1.y = rm * sin(phi); Rm1.z = rm * 0;

    dThetaM = 2 * pi * (dT1 / Tm);
    thetaM += dThetaM;

    dThetaE = 2 * pi * (dT1 / Te);
    thetaE += dThetaE;

    Rm2.x = rm * cos(thetaM); Rm2.y = rm * sin(thetaM); Rm2.z = 0;

    //Trajectory Solver
    vector V12[2];
    lambert_battin(Re1, Rm2, dT1, muSun, 0, V12);

    //Velocity Vectors
    vEarth.x = Ve * 0; vEarth.y = Ve * 1; vEarth.z = Ve * 0;
    vMars.x = Vm * -1*sin(thetaM); vMars.y = Vm * cos(thetaM); vMars.z = Vm * 0;

    //Vinf vectors
    vector VinfE1 = vinf(V12[0], vEarth);
    vector VinfM2 = vinf(V12[1], vMars);

// Mars to Mars
    dThetaM = 2 * pi * (dT2 / Tm);
    thetaM += dThetaM;

    dThetaE = 2 * pi * (dT2 / Te);
    thetaE += dThetaE;

    Rm3.x = rm * cos(thetaM); Rm3.y = rm * sin(thetaM); Rm3.z = rm * 0;

    //Trajectory Solver
    vector V34[2];
    lambert_battin(Rm2, Rm3, dT2, muSun, 0, V34);

    //Vinf Vector leaving
    vector VinfM3 = vinf(V34[0], vMars);

    vMars.x = Vm * -1*sin(thetaM); vMars.y = Vm * cos(thetaM); vMars.z = Vm * 0;

    //Vinf vector arriving
    vector VinfM4 = vinf(V34[1], vMars);

// Mars to Earth
    dThetaM = 2 * pi * (dT3 / Tm);
    thetaM += dThetaM;

    dThetaE = 2 * pi * (dT3 / Te);
    thetaE += dThetaE;

    Re4.x = re * cos(thetaE); Re4.y = re * sin(thetaE); Re4.z = re * 0;

    //Trajectory Solver
    vector V56[2];
    lambert_battin(Rm3, Re4, dT3, muSun, 0, V56);

    vEarth.x = Ve * -1*sin(thetaE); vEarth.y = Ve * cos(thetaE); vEarth.z = Ve * 0;

    //Vinf vectors
    vector VinfM5 = vinf(V56[0], vMars);
    vector VinfE6 = vinf(V56[1], vEarth);

// Earth to Earth
    dThetaM = 2 * pi * (dT4 / Tm);
    thetaM += dThetaM;

    dThetaE = 2 * pi * (dT4 / Te);
    thetaE += dThetaE;

    Re5.x = re * cos(thetaE); Re5.y = re * sin(thetaE); Re5.z = re * 0;

    //Trajectory Solver
    vector V78[2];
    lambert_battin(Re4, Re5, dT4, muSun, 0, V78);

    vector VinfE7 = vinf(V78[0], vEarth);

    vEarth.x = Ve * -1*sin(thetaE); vEarth.y = Ve * cos(thetaE); vEarth.z = Ve * 0;

    vector VinfE8 = vinf(V78[1], vEarth);

// Totals
    dV1 = fabs(norm(VinfM3) - norm(VinfM2));
    dV2 = fabs(norm(VinfM5) - norm(VinfM4));
    dV3 = fabs(norm(VinfE7) - norm(VinfE6));
    dV4 = fabs(norm(VinfE8) - norm(VinfE1));
    dV = dV1 + dV2 + dV3 + dV4;

    return dV;
}*/

// USE THIS ONE
double cycleMulti(double dT1, double dT2, double dT3, double phi, bool outFlag, bool outFlag2) {
    double re, rm, muSun, Ve, Vm, Te, Tm, SynodicT, dT4, dThetaM, thetaM, dThetaE, thetaE, dT, dV,
            dV1, dV2, dV3, dV4;
    vector Re1, Rm1, Rm2, Rm3, Re4, Re5, vEarth, vMars;

    re = 1.495979e8;
    rm = 2.279483e8;

    muSun = 1.32712440e11;

    //Assume CRP Model
    Ve = sqrt(muSun/re);
    Vm = sqrt(muSun/rm);

    Te = 2*pi/sqrt(muSun)*pow(re, 3.0/2.0);
    Tm = 2*pi/sqrt(muSun)*pow(rm, 3.0/2.0);

    SynodicT = 1/(fabs(1/Te - 1/Tm));

    //Already check for negative
    dT4 = SynodicT*2 - (dT1 + dT2 + dT3);

    thetaE = 0;
    thetaM = phi;

    FILE *outFileR;
    FILE *outFileV;
    if (outFlag2) {
        outFileR = fopen("SingleCycleR.csv", "w");
        outFileV = fopen("SingleCycleV.csv", "w");
    }

// Earth to Mars
    //Position Vectors
    Re1.x = re * 1; Re1.y = re * 0; Re1.z = re * 0;
    Rm1.x = rm * cos(phi); Rm1.y = rm * sin(phi); Rm1.z = rm * 0;

    dThetaM = 2 * pi * (dT1 / Tm);
    thetaM += dThetaM;

    dThetaE = 2 * pi * (dT1 / Te);
    thetaE += dThetaE;

    Rm2.x = rm * cos(thetaM); Rm2.y = rm * sin(thetaM); Rm2.z = 0;

    //Trajectory Solver
    vector V12[2];
    lambert_multi(Re1, Rm2, dT1, 0, muSun, V12);

    //Velocity Vectors
    vEarth.x = Ve * 0; vEarth.y = Ve * 1; vEarth.z = Ve * 0;
    vMars.x = Vm * -1*sin(thetaM); vMars.y = Vm * cos(thetaM); vMars.z = Vm * 0;

    //Vinf vectors
    vector VinfE1 = vinf(V12[0], vEarth);
    vector VinfM2 = vinf(V12[1], vMars);

    if (outFlag2) {
        fprintf(outFileR, "%f, %f, %f\n", Re1.x, Re1.y, Re1.z);
    }

// Mars to Mars
    dThetaM = 2 * pi * (dT2 / Tm);
    thetaM += dThetaM;

    dThetaE = 2 * pi * (dT2 / Te);
    thetaE += dThetaE;

    Rm3.x = rm * cos(thetaM); Rm3.y = rm * sin(thetaM); Rm3.z = rm * 0;

    //Trajectory Solver
    vector V34[2];
    vector V34pri[2];
    lambert_multi(Rm2, Rm3, dT2, -1, muSun, V34); //-1
    lambert_multi(Rm2, Rm3, dT2, 1, muSun, V34pri); //1

    //Candidate Vinf
    vector VinfM31 = vinf(V34[0], vMars);
    vector VinfM32 = vinf(V34pri[0], vMars);

    vMars.x = Vm * -1*sin(thetaM); vMars.y = Vm * cos(thetaM); vMars.z = Vm * 0;

    //Candidate Vinf
    vector VinfM41 = vinf(V34[1], vMars);
    vector VinfM42 = vinf(V34pri[1], vMars);

    if (outFlag2) {
        fprintf(outFileR, "%f, %f, %f\n", Rm2.x, Rm2.y, Rm2.z);
    }

// Mars to Earth
    dThetaM = 2 * pi * (dT3 / Tm);
    thetaM += dThetaM;

    dThetaE = 2 * pi * (dT3 / Te);
    thetaE += dThetaE;

    Re4.x = re * cos(thetaE); Re4.y = re * sin(thetaE); Re4.z = re * 0;

    //Trajectory Solver
    vector V56[2];
    lambert_multi(Rm3, Re4, dT3, 0, muSun, V56);

    vEarth.x = Ve * -1*sin(thetaE); vEarth.y = Ve * cos(thetaE); vEarth.z = Ve * 0;

    //Vinf vectors
    vector VinfM5 = vinf(V56[0], vMars);
    vector VinfE6 = vinf(V56[1], vEarth);

    if (outFlag2) {
        fprintf(outFileR, "%f, %f, %f\n", Rm3.x, Rm3.y, Rm3.z);
    }

// Earth to Earth
    dThetaM = 2 * pi * (dT4 / Tm);
    thetaM += dThetaM;

    dThetaE = 2 * pi * (dT4 / Te);
    thetaE += dThetaE;

    double test = thetaM - thetaE;

    Re5.x = re * cos(thetaE); Re5.y = re * sin(thetaE); Re5.z = re * 0;

    //Trajectory Solver
    vector V78[2];
    vector V78pri[2];

    lambert_multi(Re4, Re5, dT4, 1, muSun, V78); // 1 
    lambert_multi(Re4, Re5, dT4, -1, muSun, V78pri); // -1

    vector VinfE71 = vinf(V78[0], vEarth);
    vector VinfE72 = vinf(V78pri[0], vEarth);

    vEarth.x = Ve * -1*sin(thetaE); vEarth.y = Ve * cos(thetaE); vEarth.z = Ve * 0;

    vector VinfE81 = vinf(V78[1], vEarth);
    vector VinfE82 = vinf(V78pri[1], vEarth);

    if (outFlag2) {
        fprintf(outFileR, "%f, %f, %f\n", Re4.x, Re4.y, Re4.z);
    }

// Totals
    // pick Vinf M3 & M4
    // Start of orbit
    double dV11 = fabs(norm(VinfM31) - norm(VinfM2));  // orbit 1
    double dV12 = fabs(norm(VinfM32) - norm(VinfM2));  // orbit 2

    // end of orbit
    double dV21 = fabs(norm(VinfM5) - norm(VinfM41));  // orbit 1
    double dV22 = fabs(norm(VinfM5) - norm(VinfM42));  // orbit 2

    // Combos
    // dV11, dV21
    // dV12, dV22

    double test1 = dV11 + dV21;
    double test2 = dV12 + dV22;

    if (test1 >= test2) {
        dV1 = dV12;
        dV2 = dV22;
        if (outFlag2) {
            fprintf(outFileV, "%f, %f, %f\n", V12[0].x, V12[0].y, V12[0].z);
            fprintf(outFileV, "%f, %f, %f\n", V34pri[0].x, V34pri[0].y, V34pri[0].z);
        }
    } else {
        dV1 = dV11;
        dV2 = dV21;
        if (outFlag2) {
            fprintf(outFileV, "%f, %f, %f\n", V12[0].x, V12[0].y, V12[0].z);
            fprintf(outFileV, "%f, %f, %f\n", V34[0].x, V34[0].y, V34[0].z);
        }
    } 

    // pick Vinf E7 & E8
    double dV31 = fabs(norm(VinfE71) - norm(VinfE6));
    double dV32 = fabs(norm(VinfE72) - norm(VinfE6));

    double dV41 = fabs(norm(VinfE81) - norm(VinfE1));
    double dV42 = fabs(norm(VinfE82) - norm(VinfE1));

    double test3 = dV31 + dV41;
    double test4 = dV32 + dV42;

    if (test3 >= test4) {
        dV3 = dV32;
        dV4 = dV42;
        if (outFlag2) {
            fprintf(outFileV, "%f, %f, %f\n", V56[0].x, V56[0].y, V56[0].z);
            fprintf(outFileV, "%f, %f, %f\n", V78pri[0].x, V78pri[0].y, V78pri[0].z);
        }
    } else {
        dV3 = dV31;
        dV4 = dV41;
        if (outFlag2) {
            fprintf(outFileV, "%f, %f, %f\n", V56[0].x, V56[0].y, V56[0].z);
            fprintf(outFileV, "%f, %f, %f\n", V78[0].x, V78[0].y, V78[0].z);
        }
    }

    dV = dV1 + dV2 + dV3 + dV4;

    if (outFlag) {
        FILE *outfile = fopen("Vinf.csv", "a");
        fprintf(outfile, "%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f\n", norm(VinfE1), norm(VinfM2), norm(VinfM31), norm(VinfM32), 
            norm(VinfM41), norm(VinfM42), norm(VinfM5), norm(VinfE6), norm(VinfE71), norm(VinfE72), norm(VinfE81), norm(VinfE82));
        fclose(outfile);
    }

    return dV;
}

double cycleS1L1(double dT1, double phi) {
    double re, rm, muSun, Ve, Vm, Te, Tm, SynodicT, dT4, dThetaM, thetaM, dThetaE, thetaE, dT, dV,
            dV1, dV2, dV3, dV4;
    vector Re1, Rm1, Rm2, Rm3, Re2, Re3, Re5, vEarth, vMars, vEarth2, vEarth3;

    re = 1.495979e8;
    rm = 2.279483e8;

    muSun = 1.32712440e11;

    //Assume CRP Model
    Ve = sqrt(muSun/re);
    Vm = sqrt(muSun/rm);

    Te = 2*pi/sqrt(muSun)*pow(re, 3.0/2.0);
    Tm = 2*pi/sqrt(muSun)*pow(rm, 3.0/2.0);

    SynodicT = 1/(fabs(1/Te - 1/Tm));

    //Already check for negative
    //double dT2 = SynodicT*2 - dT1;
    double dT2 = (4 + 2.0/7.0)*(365.25*24.0*3600.0) - dT1;

    thetaE = 0;
    thetaM = phi;

// Earth to Earth 1
    //Position Vectors
    Re1.x = re * 1; Re1.y = re * 0; Re1.z = re * 0;
    Rm1.x = rm * cos(phi); Rm1.y = rm * sin(phi); Rm1.z = rm * 0;

    dThetaM = 2 * pi * (dT1 / Tm);
    thetaM += dThetaM;

    dThetaE = 2 * pi * (dT1 / Te);
    thetaE += dThetaE;

    Re2.x = re * cos(thetaE); Re2.y = re * sin(thetaE); Re2.z = 0;

    //Trajectory Solver
    vector V12[2];
    lambert_multi(Re1, Re2, dT1, -1, muSun, V12); // S1 -> Short
    

    //Velocity Vectors
    vEarth.x = Ve * 0; vEarth.y = Ve * 1; vEarth.z = Ve * 0;
    vEarth2.x = Ve * -1*sin(thetaE); vEarth2.y = Ve * cos(thetaE); vEarth2.z = Ve * 0;

    //Vinf vectors
    vector VinfE1 = vinf(V12[0], vEarth);
    vector VinfE2 = vinf(V12[1], vEarth2);

// Earth to Earth 2
    dThetaM = 2 * pi * (dT2 / Tm);
    thetaM += dThetaM;

    dThetaE = 2 * pi * (dT2 / Te);
    thetaE += dThetaE;

    Re3.x = re * cos(thetaE); Re3.y = re * sin(thetaE); Re3.z = re * 0;

    //Trajectory Solver
    vector V34[2];
    lambert_multi(Re2, Re3, dT2, 1, muSun, V34); // L1 -> long

    //Candidate Vinf
    vector VinfE3 = vinf(V34[0], vEarth2);

    vEarth3.x = Ve * -1*sin(thetaE); vEarth3.y = Ve * cos(thetaE); vEarth3.z = Ve * 0;

    //Candidate Vinf
    vector VinfE4 = vinf(V34[1], vEarth3);

// Totals
    dV1 = fabs(norm(VinfE2) - norm(VinfE3));
    dV2 = fabs(norm(VinfE4) - norm(VinfE1));

    dV = dV1 + dV2;

    FILE *outfile = fopen("VinfS1L1.csv", "a");
    fprintf(outfile, "%f, %f, %f, %f, \n", norm(VinfE1), norm(VinfE2), norm(VinfE3), norm(VinfE4));
    fclose(outfile);

    return dV;
}
