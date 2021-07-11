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

    Te = 2*PI/sqrt(muSun)*pow(re, 3.0/2.0);
    Tm = 2*PI/sqrt(muSun)*pow(rm, 3.0/2.0);

    SynodicT = 1/(fabs(1/Te - 1/Tm));

    //Already check for negative
    dT4 = SynodicT*2 - (dT1 + dT2 + dT3);

    thetaE = 0;
    thetaM = phi;


// Earth to Mars
    //Position Vectors
    Re1.x = re * 1; Re1.y = re * 0; Re1.z = re * 0;
    Rm1.x = rm * cos(phi); Rm1.y = rm * sin(phi); Rm1.z = rm * 0;

    dThetaM = 2 * PI * (dT1 / Tm);
    thetaM += dThetaM;

    dThetaE = 2 * PI * (dT1 / Te);
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
    dThetaM = 2 * PI * (dT2 / Tm);
    thetaM += dThetaM;

    dThetaE = 2 * PI * (dT2 / Te);
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
    dThetaM = 2 * PI * (dT3 / Tm);
    thetaM += dThetaM;

    dThetaE = 2 * PI * (dT3 / Te);
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
    dThetaM = 2 * PI * (dT4 / Tm);
    thetaM += dThetaM;

    dThetaE = 2 * PI * (dT4 / Te);
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

// States
// 1 -> Leave Earth
// 2 -> Arrive at Mars/Leave Mars for Mars
// 3 -> Arrive at Mars/Leave Mars for Earth
// 4 -> Arrive at Earth/Leave Earth for Earth
// 5 -> Arrive at Earth

// Vinf
// 1 -> Earth Departure
// 2 -> Mars Arrival
// 3 -> Mars Departure
// 4 -> Mars Arrival
// 5 -> Mars Departure
// 6 -> Earth Arrival
// 7 -> Earth Departure
// 8 -> Earth Arrival

// USE THIS ONE
double cycleMulti(double dT1, double dT2, double dT3, double phi, bool flag1, bool flag2, vec VinfArr[8]) {
    vector vEarth;
    vector vMars;

    double re = 1.495979e8;
    double rm = 2.279483e8;

    double muSun = 1.32712440e11;

    Orbit* Earth = new Orbit(muSun, re, 0, 0, 0, 0, 0);
    Orbit* Mars = new Orbit(muSun, rm, 0, 0, 0, 0, phi);

    //Assume CRP Model
    double Te = 2*pi/sqrt(muSun)*pow(re, 3.0/2.0);
    double Tm = 2*pi/sqrt(muSun)*pow(rm, 3.0/2.0);

    double SynodicT = 1/(fabs(1/Te - 1/Tm));

    //Already check for negative
    double dT4 = SynodicT*2 - (dT1 + dT2 + dT3);

// Earth to Mars
    //Position Vectors
    vec* state1E = Earth->getState();
    vec* state1M = Mars->getState();

    vec Re1 = state1E[0];
    vec Ve1 = state1E[1];

    Earth->Kepler(dT1);
    Mars->Kepler(dT1);

    vec* state2E = Earth->getState();
    vec* state2M = Mars->getState();

    vec Rm2 = state2M[0];
    vec Vm2 = state2M[1];

    //Trajectory Solver
    vec V12[2];
    lambert_multi(Re1, Rm2, dT1, 0, muSun, V12);

    //Velocity Vectors
    vEarth.x = Ve1(0); vEarth.y = Ve1(1); vEarth.z = Ve1(2);
    vMars.x = Vm2(0); vMars.y = Vm2(1); vMars.z = Vm2(2);

    //Vinf vectors
    vec VinfE1 = V12[0] - Ve1;
    vec VinfM2 = V12[1] - Vm2;

// Mars to Mars
    Earth->Kepler(dT2);
    Mars->Kepler(dT2);

    vec* state3E = Earth->getState();
    vec* state3M = Mars->getState();

    vec Rm3 = state3M[0];
    vec Vm3 = state3M[1];

    //Trajectory Solver
    vec V34[2];
    vec V34pri[2];
    lambert_multi(Rm2, Rm3, dT2, -1, muSun, V34); //-1
    lambert_multi(Rm2, Rm3, dT2, 1, muSun, V34pri); //1

    //Candidate Vinf
    vec VinfM31 = V34[0] - Vm2;
    vec VinfM32 = V34pri[0] - Vm2;

    vMars.x = Vm3(0); vMars.y = Vm3(1); vMars.z = Vm3(2);

    //Candidate Vinf
    vec VinfM41 = V34[1] - Vm3;
    vec VinfM42 = V34pri[1] - Vm3;

// Mars to Earth
    Earth->Kepler(dT3);
    Mars->Kepler(dT3);

    vec* state4E = Earth->getState();
    vec* state4M = Mars->getState();

    vec Re4 = state4E[0];
    vec Ve4 = state4E[1];

    //Trajectory Solver
    vec V56[2];
    lambert_multi(Rm3, Re4, dT3, 0, muSun, V56);

    vEarth.x = Ve4(0); vEarth.y = Ve4(1); vEarth.z = Ve4(2);

    //Vinf vectors
    vec VinfM5 = V56[0] - Vm3;
    vec VinfE6 = V56[1] - Ve4;

// Earth to Earth
    Earth->Kepler(dT4);
    Mars->Kepler(dT4);

    vec* state5E = Earth->getState();
    vec* state5M = Mars->getState();

    vec Re5 = state5E[0];
    vec Ve5 = state5E[1];

    //Trajectory Solver
    vec V78[2];
    vec V78pri[2];

    lambert_multi(Re4, Re5, dT4, 1, muSun, V78); // 1 
    lambert_multi(Re4, Re5, dT4, -1, muSun, V78pri); // -1


    vec VinfE71 = V78[0] -  Ve4;
    vec VinfE72 = V78pri[0] - Ve4;

    vEarth.x = state5E[1](0); vEarth.y = state5E[1](1); vEarth.z = state5E[1](2);

    vec VinfE81 = V78[1] - Ve5;
    vec VinfE82 = V78pri[1] - Ve5;

// Totals
    // PIck Vinf M3 & M4
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

    double dV1;
    double dV2;
    bool check1 = false;
    if (test1 >= test2) {
        dV1 = dV12;
        dV2 = dV22;
        check1 = true;
    } else {
        dV1 = dV11;
        dV2 = dV21;
    } 

    // PIck Vinf E7 & E8
    double dV31 = fabs(norm(VinfE71) - norm(VinfE6));
    double dV32 = fabs(norm(VinfE72) - norm(VinfE6));

    double dV41 = fabs(norm(VinfE81) - norm(VinfE1));
    double dV42 = fabs(norm(VinfE82) - norm(VinfE1));

    double test3 = dV31 + dV41;
    double test4 = dV32 + dV42;

    double dV3;
    double dV4;
    bool check2 = false;
    if (test3 >= test4) {
        dV3 = dV32;
        dV4 = dV42;
        check2 = true;
    } else {
        dV3 = dV31;
        dV4 = dV41;
    }
    
    if (flag1) {
        VinfE1.print("E1: ");
        VinfM2.print("M2: ");
        if (check1) {
            VinfM32.print("M3: ");
            VinfM42.print("M4: ");
        } else {
            VinfM31.print("M3: ");
            VinfM41.print("M4: "); 
        }

        VinfM5.print("M5: ");
        VinfE6.print("E6: ");
        if (check2) {
            VinfE72.print("E7: ");
            VinfE82.print("E8: ");
        } else {
            VinfE71.print("E7: ");
            VinfE81.print("E8: ");
        }
    }

    if (flag2) {
        VinfArr[0] = VinfE1;
        VinfArr[1] = VinfM2;
        if (check1) {
            VinfArr[2] = VinfM32;
            VinfArr[3] = VinfM42;
        } else {
            VinfArr[2] = VinfM31;
            VinfArr[3] = VinfM41; 
        }

        VinfArr[4] = VinfM5;
        VinfArr[5] = VinfE6;
        if (check2) {
            VinfArr[6] = VinfE72;
            VinfArr[7] = VinfE82;
        } else {
            VinfArr[6] = VinfE71;
            VinfArr[7] = VinfE81;
        }
    }

    delete Earth;
    delete Mars;

    double dV = dV1 + dV2 + dV3 + dV4;

    return dV;
}

/*double cycleS1L1(double dT1, double phi) {
    double re, rm, muSun, Ve, Vm, Te, Tm, SynodicT, dT4, dThetaM, thetaM, dThetaE, thetaE, dT, dV,
            dV1, dV2, dV3, dV4;
    vector Re1, Rm1, Rm2, Rm3, Re2, Re3, Re5, vEarth, vMars, vEarth2, vEarth3;

    re = 1.495979e8;
    rm = 2.279483e8;

    muSun = 1.32712440e11;

    //Assume CRP Model
    Ve = sqrt(muSun/re);
    Vm = sqrt(muSun/rm);

    Te = 2*PI/sqrt(muSun)*pow(re, 3.0/2.0);
    Tm = 2*PI/sqrt(muSun)*pow(rm, 3.0/2.0);

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

    dThetaM = 2 * PI * (dT1 / Tm);
    thetaM += dThetaM;

    dThetaE = 2 * PI * (dT1 / Te);
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
    dThetaM = 2 * PI * (dT2 / Tm);
    thetaM += dThetaM;

    dThetaE = 2 * PI * (dT2 / Te);
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
}*/
