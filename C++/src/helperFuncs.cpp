//
// Created by Daniel Owen on 2019-05-15.
//

#include "include/helperFuncs.h"
#include <cmath>
#include <stdio.h>

double norm(vector v) {
    return sqrt(pow(v.x, 2) + pow(v.y, 2) + pow(v.z, 2));
}

vector vinf(vector v, vector vPlanet) {
    vector vInf;

    vInf.x = v.x - vPlanet.x;
    vInf.y = v.y - vPlanet.y;
    vInf.z = v.z - vPlanet.z;

    return vInf;
}

vector cross(vector a, vector b) {
    vector c;

    c.x = a.y * b.z - a.z * b.y;
    c.y = a.z * b.x - a.x * b.z;
    c.z = a.x * b.y - a.y * b.x;

    return c;
}

double dot(vector a, vector b) {
    double c = 0;

    c += a.x*b.x;
    c += a.y*b.y;
    c += a.z*b.z;

    return c;
}

void lambert(vector R1, vector R2, double dT, double mu, int k, vector V[2]) {

    double r1 = norm(R1);
    double r2 = norm(R2);

    vector c = cross(R1, R2);

    double dTheta;

    if (k > 0) {
        if (c.z > 0) {
            dTheta = acos(dot(R1, R2) / (r1 * r2));
        } else {
            dTheta = 2 * pi - acos(dot(R1, R2) / (r1 * r2));
        }
    } else {
        if (c.z > 0) {
            dTheta = 2 * pi - acos(dot(R1, R2) / (r1 * r2));
        } else {
            dTheta = acos(dot(R1, R2) / (r1 * r2));
        }
    }
    double A = sin(dTheta)*sqrt(r1*r2/(1 - cos(dTheta)));

    double Znew, Zold, dZ, C, S, y, F, FPri;
    Zold = pi;
    int count = 0;
    bool biSect = 0;
    do {
        C = 0.5 - Zold/24 + pow(Zold, 2)/720;
        S = 1.0/6.0 - Zold/120 + pow(Zold, 2)/5040;
        y = r1 + r2 + A*(Zold*S - 1)/sqrt(C);

        F = pow(y/C, 3.0/2.0)*S + A*sqrt(y) - sqrt(mu)*dT;
        FPri = pow(y/C, 3.0/2.0)*(1/(2*Zold)*(C - 3*S/(2*C)) + 3*pow(S, 2)/(4*C)) + A/8*(3*S/C*sqrt(y) + A*sqrt(C/y));

        Znew = Zold - F/FPri;
        dZ = fabs(Znew - Zold);
        Zold = Znew;
        count++;
        if (count > 1000) {
            biSect = 1;
        }
    } while (dZ > 0.00001 && count < 1001);

    double Zl, Zu, Z;

    if (biSect) {
        Zl = -4*pi;
        Zu = 4*pow(pi, 2);
        Z = 0;

        while (fabs(F) >= 0.001) {
            C = 0.5 - Z/24 + pow(Z, 2)/720 - pow(Z, 3)/40320;
            S = 1.0/6.0 - Z/120 + pow(Z, 2)/5040 - pow(Z, 3)/362880;
            y = r1 + r2 + A*(Z*S - 1)/sqrt(C);

            F = pow(y/C, 3.0/2.0)*S + A*sqrt(y) - sqrt(mu)*dT;

            if (F > 0) {
                Zu = Z;
            } else {
                Zl = Z;
            }

            Z = (Zl + Zu)/2;
        }
    }

    y = r1 + r2 + A*(Zold*S - 1)/sqrt(C);

    double f, g, gDot;

    f = 1 - y/r1;
    g = A*sqrt(y/mu);
    gDot = 1 - y/r2;

    vector v1;
    v1.x = 1/g*(R2.x - f*R1.x); v1.y = 1/g*(R2.y - f*R1.y); v1.z = 1/g*(R2.z - f*R1.z);

    vector v2;
    v2.x = 1/g*(gDot*R2.x - R1.x); v2.y = 1/g*(gDot*R2.y - R1.y); v2.z = 1/g*(gDot*R2.z - R1.z);

    V[0] = v1;
    V[1] = v2;
}
