//
// Created by Daniel Owen on 2019-06-21.
//

#include "include/Lambert_Battin_Multi.h"

double eps_m = 10e-15;

using std::abs;

void lambert_battin_multi(vector R1, vector R2, double dT, double mu, double dir, int N, vector V[2]) {
    double r1, r2, kk, theta, c, s, L, r0p, m, n, v, l, xL, xR, d;


    r1 = norm(R1);
    r2 = norm(R2);

    d = norm(vinf(R1, R2));

    kk = cross(R1, R2).z;

    theta = acos(dot(R1, R2)/r1/r2);

    if ((dir == 0 && kk < 0) || (dir > 0 && kk >= 0)) {
        theta = 2 * pi - theta;
    }

    c = sqrt(pow(r1, 2) + pow(r2, 2) - 2*r1*r2*cos(theta));
    s = 0.5*(r1 + r2 + c);                              // Page 304

    L = sqrt((s-c)/s);                                  // 7.122
    if (theta > pi)
        L = -L;

    r0p = 0.25*s*pow((1+L), 2);                          // 6.77 -> worked out proof

    m = mu*pow(dT, 2)/(8*pow(r0p, 3));                   // 7.89

    n = (r1 + r2)/2;
    v = acos(sqrt(r1*r2)*cos(theta/2)/n);
    l = pow(tan(v/2), 2);

    double yR, yL;

    xR = revSucSub(N, m, l, l, yR);
    xL = sucSub(N, m, l, l, yL);

    printf("xR is: %f\n", xR);
    printf("xL is: %f\n", xL);

    double EL = 2*atan(sqrt(xL));
    double ER = 2*atan(sqrt(xR));

    printf("ER is: %f\n", ER);
    printf("EL is: %f\n", EL);


    // Check these
    //double aL = (-1*sqrt(4*pow(n, 2) + pow(d, 2))* cos(EL) + 2*n)/(2*pow(sin(EL), 2));
    //double aR = (-1*sqrt(4*pow(n, 2) + pow(d, 2))* cos(ER) + 2*n)/(2*pow(sin(ER), 2));

    double aL = mu*pow(dT, 2)/16.0/pow(r0p, 2)/xL/pow(yL, 2);
    double aR = mu*pow(dT, 2)/16.0/pow(r0p, 2)/xR/pow(yR, 2);
    double aR2 = m*s*pow((1+L), 2)/8/xR/pow(yR, 2);
    //          mu*pow(dT, 2)/16.0/pow(r0p, 2)/x /pow(y , 2);

    printf("aR is: %f\n", aR);
    printf("aR2 is: %f\n", aR2);
    printf("aL is: %f\n", aL);

    printf("yR is: %f\n", yR);
    printf("yL is: %f\n", yL);

    double b, amin, tmin, ae, dE, f, g, gdot;

    // Still need to figure out what to do with xL and xR
    if (xL > 0) {
        f = 1 - aL / r1 * (1 - cos(EL));
        g = dT - sqrt(pow(aL, 3) / mu) * (EL - sin(EL));
        gdot = 1 - aL / r2 * (1 - cos(EL));
    } else {
        printf("Hyperbolic orbit");
    }

    /*b = 2 * asin(sqrt(0.5 * (s - c) / aL));
    if (theta > pi) {
        b = -b;
    }

    printf("b is: %f\n", b);
    printf("arg is: %f\n", 0.5 * (s - c) / aR);

    amin = 0.5 * s;
    tmin = sqrt(pow(amin, 3) / mu) * (pi - b + sin(b));
    ae = 2 * asin(sqrt(0.5 * s / aL));

    //printf("amin is: %f\n", amin);
    //printf("tmin is: %f\n", tmin);
    //printf("ae is: %f\n", ae);

    if (dT > tmin) {
        ae = 2 * pi - ae;
    }

    dE = ae - b;
    f = 1 - aL / r1 * (1 - cos(dE));
    g = dT - sqrt(pow(aL, 3) / mu) * (dE - sin(dE));
    gdot = 1 - aL / r2 * (1 - cos(dE));
*/
    printf("f is: %f\n", f);
    printf("g is: %f\n", g);
    printf("gdot is: %f\n", gdot);

    vector V1, V2;

    V1.x = (R2.x - f*R1.x)/g; V1.y = (R2.y - f*R1.y)/g; V1.z = (R2.z - f*R1.z)/g;
    V2.x = (gdot*R2.x - R1.x)/g; V2.y = (gdot*R2.y - R1.y)/g; V2.z = (gdot*R2.z - R1.z)/g;

    V[0] = V1;
    V[1] = V2;
}


// Solving for xL
double sucSub(int N, double m, double l, double x0, double y) {
    double x, yold, E, rhs;

    x = x0;
    y = 0;
    yold = 1;

    for (int i = 0; i < 100; i++) {
        if (abs(yold - y) > eps_m) {
            break;
        }
        E = 2*atan(sqrt(x));

        yold = y;

        rhs = m*(N*pi + E - sin(E))/(4*pow(tan(E/2), 3));

        y = ySolve(rhs);

        if (y > sqrt(m/l)) {
            x = revSucSub(N, m, l, x0, 0)*2;
            continue;
        }

        x = (sqrt((pow(l, 2) - 2*l + 1) * pow(y, 2) + 4*m) - (l + 1)*y)/(2*y);
    }

    return x;
}

// Solving for xR
double revSucSub(int N, double m, double l, double x0, double y1) {
    double x, y1old, y2, q, E0, h, hpri, Enew, E;

    x = x0;
    y1 = 0;
    y1old = 1;

    for (int i = 0; i < 100; i++) {
        if (abs(y1old - y1) < eps_m) {
            break;
        }

        y1old = y1;
        y1 = sqrt(m / ((l + x) * (1 + x)));
        y2 = y1;

        if (y1 < 1) {
            x = sucSub(N, m, l, x0, 0)/2;
            continue;
        }

        E0 = pi;                                                    //Guess
        q = 4 / m * pow(y2, 3) - pow(y2, 2);
        h = (N * pi + E0 - sin(E0)) / pow(tan(E0 / 2), 3) - q;

        if (h < 0) {                                               //Make E0 so that h(E0) > 0
            while (h < 0) {
                E0 = E0 / 2;
                h = (N * pi + E0 - sin(E0)) / pow(tan(E0 / 2), 3) - q;
            }
        }

        Enew = E0;
        E = 0;

        for (int j = 0; j < 100; j ++) {
            if (abs(Enew - E) < eps_m) {
                break;
            }

            E = Enew;
            h = (N * pi + E - sin(E)) / pow(tan(E / 2), 3) - q;
            hpri = -1 * pow(cos(E / 2), 2) * (2 * (cos(E) - 1) * sin(E / 2) * cos(E / 2) - 3 * (sin(E) - E - N * pi)) /
                   (2 * pow(sin(E / 2), 4));

            Enew = E - h / hpri;
        }

        x = pow(tan(Enew/2), 2);
    }


    return x;
}

double ySolve(double c) {
    double num1 = 3*sqrt(3)*sqrt(27*pow(c, 2) + 4*c) + 27*c + 2;

    double term1 = pow(num1, 1.0/3.0)/pow(2, 1.0/3.0);

    double term2 = pow(term1, -1);

    return (term1 + term2 + 1)/3;
}