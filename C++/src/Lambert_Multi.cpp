//
// Created by Daniel Owen on 2/19/21
// Using lambert_multi.m from Dr. Brian Kaplinger
//

#include "include/Lambert_Multi.h"
#include <iostream>
//#include </usr/local/include/armadillo>

//using namespace arma;

void lambert_multi(vector R1, vector R2, double tof, int m, double mu, vector V[2]) {
    double r1 = norm(R1);
    double r2 = norm(R2);

    double th = acos(dot(R1, R2)/r1/r2);

    vector z; z.z = 1;
    if (dot(cross(R1, R2), z) < 0) {
        th = 2*pi - th;
    }

    int s = 0;

    if (th < pi) {
        s = 1;
    } else {
        s = -1;
    }

    s = s*sgn(tof);

    double tf = s*tof/86400;

    mlambert(R1, R2, tf, m, mu, V);
}
    
void mlambert(vector R1, vector R2, double tf, int m, double mu, vector Vel[2]) {
    double tol = 1e-14; bool bad = 0; int days = 86400;

    double r1norm = norm(R1); vector r1; r1.x = R1.x/r1norm; r1.y = R1.y/r1norm; r1.z = R1.z/r1norm;
    double V = sqrt(mu/r1norm); vector r2; r2.x = R2.x/r1norm; r2.y = R2.y/r1norm; r2.z = R2.z/r1norm;
    double T = r1norm/V;   tf = tf*days/T;

    double mr2vec = norm(r2);
    double test = dot(r1, r2);
    double dth = acos( fmax(-1.0, fmin(1.0, dot(r1, r2)/mr2vec)));

    int leftbranch = sgn(m);    int longway = sgn(tf);

    m = fabs(m);                tf = fabs(tf);


    if (longway < 0) {
        dth = 2*pi - dth;
    }

    double c = sqrt(1 + pow(mr2vec, 2) - 2*mr2vec*cos(dth));
    double s = (1 + mr2vec + c)/2;
    double amin = s/2;
    double Lambda = sqrt(mr2vec)*cos(dth/2)/s;
    vector crossprd = cross(r1, r2);
    double mcr = norm(crossprd);
    vector nrmunit; nrmunit.x = crossprd.x/mcr; nrmunit.y = crossprd.y/mcr; nrmunit.z = crossprd.z/mcr;

    double logt = log(tf);

    double inn1, inn2, x1, x2;
    if (m==0) {
        inn1 = -0.5233;
        inn2 = 0.5233;
        x1 = log(1 + inn1);
        x2 = log(1 + inn2);

    } else {
        if (leftbranch < 0) {
            inn1 = -0.5234;
            inn2 = -0.2234;
        } else {
            inn1 = 0.7234;
            inn2 = 0.5234;
        }

        x1 = tan(inn1*pi/2);
        x2 = tan(inn2*pi/2);
    }

    double xx[2] = {inn1, inn2};
    double aa[2], bbeta[2], aalfa[2], y12[2];

    for (int i = 0; i < 2; i++) {
        aa[i] = amin/(1 - pow(xx[i], 2));
        bbeta[i] = longway*2*asin(sqrt((s-c)/2/aa[i]));
        aalfa[i] = 2*acos(fmax(-1, fmin(1, xx[i])));
        y12[i] = aa[i]*sqrt(aa[i])*((aalfa[i] - sin(aalfa[i])) - (bbeta[i] - sin(bbeta[i])) + 2*pi*m);
    } 

    double y1, y2;
    if (m==0) {
        y1 = log(y12[0]) - logt;
        y2 = log(y12[1]) - logt;
    } else {
        y1 = y12[0] - tf;
        y2 = y12[1] - tf;
    }

    double err = 10000; int iterations = 0; double xnew = 0;

    double x, a, beta, alfa, tof, ynew;
    while (err > tol) {
        iterations = iterations + 1;
        xnew = (x1*y2 - y1*x2)/(y2-y1);
        
        if (m==0) {
            x = exp(xnew) - 1;
        } else {
            x = atan(xnew)*2/(pi);
        }

        a = amin/(1 - pow(x, 2));

        if (x < 1) {
            beta = longway*2*asin(sqrt((s-c)/2/a));
            alfa = 2*acos(fmax(-1, fmin(1, x)));
        } else {
            alfa = 2*acosh(x);
            beta = longway*2*asinh(sqrt((s-c)/(-2*a)));
        }

        if (a > 0) {
            tof = a*sqrt(a)*((alfa - sin(alfa)) - (beta - sin(beta)) + 2*pi*m);
        } else {
            tof = -a*sqrt(-a)*((sinh(alfa) - alfa) - (sinh(beta) - beta));
        }
        
        if (m == 0) {
            ynew = log(tof) - logt;
        } else {
            ynew = tof - tf;
        }

        x1 = x2; x2 = xnew;
        y1 = y2; y2 = ynew;

        err = fabs(x1 - xnew);

        if (iterations > 15) {
            bad = true;
            break;
        }
    }

    if (bad) {
        lambert_LancasterBlanchard(R1, R2, longway*tf*T, leftbranch*m, mu, Vel);
        return;
    }

    if (m==0) {
        x = exp(xnew) - 1;
    } else {
        x = atan(xnew)*2/(pi);
    }

    a = amin/(1 - pow(x, 2));

    double psi, eta2, eta;
    if (x < 1) {
        beta = longway*2*asin(sqrt((s - c)/2/a));
        alfa = 2*acos(fmax(-1, fmin(1, x)));
        psi = (alfa - beta)/2;
        eta2 = 2*a*pow(sin(psi), 2)/s;
        eta = sqrt(eta2);
    } else {
        beta = longway*2*asinh(sqrt((c - s)/2/a));
        alfa = 2*acosh(x);
        eta2 = -2*a*pow(sinh(psi), 2)/s;
        eta = sqrt(eta2);
    }

    vector ih; ih.x = longway*nrmunit.x; ih.y = longway*nrmunit.y; ih.z = longway*nrmunit.z;
    vector r2n; r2n.x = r2.x/mr2vec;  r2n.y = r2.y/mr2vec; r2n.z = r2.z/mr2vec;

    vector crsprd1 = cross(ih, r1);
    vector crsprd2 = cross(ih, r2n);

    double Vr1 = 1/eta/sqrt(amin)*(2*Lambda*amin - Lambda - x*eta);
    double Vt1 = sqrt(mr2vec/amin/eta2*pow(sin(dth/2), 2));
    double Vt2 = Vt1/mr2vec;
    double Vr2 = (Vt1 - Vt2)/tan(dth/2) - Vr1;

    vector V1; V1.x = (Vr1*r1.x + Vt1*crsprd1.x)*V; V1.y = (Vr1*r1.y + Vt1*crsprd1.y)*V; V1.z = (Vr1*r1.z + Vt1*crsprd1.z)*V;
    vector V2; V2.x = (Vr2*r2n.x + Vt2*crsprd2.x)*V; V2.y = (Vr2*r2n.y + Vt2*crsprd2.y)*V; V2.z = (Vr2*r2n.z + Vt2*crsprd2.z)*V;

    Vel[0] = V1;
    Vel[1] = V2;
}

void lambert_LancasterBlanchard(vector R1, vector R2, double tf, int m, double muC, vector V[2]) {
    double tol = 1e-12;
    double r1 = norm(R1);
    double r2 = norm(R2);
    vector r1unit; r1unit.x = R1.x/r1; r1unit.y = R1.y/r1; r1unit.z = R1.z/r1;
    vector r2unit; r2unit.x = R2.x/r2; r2unit.y = R2.y/r2; r2unit.z = R2.z/r2;
    vector crsprd = cross(R1, R2);
    double mcrsprd = norm(crsprd);
    vector temp; temp.x = crsprd.x/mcrsprd; temp.y = crsprd.y/mcrsprd; temp.z = crsprd.z/mcrsprd;
    vector th1unit = cross(temp, r1unit);
    vector th2unit = cross(temp, r2unit);

    double dth = acos(fmax(-1, fmin(1, (dot(R1, R2)/r1/r2))));

    int longway = sgn(tf); tf = fabs(tf);

    if (longway < 0) {
        dth = dth - 2*pi;
    }

    int leftbranch = sgn(m); m = fabs(m);

    double c = sqrt(pow(r1, 2) + pow(r2, 2) - 2*r1*r2*cos(dth));
    double s = (r1 + r2+ c)/2;
    double T = sqrt(8*muC/pow(s, 3))*tf;
    double q = sqrt(r1*r2)/s*cos(dth/2);

    double Ts[4];
    LancasterBlanchard(0, q, m, Ts);
    double T0 = Ts[0];

    double Td = T0 - T;
    double phr = fmod((2*atan2(1 - pow(q, 2), 2*q)), (2*pi));

    double x01, x0, x02, W, x03, lambda, xMpi, xM0, xM, Tp, Tpp, Tppp, xMp, TM, TmTM, T0mTM, x, x00, Tx, xp;
    if (m==0) {
        x01 = T0*Td/4/T;

        if (Td > 0) {
            x0 = x01;
        } else {
            x01 = Td/(4 - Td);
            x02 = -sqrt(-Td/(T + T0/2));
            W = x01 + 1.7*sqrt(2 - phr/(pi));

            if ( W >= 0) {
                x03 = x01;
            } else {
                x03 = x01 + pow(-W, 1/16)*(x02 - x01);
            }

            lambda = 1 + x03*(1 + x01)/2 - 0.03*pow(x03, 2)*sqrt(1 + x01);
            x0 = lambda*x03;
        }

        if (x0 < -1) {  // shouldn't happen
            //exitflag = -1;
            //return;
        }
    } else {
        xMpi = 4/(3*pi*(2*m + 1));

        if (phr < pi) {
            xM0 = xMpi*pow((phr/(pi)), 1.0/8.0);
        } else if (phr > pi) {
            xM0 = xMpi*(2 - pow((2 - phr/(pi)), 1.0/8.0));
        } else {
            xM0 = 0;
        }

        xM = xM0;  Tp = 10000000; int n = 0;

        while (fabs(Tp) > tol) {
            n++;
            LancasterBlanchard(xM, q, m, Ts);
            Tp = Ts[1]; Tpp = Ts[2]; Tppp = Ts[3];

            xMp = xM;
            xM = xM - 2*Tp*Tpp/(2*pow(Tpp, 2) - Tp*Tppp);

            if (n%7) {
                xM = (xMp + xM)/2;
            }

            if (n > 25) {
                //exitflag = -2;
                //return;
            }
        }

        if ((xM < -1) || (xM > 1)) { //shouldnt happen
            //exitflag = -1;
            //return;
        }

        LancasterBlanchard(xM, q, m, Ts);
        TM = Ts[0];

        if (TM > T) {
            //exitflag = -1;
            //return;
        }

        TmTM = T - TM;   T0mTM = T0 - TM;
        LancasterBlanchard(xM, q, m, Ts);
        Tp = Ts[1]; Tpp = Ts[2];

        if (leftbranch > 0) {
            x = sqrt(TmTM/(Tpp/2 + TmTM/pow((1 - xM), 2)));
            W = xM + x;
            W = 4*W/(4 + TmTM) + pow((1 - W), 2);
            x0 = x*(1 - (1 + m + (dth - 1/2))/(1 + 0.15*m)*x*(W/2 + 0.03*x*sqrt(W))) + xM;

            if (x0 > 1) {
                //exitflag = -1;
                //return;
            }
        } else {
            if (Td > 0) {
                x0 = xM - sqrt(TM/(Tpp/2 - TmTM*(Tpp/2/T0mTM - 1/pow(xM, 2))));
            } else {
                x00 = Td/(4 - Td);
                W = x00 + 1.7*sqrt(2*(1 - phr));

                if (W >= 0) {
                    x03 = x00;
                } else {
                    x03 = x00 - sqrt(pow(-W, 1.0/8.0))*(x00 + sqrt(-Td/(1.5*T0 - Td)));
                }

                W = 4/(4 - Td);
                lambda = (1 + (1 + m + 0.24*(dth - 1.0/2.0))/(1 + 0.15*m)*x03*(W/2 - 0.03*x03*sqrt(W)));
                x0 = x03*lambda;
            }

            if (x0 < -1) {
                //exitflag = -1;
                //return;
            }
        }
    }

    x = x0;   Tx = 100000; int iterations = 0;

    while (fabs(Tx) > tol) {
        iterations++;
        LancasterBlanchard(x, q, m, Ts);
        Tx = Ts[0]; Tp = Ts[1]; Tpp = Ts[2];

        Tx = Tx - T;
        xp = x;
        x = x - 2*Tx*Tp/(2*pow(Tp, 2) - Tx*Tpp);

        if (iterations%7) {
            x = (xp+x)/2;
        }

        if (iterations > 25) {
            //exitflag = -2;
            //return;
        }
    }

    double gamma = sqrt(muC*s/2);
    double sigma, rho, z;
    if (c == 0) {
        sigma = 1;
        rho = 0;
        z = fabs(x);
    } else {
        sigma = 2*sqrt(r1*r2/pow(x, 2))*sin(dth/2);
        rho = (r1 - r2)/c;
        z = sqrt(1 + pow(1, 2)*(pow(x, 2) - 1));
    }

    double Vr1 = gamma*((q*z - x) - rho*(q*z + x))/r1;
    double Vr2 = -gamma*((q*z - x) + rho*(q*z + x))/r2;
    
    double Vtan1 = sigma*gamma*(z + q*x)/r1;
    double Vtan2 = sigma*gamma*(z + q*x)/r2;

    vector V1, V2;

    V1.x = Vr1*r1unit.x + Vtan1*th1unit.x;
    V1.y = Vr1*r1unit.y + Vtan1*th1unit.y;
    V1.z = Vr1*r1unit.z + Vtan1*th1unit.z;

    V2.x = Vr2*r2unit.x + Vtan2*th2unit.x;
    V2.y = Vr2*r2unit.y + Vtan2*th2unit.y;
    V2.z = Vr2*r2unit.z + Vtan2*th2unit.z;

    V[0] = V1;
    V[1] = V2;

    //exitflag = 1;
}

void LancasterBlanchard(double x, double q, int m, double Ts[4]) {
    if (x < -1) { //Should be impossible
        x = fabs(x) - 2;
    } else if (x == -1) {
        x = x + 1e-16;
    }

    double E = pow(x, 2) - 1;
    double T, Tp, Tpp, Tppp;

    if (x==1) {
        T = 4.0/3.0*(1 - pow(q, 3));
        Tp = 4.0/5.0*(pow(q, 5) - 1);
        Tpp = Tp + 120.0/70.0*(1 - pow(q, 7));
        Tppp = 3*(Tpp - Tp) + 2400.0/1080.0*(pow(q, 9));
    } else if (fabs(x - 1) < 1e-2) {
        double sig1, sig2, dsigdx1, dsigdx2, d2sigdx21, d2sigdx22, d3sigdx31, d3sigdx32;
        
        double sigs[4], sigs2[4];
        sigmax(-E, sigs);
        sig1 = sigs[0]; dsigdx1 = sigs[1]; d2sigdx21 = sigs[2]; d3sigdx31 = sigs[3];

        sigmax(-E*pow(q, 2), sigs2);
        sig2 = sigs2[0]; dsigdx2 = sigs2[1]; d2sigdx22 = sigs2[2]; d3sigdx32 = sigs2[3];

        T = sig1 - pow(q, 3)*sig2;
        Tp = 2*x*(pow(q, 5)*dsigdx2 - dsigdx1);
        Tpp = Tp/x + 4*pow(x, 2)*(d2sigdx21 - pow(q, 7)*d2sigdx22);
        Tppp = 3*(Tpp - Tp/x)/x + 8*pow(x, 2)*(pow(q, 9)*d3sigdx32 - d3sigdx31);
    } else {
        double y = sqrt(abs(E));
        double z = sqrt(1 + pow(q, 2)*E);
        double f = y*(z - q*x);
        double g = x*z - q*E;

        double d;
        if (E < 0) {
            d = atan2(f, g) + pi*m;
        } else if (E == 0) {
            d = 0;
        } else {
            d = log(fmax(0, f+g));
        }

        T = 2*(x - q*z - d/y)/E;
        Tp = (4 - 4*pow(q, 3)*x/z - 3*x*T)/E;
        Tpp = (-4*pow(q, 3)/z * (1 - pow(q, 2)*pow(x, 2)/pow(z, 2)) - 3*T - 3*x*Tp)/E;
        Tppp = (4*pow(q, 3)/pow(z, 2)*((1 - pow(q, 2)*pow(x, 2)/pow(z, 2)) + 2*pow(q, 2)*x/pow(z, 2)*(z - x)) - 8*Tp - 7*x*Tpp)/E;
    }

    Ts[0] = T; Ts[1] = Tp; Ts[2] = Tpp; Ts[3] = Tppp;
}

double an[25] = {4.000000000000000e-001, 2.142857142857143e-001, 4.629629629629630e-002, 6.628787878787879e-003, 7.211538461538461e-004, 
                 6.365740740740740e-005, 4.741479925303455e-006, 3.059406328320802e-007, 1.742836409255060e-008, 8.892477331109578e-010, 
                 4.110111531986532e-011, 1.736709384841458e-012, 6.759767240041426e-014, 2.439123386614026e-015, 8.203411614538007e-017,
                 2.583771576869575e-018, 7.652331327976716e-020, 2.138860629743989e-021, 5.659959451165552e-023, 1.422104833817366e-024,
                 3.401398483272306e-026, 7.762544304774155e-028, 1.693916882090479e-029, 3.541295006766860e-031, 7.105336187804402e-033};

void sigmax(double y, double sigs[4]) {
    double sig, dsigdx, d2sigdx2, d3sigdx3;
    double count = 0;   double count2 = 0;   double count3 = 0;   double count4 = 0;
    int var = 0;
    for (int i = 0; i < 25; i++) {
        var = i + 1;
        count += pow(y, var)*an[i];
        count2 += var*pow(y, var - 1)*an[i];
        count3 += var*(var - 1)*pow(y, var - 2)*an[i];
        count4 += var*(var - 1)*(var - 2)*pow(y, var - 3)*an[i];
    }

    sig = 4.0/3.0 + count;
    dsigdx = count2;
    d2sigdx2 = count3;
    d3sigdx3 = count4;


    sigs[0] = sig; sigs[1] = dsigdx; sigs[2] = d2sigdx2; sigs[3] = d3sigdx3;
}