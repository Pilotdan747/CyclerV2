#include "include/Orbit.h"

Orbit::Orbit(double mu) {
    a = 1;
    e = 0;
    inc = 0;
    RAAN = 0;
    argPeri = 0;
    TA = 0;
    this->mu = mu;
}

Orbit::Orbit(double muIn, vec R, vec V) {
    mu = muIn;
    state2OE(R, V);
}

Orbit::Orbit(double muIn, double a, double e, double inc, double RAAN, double argPeri, double TA) {
    mu = muIn;
    this->a = a;
    this->e = e;
    this->inc = inc;
    this->RAAN = RAAN;
    this->argPeri = argPeri;
    this->TA = TA;
}

//Orbit::~Orbit() {
    //delete state;
//}

void Orbit::Kepler(double t) {
    double n = sqrt(mu/pow(a, 3.0));
    double E0 = atan2(sqrt(1 - pow(e, 2.0))*sin(TA), e + cos(TA));
    double M0 = E0 - e*sin(E0);

    double M = M0 + n*t;

    double Enew = 0;
    if (M > pi || (M < 0 && M > -pi))
        Enew = M - e;
    else
        Enew = M + e;


    double E = Enew + 1;
    while (abs(Enew - E) > 1e-12) {
        E = Enew;
        Enew = E + (M - E + e*sin(E))/(1 - e*cos(E));
    }

    double s = sin(E)*sqrt(1 - pow(e, 2.0))/(1 - e*cos(E));
    double c = (cos(E) - e)/(1 - e*cos(E));

    double TA2 = atan2(s, c);
    TA = TA2;
}

void Orbit::state2OE(vec R, vec V) {
    double parabTol = 1e-5;
    vec k = {0, 0, 1};
    vec j = {0, 1, 0};
    vec i = {1, 0, 0};

    vec h = cross(R, V);
    vec n = cross(k, h);

    double E = pow(norm(V), 2.0)/2 - mu/norm(R);

    vec ee = 1/mu*((pow(norm(V), 2.0) - mu/norm(R))*R - dot(R, V)*V);
    e = norm(ee);

    double p = 0;

    if (fabs(E) < parabTol) {
        a = NULL;
    } else {
        a = -mu/(2*E);
    }

    inc = acos(dot(k, h)/norm(h));

    RAAN = acos(dot(i, n)/norm(n)); 
    if (dot(n, j) < 0) 
        RAAN = 2*pi - RAAN; 

    argPeri = acos(dot(n, ee)/(norm(n)*norm(ee))); 
    if (dot(ee, k) < 0) 
        argPeri = 2*pi - argPeri; 

    TA = acos(dot(ee, R)/(norm(ee)*norm(R))); 
    if (dot(R, V) < 0) 
        TA = 2*pi - TA; 
}

vec* Orbit::getState() {
    double h = sqrt(mu*a*(1-pow(e, 2.0)));
    double p = pow(h, 2.0)/mu;
    double r = p/(1 + e*cos(TA));
    double v = sqrt(mu/p);

    vec rPQW = {r*cos(TA), r*sin(TA), 0};
    vec vPQW = {v*-sin(TA), v*(e + cos(TA)), 0};

    mat C1 = {{cos(-argPeri), sin(-argPeri), 0}, {-sin(-argPeri), cos(-argPeri), 0}, {0, 0, 1}};
    mat C2 = {{1, 0, 0}, {0, cos(-inc), sin(-inc)}, {0, -sin(-inc), cos(-inc)}};
    mat C3 = {{cos(-RAAN), sin(-RAAN), 0}, {-sin(-RAAN), cos(-RAAN), 0}, {0, 0, 1}};

    mat C = C3*C2*C1;

    vec R = C*rPQW;
    vec V = C*vPQW;

    state[0] = R;
    state[1] = V;
    return state;
}