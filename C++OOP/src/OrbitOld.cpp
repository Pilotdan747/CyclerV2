#include "include/Orbit.h"
#include </usr/local/include/armadillo>

using namespace arma;

Orbit::Orbit(double mu) {
    a = 1;
    e = 0;
    inc = 0;
    RAAN = 0;
    argPeri = 0;
    TA = 0;
    this->mu = mu;
}

Orbit::Orbit(double muIn, vector R, vector V) {
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

void Orbit::Kepler(double t) {
    double n = sqrt(mu/pow(a, 3.0));
    double E0 = atan2(sqrt(1 - pow(e, 2.0))*sin(TA), e + cos(TA));
    double M0 = E0 - e*sin(E0);

    double M = M0 + n*t;

    double Enew = 0;
    if (M > PI || (M < 0 && M > -PI))
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

void Orbit::state2OE(vector R, vector V) {
    double parabTol = 1e-5;
    vec kk(0, 0, 1);
    vector k; k.x = 0; k.y = 0; k.z = 1;
    vector j; j.x = 0; j.y = 1; j.z = 0;
    vector i; i.x = 1; i.y = 0; i.z = 0;

    vector h = cross(R, V);
    vector n = cross(k, h);
    double E = pow(norm(V), 2.0)/2 - mu/norm(R);

    vector ee;
    ee.x = 1/mu*((pow(norm(V), 2.0) - mu/norm(R))*R.x - dot(R, V)*V.x);
    ee.y = 1/mu*((pow(norm(V), 2.0) - mu/norm(R))*R.y - dot(R, V)*V.y);
    ee.z = 1/mu*((pow(norm(V), 2.0) - mu/norm(R))*R.z - dot(R, V)*V.z);

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
        RAAN = 2*PI - RAAN; 

    argPeri = acos(dot(n, ee)/(norm(n)*norm(ee))); 
    if (dot(ee, k) < 0) 
        argPeri = 2*PI - argPeri; 

    TA = acos(dot(ee, R)/(norm(ee)*norm(R))); 
    if (dot(R, V) < 0) 
        TA = 2*PI - TA; 
}

vector* Orbit::getState() {
    double h = sqrt(mu*a*(1-pow(e, 2.0)));
    double p = pow(h, 2.0)/mu;
    double r = p/(1 + e*cos(TA));
    double v = sqrt(mu/p);

    vector rPQW;
    rPQW.x = r*cos(TA); 
    rPQW.y = r*sin(TA); 
    rPQW.z = 0;

    vector vPQW;
    vPQW.x = v*-sin(TA); 
    vPQW.y = v*(e + cos(TA)); 
    vPQW.z = 0;

    double C1[3][3] = {{cos(-argPeri), sin(-argPeri), 0}, {-sin(-argPeri), cos(-argPeri), 0}, {0, 0, 1}};
    double C2[3][3] = {{1, 0, 0}, {0, cos(-inc), sin(-inc)}, {0, -sin(-inc), cos(-inc)}};
    double C3[3][3] = {{cos(-RAAN), sin(-RAAN), 0}, {-sin(-RAAN), cos(-RAAN), 0}, {0, 0, 1}};

    double C[3][3];
    multiply(C3, C2, C1, C);

    vector R = multiplyVec(C, rPQW);
    vector V = multiplyVec(C, vPQW);

    state[0] = R;
    state[1] = V;
    return state;
}

vector Orbit::multiplyVec(double A[3][3], vector x) {
	vector result;

    result.x = A[0][0]*x.x + A[0][1]*x.y + A[0][2]*x.z;
    result.y = A[1][0]*x.x + A[1][1]*x.y + A[1][2]*x.z;
    result.z = A[2][0]*x.x + A[2][1]*x.y + A[2][2]*x.z;

    return result;
}
         
void Orbit::multiply(double A[3][3], double B[3][3], double C[3][3], double result[3][3]) {
	double inter[3][3];
    int length = 3;

    for (int i = 0; i < length; i++) {
		for (int j = 0; j < length; j++) {
			inter[i][j] = 0;
			for (int k = 0; k < length; k++) {
				inter[i][j] += A[i][k] * B[k][j];
			}
		}
	}

	for (int i = 0; i < length; i++) {
		for (int j = 0; j < length; j++) {
			result[i][j] = 0;
			for (int k = 0; k < length; k++) {
				result[i][j] += inter[i][k] * C[k][j];
			}
		}
	}
}

int main() {
    vector R = {6678.0, 0, 0};
    vector V = {0, 7.725835, 0};

    //Orbit* ob1 = new Orbit(3.986e5, 6378+300, 0, 0, 0, 0, 0);
    Orbit* ob1 = new Orbit(3.986e5, R, V);

    ob1->Kepler(2*PI/sqrt(3.986e5)*pow(6378+300, 3.0/2.0));
    vector *state = ob1->getState();

    vector r = state[0];
    vector v = state[1];

    printf("R = [%f, %f, %f]\n", r.x, r.y, r.z);
    printf("V = [%f, %f, %f]\n", v.x, v.y, v.z);

    return 0;
}