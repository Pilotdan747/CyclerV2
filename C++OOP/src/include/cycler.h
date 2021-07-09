//
// Created by Daniel Owen on 2019-05-17.
//

#ifndef C___CYCLER_H
#define C___CYCLER_H

#include "helperFuncs.h"
#include "math.h"
#include <cmath>
#include "Lambert_Multi.h"
#include "Orbit.h"

double cycle(double dT1, double dT2, double dT3, double phi);
double cycleMulti(double dT1, double dT2, double dT3, double phi, bool flag1 = false, bool flag2 = false, vec VinfArr[8] = NULL);
double cycleS1L1(double dT1, double phi);

#endif //C___CYCLER_H
