//
// Created by Daniel Owen on 2019-05-17.
//

#ifndef C___CYCLER_H
#define C___CYCLER_H

#include "helperFuncs.h"
#include "math.h"
#include <cmath>
#include "Lambert_Battin.h"
#include "Lambert_Battin_Multi.h"
#include "Lambert_Multi.h"


#define pi 4*atan(1)

double cycle(double dT1, double dT2, double dT3, double phi);
double cycleMulti(double dT1, double dT2, double dT3, double phi, bool outFlag = false, bool outFlag2 = false);
double cycleS1L1(double dT1, double phi);

#endif //C___CYCLER_H
