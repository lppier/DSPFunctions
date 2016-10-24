//
// Created by pier on 9/21/16.
//

#ifndef DSPFUNCTIONS_MAPPINGFUNCTIONS_H
#define DSPFUNCTIONS_MAPPINGFUNCTIONS_H

#include "DSPFunctions.hpp"

typedef struct double3
{
    double x;
    double y;
    double z;
} double3;

double3 geodetic2geocentric(double lat, double lon, double alt);
double3 geocentric2lv(double xp, double yp, double zp, double ref_lat, double ref_lon, double height);
#endif //DSPFUNCTIONS_MAPPINGFUNCTIONS_H
