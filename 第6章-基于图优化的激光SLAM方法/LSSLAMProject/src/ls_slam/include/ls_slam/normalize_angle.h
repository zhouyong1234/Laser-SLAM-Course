#ifndef NORMALIZE_ANGLE_H
#define NORMALIZE_ANGLE_H

#include <cmath>
#include "ceres/ceres.h"

template<typename T>
inline T NormalizeAngle(const T& angle_radians)
{
    T two_pi(2.0 * M_PI);
    return angle_radians - two_pi * ceres::floor((angle_radians + T(M_PI)) / two_pi);
}

#endif