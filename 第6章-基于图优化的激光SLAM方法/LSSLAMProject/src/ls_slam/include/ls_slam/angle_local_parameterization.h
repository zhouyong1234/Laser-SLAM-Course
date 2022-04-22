#ifndef ANGLE_LOCAL_PARAMETERIZATION_H
#define ANGLE_LOCAL_PARAMETERIZATION_H

#include "ceres/local_parameterization.h"
#include "normalize_angle.h"

class AngleLocalParameterization{
public:
    template<typename T>
    bool operator()(const T* theta_radians, const T* delta_theta_radians, T* theta_radians_plus_delta) const {
        *theta_radians_plus_delta = NormalizeAngle(*theta_radians + *delta_theta_radians);
        return true;
    }

    static ceres::LocalParameterization* Create() {
        return (new ceres::AutoDiffLocalParameterization<AngleLocalParameterization, 1, 1>);
    }
};

#endif