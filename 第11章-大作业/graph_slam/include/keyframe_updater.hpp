#ifndef KEYFRAME_UPDATER_HPP
#define KEYFRAME_UPDATER_HPP

#include <ros/ros.h>
#include <Eigen/Dense>

#include <keyframe.hpp>

class KeyFrameUpdater
{
public:
    KeyFrameUpdater() : is_first(true), prev_keypose(Eigen::Isometry2d::Identity())
    {
        keyframe_delta_trans = 0.5;
        keyframe_delta_angle = 1.0;
        accum_distance = 0.0;
    }
    // ~KeyFrameUpdater() {}

    bool update(const Eigen::Isometry2d& pose)
    {
        if(is_first)
        {
            is_first = false;
            prev_keypose = pose;
            return true;
        }

        Eigen::Isometry2d delta = prev_keypose.inverse() * pose;
        double dx = delta.translation().norm();
        double da = atan2(delta.matrix().coeff(1,0), delta.matrix().coeff(0,0));

        // std::cout << delta.matrix() << std::endl;

        // std::cout << "dx: " << dx << " da: " << da << std::endl;

        if(dx < keyframe_delta_trans && abs(da) < keyframe_delta_angle)
        {
            return false;
        }

        accum_distance += dx;
        prev_keypose = pose;
        return true;
    }

    double get_accum_distance() const {
        return accum_distance;
    }
private:
    double keyframe_delta_trans;
    double keyframe_delta_angle;

    bool is_first;
    double accum_distance;
    Eigen::Isometry2d prev_keypose;
};


#endif