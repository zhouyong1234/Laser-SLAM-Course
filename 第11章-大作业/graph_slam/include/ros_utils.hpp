#ifndef ROS_UTILS_HPP
#define ROS_UTILS_HPP

#include <Eigen/Dense>
#include <Eigen/Core>
#include <ros/ros.h>
#include <nav_msgs/Odometry.h>

static Eigen::Isometry2d odom2isometry(const nav_msgs::OdometryConstPtr& odom_msg)
{
    const auto& orientation = odom_msg->pose.pose.orientation;
    const auto& position = odom_msg->pose.pose.position;

    Eigen::Quaterniond quat;
    quat.w() = orientation.w;
    quat.x() = orientation.x;
    quat.y() = orientation.y;
    quat.z() = orientation.z;

    Eigen::Isometry2d isometry = Eigen::Isometry2d::Identity();
    Eigen::Vector3d euler = quat.toRotationMatrix().eulerAngles(2,1,0);

    // std::cout << euler  << std::endl;
    // std::cout << quat.toRotationMatrix() << std::endl;

    isometry.linear() << cos(euler[0]), -sin(euler[0]),
                         sin(euler[0]), cos(euler[0]);

    // std::cout << isometry.linear() << std::endl; 

    isometry.translation() = Eigen::Vector2d(position.x, position.y);
    return isometry;
}

#endif