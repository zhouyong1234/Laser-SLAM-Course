#ifndef KEYFRAME_HPP
#define KEYFRAME_HPP

#include <ros/ros.h>
#include <sensor_msgs/LaserScan.h>
#include <pcl/point_types.h>
#include <pcl/point_cloud.h>
#include <graph_slam.hpp>

struct KeyFrame
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    using PointT = pcl::PointXYZI;
    using Ptr = std::shared_ptr<KeyFrame>;

    KeyFrame(const ros::Time& stamp, const Eigen::Isometry2d& odom, double accum_distance, const pcl::PointCloud<PointT>::ConstPtr& cloud);

    long id() const;

    SE2 estimate() const;
public:
    ros::Time stamp;
    Eigen::Isometry2d odom;
    double accum_distance;
    pcl::PointCloud<PointT>::ConstPtr cloud;
    VertexSE2* node;
};




#endif