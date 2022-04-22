#include <keyframe.hpp>

KeyFrame::KeyFrame(const ros::Time& stamp, const Eigen::Isometry2d& odom, double accum_distance, const pcl::PointCloud<PointT>::ConstPtr& cloud) : stamp(stamp), odom(odom), accum_distance(accum_distance), cloud(cloud), node(nullptr) {}


long KeyFrame::id() const {
    return node->id();
}

SE2 KeyFrame::estimate() const {
    return node->estimate();
}