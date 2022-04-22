#include "gaussian_newton_method.h"

#include "ros/ros.h"
#include "tf/transform_listener.h"
#include "sensor_msgs/LaserScan.h"
#include "nav_msgs/Path.h"
#include "nav_msgs/Odometry.h"
#include "geometry_msgs/Pose.h"
#include <iostream>
#include <chrono>
#include <string>

using namespace std::chrono;

double GN_NormalizationAngle(double angle);

class FrontEnd
{
public:
    FrontEnd()
    {
        m_laserscanSub = m_nh.subscribe("scan",5,&FrontEnd::rosLaserScanCallback,this);

        m_robotPoseSub = m_nh.subscribe("robot_pose", 5, &FrontEnd::robotPoseCallback, this);

        m_odomPub = m_nh.advertise<nav_msgs::Path>("odom_path",1,true);

        m_robotPosePub = m_nh.advertise<nav_msgs::Path>("robot_pose_path",1,true);

        m_gaussianNewtonPub = m_nh.advertise<nav_msgs::Path>("gaussian_newton_path",1,true);

        m_odomBeforePub = m_nh.advertise<nav_msgs::Odometry>("odom_before",1,true);

        m_odomAfterPub = m_nh.advertise<nav_msgs::Odometry>("odom_after",1,true);
    }

    //单纯的数据类型转换，不进行坐标系转换．
    void ConvertChampionLaserScanToEigenPointCloud(const sensor_msgs::LaserScanConstPtr& msg,
                                                   std::vector<Eigen::Vector2d>& eigen_pts)
    {
        eigen_pts.clear();
        for(int i = 0; i < msg->ranges.size();i++)
        {
            if(msg->ranges[i] < msg->range_min || msg->ranges[i] > msg->range_max)
                continue;

            double angle = msg->angle_min + msg->angle_increment * i;

            double lx = msg->ranges[i] * std::cos(angle);
            double ly = msg->ranges[i] * std::sin(angle);

            if(std::isnan(lx) || std::isinf(ly)||
               std::isnan(ly) || std::isinf(ly))
                continue;

            eigen_pts.push_back(Eigen::Vector2d(lx,ly));
        }
    }

    void PublishPath(ros::Publisher& puber,
                     std::vector<Eigen::Vector3d>& path)
    {
        nav_msgs::Path path_msg;
        path_msg.header.stamp = ros::Time::now();
        path_msg.header.frame_id = "/odom";

        geometry_msgs::PoseStamped pose;
        pose.header.stamp = ros::Time::now();
        pose.header.frame_id = "/odom";
        for(int i = 0; i < path.size();i++)
        {
            Eigen::Vector3d traj_node = path[i];
            pose.pose.position.x = traj_node(0);
            pose.pose.position.y = traj_node(1);
            pose.pose.orientation = tf::createQuaternionMsgFromYaw(traj_node(2));
            path_msg.poses.push_back(pose);
        }

        puber.publish(path_msg);
    }

    void rosLaserScanCallback(const sensor_msgs::LaserScanConstPtr& msg)
    {

        
        static bool isFirstFrame = true;
        Eigen::Vector3d nowPose;
        // std::cout << msg->header.stamp << std::endl;
        if(getOdomPose(msg->header.stamp,nowPose) == false)
        {
            std::cout <<"Failed to get Odom Pose"<<std::endl;
            return ;
        }

        

        if(isFirstFrame == true)
        {
            std::cout <<"First Frame"<<std::endl;
            isFirstFrame = false;

            m_prevLaserPose = nowPose;
            ConvertChampionLaserScanToEigenPointCloud(msg,m_prevPts);

            m_odomPath.push_back(nowPose);
            m_gaussianNewtonPath.push_back(nowPose);

            return ;
        }

        auto pre_odom_pose = m_odomPath.back();
        double delta_dist2 = std::pow(nowPose(0) - pre_odom_pose(0),2) + std::pow(nowPose(1) - pre_odom_pose(1),2);
        double delta_angle = std::fabs(tfNormalizeAngle(nowPose(2) - pre_odom_pose(2)));

        // std::cout << "delta_dist2: " << delta_dist2 << std::endl;
        // std::cout << "delta_angle: " << delta_angle << std::endl;

        // if(delta_angle < tfRadians(5.0))
        // {
        //     return ;
        // }
        //数据类型转换．
        std::vector<Eigen::Vector2d> nowPts;
        ConvertChampionLaserScanToEigenPointCloud(msg,nowPts);

        time_point<system_clock> start = system_clock::now();

        //生成地图
        map_t* map = CreateMapFromLaserPoints(m_prevLaserPose,m_prevPts,0.1);

        time_point<system_clock> end = system_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        // std::cout << "Elapsed time: " << elapsed.count() << "s" << std::endl;

        //进行优化．
        //初始解为上一帧激光位姿+运动增量
        Eigen::Vector3d deltaPose = nowPose - m_odomPath.back();
        deltaPose(2) = GN_NormalizationAngle(deltaPose(2));

        Eigen::Matrix3d R_laser;
        double theta = m_prevLaserPose(2);
        R_laser << cos(theta), -sin(theta), 0, 
                   sin(theta),  cos(theta), 0,
                        0,          0,      1;

        // std::cout << "R_laser: " << std::endl << R_laser << std::endl;

        Eigen::Matrix3d R_odom;
        theta = m_odomPath.back()(2);
        R_odom << cos(theta), -sin(theta), 0, 
                  sin(theta),  cos(theta), 0,
                       0,          0,      1;

        // std::cout << "R_odom: " << std::endl << R_odom << std::endl;

        // std::cout << "R_laser - R_odom: " << std::endl << R_laser * R_odom.transpose() << std::endl;

        // std::cout << "deltaPose: " << deltaPose << std::endl;

        Eigen::Vector3d finalPose = m_prevLaserPose + R_laser * R_odom.transpose() * deltaPose;
        finalPose(2) = GN_NormalizationAngle(finalPose(2));

        // auto t = std::chrono::system_clock::now();
        // std::time_t tt = std::chrono::system_clock::to_time_t(t);
        // std::string stt = std::ctime(&tt);
        // std::cout << stt;

    
        // std::cout << "Init Pose: " << finalPose.transpose() << std::endl;
        GaussianNewtonOptimization(map,finalPose,nowPts);

        //更新数据．
        m_prevLaserPose = finalPose;
        m_prevPts = nowPts;

        // std::cout <<"Final Pose: "<<finalPose.transpose()<<std::endl<< std::endl;

        
        //释放地图
        map_free(map);

        //保存路径．
        m_odomPath.push_back(nowPose);
        m_gaussianNewtonPath.push_back(finalPose);

        PublishPath(m_odomPub,m_odomPath);
        PublishPath(m_gaussianNewtonPub,m_gaussianNewtonPath);

        publishOdom(m_odomBeforePub, nowPose);
        publishOdom(m_odomAfterPub, finalPose);

        
    }

    void robotPoseCallback(const geometry_msgs::PoseConstPtr &msg)
    {
        // std::cout << "robot pose callback..." << std::endl;

        Eigen::Vector3d robot_pose;
        tf::Quaternion q(msg->orientation.x, msg->orientation.y, msg->orientation.z, msg->orientation.w);
        tf::Matrix3x3 m(q);
        double roll, pitch, yaw;
        m.getRPY(roll, pitch, yaw);
        robot_pose << msg->position.x, msg->position.y, yaw;

        m_robotPosePath.push_back(robot_pose);
        PublishPath(m_robotPosePub, m_robotPosePath);
    }

    void publishOdom(ros::Publisher& puber, Eigen::Vector3d& pose)
    {
        Eigen::Quaterniond quaternion;
        quaternion = Eigen::AngleAxisd(pose(2), Eigen::Vector3d::UnitZ()) * Eigen::AngleAxisd(0, Eigen::Vector3d::UnitY()) * Eigen::AngleAxisd(0, Eigen::Vector3d::UnitX());

        nav_msgs::Odometry odom_msg;
        odom_msg.header.stamp = ros::Time::now();
        odom_msg.header.frame_id = "/odom";
        odom_msg.pose.pose.position.x = pose(0);
        odom_msg.pose.pose.position.y = pose(1);
        odom_msg.pose.pose.orientation.x = quaternion.x();
        odom_msg.pose.pose.orientation.y = quaternion.y();
        odom_msg.pose.pose.orientation.z = quaternion.z();
        odom_msg.pose.pose.orientation.w = quaternion.w();
        puber.publish(odom_msg);
    }

    bool getOdomPose(ros::Time t,
                     Eigen::Vector3d& pose)
    {
        // Get the robot's pose
        // tf::Stamped<tf::Pose> ident (tf::Transform(tf::createQuaternionFromRPY(0,0,0),
        //                                            tf::Vector3(0,0,0)), t, "/base_footprint");
        // tf::Stamped<tf::Transform> odom_pose;

        tf::StampedTransform transform;
        transform.setIdentity();

        try
        {
            m_tfListener.waitForTransform("/odom", "/base_link", t, ros::Duration(5.0));
            m_tfListener.lookupTransform("/odom", "/base_link", t, transform);
        }
        catch(tf::TransformException &e)
        {
            ROS_WARN("Failed to compute odom pose, skipping scan (%s)", e.what());
            return false;
        }


        // std::cout << transform << std::endl;
        

        // try
        // {
        //     m_tfListener.transformPose("/odom", ident, odom_pose);
        // }
        // catch(tf::TransformException e)
        // {
        //     ROS_WARN("Failed to compute odom pose, skipping scan (%s)", e.what());
        //     return false;
        // }

        // double yaw = tf::getYaw(odom_pose.getRotation());
        // pose << odom_pose.getOrigin().x(),odom_pose.getOrigin().y(),yaw;

        double yaw = tf::getYaw(transform.getRotation());
        pose << transform.getOrigin().x(), transform.getOrigin().y(), yaw;

        return true;
    }

    ros::NodeHandle m_nh;

    Eigen::Vector3d m_prevLaserPose;

    std::vector<Eigen::Vector2d> m_prevPts;

    std::vector<Eigen::Vector3d> m_odomPath;
    std::vector<Eigen::Vector3d> m_gaussianNewtonPath;
    std::vector<Eigen::Vector3d> m_robotPosePath;


    tf::TransformListener m_tfListener;
    ros::Subscriber m_laserscanSub;
    ros::Subscriber m_robotPoseSub;
    ros::Publisher m_odomPub;
    ros::Publisher m_gaussianNewtonPub;
    ros::Publisher m_robotPosePub;
    ros::Publisher m_odomBeforePub;
    ros::Publisher m_odomAfterPub;
};


int main(int argc, char** argv)
{
    ros::init(argc, argv, "front_end");

    FrontEnd front_end;

    ros::spin();

    return (0);
}

