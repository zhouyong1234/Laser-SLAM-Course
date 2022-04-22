#include <graph_slam.hpp>
#include <ros_utils.hpp>
#include <keyframe_updater.hpp>
#include <keyframe.hpp>
#include <graph_slam.hpp>
#include <loop_detector.hpp>
#include <iostream>
#include <memory>
#include <boost/format.hpp>
#include <boost/thread.hpp>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>

#include <message_filters/subscriber.h>
#include <message_filters/time_synchronizer.h>
#include <message_filters/sync_policies/approximate_time.h>

#include <nav_msgs/Odometry.h>
#include <sensor_msgs/LaserScan.h>
#include <sensor_msgs/PointCloud2.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl_conversions/pcl_conversions.h>
#include <laser_geometry/laser_geometry.h>
#include <tf/transform_listener.h>



class BackEnd
{
public:
    typedef pcl::PointXYZI PointT;
    typedef message_filters::sync_policies::ApproximateTime<nav_msgs::Odometry, sensor_msgs::PointCloud2> ApproxSyncPolicy;

    BackEnd(ros::NodeHandle& nh)
    {

        trans_odom2map.setIdentity();
        max_keyframes_per_update = 10;

        // tf_listener.setExtrapolationLimit(ros::Duration(0.1));

        graph_slam = std::make_shared<GraphSlam>();
        keyframe_updater = std::make_shared<KeyFrameUpdater>();
        loop_detector = std::make_shared<LoopDetector>();

        // graph_slam.reset(new GraphSlam("lm_var"));
        // keyframe_updater.reset(new KeyFrameUpdater());
        // loop_detector.reset(new LoopDetector());


        point_cloud_pub = nh.advertise<sensor_msgs::PointCloud2>("/cloud", 10, false);

        // scan -> point_cloud
        scan_sub = nh.subscribe<sensor_msgs::LaserScan>("/scan", 10, &BackEnd::scan_callback, this);

        odom_sub = std::make_shared<message_filters::Subscriber<nav_msgs::Odometry> >(nh, "/odom_after", 10);
        cloud_sub = std::make_shared<message_filters::Subscriber<sensor_msgs::PointCloud2> >(nh, "/cloud", 10);
        sync = std::make_shared<message_filters::Synchronizer<ApproxSyncPolicy> >(ApproxSyncPolicy(10), *odom_sub, *cloud_sub);
        sync->registerCallback(boost::bind(&BackEnd::cloud_callback, this, _1, _2));

        // odom_sub.reset(new message_filters::Subscriber<nav_msgs::Odometry>(nh, "/odom", 256));
        // cloud_sub.reset(new message_filters::Subscriber<sensor_msgs::PointCloud2>(nh, "/cloud", 32));
        // sync.reset(new message_filters::Synchronizer<ApproxSyncPolicy>(ApproxSyncPolicy(32), *odom_sub, *cloud_sub));
        // sync->registerCallback(boost::bind(&BackEnd::cloud_callback, this, _1, _2));

        double graph_update_interval = 3.0;

        optimization_timer = nh.createWallTimer(ros::WallDuration(graph_update_interval), &BackEnd::optimization_timer_callback, this);

    }

private:

    void scan_callback(const sensor_msgs::LaserScan::ConstPtr& scan)
    {
        sensor_msgs::PointCloud2 cloud;
        projector.transformLaserScanToPointCloud("laser", *scan, cloud, tf_listener);

        point_cloud_pub.publish(cloud);

    }

    void cloud_callback(const nav_msgs::OdometryConstPtr& odom_msg, const sensor_msgs::PointCloud2::ConstPtr& cloud_msg)
    {
        const ros::Time& stamp = cloud_msg->header.stamp;
        Eigen::Isometry2d odom = odom2isometry(odom_msg);

        pcl::PointCloud<PointT>::Ptr cloud(new pcl::PointCloud<PointT>());
        pcl::fromROSMsg(*cloud_msg, *cloud);

        // std::cout << "cloud callback..." << std::endl;

        if(!keyframe_updater->update(odom))
        {
            // std::cout << "keyframe not updated..." << std::endl;
            return;
        }

        double accum_d = keyframe_updater->get_accum_distance();

        KeyFrame::Ptr keyframe(new KeyFrame(stamp, odom, accum_d, cloud));

        keyframe_queue.push_back(keyframe);

        // std::cout << keyframe_queue.size() << std::endl;
    }

    /**
     * @brief this methods adds all the data in the queues to the pose graph, adn then optimized the pose graph
     * 
     * @param event 
     */
    void optimization_timer_callback(const ros::WallTimerEvent& event)
    {

        // std::cout << "optimization callback" << std::endl;
        // add keyframes to the pose graph

        bool keyframe_updated = flush_keyframe_queue();
        
        if(!keyframe_updated)
        {
            return;
        }

        // std::cout << "loop detection" << std::endl;

        // loop detection
        std::vector<Loop::Ptr> loops = loop_detector->detect(keyframes, new_keyframes);

        // std::cout << loops.size() << std::endl;
        // std::cout << keyframes.size() << std::endl;

        for(const auto& loop : loops)
        {

            matrixAsTransform(loop->relative_pose, transform);

            double rel_x = transform.getOrigin().x();
            double rel_y = transform.getOrigin().y();
            double rel_theta = transform.getRotation().getAngle();

            const SE2& rel_pose = SE2(rel_x, rel_y, rel_theta);

            Eigen::Matrix3d information = Eigen::Matrix3d::Identity();
            information.topLeftCorner(2,2).array() /= 0.5;
            information.bottomRightCorner(1,1).array() /= 0.1;

            auto edge = graph_slam->add_se2_edge(loop->key1->node, loop->key2->node, rel_pose, information);
        }

        std::copy(new_keyframes.begin(), new_keyframes.end(), std::back_inserter(keyframes));
        new_keyframes.clear();

        std::cout << "begin to optimize..." << std::endl;

        // optimize the pose graph
        int num_iterations = 1024;
        graph_slam->optimize(num_iterations);

        std::cout << "optimize finished..." << std::endl;

        // publish tf
        const auto& keyframe = keyframes.back();

        // std::cout << "keyframe node estimate: " << std::endl << keyframe->node->estimate().rotation().toRotationMatrix() << std::endl;
        // std::cout << "keyframe odom: " << std::endl << keyframe->odom.translation() << std::endl;

        Eigen::Isometry2d trans_estimate = Eigen::Isometry2d::Identity();

        trans_estimate.linear() = keyframe->node->estimate().rotation().toRotationMatrix();
        trans_estimate.translation() = keyframe->node->estimate().translation();

        // Eigen::Isometry2d trans = trans_estimate * keyframe->odom.inverse();

        std::cout << "trans_estimate: " << std::endl << trans_estimate.matrix() << std::endl;
        std::cout << "keyframe_odom: " << std::endl << keyframe->odom.matrix() << std::endl;

        Eigen::Matrix3d trans = trans_estimate.matrix() * keyframe->odom.matrix().inverse();

        std::cout << "trans: " << std::endl << trans << std::endl;


        // trans_odom2map = trans.matrix().cast<double>();

        // std::cout << "trans_odom2map: " << std::endl << trans_odom2map << std::endl;

    }

    /**
     * @brief adds all the keyframes to the pose graph
     * 
     * @return true 
     * @return false 
     */
    bool flush_keyframe_queue()
    {
        if(keyframe_queue.empty())
        {
            return false;
        }

        Eigen::Isometry2d odom2map(trans_odom2map.cast<double>());

        int num_processed = 0;
        for(int i = 0; i < std::min<int>(keyframe_queue.size(), max_keyframes_per_update); i++)
        {
            num_processed = i;
            const auto& keyframe = keyframe_queue[i];
            new_keyframes.push_back(keyframe);

            // add pose node
            Eigen::Isometry2d odom = odom2map * keyframe->odom;

            double x = odom.translation()[0];
            double y = odom.translation()[1];
            double theta = atan2(odom.matrix().coeff(1,0), odom.matrix().coeff(0,0));

            const SE2& t = SE2(x,y,theta);
            keyframe->node = graph_slam->add_se2_node(t);

            // fix the first node
            if(keyframes.empty() && new_keyframes.size() == 1)
            {
                Eigen::Matrix3d inf = Eigen::Matrix3d::Identity();
                anchor_node = graph_slam->add_se2_node(SE2(0,0,0));
                anchor_node->setFixed(true);
                anchor_edge = graph_slam->add_se2_edge(anchor_node, keyframe->node, SE2(0,0,0), inf);
            }

            if(i == 0 && keyframes.empty())
            {
                continue;
            }

            // add edge between consecutive keyframes

            const auto& prev_keyframe = i == 0 ? keyframes.back() : keyframe_queue[i-1];

            // relative_pose between two frames
            Eigen::Isometry2d relative_pose = keyframe->odom.inverse() * prev_keyframe->odom;

            double rel_x = relative_pose.translation()[0];
            double rel_y = relative_pose.translation()[1];
            double rel_theta = atan2(relative_pose.matrix().coeff(1,0), relative_pose.matrix().coeff(0,0));
            
            const SE2& rel_pose = SE2(rel_x, rel_y, rel_theta);

            Eigen::Matrix3d information = Eigen::Matrix3d::Identity();
            information.topLeftCorner(2,2).array() /= 0.5;
            information.bottomRightCorner(1,1).array() /= 0.1;

            EdgeSE2* edge = graph_slam->add_se2_edge(keyframe->node, prev_keyframe->node, rel_pose, information);

        }

        keyframe_queue.erase(keyframe_queue.begin(), keyframe_queue.begin() + num_processed + 1);

        return true;

    }

private:
    
    ros::WallTimer optimization_timer;

    std::shared_ptr<message_filters::Subscriber<nav_msgs::Odometry>> odom_sub;
    std::shared_ptr<message_filters::Subscriber<sensor_msgs::PointCloud2>> cloud_sub;
    std::shared_ptr<message_filters::Synchronizer<ApproxSyncPolicy>> sync;
    
    // std::unique_ptr<message_filters::Subscriber<nav_msgs::Odometry>> odom_sub;
    // std::unique_ptr<message_filters::Subscriber<sensor_msgs::PointCloud2>> cloud_sub;
    // std::unique_ptr<message_filters::Synchronizer<ApproxSyncPolicy>> sync;

    std::shared_ptr<GraphSlam> graph_slam;
    std::shared_ptr<KeyFrameUpdater> keyframe_updater;
    std::shared_ptr<LoopDetector> loop_detector;

    // std::unique_ptr<GraphSlam> graph_slam;
    // std::unique_ptr<KeyFrameUpdater> keyframe_updater;
    // std::unique_ptr<LoopDetector> loop_detector;

    VertexSE2* anchor_node;
    EdgeSE2* anchor_edge;

    std::deque<KeyFrame::Ptr> keyframe_queue;
    std::deque<KeyFrame::Ptr> new_keyframes;
    std::vector<KeyFrame::Ptr> keyframes;

    Eigen::Matrix3d trans_odom2map;

    int max_keyframes_per_update;

    ros::Subscriber scan_sub;
    ros::Publisher point_cloud_pub;

    laser_geometry::LaserProjection projector;
    
    tf::TransformListener tf_listener;

    tf::Transform transform;

private:
    inline void matrixAsTransform(const Eigen::Matrix4f& out_mat, tf::Transform& bt)
    {
        double mv[12];
        mv[0] = out_mat(0,0);
        mv[4] = out_mat(0,1);
        mv[8] = out_mat(0,2);
        mv[1] = out_mat(1,0);
        mv[5] = out_mat(1,1);
        mv[9] = out_mat(1,2);
        mv[2] = out_mat(2,0);
        mv[6] = out_mat(2,1);
        mv[10] = out_mat(2,2);
        tf::Matrix3x3 basis;
        basis.setFromOpenGLSubMatrix(mv);
        tf::Vector3 origin(out_mat(0,3), out_mat(1,3), out_mat(2,3));
        bt = tf::Transform(basis, origin);
    }

};

int main(int argc, char** argv)
{
    ros::init(argc, argv, "back_end");
    ros::NodeHandle nh;

    BackEnd back_end(nh);

    ros::spin();

    return (0);
}
