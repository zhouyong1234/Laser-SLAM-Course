#ifndef LOOP_DETECTOR_HPP
#define LOOP_DETECTOR_HPP

#include <boost/format.hpp>
#include <keyframe.hpp>
#include <graph_slam.hpp>
#include <pcl/registration/registration.h>
#include <pcl/registration/icp.h>
#include <tf/transform_listener.h>

struct Loop
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    using Ptr = std::shared_ptr<Loop>;

    Loop(const KeyFrame::Ptr& key1, const KeyFrame::Ptr& key2, const Eigen::Matrix4f& relative_pose) : key1(key1), key2(key2), relative_pose(relative_pose) {}

public:
    KeyFrame::Ptr key1;
    KeyFrame::Ptr key2;
    Eigen::Matrix4f relative_pose;
};

class LoopDetector
{
public:

    typedef pcl::PointXYZI PointT;

    LoopDetector()
    {
        last_edge_accum_distance = 0.0;
        distance_from_last_edge_thresh = 5.0;
        accum_distance_thresh = 8.0;
        distance_thresh = 5.0;

        fitness_score_max_range = std::numeric_limits<double>::max();
        fitness_score_thresh = 0.5;

        registration.reset(new pcl::IterativeClosestPoint<PointT, PointT>());
        registration->setTransformationEpsilon(0.01);
        registration->setMaximumIterations(64);
        registration->setMaxCorrespondenceDistance(2.5);
    }
    // ~LoopDetector() {}

    std::vector<Loop::Ptr> detect(const std::vector<KeyFrame::Ptr>& keyframes, const std::deque<KeyFrame::Ptr>& new_keyframes)
    {
        std::vector<Loop::Ptr> detected_loops;
        for(const auto& new_keyframe : new_keyframes)
        {
            auto candidates = find_candidates(keyframes, new_keyframe);

            auto loop = matching(candidates, new_keyframe);
            if(loop)
            {
                detected_loops.push_back(loop);
            }

        }

        return detected_loops;
    }

    std::vector<KeyFrame::Ptr> find_candidates(const std::vector<KeyFrame::Ptr>& keyframes, const KeyFrame::Ptr& new_keyframe) const
    {
        if(new_keyframe->accum_distance - last_edge_accum_distance < distance_from_last_edge_thresh)
        {
            return std::vector<KeyFrame::Ptr>();
        }

        std::vector<KeyFrame::Ptr> candidates;
        candidates.reserve(32);

        for(const auto& k : keyframes)
        {
            if(new_keyframe->accum_distance - k->accum_distance < accum_distance_thresh)
            {
                continue;
            }

            const auto& pos1 = k->node->estimate().translation();
            const auto& pos2 = new_keyframe->node->estimate().translation();

            double dist = (pos1.head<2>() - pos2.head<2>()).norm();

            if(dist > distance_thresh)
            {
                continue;
            }

            candidates.push_back(k);
        }

        return candidates;
    }

    Loop::Ptr matching(const std::vector<KeyFrame::Ptr>& candidate_keyframes, const KeyFrame::Ptr& new_keyframe)
    {
        if(candidate_keyframes.empty())
        {
            return nullptr;
        }

        registration->setInputTarget(new_keyframe->cloud);
        double best_score = std::numeric_limits<double>::max();
        KeyFrame::Ptr best_matched;
        Eigen::Matrix4f relative_pose;

        auto t1 = ros::Time::now();

        pcl::PointCloud<PointT>::Ptr aligned(new pcl::PointCloud<PointT>());

        for(const auto& candidate : candidate_keyframes)
        {
            registration->setInputSource(candidate->cloud);
            registration->align(*aligned);

            double score = registration->getFitnessScore(fitness_score_max_range);

            if(!registration->hasConverged() || score > best_score)
            {
                continue;
            }

            best_score = score;
            best_matched = candidate;
            relative_pose = registration->getFinalTransformation();
        }

        auto t2 = ros::Time::now();

        std::cout << "best_score: " << boost::format("%.3f") % best_score << "  time: " << boost::format("%.3f") % (t2 - t1).toSec() << "[sec]" << std::endl;

        if(best_score > fitness_score_thresh)
        {
            std::cout << "loop not found..." << std::endl;
            return nullptr;
        }

        std::cout << "loop found!!!" << std::endl;
        last_edge_accum_distance = new_keyframe->accum_distance;

        return std::make_shared<Loop>(new_keyframe, best_matched, relative_pose);
        
    }

private:
    /* data */
    ros::NodeHandle nh;
    double last_edge_accum_distance;
    double distance_from_last_edge_thresh;
    double accum_distance_thresh;
    double distance_thresh;
    double fitness_score_max_range;
    double fitness_score_thresh;

    pcl::Registration<PointT, PointT>::Ptr registration;
};



#endif