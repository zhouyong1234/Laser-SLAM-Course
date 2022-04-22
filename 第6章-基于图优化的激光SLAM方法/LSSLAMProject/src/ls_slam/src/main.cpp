#include <ros/ros.h>
#include <visualization_msgs/MarkerArray.h>
#include <gaussian_newton.h>
#include <readfile.h>
#include <iostream>
#include <g2o/core/base_vertex.h>
#include <g2o/core/base_binary_edge.h>
// #include <g2o/core/base_unary_edge.h>
// #include <g2o/core/base_fixed_sized_edge.h>
#include <g2o/core/block_solver.h>
#include <g2o/core/optimization_algorithm_levenberg.h>
// #include <g2o/core/optimization_algorithm_gauss_newton.h>
#include <g2o/core/optimization_algorithm_dogleg.h>
#include <g2o/solvers/dense/linear_solver_dense.h>
#include <g2o/solvers/eigen/linear_solver_eigen.h>
#include <g2o/core/robust_kernel_impl.h>
#include <g2o/core/sparse_optimizer.h>
#include <g2o/core/factory.h>
#include <g2o/core/hyper_graph.h>

#include <Eigen/Core>
#include <cmath>

#include <ceres/ceres.h>
#include "glog/logging.h"
#include "gflags/gflags.h"
#include "angle_local_parameterization.h"
#include "pose_graph_2d_error_term.h"

const double GN_PI = 3.1415926;

// double GN_NormalizationAngle(double angle)
// {
//     if(angle > GN_PI)
//         angle -= 2*GN_PI;
//     else if(angle < - GN_PI)
//         angle += 2*GN_PI;
//     return angle;
// }

Eigen::Matrix3d PosetoTrans(Eigen::Vector3d x)
{
    Eigen::Matrix3d T ;
    T << cos(x(2)),-sin(x(2)),x(0),
         sin(x(2)),cos(x(2)),x(1),
         0.0,0.0,1;
         return T;
}


class PoseVertex: public g2o::BaseVertex<3,Eigen::Vector3d>
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
     PoseVertex(){}
    PoseVertex( Eigen::Vector3d org ){
        _org = org;
    }
    virtual void setToOriginImpl()
    {
        _estimate = _org;
    }
    virtual void oplusImpl(const double* update)
    {
        _estimate +=Eigen::Vector3d(update);
        _estimate(2) = GN_NormalizationAngle(_estimate(2));
    }
    virtual bool read(std::istream& in) {}//return true;
    virtual bool write(std::ostream& out) const {}
private:
   Eigen::Vector3d _org;   
};


class  PoseEdge: public g2o::BaseBinaryEdge<3, Eigen::Vector3d  ,PoseVertex  ,PoseVertex>
{
public:
EIGEN_MAKE_ALIGNED_OPERATOR_NEW  
    PoseEdge(){}
    void computeError()
    {
        const PoseVertex* v0 = static_cast<const PoseVertex*>(_vertices[0]);
        const PoseVertex* v1 = static_cast<const PoseVertex*>(_vertices[1]);
        const Eigen::Vector3d xi = v0->estimate();
        const Eigen::Vector3d xj = v1->estimate();
        const Eigen::Vector3d z = _measurement;
        Eigen::Matrix3d Ti = PosetoTrans(xi);
        Eigen::Matrix3d Tj = PosetoTrans(xj);
        Eigen::Matrix3d Tij = Ti.transpose()*Tj;
        _Ri  = Ti.block(0,0,2,2);
        _Rj = Tj.block(0,0,2,2);
        _Rij = Tij.block(0,0,2,2);
        _Ri_t = Ti.transpose().block(0,0,2,2);
        _Rj_t = Tj.transpose().block(0,0,2,2);
        _Rij_t = Tij.transpose().block(0,0,2,2);
        _ti =Eigen::Vector2d(xi(0),xi(1));
        _tj =Eigen::Vector2d(xj(0),xj(1));
        _tij=Eigen::Vector2d(z(0),z(1));
        _d_Ri_T = Ti.transpose().block(0,0,2,2);
        //  _Ri << cos(xi(2)),-sin(xi(2)),
        //         sin(xi(2)),cos(xi(2));
        //  _Rj << cos(xj(2)),-sin(xj(2)),
        //         sin(xj(2)),cos(xj(2));
        //  _Rij << cos(z(2)),-sin(z(2)),
        //         sin(z(2)),cos(z(2));  
        // _detc_Ri_T <<   -sin(xi(2)),cos(xi(2)),
        //        -cos(xi(2)),-sin(xi(2)); 
        //  _ti << xi(0),xi(1);
        //  _tj << xj(0),xj(1);
        //  _tij << z(0),z(1);
        Eigen::Vector2d error = _Rij_t*(_Ri_t*(_tj-_ti)-_tij);
        double col =GN_NormalizationAngle(xj(2)-xi(2)-z(2));
        _error  = Eigen::Vector3d(error(0),error(1),col);
    }
    void linearizeOplus()   
    {
        _e1 = _Rij_t*_Ri_t;
        _e2 = _Rij_t*_d_Ri_T*(_tj-_ti);
        // _jacobianOplusXi.block(0,0,3,3) = -Eigen::Matrix3d().setIdentity(); 
        // _jacobianOplusXi.block(0,0,2,2) = _e1.block(0,0,2,2);
        // _jacobianOplusXi.block(0,2,2,1) = _e2.block(0,0,2,1);
        // _jacobianOplusXj.block(0,0,3,3) = Eigen::Matrix3d().setIdentity(); 
        // _jacobianOplusXj.block(0,0,2,2) = _e1.block(0,0,2,2);
        _jacobianOplusXi.block(0,0,2,2) = -_e1;//-_Rij_t*_Ri_t;
        _jacobianOplusXi.block(0,2,2,1) = _e2;
        _jacobianOplusXi.block(2,0,1,3)  = Eigen::Vector3d(0.0,0.0,-1.0).transpose();
        _jacobianOplusXj.block(0,0,2,2) = _e1;
        _jacobianOplusXj.block(0,2,2,1) =Eigen::Vector2d(0.0,0.0);
        _jacobianOplusXj.block(2,0,1,3) = Eigen::Vector3d(0.0,0.0,1.0).transpose();
    }
    virtual bool read(std::istream& in)   {}
    virtual bool write(std::ostream& out) const {}
private:
         Eigen::Matrix2d _Ri ,_Ri_t,_Rj,_Rj_t,_Rij ,_Rij_t ,_d_Ri_T,_e1;
         Eigen::Vector2d _ti,_tj,_tij,_e2;
};


void PublishGraphForVisulization(ros::Publisher* pub, std::map<int, Pose2d>* poses, std::vector<Constraint2d>& constraints, int color = 0)
{

    CHECK(poses != NULL);
    if(constraints.empty())
    {
        LOG(INFO) << "No constraints, no problem to optimize.";
        return;
    }

    visualization_msgs::MarkerArray marray;

    //point--red
    visualization_msgs::Marker m;
    m.header.frame_id = "map";
    m.header.stamp = ros::Time::now();
    m.id = 0;
    m.ns = "ls-slam";
    m.type = visualization_msgs::Marker::SPHERE;
    m.pose.position.x = 0.0;
    m.pose.position.y = 0.0;
    m.pose.position.z = 0.0;
    m.scale.x = 0.1;
    m.scale.y = 0.1;
    m.scale.z = 0.1;

    if(color == 0)
    {
        m.color.r = 1.0;
        m.color.g = 0.0;
        m.color.b = 0.0;
    }
    else
    {
        m.color.r = 0.0;
        m.color.g = 1.0;
        m.color.b = 0.0;
    }

    m.color.a = 1.0;
    m.lifetime = ros::Duration(0);


    //linear--blue
    visualization_msgs::Marker edge;
    edge.header.frame_id = "map";
    edge.header.stamp = ros::Time::now();
    edge.action = visualization_msgs::Marker::ADD;
    edge.ns = "karto";
    edge.id = 0;
    edge.type = visualization_msgs::Marker::LINE_STRIP;
    edge.scale.x = 0.1;
    edge.scale.y = 0.1;
    edge.scale.z = 0.1;

    if(color == 0)
    {
        edge.color.r = 0.0;
        edge.color.g = 0.0;
        edge.color.b = 1.0;
    }
    else
    {
        edge.color.r = 1.0;
        edge.color.g = 0.0;
        edge.color.b = 1.0;
    }
    edge.color.a = 1.0;

    m.action = visualization_msgs::Marker::ADD;
    uint id = 0;

    // add poses
    for(std::map<int, Pose2d>::const_iterator pose_iter = poses->begin(); pose_iter != poses->end(); ++pose_iter)
    {
        const std::map<int, Pose2d>::value_type& pair = *pose_iter;
        m.id = id;
        m.pose.position.x = pair.second.x;
        m.pose.position.y = pair.second.y;
        marray.markers.push_back(visualization_msgs::Marker(m));
        id++;
    }


    // add constraints
    for(std::vector<Constraint2d>::const_iterator constraint_iter = constraints.begin(); constraint_iter != constraints.end(); ++constraint_iter)
    {
        const Constraint2d& constraint = *constraint_iter;
        edge.points.clear();

        std::map<int, Pose2d>::iterator pose_begin_iter = poses->find(constraint.id_begin);


        std::map<int, Pose2d>::iterator pose_end_iter = poses->find(constraint.id_end);

        geometry_msgs::Point p;
        p.x = pose_begin_iter->second.x;
        p.y = pose_begin_iter->second.y;
        edge.points.push_back(p);

        p.x = pose_end_iter->second.x;
        p.y = pose_end_iter->second.y;
        edge.points.push_back(p);
        edge.id = id;

        marray.markers.push_back(visualization_msgs::Marker(edge));
        id++;
    }

    pub->publish(marray);

}


//for visual
void PublishGraphForVisulization(ros::Publisher* pub,
                                 std::vector<Eigen::Vector3d>& Vertexs,
                                 std::vector<Edge>& Edges,
                                 int color = 0)
{
    visualization_msgs::MarkerArray marray;

    //point--red
    visualization_msgs::Marker m;
    m.header.frame_id = "map";
    m.header.stamp = ros::Time::now();
    m.id = 0;
    m.ns = "ls-slam";
    m.type = visualization_msgs::Marker::SPHERE;
    m.pose.position.x = 0.0;
    m.pose.position.y = 0.0;
    m.pose.position.z = 0.0;
    m.scale.x = 0.1;
    m.scale.y = 0.1;
    m.scale.z = 0.1;

    if(color == 0)
    {
        m.color.r = 1.0;
        m.color.g = 0.0;
        m.color.b = 0.0;
    }
    else
    {
        m.color.r = 0.0;
        m.color.g = 1.0;
        m.color.b = 0.0;
    }

    m.color.a = 1.0;
    m.lifetime = ros::Duration(0);

    //linear--blue
    visualization_msgs::Marker edge;
    edge.header.frame_id = "map";
    edge.header.stamp = ros::Time::now();
    edge.action = visualization_msgs::Marker::ADD;
    edge.ns = "karto";
    edge.id = 0;
    edge.type = visualization_msgs::Marker::LINE_STRIP;
    edge.scale.x = 0.1;
    edge.scale.y = 0.1;
    edge.scale.z = 0.1;

    if(color == 0)
    {
        edge.color.r = 0.0;
        edge.color.g = 0.0;
        edge.color.b = 1.0;
    }
    else
    {
        edge.color.r = 1.0;
        edge.color.g = 0.0;
        edge.color.b = 1.0;
    }
    edge.color.a = 1.0;

    m.action = visualization_msgs::Marker::ADD;
    uint id = 0;

    //加入节点
    for (uint i=0; i<Vertexs.size(); i++)
    {
        m.id = id;
        m.pose.position.x = Vertexs[i](0);
        m.pose.position.y = Vertexs[i](1);
        marray.markers.push_back(visualization_msgs::Marker(m));
        id++;
    }

    //加入边
    for(int i = 0; i < Edges.size();i++)
    {
        Edge tmpEdge = Edges[i];
        edge.points.clear();

        geometry_msgs::Point p;
        p.x = Vertexs[tmpEdge.xi](0);
        p.y = Vertexs[tmpEdge.xi](1);
        edge.points.push_back(p);

        p.x = Vertexs[tmpEdge.xj](0);
        p.y = Vertexs[tmpEdge.xj](1);
        edge.points.push_back(p);
        edge.id = id;

        marray.markers.push_back(visualization_msgs::Marker(edge));
        id++;
    }

    pub->publish(marray);
}

bool CeresOptimize(const int maxIteration)
{
    for(int i = 0; i < maxIteration; i++)
    {
        ceres::Problem problem;
    }
}


void BuildOptimizationProblem(const std::vector<Constraint2d>& constraints, std::map<int, Pose2d>* poses, ceres::Problem* problem)
{
    CHECK(poses != NULL);
    CHECK(problem != NULL);
    if(constraints.empty())
    {
        LOG(INFO) << "No constraints, no problem to optimize.";
        return;
    }

    ceres::LossFunction* loss_function = NULL;
    // 更新角度 yaw_new = yaw + delta_yaw
    ceres::LocalParameterization* angle_local_parameterization = AngleLocalParameterization::Create();

    for(std::vector<Constraint2d>::const_iterator constraints_iter = constraints.begin(); constraints_iter != constraints.end(); ++ constraints_iter)
    {
        const Constraint2d& constraint = *constraints_iter;

        std::map<int, Pose2d>::iterator pose_begin_iter = poses->find(constraint.id_begin);
        CHECK(pose_begin_iter != poses->end())
            << "Pose with ID: " << constraint.id_begin << " not found.";
        std::map<int, Pose2d>::iterator pose_end_iter = poses->find(constraint.id_end);
        CHECK(pose_end_iter != poses->end())
            << "Pose with ID: " << constraint.id_end << " not found.";

        const Eigen::Matrix3d sqrt_information = constraint.information.llt().matrixL();

        ceres::CostFunction* cost_function = PoseGraph2dErrorTerm::Create(constraint.x, constraint.y, constraint.yaw_radians, sqrt_information);

        problem->AddResidualBlock(cost_function, loss_function, &pose_begin_iter->second.x, &pose_begin_iter->second.y, &pose_begin_iter->second.yaw_radians, &pose_end_iter->second.x, &pose_end_iter->second.y, &pose_end_iter->second.yaw_radians);

        problem->SetParameterization(&pose_begin_iter->second.yaw_radians, angle_local_parameterization);
        problem->SetParameterization(&pose_end_iter->second.yaw_radians, angle_local_parameterization);
    }

    std::map<int, Pose2d>::iterator pose_start_iter = poses->begin();

    CHECK(pose_start_iter != poses->end()) << "There are no poses.";

    // 固定初始位姿
    problem->SetParameterBlockConstant(&pose_start_iter->second.x);
    problem->SetParameterBlockConstant(&pose_start_iter->second.y);
    problem->SetParameterBlockConstant(&pose_start_iter->second.yaw_radians);
}

bool SolveOptimizationProblem(ceres::Problem* problem)
{
    CHECK(problem != NULL);

    ceres::Solver::Options options;
    options.max_num_iterations = 100;
    options.linear_solver_type = ceres::SPARSE_NORMAL_CHOLESKY;

    ceres::Solver::Summary summary;
    ceres::Solve(options, problem, &summary);

    // std::cout << summary.FullReport() << "\n";
    std::cout << "finish" << std::endl;

    return summary.IsSolutionUsable();
}


bool G2oOptimize(std::vector<Eigen::Vector3d>& Vertexs, std::vector<Edge>& Edges, const int maxIteration)
{
    typedef g2o::BlockSolver< g2o::BlockSolverTraits<3,3> >Block;
    Block::LinearSolverType* linearSolverType = new g2o::LinearSolverEigen<Block::PoseMatrixType>();
    Block* solver_ptr = new Block(std::unique_ptr<Block::LinearSolverType>(linearSolverType));

    g2o::OptimizationAlgorithmLevenberg* solver = new g2o::OptimizationAlgorithmLevenberg(std::unique_ptr<Block>(solver_ptr));
    g2o::SparseOptimizer optimizer;
    optimizer.setAlgorithm(solver);
    // optimizer.setVerbose(true);
    std::vector< PoseVertex* > vertice;
    for(int i = 0; i < Vertexs.size(); i++)
    {
        PoseVertex* v = new PoseVertex(Vertexs[i]);
        if(i == 0)
        v->setFixed(true);
        v->setEstimate(Vertexs[i]);
        v->setId(i);
        optimizer.addVertex(v);
        vertice.push_back(v);
    }

    for(int i = 0; i < Edges.size(); i++)
    {
        PoseEdge* e = new PoseEdge();
        e->setId(i);
        e->setVertex(0, vertice[Edges[i].xi]);
        e->setVertex(1, vertice[Edges[i].xj]);
        e->setMeasurement(Edges[i].measurement);
        e->setInformation(Edges[i].infoMatrix);
        optimizer.addEdge(e);
    }

    optimizer.initializeOptimization();
    optimizer.optimize(maxIteration);

    for(int i = 0; i < vertice.size(); i++)
        Vertexs[i] = vertice[i]->estimate();
    std::cout << "finish" << std::endl;
    return true;
}



int main(int argc, char **argv)
{
    ros::init(argc, argv, "ls_slam");

    ros::NodeHandle nodeHandle;

    // beforeGraph
    ros::Publisher beforeGraphPub,afterGraphPub;
    beforeGraphPub = nodeHandle.advertise<visualization_msgs::MarkerArray>("beforePoseGraph",1,true);
    afterGraphPub  = nodeHandle.advertise<visualization_msgs::MarkerArray>("afterPoseGraph",1,true);


    std::string VertexPath = "/home/touchair/Downloads/HW6/LSSLAMProject/src/ls_slam/data/test_quadrat-v.dat";
    std::string EdgePath = "/home/touchair/Downloads/HW6/LSSLAMProject/src/ls_slam/data/test_quadrat-e.dat";

    std::string originalPath = "/home/touchair/Downloads/HW6/LSSLAMProject/src/ls_slam/data/poses_original.txt";
    std::string optimizedPath = "/home/touchair/Downloads/HW6/LSSLAMProject/src/ls_slam/data/poses_optimized.txt";


    /***************************************************************************************************
     * Ceres Optimization
    ****************************************************************************************************/
    #if 1
    std::map<int, Pose2d> poses;
    std::vector<Constraint2d> constraints;

    ReadPoses(VertexPath, poses);
    ReadConstraints(EdgePath, constraints);

    PublishGraphForVisulization(&beforeGraphPub, &poses, constraints);


    double initError = ComputeError(poses, constraints);
    std::cout << "initError: " << initError << std::endl;


    // std::cout << "Number of poses: " << poses.size() << '\n';
    // std::cout << "Number of constraints: " << constraints.size() << '\n';

    if(!OutputPoses(originalPath, poses))
    {
        std:cerr << "Error writing file..." << "\n";
    }

    ceres::Problem problem;
    BuildOptimizationProblem(constraints, &poses, &problem);

    CHECK(SolveOptimizationProblem(&problem))
        << "The solve was not successful, exiting.";

    if(!OutputPoses(optimizedPath, poses))
    {
        std::cerr << "Error outputting to poses_original.txt";
    }

    PublishGraphForVisulization(&afterGraphPub, &poses, constraints, 1);

    double finalError = ComputeError(poses, constraints);
    std::cout << "finalError: " << finalError << std::endl;

    #else

    std::vector<Eigen::Vector3d> Vertexs;
    std::vector<Edge> Edges;

    ReadVertexInformation(VertexPath,Vertexs);
    ReadEdgesInformation(EdgePath,Edges);

    // for(int i = 0; i < Edges.size(); i++)
    // {
    //     std::cout << Edges[i].infoMatrix << std::endl;
    // }

    PublishGraphForVisulization(&beforeGraphPub,
                                Vertexs,
                                Edges);

    double initError = ComputeError(Vertexs,Edges);
    std::cout <<"initError: "<<initError<<std::endl;

    int maxIteration = 100;
    double epsilon = 1e-4;

    #if 0
    /***************************************************************************************************
     * Gauss Newton Optimization
    ****************************************************************************************************/
    for(int i = 0; i < maxIteration;i++)
    {
        std::cout <<"Iterations:"<<i<<std::endl;
        Eigen::VectorXd dx = LinearizeAndSolve(Vertexs,Edges);

        // std::cout << dx << std::endl;

        //进行更新
        //TODO--Start
        // std::cout << "111111111111111111" << std::endl;
        for(int j = 0; j < Vertexs.size(); j++)
        {
            Vertexs[j] += Eigen::Vector3d(dx(3 * j), dx(3 * j + 1), dx(3 * j + 2));
        }
        //TODO--End

        double maxError = -1;
        for(int k = 0; k < 3 * Vertexs.size();k++)
        {
            if(maxError < std::fabs(dx(k)))
            {
                maxError = std::fabs(dx(k));
            }
        }

        if(maxError < epsilon)
            break;
    }
    #else
    /***************************************************************************************************
     * G2o Optimization
    ****************************************************************************************************/
    if(! G2oOptimize(Vertexs, Edges, maxIteration))
    {
        std::cout << "error" << std::endl;
        return 0;
    }
    #endif


    double finalError  = ComputeError(Vertexs,Edges);

    std::cout <<"FinalError:"<<finalError<<std::endl;

    PublishGraphForVisulization(&afterGraphPub,
                                Vertexs,
                                Edges,1);

    #endif

    ros::spin();

    return 0;
}




