#ifndef GRAPH_SLAM_HPP
#define GRAPH_SLAM_HPP

// #include <ros/time.h>
#include <g2o/core/hyper_graph.h>
#include <g2o/core/base_vertex.h>
#include <g2o/core/base_binary_edge.h>
#include <g2o/core/sparse_optimizer.h>
#include "se2.h"
#include "vertex_se2.hpp"
#include "edge_se2.hpp"


#include <g2o/solvers/dense/linear_solver_dense.h>
#include <g2o/core/optimization_algorithm_levenberg.h>
#include <g2o/stuff/macros.h>
#include <g2o/core/factory.h>
#include <g2o/core/block_solver.h>
#include <g2o/core/linear_solver.h>
#include <g2o/core/optimization_algorithm_factory.h>
#include <g2o/solvers/pcg/linear_solver_pcg.h>
#include <g2o/solvers/eigen/linear_solver_eigen.h>
#include <g2o/core/optimization_algorithm_dogleg.h>
#include <g2o/core/robust_kernel_impl.h>



#include <Eigen/Core>
#include <cmath>
#include <iostream>


class GraphSlam
{

public:
    g2o::SparseOptimizer optimizer;

public:
    GraphSlam()
    {

        typedef g2o::BlockSolver< g2o::BlockSolverTraits<3,3> > Block;
        Block::LinearSolverType* linearSolver = new g2o::LinearSolverDense<Block::PoseMatrixType>();
        Block* solver_ptr = new Block( std::unique_ptr<Block::LinearSolverType>(linearSolver) );
        g2o::OptimizationAlgorithmLevenberg* solver = new g2o::OptimizationAlgorithmLevenberg(std::unique_ptr<Block>(solver_ptr) );

        // optimizer.setVerbose(true);
        optimizer.setAlgorithm(solver);


    }
    ~GraphSlam()
    {

    }

    // VertexSE2* add_se2_node(const SE2& pose);
    // EdgeSE2* add_se2_edge(VertexSE2* v1, VertexSE2* v2, const SE2& relative_pose, const Eigen::Matrix3d& information_matrix);
    // int optimize(int num_iterations);

    int optimize(int num_iterations)
    {
        // optimizer.save("before.g2o");

        optimizer.initializeOptimization();
        optimizer.optimize(num_iterations);

        // optimizer.save("after.g2o");

    }



    VertexSE2* add_se2_node(const SE2& pose)
    {
        VertexSE2* vertex = new VertexSE2();
        vertex->setId(static_cast<int>(optimizer.vertices().size()));
        vertex->setEstimate(pose);
        optimizer.addVertex(vertex);
        return vertex;
    }

    EdgeSE2* add_se2_edge(VertexSE2* v1, VertexSE2* v2, const SE2& relative_pose, const Eigen::Matrix3d& information_matrix)
    {
        EdgeSE2* edge = new EdgeSE2();
        edge->setMeasurement(relative_pose);
        edge->setInformation(information_matrix);
        edge->vertices()[0] = v1;
        edge->vertices()[1] = v2;
        optimizer.addEdge(edge);
        return edge;
    }

};



#endif