#ifndef GAUSSIAN_NEWTON_H
#define GAUSSIAN_NEWTON_H

#include <vector>
#include <cmath>
#include <map>
#include <eigen3/Eigen/Core>
// #include "readfile.h"

using namespace std;

typedef struct edge
{
  int xi,xj;
  Eigen::Vector3d measurement;
  Eigen::Matrix3d infoMatrix;
}Edge;


struct Pose2d {
    double x;
    double y;
    double yaw_radians;

    static std::string name() {
        return "VERTEX2";
    }
};


struct Constraint2d {
    int id_begin;
    int id_end;

    double x;
    double y;
    double yaw_radians;

    Eigen::Matrix3d information;

    static std::string name() {
        return "EDGE2";
    }
};


Eigen::VectorXd  LinearizeAndSolve(std::vector<Eigen::Vector3d>& Vertexs,
                                   std::vector<Edge>& Edges);

double ComputeError(std::vector<Eigen::Vector3d>& Vertexs,
                    std::vector<Edge>& Edges);

double ComputeError(std::map<int, Pose2d>& poses, std::vector<Constraint2d>& constraints);


double GN_NormalizationAngle(double angle);



#endif
