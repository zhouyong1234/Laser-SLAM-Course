#ifndef READFILE_H
#define READFILE_H

#include <vector>
#include <eigen3/Eigen/Core>
#include <iostream>
#include <map>

#include "gaussian_newton.h"

using namespace std;




void ReadVertexInformation(const std::string path,std::vector<Eigen::Vector3d>& nodes);
void ReadEdgesInformation(const std::string path,std::vector<Edge>& edges);


void ReadPoses(const std::string path, std::map<int, Pose2d>& poses);
void ReadConstraints(const std::string path, std::vector<Constraint2d>& constraints);

bool OutputPoses(const std::string filename, const std::map<int, Pose2d>& poses);


#endif
