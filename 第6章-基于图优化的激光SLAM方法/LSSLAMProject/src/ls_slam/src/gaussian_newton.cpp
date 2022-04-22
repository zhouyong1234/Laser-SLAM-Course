#include "gaussian_newton.h"
#include <eigen3/Eigen/Jacobi>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Householder>
#include <eigen3/Eigen/Cholesky>
#include <eigen3/Eigen/LU>

#include <iostream>
#include <chrono>
#include <string>


const double GN_PI = 3.1415926;

using namespace std::chrono;


//进行角度正则化．
double GN_NormalizationAngle(double angle)
{
    if(angle > GN_PI)
        angle -= 2*GN_PI;
    else if(angle < -GN_PI)
        angle += 2*GN_PI;

    return angle;
}


//位姿-->转换矩阵
Eigen::Matrix3d PoseToTrans(Eigen::Vector3d x)
{
    Eigen::Matrix3d trans;
    trans << cos(x(2)),-sin(x(2)),x(0),
             sin(x(2)), cos(x(2)),x(1),
                     0,         0,    1;

    return trans;
}


//转换矩阵－－＞位姿
Eigen::Vector3d TransToPose(Eigen::Matrix3d trans)
{
    Eigen::Vector3d pose;
    pose(0) = trans(0,2);
    pose(1) = trans(1,2);
    pose(2) = atan2(trans(1,0),trans(0,0));

    return pose;
}

double ComputeError(std::map<int, Pose2d>& poses, std::vector<Constraint2d>& constraints)
{
    double sumError = 0;
    for(std::vector<Constraint2d>::const_iterator constraint_iter = constraints.begin(); constraint_iter != constraints.end(); ++constraint_iter)
    {
        const Constraint2d& constraint = *constraint_iter;
        std::map<int, Pose2d>::iterator pose_begin_ptr = poses.find(constraint.id_begin);
        std::map<int, Pose2d>::iterator pose_end_ptr = poses.find(constraint.id_end);

        Eigen::Vector3d xi = Eigen::Vector3d(pose_begin_ptr->second.x, pose_begin_ptr->second.y, pose_begin_ptr->second.yaw_radians);
        Eigen::Vector3d xj = Eigen::Vector3d(pose_end_ptr->second.x, pose_end_ptr->second.y, pose_end_ptr->second.yaw_radians);
        Eigen::Vector3d z = Eigen::Vector3d(constraint.x, constraint.y, constraint.yaw_radians);
        Eigen::Matrix3d infoMatrix = constraint.information;

        Eigen::Matrix3d Xi = PoseToTrans(xi);
        Eigen::Matrix3d Xj = PoseToTrans(xj);
        Eigen::Matrix3d Z = PoseToTrans(z);

        Eigen::Matrix3d Ei = Z.inverse() * Xi.inverse() * Xj;
        Eigen::Vector3d ei = TransToPose(Ei);

        sumError += ei.transpose() * infoMatrix * ei;
    }
    return sumError;
}

//计算整个pose-graph的误差
double ComputeError(std::vector<Eigen::Vector3d>& Vertexs,
                    std::vector<Edge>& Edges)
{
    double sumError = 0;
    for(int i = 0; i < Edges.size();i++)
    {
        Edge tmpEdge = Edges[i];
        Eigen::Vector3d xi = Vertexs[tmpEdge.xi];
        Eigen::Vector3d xj = Vertexs[tmpEdge.xj];
        Eigen::Vector3d z = tmpEdge.measurement;
        Eigen::Matrix3d infoMatrix = tmpEdge.infoMatrix;

        Eigen::Matrix3d Xi = PoseToTrans(xi);
        Eigen::Matrix3d Xj = PoseToTrans(xj);
        Eigen::Matrix3d Z  = PoseToTrans(z);

        Eigen::Matrix3d Ei = Z.inverse() *  Xi.inverse() * Xj;

        Eigen::Vector3d ei = TransToPose(Ei);


        sumError += ei.transpose() * infoMatrix * ei;
    }
    return sumError;
}


/**
 * @brief CalcJacobianAndError
 *         计算jacobian矩阵和error
 * @param xi    fromIdx
 * @param xj    toIdx
 * @param z     观测值:xj相对于xi的坐标
 * @param ei    计算的误差
 * @param Ai    相对于xi的Jacobian矩阵
 * @param Bi    相对于xj的Jacobian矩阵
 */
void CalcJacobianAndError(Eigen::Vector3d xi,Eigen::Vector3d xj,Eigen::Vector3d z,
                          Eigen::Vector3d& ei,Eigen::Matrix3d& Ai,Eigen::Matrix3d& Bi)
{
    //TODO--Start
    Eigen::Matrix3d Xi = PoseToTrans(xi);
    Eigen::Matrix3d Xj = PoseToTrans(xj);
    Eigen::Matrix3d Z = PoseToTrans(z);
    Eigen::Matrix2d Rij = Z.block(0, 0, 2, 2);
    Eigen::Matrix2d Ri = Xi.block(0, 0, 2, 2);
    Eigen::Matrix2d Rijt = Rij.transpose();
    Eigen::Matrix2d Rit = Ri.transpose();
    Eigen::Vector2d ti;
    Eigen::Vector2d tj;
    Eigen::Vector2d tij;
    ti << xi(0), xi(1);
    tj << xj(0), xj(1);
    tij << z(0), z(1);

    Eigen::Matrix2d Rith;
    Rith << -sin(xi(2)), cos(xi(2)),
           -cos(xi(2)), -sin(xi(2));

    Ai = Eigen::Matrix3d::Zero();
    Bi = Eigen::Matrix3d::Zero();
    Ai.block(0, 0, 2, 2) = -Rijt * Rit;
    Ai.block(0, 2, 2, 1) = Rijt * Rith * (tj - ti);
    Ai(2, 2) = -1;
    Bi.block(0, 0, 2, 2) = Rijt * Rit;
    Bi(2, 2) = 1;

    Eigen::Vector2d error = Rijt * (Rit * (tj - ti) - tij);
    double angle = GN_NormalizationAngle(xj(2) - xi(2) - z(2));
    ei << error(0), error(1), angle;

    // std::cout << "Ai: " << "\n" << Ai << std::endl;
    // std::cout << "Bi: " << "\n" << Bi << std::endl;
    // std::cout << "ei: " << "\n" << ei << std::endl;
    //TODO--end
}

/**
 * @brief LinearizeAndSolve
 *        高斯牛顿方法的一次迭代．
 * @param Vertexs   图中的所有节点
 * @param Edges     图中的所有边
 * @return          位姿的增量
 */
Eigen::VectorXd  LinearizeAndSolve(std::vector<Eigen::Vector3d>& Vertexs,
                                   std::vector<Edge>& Edges)
{
    //申请内存
    Eigen::MatrixXd H(Vertexs.size() * 3,Vertexs.size() * 3);
    Eigen::VectorXd b(Vertexs.size() * 3);

    H.setZero();
    b.setZero();

    //固定第一帧
    Eigen::Matrix3d I;
    I.setIdentity();
    H.block(0,0,3,3) += I;

    //构造H矩阵　＆ b向量
    for(int i = 0; i < Edges.size();i++)
    {
        //提取信息
        Edge tmpEdge = Edges[i];
        Eigen::Vector3d xi = Vertexs[tmpEdge.xi];
        Eigen::Vector3d xj = Vertexs[tmpEdge.xj];
        Eigen::Vector3d z = tmpEdge.measurement;
        Eigen::Matrix3d infoMatrix = tmpEdge.infoMatrix;

        //计算误差和对应的Jacobian
        Eigen::Vector3d ei;
        Eigen::Matrix3d Ai;
        Eigen::Matrix3d Bi;

        

        CalcJacobianAndError(xi,xj,z,ei,Ai,Bi);

        

        //TODO--Start

        int x = tmpEdge.xi;
        int y = tmpEdge.xj;


        H.block(3 * x, 3 * x, 3, 3) += Ai.transpose() * infoMatrix * Ai;
        H.block(3 * x, 3 * y, 3, 3) += Ai.transpose() * infoMatrix * Bi;
        H.block(3 * y, 3 * x, 3, 3) += Bi.transpose() * infoMatrix * Ai;
        H.block(3 * y, 3 * y, 3, 3) += Bi.transpose() * infoMatrix * Bi;

        // std::cout << "x: " << x << " y: " << y << std::endl;
        // std::cout << "H: " << std::endl << H << std::endl;

        // std::cout << "AA: " << std::endl << H.block(3 * x, 3 * x, 3, 3) << std::endl;
        // std::cout << "AB: " << std::endl << H.block(3 * x, 3 * y, 3, 3) << std::endl;
        // std::cout << "BA: " << std::endl << H.block(3 * y, 3 * x, 3, 3) << std::endl;
        // std::cout << "BB: " << std::endl << H.block(3 * y, 3 * y, 3, 3) << std::endl;

        Eigen::Vector3d bi, bj;
        bi = (ei.transpose() * infoMatrix * Ai).transpose();
        bj = (ei.transpose() * infoMatrix * Bi).transpose();

        b(3 * tmpEdge.xi) += bi(0);
        b(3 * tmpEdge.xi+1) += bi(1);
        b(3 * tmpEdge.xi+2) += bi(2);
        b(3 * tmpEdge.xj) += bj(0);
        b(3 * tmpEdge.xj+1) += bj(1);
        b(3 * tmpEdge.xj+2) += bj(2);
        // std::cout << "b: " << std::endl << b << std::endl;
        //TODO--End
    }

    // std::cout << "H: " << std::endl << H << std::endl;
    // std::cout << "b: " << std::endl << b << std::endl;

    //求解
    Eigen::VectorXd dx(Vertexs.size() * 3);

    //TODO--Start
    // dx = -H.colPivHouseholderQr().solve(b);
    time_point<system_clock> start = system_clock::now();
    
    dx = (H.transpose() * H).llt().solve(H.transpose() * (-b));

    time_point<system_clock> end = system_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        std::cout << "Elapsed time: " << elapsed.count() << "s" << std::endl;

    // std::cout << "dx: " << std::endl << dx << std::endl;
    //TODO-End

    return dx;
}











