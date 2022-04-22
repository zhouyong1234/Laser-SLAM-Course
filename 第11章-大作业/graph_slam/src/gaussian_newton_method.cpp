#include <map.h>
#include "gaussian_newton_method.h"

const double GN_PI = 3.1415926;

//进行角度正则化．
double GN_NormalizationAngle(double angle)
{
    if(angle > GN_PI)
        angle -= 2*GN_PI;
    else if(angle < -GN_PI)
        angle += 2*GN_PI;

    return angle;
}

Eigen::Matrix3d GN_V2T(Eigen::Vector3d vec)
{
    Eigen::Matrix3d T;
    T  << cos(vec(2)),-sin(vec(2)),vec(0),
            sin(vec(2)), cos(vec(2)),vec(1),
            0,           0,     1;

    return T;
}

//对某一个点进行转换．
Eigen::Vector2d GN_TransPoint(Eigen::Vector2d pt,Eigen::Matrix3d T)
{
    Eigen::Vector3d tmp_pt(pt(0),pt(1),1);
    tmp_pt = T * tmp_pt;
    return Eigen::Vector2d(tmp_pt(0),tmp_pt(1));
}



//用激光雷达数据创建势场．
map_t* CreateMapFromLaserPoints(Eigen::Vector3d map_origin_pt,
                                std::vector<Eigen::Vector2d> laser_pts,
                                double resolution)
{
    map_t* map = map_alloc();

    map->origin_x = map_origin_pt(0);
    map->origin_y = map_origin_pt(1);
    map->resolution = resolution;

    //固定大小的地图，必要时可以扩大．
    map->size_x = 6000;
    map->size_y = 6000;

    map->cells = (map_cell_t*)malloc(sizeof(map_cell_t)*map->size_x*map->size_y);

    //高斯平滑的sigma－－固定死
    map->likelihood_sigma = 0.5;

    Eigen::Matrix3d Trans = GN_V2T(map_origin_pt);

    //设置障碍物
    for(int i = 0; i < laser_pts.size();i++)
    {
        Eigen::Vector2d tmp_pt = GN_TransPoint(laser_pts[i],Trans);

        int cell_x,cell_y;
        cell_x = MAP_GXWX(map,tmp_pt(0));
        cell_y = MAP_GYWY(map,tmp_pt(1));

        map->cells[MAP_INDEX(map,cell_x,cell_y)].occ_state = CELL_STATUS_OCC;
    }

    //进行障碍物的膨胀--最大距离固定死．
    map_update_cspace(map,0.5);

    return map;
}


/**
 * @brief InterpMapValueWithDerivatives
 * 在地图上的进行插值，得到coords处的势场值和对应的关于位置的梯度．
 * 返回值为Eigen::Vector3d ans
 * ans(0)表示市场值
 * ans(1:2)表示梯度
 * @param map
 * @param coords
 * @return
 */
Eigen::Vector3d InterpMapValueWithDerivatives(map_t* map,Eigen::Vector2d& coords)
{
    Eigen::Vector3d ans;
    //TODO
    int x0 = floor((coords(0) - map->origin_x) / map->resolution + 0.5) + map->size_x / 2;
    int y0 = floor((coords(1) - map->origin_y) / map->resolution + 0.5) + map->size_y / 2;

    double u = (coords(0) - map->origin_x) / map->resolution + 0.5 + double(map->size_x / 2) - (double)x0;
    double v = (coords(1) - map->origin_y) / map->resolution + 0.5 + double(map->size_y / 2) - (double)y0;

    double P1, P2, P3, P4;
    P1 = map->cells[MAP_INDEX(map, x0, y0)].score;
    P2 = map->cells[MAP_INDEX(map, x0+1, y0)].score;
    P3 = map->cells[MAP_INDEX(map, x0+1, y0+1)].score;
    P4 = map->cells[MAP_INDEX(map, x0, y0+1)].score;

    ans(0) = u * v * P3 + (1 - u) * v * P4 + (1 - u) * (1 - v) * P1 + u * (1 - v) * P2;
    ans(1) = (v * (P3 - P4) + (1 - v) * (P2 - P1)) / map->resolution;
    ans(2) = (u * (P3 - P2) + (1 - u) * (P4 - P1)) / map->resolution;
    //END OF TODO

    return ans;
}


/**
 * @brief ComputeCompleteHessianAndb
 * 计算H*dx = b中的H和b
 * @param map
 * @param now_pose
 * @param laser_pts
 * @param H
 * @param b
 */
double ComputeHessianAndb(map_t* map, Eigen::Vector3d now_pose,
                        std::vector<Eigen::Vector2d>& laser_pts,
                        Eigen::Matrix3d& H, Eigen::Vector3d& b)
{
    H = Eigen::Matrix3d::Zero();
    b = Eigen::Vector3d::Zero();

    //TODO
    Eigen::Matrix3d T;
    Eigen::Vector3d ST;
    Eigen::Vector2d ST2D;
    Eigen::Vector3d ans;
    Eigen::Vector2d DMST;
    Eigen::Matrix<double, 2, 3> DST;
    T << cos(now_pose(2)), -sin(now_pose(2)), now_pose(0),
         sin(now_pose(2)),  cos(now_pose(2)), now_pose(1),
         0,                 0,                1;

    Eigen::Vector3d laser_pose;
    double error = 0.0;
    for(int i = 0; i < laser_pts.size(); i++)
    {
        laser_pose << laser_pts[i](0), laser_pts[i](1), 1;
        ST = T * laser_pose;
        ST2D << ST(0), ST(1);
        ans = InterpMapValueWithDerivatives(map, ST2D);
        DST << 1, 0, -sin(now_pose(2)) * laser_pose(0) - cos(now_pose(2)) * laser_pose(1),
               0, 1,  cos(now_pose(2)) * laser_pose(0) - sin(now_pose(2)) * laser_pose(1);
        DMST << ans(1), ans(2);
        double J1 = ans(1);
        double J2 = ans(2);
        double J3 = ans(1) * DST(0, 2) + ans(2) * DST(1, 2);
        
        H(0, 0) += J1 * J1;
        H(0, 1) += J1 * J2;
        H(0, 2) += J1 * J3;
        H(1, 0) += H(0, 1);
        H(1, 1) += J2 * J2;
        H(1, 2) += J2 * J3;
        H(2, 0) += H(0, 2);
        H(2, 1) += H(1, 2);
        H(2, 2) += J3 * J3;

        double FX = 1 - ans(0);
        Eigen::Vector3d JFX;
        JFX << J1 * FX, J2 * FX, J3 * FX;
        b += JFX;
        error += FX * FX;
    }

    return error;
    //END OF TODO
}


/**
 * @brief GaussianNewtonOptimization
 * 进行高斯牛顿优化．
 * @param map
 * @param init_pose
 * @param laser_pts
 */
void GaussianNewtonOptimization(map_t*map,Eigen::Vector3d& init_pose,std::vector<Eigen::Vector2d>& laser_pts)
{
    int maxIteration = 20;
    Eigen::Vector3d now_pose = init_pose;

    double lasterror = 0.0;
    for(int i = 0; i < maxIteration;i++)
    {
        //TODO
        Eigen::Matrix3d H;
        Eigen::Vector3d b;
        double error = ComputeHessianAndb(map, now_pose, laser_pts, H, b);
        // Eigen::Vector3d dT = H.colPivHouseholderQr().solve(b);
        Eigen::Vector3d dT = H.ldlt().solve(b);
        
        if(lasterror > 0 && error > lasterror)
            break;
        
        if(std::isnan(dT[0]))
            break;
        lasterror = error;
        now_pose += dT;
        if(error < 10e-5)
            break;

        //END OF TODO
    }
    init_pose = now_pose;

}
