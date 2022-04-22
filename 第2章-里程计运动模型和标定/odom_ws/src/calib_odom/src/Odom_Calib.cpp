#include "../include/calib_odom/Odom_Calib.hpp"


//设置数据长度,即多少数据计算一次
void OdomCalib::Set_data_len(int len)
{
    data_len = len;
    A.conservativeResize(len*3,9);
    b.conservativeResize(len*3);
    A.setZero();
    b.setZero();
}


/*
输入:里程计和激光数据

TODO:
构建最小二乘需要的超定方程组
Ax = b

*/
bool OdomCalib::Add_Data(Eigen::Vector3d Odom,Eigen::Vector3d scan)
{

    if(now_len<INT_MAX)
    {
        //TODO: 构建超定方程组
        A(3 * now_len, 0) = Odom(0);
        A(3 * now_len, 1) = Odom(1);
        A(3 * now_len, 2) = Odom(2);
        A(3 * now_len + 1, 3) = Odom(0);
        A(3 * now_len + 1, 4) = Odom(1);
        A(3 * now_len + 1, 5) = Odom(2);
        A(3 * now_len + 2, 6) = Odom(0);
        A(3 * now_len + 2, 7) = Odom(1);
        A(3 * now_len + 2, 8) = Odom(2);
        b(3 * now_len) = scan(0);
        b(3 * now_len + 1) = scan(1);
        b(3 * now_len + 2) = scan(2);
        //end of TODO
        now_len++;
        return true;
    }
    else
    {
        return false;
    }
}

/*
 * TODO:
 * 求解线性最小二乘Ax=b
 * 返回得到的矫正矩阵
*/
Eigen::Matrix3d OdomCalib::Solve()
{
    Eigen::Matrix3d correct_matrix;

    //TODO: 求解线性最小二乘
    Eigen::VectorXd x(9);
    x = A.colPivHouseholderQr().solve(b);
    //end of TODO
    correct_matrix(0, 0) = x(0);
    correct_matrix(0, 1) = x(1);
    correct_matrix(0, 2) = x(2);
    correct_matrix(1, 0) = x(3);
    correct_matrix(1, 1) = x(4);
    correct_matrix(1, 2) = x(5);
    correct_matrix(2, 0) = x(6);
    correct_matrix(2, 1) = x(7);
    correct_matrix(2, 2) = x(8);
    return correct_matrix;
}

/* 用于判断数据是否满
 * 数据满即可以进行最小二乘计算
*/
bool OdomCalib::is_full()
{
    if(now_len%data_len==0&&now_len>=1)
    {
        now_len = data_len;
        return true;
    }
    else
        return false;
}

/*
 * 数据清零
*/
void OdomCalib::set_data_zero()
{
    A.setZero();
    b.setZero();
}
