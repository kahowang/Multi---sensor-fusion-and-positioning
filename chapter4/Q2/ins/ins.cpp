#include <vector>
#include <Eigen/Core>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <Eigen/Geometry>
#include "Geocentric/LocalCartesian.hpp"
#include <deque>
#include <random>
#include <yaml-cpp/yaml.h>

#define D2R 0.017453292519943295
struct State
{
    Eigen::Vector3d p;
    Eigen::Vector3d v;
    Eigen::Quaterniond q;
    Eigen::Vector3d bg;
    Eigen::Vector3d ba;
};
struct ErrorState
{
    Eigen::Matrix<double, 15, 1> x;         
    Eigen::Matrix<double, 15, 15> p;       
};
struct IMUData
{
    double time;
    Eigen::Vector3d acc;
    Eigen::Vector3d gyro;
};
struct GNSSData
{
    double time;
    Eigen::Vector3d lla;
    Eigen::Vector3d v;
};
struct PoseData
{
    double time;
    State state;
};

std::deque<GNSSData> gnss_buff;
std::deque<IMUData> imu_buff;
std::deque<PoseData> gt_buff;
std::deque<IMUData> current_imu;
GNSSData current_gnss;
PoseData current_gt;
PoseData current_pose;
ErrorState current_error_state;
GeographicLib::LocalCartesian geo_converter(32, 120, 0);
std::ofstream gt_ofs;
std::ofstream pose_ofs;
std::ofstream sv_ofs;
double gyro_noise = 1e-6;
double acc_noise = 1e-5;
double dp_noise = 1e-3;
std::vector<double> init_noise;
int FGsize = 20;
double time_interval = 10;
double end_time = 20;
double T = 0.1;
struct FG
{
    double time;
    Eigen::Matrix<double, 15, 15> F;                             //状态转移矩阵
    Eigen::Matrix<double, 3, 15> G;                               //观测转移矩阵 
    std::vector<Eigen::Matrix<double, 3, 1>> Y;       // 观测值
};
std::vector<FG> FGs;                                                       // 当前时刻的转移矩阵
Eigen::Matrix<double, 15, 15> Ft = Eigen::Matrix<double, 15, 15>::Zero();
Eigen::Matrix<double, 3, 15> Gt = Eigen::Matrix<double, 3, 15>::Zero();
Eigen::Matrix<double, 3, 1> Y = Eigen::Matrix<double, 3, 1>::Zero();
bool correct = true;

//   读取仿真数据
bool ReadData(const std::vector<std::string> &path)
{
    gnss_buff.clear();
    imu_buff.clear();
    gt_buff.clear();
    std::vector<std::ifstream> reads;
    // int count = 0;
    for (int i = 0; i < path.size(); ++i)
    {
        reads.push_back(std::ifstream(path[i]));
    }
    for (int i = 0; i < path.size(); ++i)
    {
        std::string strs;
        std::getline(reads[i], strs);
    }
    std::string strs;
    while (std::getline(reads[0], strs))
    {
        double time = std::stod(strs);
        std::getline(reads[1], strs);
        std::string temp = "";
        std::vector<double> acc;
        for (int i = 0; i < strs.size(); ++i)
        {
            if (strs[i] == ',')
            {
                acc.push_back(std::stod(temp));
                temp = "";
            }
            else
            {
                temp = temp + strs[i];
            }
        }
        acc.push_back(std::stod(temp));

        std::getline(reads[2], strs);
        temp = "";
        std::vector<double> gyro;
        for (int i = 0; i < strs.size(); ++i)
        {
            if (strs[i] == ',')
            {
                gyro.push_back(std::stod(temp));
                temp = "";
            }
            else
            {
                temp = temp + strs[i];
            }
        }
        gyro.push_back(std::stod(temp));
        IMUData imu;
        imu.time = time;
        imu.acc = Eigen::Vector3d(acc[0], acc[1], acc[2]);
        imu.gyro = Eigen::Vector3d(gyro[0] * D2R, gyro[1] * D2R, gyro[2] * D2R);
        imu_buff.push_back(imu);

        std::getline(reads[5], strs);
        temp = "";
        std::vector<double> ref_pos;
        for (int i = 0; i < strs.size(); ++i)
        {
            if (strs[i] == ',')
            {
                ref_pos.push_back(std::stod(temp));
                temp = "";
            }
            else
            {
                temp = temp + strs[i];
            }
        }
        ref_pos.push_back(std::stod(temp));

        std::getline(reads[6], strs);
        temp = "";
        std::vector<double> ref_vel;
        for (int i = 0; i < strs.size(); ++i)
        {
            if (strs[i] == ',')
            {
                ref_vel.push_back(std::stod(temp));
                temp = "";
            }
            else
            {
                temp = temp + strs[i];
            }
        }
        ref_vel.push_back(std::stod(temp));

        std::getline(reads[7], strs);
        temp = "";
        std::vector<double> ref_att_quat;
        for (int i = 0; i < strs.size(); ++i)
        {
            if (strs[i] == ',')
            {
                ref_att_quat.push_back(std::stod(temp));
                temp = "";
            }
            else
            {
                temp = temp + strs[i];
            }
        }
        ref_att_quat.push_back(std::stod(temp));

        Eigen::Quaterniond q = Eigen::AngleAxisd(90 * D2R, Eigen::Vector3d::UnitZ()) *
                               Eigen::AngleAxisd(0, Eigen::Vector3d::UnitY()) *
                               Eigen::AngleAxisd(180 * D2R, Eigen::Vector3d::UnitX());
        q = q.inverse();

        PoseData pose;
        pose.time = time;
        double geo_x, geo_y, geo_z;
        geo_converter.Forward(ref_pos[0], ref_pos[1], ref_pos[2], geo_x, geo_y, geo_z);
        pose.state.p = Eigen::Vector3d(geo_x, geo_y, geo_z);
        // pose.state.p = q * Eigen::Vector3d(ref_pos[0], ref_pos[1], ref_pos[2]);
        pose.state.v = q * Eigen::Vector3d(ref_vel[0], ref_vel[1], ref_vel[2]);
        pose.state.q = q * Eigen::Quaterniond(ref_att_quat[0], ref_att_quat[1], ref_att_quat[2], ref_att_quat[3]);
        pose.state.q.normalize();
        pose.state.bg = Eigen::Vector3d(0, 0, 0);
        pose.state.ba = Eigen::Vector3d(0, 0, 0);
        gt_buff.push_back(pose);
    }
    while (std::getline(reads[3], strs))
    {
        double time = std::stod(strs);
        std::getline(reads[4], strs);
        std::string temp = "";
        std::vector<double> gps;
        for (int i = 0; i < strs.size(); ++i)
        {
            if (strs[i] == ',')
            {
                gps.push_back(std::stod(temp));
                temp = "";
            }
            else
            {
                temp = temp + strs[i];
            }
        }
        gps.push_back(std::stod(temp));
        GNSSData gnss;
        gnss.time = time;
        gnss.lla = Eigen::Vector3d(gps[0], gps[1], gps[2]);
        // 北东地   ->   东北天
        Eigen::Quaterniond q = Eigen::AngleAxisd(90 * D2R, Eigen::Vector3d::UnitZ()) *
                               Eigen::AngleAxisd(0, Eigen::Vector3d::UnitY()) *
                               Eigen::AngleAxisd(180 * D2R, Eigen::Vector3d::UnitX());
        q = q.inverse();
        gnss.v = q * Eigen::Vector3d(gps[3], gps[4], gps[5]);
        gnss_buff.push_back(gnss);
    }
}

// 同步 GPS 和 IMU 数据
bool SyncData(bool inited)
{
    if (gnss_buff.empty())
    {
        return false;
    }
    current_gnss = gnss_buff.front();   //读取当前 gnss 数据的第一个时间
    double sync_time = current_gnss.time;    // 同步时间以当前 gnss 时间为准

    while (gt_buff.size() > 1)
    {
        if (gt_buff[1].time < sync_time)
        {
            gt_buff.pop_front();      // 删除 gt.buff  第一个值
        }
        else
        {
            break;
        }
    }

    //   同步当前的 groundtruth
    if (gt_buff.size() > 1)
    {
        PoseData front_data = gt_buff.at(0);
        PoseData back_data = gt_buff.at(1);
        double front_scale = (back_data.time - sync_time) / (back_data.time - front_data.time);
        double back_scale = (sync_time - front_data.time) / (back_data.time - front_data.time);
        current_gt.time = sync_time;
        current_gt.state.p = front_data.state.p * front_scale + back_data.state.p * back_scale;
        current_gt.state.v = front_data.state.v * front_scale + back_data.state.v * back_scale;
        current_gt.state.q = front_data.state.q.slerp(front_scale, back_data.state.q);
        current_gt.state.bg = front_data.state.bg * front_scale + back_data.state.bg * back_scale;
        current_gt.state.ba = front_data.state.ba * front_scale + back_data.state.ba * back_scale;
    }
    else
    {
        return false;
    }


// 同步当前 imu的数据
    while (!inited && imu_buff.size() > 1)
    {
        if (imu_buff[1].time < sync_time)
        {
            imu_buff.pop_front();
        }
        else
        {
            break;
        }
    }

    if (imu_buff.size() > 1)
    {
        if (!inited)
        {
            current_imu.clear();
            IMUData front_data = imu_buff.at(0);
            IMUData back_data = imu_buff.at(1);
            IMUData synced_data;

            double front_scale = (back_data.time - sync_time) / (back_data.time - front_data.time);
            double back_scale = (sync_time - front_data.time) / (back_data.time - front_data.time);
            synced_data.time = sync_time;
            synced_data.acc = front_data.acc * front_scale + back_data.acc * back_scale;
            synced_data.gyro = front_data.gyro * front_scale + back_data.gyro * back_scale;
            current_imu.push_back(synced_data);
            imu_buff.pop_front();
            gnss_buff.pop_front();
            // std::cout << std::setprecision(12) << "sync_time " << sync_time
            //           << " current_imu.time " << current_imu.front().time
            //           << "  " << current_imu.back().time << std::endl;
            return true;
        }

        if (imu_buff.back().time < sync_time)
        {
            return false;
        }
        while (current_imu.size() > 1)
        {
            current_imu.pop_front();
        }
        while (imu_buff.front().time < sync_time)
        {
            IMUData temp = imu_buff.front();
            imu_buff.pop_front();
            current_imu.push_back(temp);
        }
        IMUData front_data = current_imu.back();
        IMUData back_data = imu_buff.at(0);
        IMUData synced_data;

        double front_scale = (back_data.time - sync_time) / (back_data.time - front_data.time);
        double back_scale = (sync_time - front_data.time) / (back_data.time - front_data.time);
        synced_data.time = sync_time;
        synced_data.acc = front_data.acc * front_scale + back_data.acc * back_scale;
        synced_data.gyro = front_data.gyro * front_scale + back_data.gyro * back_scale;
        current_imu.push_back(synced_data);
        gnss_buff.pop_front();
        // std::cout << std::setprecision(12) << "sync_time " << sync_time
        //           << " current_imu.time " << current_imu.front().time
        //           << "  " << current_imu.back().time << std::endl;
        return true;
    }
    else
    {
        return false;
    }
}

bool InitSensor()
{
    while (!gnss_buff.empty())
    {
        if (imu_buff.front().time > gnss_buff.front().time)
        {
            gnss_buff.pop_front();    //删掉时间不匹配的 gnss 数据
        }
        else
        {
            return true;
        }
    }
    return false;
}

bool InitPose()
{
    static bool pose_inited = false;
    if (pose_inited)
    {
        return true;
    }
    if (!SyncData(false))
    {
        return false;
    }
    current_pose.time = current_gt.time;
    current_pose.state.p = current_gt.state.p;
    current_pose.state.q = current_gt.state.q;
    current_pose.state.v = current_gt.state.v;
    current_pose.state.bg = current_gt.state.bg;
    current_pose.state.ba = current_gt.state.ba;
    current_error_state.x.setZero();
    current_error_state.p.setZero();
    current_error_state.p.block<3, 3>(0, 0) = Eigen::Matrix3d::Identity() * init_noise[0];
    current_error_state.p.block<3, 3>(3, 3) = Eigen::Matrix3d::Identity() * init_noise[1];
    current_error_state.p.block<3, 3>(6, 6) = Eigen::Matrix3d::Identity() * init_noise[2];
    current_error_state.p.block<3, 3>(9, 9) = Eigen::Matrix3d::Identity() * init_noise[3];
    current_error_state.p.block<3, 3>(12, 12) = Eigen::Matrix3d::Identity() * init_noise[4];
    pose_inited = true;
    return true;
}
bool Predict()
{
    current_pose.time = current_gt.time;
    Eigen::Vector3d pp = current_pose.state.p;
    Eigen::Vector3d vv = current_pose.state.v;
    Eigen::Quaterniond qq = current_pose.state.q;
    double w = 7.27220521664304e-05;
    Eigen::Vector3d gn(0, 0, -9.79484197226504);
    Eigen::Vector3d w_ie_n(0, w * std::cos(current_gnss.lla[0] * M_PI / 180),
                           w * std::sin(current_gnss.lla[0] * M_PI / 180));
    double rm = 6353346.18315;       // 短半轴
    double rn = 6384140.52699;        // 长半轴
    Eigen::Vector3d w_en_n(-vv[1] / (rm + current_gnss.lla[2]),
                           vv[0] / (rn + current_gnss.lla[2]),
                           vv[0] / (rn + current_gnss.lla[2]) * std::tan(current_gnss.lla[0] * M_PI / 180));
    Eigen::Vector3d w_in_n = w_ie_n + w_en_n;
    for (int i = 1; i < current_imu.size(); ++i)
    {
        double dt = current_imu[i].time - current_imu[i - 1].time;
        Eigen::Vector3d wtemp = w_in_n * dt;
        double angle = wtemp.norm();
        Eigen::Quaterniond qn(1, 0, 0, 0);
        if (angle != 0)
        {
            wtemp = wtemp / angle;
            wtemp = std::sin(angle / 2) * wtemp;
            qn = Eigen::Quaterniond(std::cos(angle / 2), wtemp[0], wtemp[1], wtemp[2]);
        }
        qn.normalize();   // 地球自转的叫角增量

        Eigen::Vector3d wb = 0.5 * current_imu[i - 1].gyro + 0.5 * current_imu[i].gyro;
        wb = wb + current_pose.state.bg;
        wb = wb * dt;
        angle = wb.norm();
        Eigen::Quaterniond qb(1, 0, 0, 0);
        if (angle != 0)
        {
            wb = wb / angle;
            wb = std::sin(angle / 2) * wb;
            qb = Eigen::Quaterniond(std::cos(angle / 2), wb[0], wb[1], wb[2]);
        }
        qb.normalize();        //  载体角增量

        Eigen::Quaterniond qq2 = qn.inverse() * qq * qb;   //旋转
        Eigen::Vector3d f1 = current_imu[i - 1].acc;
        f1 = f1 + current_pose.state.ba;
        Eigen::Vector3d f2 = current_imu[i].acc;
        f2 = f2 + current_pose.state.ba;
        Eigen::Vector3d vv2 = vv + dt * (0.5 * (qq * f1 + qq2 * f2) + gn);
        Eigen::Vector3d pp2 = pp + 0.5 * dt * (vv + vv2);
        pp = pp2;
        vv = vv2;
        qq = qq2;
    }
    current_pose.state.p = pp;       // 更新位置
    current_pose.state.v = vv;         //更新速度
    current_pose.state.q = qq;        //更新旋转

    Ft = Eigen::Matrix<double, 15, 15>::Zero();
    Ft.block<3, 3>(0, 3) = Eigen::Matrix<double, 3, 3>::Identity();
    Eigen::Matrix<double, 3, 3> temp = Eigen::Matrix<double, 3, 3>::Zero();
    Eigen::Vector3d ff = current_imu.back().acc;
    ff = qq * ff;
    temp(0, 1) = -ff[2];
    temp(0, 2) = ff[1];
    temp(1, 0) = ff[2];
    temp(1, 2) = -ff[0];
    temp(2, 0) = -ff[1];
    temp(2, 1) = ff[0];
    Ft.block<3, 3>(3, 6) = temp;
    Ft.block<3, 3>(3, 12) = qq.toRotationMatrix();
    temp.setZero();
    temp(0, 1) = w_ie_n(2);
    temp(0, 2) = -w_ie_n(1);
    temp(1, 0) = -w_ie_n(2);
    temp(2, 0) = w_ie_n(1);
    Ft.block<3, 3>(6, 6) = temp;
    Ft.block<3, 3>(6, 9) = -Ft.block<3, 3>(3, 12);
    Eigen::Matrix<double, 15, 6> Bt = Eigen::Matrix<double, 15, 6>::Zero();
    Bt.block<3, 3>(3, 3) = Ft.block<3, 3>(3, 12);
    Bt.block<3, 3>(6, 0) = Ft.block<3, 3>(6, 9);
    T = current_imu.back().time - current_imu.front().time;
    Ft = Eigen::Matrix<double, 15, 15>::Identity() + Ft * T;      // 求离散滤波器下的状态转移矩阵
    Bt = Bt * T;                                                                                            // 求观测矩阵
    Eigen::Matrix<double, 6, 1> W = Eigen::Matrix<double, 6, 1>::Zero();     // 器件噪声
    // std::random_device rd;
    // std::default_random_engine generator(rd());
    // std::normal_distribution<double> distribution(0.0, 1.0);
    // Eigen::Vector3d noise_gyro(distribution(generator), distribution(generator), distribution(generator));
    // Eigen::Vector3d noise_acc(distribution(generator), distribution(generator), distribution(generator));
    // noise_gyro = noise_gyro * gyro_noise;
    // noise_acc = noise_acc * acc_noise;
    // W.head(3) = noise_gyro;
    // W.tail(3) = noise_acc;
    current_error_state.x = Ft * current_error_state.x + Bt * W;     //   k 时刻的状态先验
    Eigen::Matrix<double, 6, 6> Q = Eigen::Matrix<double, 6, 6>::Identity();                //  惯性噪声
    Q.block<3, 3>(0, 0) = Eigen::Matrix<double, 3, 3>::Identity() * gyro_noise * gyro_noise;
    Q.block<3, 3>(3, 3) = Eigen::Matrix<double, 3, 3>::Identity() * acc_noise * acc_noise;                   
    current_error_state.p = Ft * current_error_state.p * Ft.transpose() + Bt * Q * Bt.transpose();
    return true;
}
bool Correct()
{
    double geo_x, geo_y, geo_z;
    geo_converter.Forward(current_gnss.lla(0), current_gnss.lla(1),
                          current_gnss.lla(2), geo_x, geo_y, geo_z);      // 将经纬度处理为  东北天
    Eigen::Vector3d gnss_xyz(geo_x, geo_y, geo_z);
    // std::cout << gnss_xyz.transpose() << std::endl;
    Y.block<3, 1>(0, 0) = current_pose.state.p - gnss_xyz;      //观测误差
    // Y.block<3, 1>(3, 0) = current_pose.state.v - current_gnss.v;

    Gt = Eigen::Matrix<double, 3, 15>::Zero();
    Gt.block<3, 3>(0, 0) = Eigen::Matrix<double, 3, 3>::Identity();
    // Gt.block<3, 3>(3, 3) = Eigen::Matrix<double, 3, 3>::Identity();
    Eigen::Matrix<double, 3, 3> Ct = Eigen::Matrix<double, 3, 3>::Identity();

    // Eigen::Matrix<double, 6, 1> N = Eigen::Matrix<double, 6, 1>::Zero();
    // std::random_device rd;
    // std::default_random_engine generator(rd());
    // std::normal_distribution<double> distribution(0.0, 1.0);
    // Eigen::Vector3d noise_dp(distribution(generator),distribution(generator),distribution(generator));
    // noise_dp = noise_dp * dp_noise;
    Eigen::Matrix<double, 3, 3> R = Eigen::Matrix<double, 3, 3>::Identity();    // 观测噪声， gps噪声
    R.block<3, 3>(0, 0) = Eigen::Matrix3d::Identity() * dp_noise * dp_noise;
    Eigen::Matrix<double, 15, 3> K = current_error_state.p * Gt.transpose() * (Gt * current_error_state.p * Gt.transpose() + Ct * R * Ct.transpose()).inverse();
    // update
    current_error_state.p = (Eigen::Matrix<double, 15, 15>::Identity() - K * Gt) * current_error_state.p;
    current_error_state.x = current_error_state.x + K * (Y - Gt * current_error_state.x);
    current_pose.state.p = current_pose.state.p - current_error_state.x.block<3, 1>(0, 0);
    current_pose.state.v = current_pose.state.v - current_error_state.x.block<3, 1>(3, 0);
    Eigen::Vector3d dphi_dir = current_error_state.x.block<3, 1>(6, 0);
    double dphi_norm = dphi_dir.norm();
    if (dphi_norm != 0)
    {
        dphi_dir = dphi_dir / dphi_norm;
        dphi_dir = dphi_dir * std::sin(dphi_norm / 2);
    }
    Eigen::Quaterniond temp2(std::cos(dphi_norm / 2), dphi_dir[0], dphi_dir[1], dphi_dir[2]);
    current_pose.state.q = temp2 * current_pose.state.q;     // 更新旋转矢量
    current_pose.state.q.normalize();
    current_pose.state.bg = current_pose.state.bg - current_error_state.x.block<3, 1>(9, 0);
    current_pose.state.ba = current_pose.state.ba - current_error_state.x.block<3, 1>(12, 0);
    current_error_state.x.setZero();        // 清空误差
    return true;
}

/*保存观测度分析所需的F G 和 Y*/
bool SaveFG()    
{
    if (FGs.size() > FGsize)
    {
        return true;
    }
    if (FGs.empty())
    {
        FG fg;
        fg.time = current_gt.time;
        // fg.F = Ft;
        fg.F = Ft - Eigen::Matrix<double, 15, 15>::Identity();
        // fg.F = (Ft - Eigen::Matrix<double, 15, 15>::Identity()) / T;
        fg.G = Gt;
        fg.Y.push_back(Y);
        FGs.push_back(fg);
    }
    else
    {
        if (FGs.back().Y.size() == 15)
        {
            if (current_gt.time - FGs.back().time < time_interval || FGs.size() >= FGsize)
            {
                return true;
            }
            FG fg;
            fg.time = current_gt.time;
            // fg.F = Ft;
            fg.F = Ft - Eigen::Matrix<double, 15, 15>::Identity();
            // fg.F = (Ft - Eigen::Matrix<double, 15, 15>::Identity()) / T;
            fg.G = Gt;
            fg.Y.push_back(Y);
            FGs.push_back(fg);
        }
        else
        {
            FGs.back().Y.push_back(Y);
        }
        
    }
    return true;
}
/*滤波*/
bool Filter()    
{
    Predict();
    if (correct)
    {
        Correct();
    }
    return true;
}
void SavePose(std::ofstream &save_points, PoseData &pose)
{
    Eigen::Quaterniond qtemp = pose.state.q;
    // if (qtemp.w() < 0)
    // {
    //     qtemp.coeffs() = -1.0 * qtemp.coeffs();
    // }
    double angle = std::acos(qtemp.w()) * 2;
    double sin_angle = std::sin(angle / 2);
    Eigen::Vector3d dir(0, 0, 0);
    if (sin_angle != 0)
    {
        dir(0) = qtemp.x() / sin_angle;
        dir(1) = qtemp.y() / sin_angle;
        dir(2) = qtemp.z() / sin_angle;
        dir = dir * angle;
    }
    save_points.precision(12);
    save_points << pose.time
                << "," << pose.state.p(0)
                << "," << pose.state.p(1)
                << "," << pose.state.p(2)
                << "," << pose.state.v(0)
                << "," << pose.state.v(1)
                << "," << pose.state.v(2)
                // << "," << pose.state.q.x()
                // << "," << pose.state.q.y()
                // << "," << pose.state.q.z()
                // << "," << pose.state.q.w()
                << "," << dir(0)
                << "," << dir(1)
                << "," << dir(2)
                << "," << pose.state.bg(0)
                << "," << pose.state.bg(1)
                << "," << pose.state.bg(2)
                << "," << pose.state.ba(0)
                << "," << pose.state.ba(1)
                << "," << pose.state.ba(2)
                << std::endl;
}
bool SaveData()    // 保存数据
{
    SavePose(gt_ofs, current_gt);
    SavePose(pose_ofs, current_pose);
}

int main(int argc, char **argv)
{
    std::vector<std::string> path;

    path.push_back("../../gnss-ins-sim/imu_data/data7/time.csv");
    path.push_back("../../gnss-ins-sim/imu_data/data7/accel-0.csv");
    path.push_back("../../gnss-ins-sim/imu_data/data7/gyro-0.csv");
    path.push_back("../../gnss-ins-sim/imu_data/data7/gps_time.csv");
    path.push_back("../../gnss-ins-sim/imu_data/data7/gps-0.csv");
    path.push_back("../../gnss-ins-sim/imu_data/data7/ref_pos.csv");
    path.push_back("../../gnss-ins-sim/imu_data/data7/ref_vel.csv");
    path.push_back("../../gnss-ins-sim/imu_data/data7/ref_att_quat.csv");
    ReadData(path);
    gt_ofs.open("../../gnss-ins-sim/imu_data/data7/gt.txt", std::fstream::out);
    pose_ofs.open("../../gnss-ins-sim/imu_data/data7/pose.txt", std::fstream::out);
    sv_ofs.open("../../gnss-ins-sim/imu_data/data7/sv.txt", std::fstream::out);
    FGs.clear();
    
    YAML::Node yaml_node = YAML::LoadFile("../param.yaml");
    gyro_noise = yaml_node["gyro_noise"].as<double>();
    acc_noise = yaml_node["acc_noise"].as<double>();
    dp_noise = yaml_node["dp_noise"].as<double>();
    init_noise = yaml_node["init_noise"].as<std::vector<double>>();
    FGsize = yaml_node["FGsize"].as<int>();
    end_time = yaml_node["end_time"].as<double>();
    time_interval = yaml_node["time_interval"].as<double>();
    correct = yaml_node["correct"].as<bool>();

    if (!InitSensor())
    {
        std::cerr << "InitSensor Error!!!" << std::endl;
        return -1;
    }
    if (!InitPose())
    {
        std::cerr << "InitPose Error!!!" << std::endl;
        return -1;
    }
    SaveData();
    while (SyncData(true))
    {
        Filter();
        SaveData();
        SaveFG();
        if (current_gt.time > end_time)
        {
            break;
        }
    }

    Eigen::MatrixXd Qso(3*15*FGs.size(), 15);
    Eigen::MatrixXd Ys(3*15*FGs.size(), 1);
    Eigen::Matrix<double, 15, 15> Fnn = Eigen::Matrix<double, 15, 15>::Identity();
    for (int i = 0; i < FGs.size(); ++i)
    {
        Eigen::Matrix<double, 15, 15> Fn = Eigen::Matrix<double, 15, 15>::Identity();
        for (int j = 0; j < 15; j++)
        {
            if (j > 0)
            {
                Fn = Fn * FGs[i].F;
            }
            Ys.block<3, 1>(i*3*15+3*j, 0) = FGs[i].Y[j];
            Qso.block<3, 15>(i*3*15+3*j, 0) = FGs[i].G * Fn * Fnn;
        }
        // Fnn = Fnn * Fn;
    }

    Eigen::JacobiSVD<Eigen::MatrixXd> svd(Qso, Eigen::ComputeFullU | Eigen::ComputeFullV);
    std::cout << Qso.rows() << ", " << Qso.cols() << std::endl;
    // std::cout << svd.singularValues() << std::endl;
    for (int i = 0; i < 15; ++i)
    {
        double temp = (svd.matrixU().row(i) * Ys)[0] / svd.singularValues()[i];
        Eigen::MatrixXd Xi = temp * svd.matrixV().col(i);
        // std::cout << Xi.transpose() << std::endl;
        Eigen::MatrixXd::Index maxRow, maxCol;
        Xi = Xi.cwiseAbs();   //  Xi的绝对值
        double maxvalue = Xi.maxCoeff(&maxRow, &maxCol);
        std::cout << svd.singularValues()(i) / svd.singularValues()(0) << "," << maxRow << std::endl;
        sv_ofs << svd.singularValues()(i) / svd.singularValues()(0) << "," << maxRow << std::endl;
    }
    return 0;
}