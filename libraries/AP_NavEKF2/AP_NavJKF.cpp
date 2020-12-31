#include "AP_NavJKF.h"

AP_NavJKF::AP_NavJKF(const AP_AHRS *ahrs) : 
    _ahrs(ahrs)
{

}


AP_NavJKF::~AP_NavJKF()
{
    memset(X, 0, sizeof(X));
    memset(P, 0, sizeof(P));
    memset(Q, 0, sizeof(Q));
    memset(R, 0, sizeof(R));
}

//过程噪声
#define VELNE_P_NSE_DEFAULT     1e-8f
#define VELD_P_NSE_DEFAULT     1e-8f

#define ABIAS_NE_P_NSE_DEFAULT     1e-8f
#define ABIAS_D_P_NSE_DEFAULT     1e-8f

//测量噪声
#define VELNE_M_NSE_DEFAULT     0.5f
#define VELD_M_NSE_DEFAULT      0.7f

void AP_NavJKF::sys_init()
{
    printf("NavJKF:: init\n");

    // 系统状态初始化
    X[0] = 0;
    X[1] = 0;
    X[2] = 0;
    X[3] = 0.01;
    X[4] = 0.01;
    X[5] = 0.01;

    // 协方差矩阵初始化
    P[0][0] = VELNE_M_NSE_DEFAULT;
    P[1][1] = VELNE_M_NSE_DEFAULT;
    P[2][2] = VELD_M_NSE_DEFAULT;
    P[3][3] = 0.01;
    P[4][4] = 0.01;
    P[5][5] = 0.01;

    // 过程噪声矩阵初始化
    // velocity variance
    Q[0][0] = VELNE_P_NSE_DEFAULT;
    Q[1][1] = VELNE_P_NSE_DEFAULT;
    Q[2][2] = VELD_P_NSE_DEFAULT;
    // acc bias variance
    Q[3][3] = ABIAS_NE_P_NSE_DEFAULT;
    Q[4][4] = ABIAS_NE_P_NSE_DEFAULT;
    Q[5][5] = ABIAS_D_P_NSE_DEFAULT;

    // 观测噪声矩阵初始化
    // measurement noise
    R[0][0]= VELNE_M_NSE_DEFAULT;
    R[1][1]= VELNE_M_NSE_DEFAULT;
    R[2][2]= VELD_M_NSE_DEFAULT; 
}

#define earthRate 0.000072921f // earth rotation rate (rad/sec)

// calculate the NED earth spin vector in rad/sec
void AP_NavJKF::calcEarthRateNED(Vector3f &omega, int32_t latitude) const
{
    float lat_rad = radians(latitude*1.0e-7f);
    omega.x  = earthRate*cosf(lat_rad);
    omega.y  = 0;
    omega.z  = -earthRate*sinf(lat_rad);
}

void AP_NavJKF::sys_predict(float delta_t, Vector3f &acc, Matrix3f &dcm)
{
    float velX = X[0];
    float velY = X[1];
    float velZ = X[2];
    float a_bx = X[3];
    float a_by = X[4];
    float a_bz = X[5];
    float acc_x = acc.x;
    float acc_y = acc.y;
    float acc_z = acc.z;

    Vector3f acc_data = {acc_x-a_bx, acc_y-a_by, acc_z-a_bz};
    Vector3f acc_truth = dcm*acc_data;

    // kf_state_prediction.cpp
    float X_next[6];
    #if 0
    X_next[0] = velX + delta_t*(-a_bx + acc_x);
    X_next[1] = velY + delta_t*(-a_by + acc_y);
    X_next[2] = velZ + delta_t*(-a_bz + GRAVITY_MSS + acc_z);
    #else
    X_next[0] = velX + delta_t*(acc_truth.x);
    X_next[1] = velY + delta_t*(acc_truth.y);
    X_next[2] = velZ + delta_t*(acc_truth.z + GRAVITY_MSS);
    #endif
    X_next[3] = a_bx;
    X_next[4] = a_by;
    X_next[5] = a_bz;

    // kf_covariance_prediction.cpp
    float P_item0 = P[0][3] - P[3][3]*delta_t;
    float P_item1 = -P[3][4]*delta_t;
    float P_item2 = P[0][4] + P_item1;
    float P_item3 = -P[3][5]*delta_t;
    float P_item4 = P[0][5] + P_item3;
    float P_item5 = P[1][3] + P_item1;
    float P_item6 = P[1][4] - P[4][4]*delta_t;
    float P_item7 = -P[4][5]*delta_t;
    float P_item8 = P[1][5] + P_item7;
    float P_item9 = P[2][3] + P_item3;
    float P_item10 = P[2][4] + P_item7;
    float P_item11 = P[2][5] - P[5][5]*delta_t;
    
    
    float P_next[6][6];
    P_next[0][0] = P[0][0] - P[0][3]*delta_t - P_item0*delta_t + Q[0][0]*delta_t;
    P_next[1][0] = P[0][1] - P[0][4]*delta_t - P_item5*delta_t + Q[1][0]*delta_t;
    P_next[2][0] = P[0][2] - P[0][5]*delta_t - P_item9*delta_t + Q[2][0]*delta_t;
    P_next[3][0] = P_item0 + Q[3][0]*delta_t;
    P_next[4][0] = P_item2 + Q[4][0]*delta_t;
    P_next[5][0] = P_item4 + Q[5][0]*delta_t;
    P_next[0][1] = P[0][1] - P[1][3]*delta_t - P_item2*delta_t + Q[0][1]*delta_t;
    P_next[1][1] = P[1][1] - P[1][4]*delta_t - P_item6*delta_t + Q[1][1]*delta_t;
    P_next[2][1] = P[1][2] - P[1][5]*delta_t - P_item10*delta_t + Q[2][1]*delta_t;
    P_next[3][1] = P_item5 + Q[3][1]*delta_t;
    P_next[4][1] = P_item6 + Q[4][1]*delta_t;
    P_next[5][1] = P_item8 + Q[5][1]*delta_t;
    P_next[0][2] = P[0][2] - P[2][3]*delta_t - P_item4*delta_t + Q[0][2]*delta_t;
    P_next[1][2] = P[1][2] - P[2][4]*delta_t - P_item8*delta_t + Q[1][2]*delta_t;
    P_next[2][2] = P[2][2] - P[2][5]*delta_t - P_item11*delta_t + Q[2][2]*delta_t;
    P_next[3][2] = P_item9 + Q[3][2]*delta_t;
    P_next[4][2] = P_item10 + Q[4][2]*delta_t;
    P_next[5][2] = P_item11 + Q[5][2]*delta_t;
    P_next[0][3] = P_item0 + Q[0][3]*delta_t;
    P_next[1][3] = P_item5 + Q[1][3]*delta_t;
    P_next[2][3] = P_item9 + Q[2][3]*delta_t;
    P_next[3][3] = P[3][3] + Q[3][3]*delta_t;
    P_next[4][3] = P[3][4] + Q[4][3]*delta_t;
    P_next[5][3] = P[3][5] + Q[5][3]*delta_t;
    P_next[0][4] = P_item2 + Q[0][4]*delta_t;
    P_next[1][4] = P_item6 + Q[1][4]*delta_t;
    P_next[2][4] = P_item10 + Q[2][4]*delta_t;
    P_next[3][4] = P[3][4] + Q[3][4]*delta_t;
    P_next[4][4] = P[4][4] + Q[4][4]*delta_t;
    P_next[5][4] = P[4][5] + Q[5][4]*delta_t;
    P_next[0][5] = P_item4 + Q[0][5]*delta_t;
    P_next[1][5] = P_item8 + Q[1][5]*delta_t;
    P_next[2][5] = P_item11 + Q[2][5]*delta_t;
    P_next[3][5] = P[3][5] + Q[3][5]*delta_t;
    P_next[4][5] = P[4][5] + Q[4][5]*delta_t;
    P_next[5][5] = P[5][5] + Q[5][5]*delta_t;
    
    memcpy(X, X_next, sizeof(X_next));
    memcpy(P, P_next, sizeof(P_next));
}
void AP_NavJKF::sys_update(Vector3f &gps, Matrix3f &dcm){
    float velX = X[0];
    float velY = X[1];
    float velZ = X[2];
    float a_bx = X[3];
    float a_by = X[4];
    float a_bz = X[5];
    float gpsVelX = gps.x;
    float gpsVelY = gps.y;
    float gpsVelZ = gps.z;
    float R_meas_x = R[0][0];
    float R_meas_y = R[1][1];
    float R_meas_z = R[2][2];
    
    // kf_SP.cpp
    float SP[3][3];
    SP[0][0] = P[0][0] + R_meas_x;
    SP[1][0] = P[0][1];
    SP[2][0] = P[0][2];
    SP[0][1] = P[0][1];
    SP[1][1] = P[1][1] + R_meas_y;
    SP[2][1] = P[1][2];
    SP[0][2] = P[0][2];
    SP[1][2] = P[1][2];
    SP[2][2] = P[2][2] + R_meas_z;
    
    // kf_SPINV.cpp
    float SPINV_item0 = SP[1][1]*SP[2][2];
    float SPINV_item1 = SP[1][2]*SP[2][1];
    float SPINV_item2 = SP[0][1]*SP[1][2];
    float SPINV_item3 = SP[0][2]*SP[2][1];
    float SPINV_item4 = SP[0][1]*SP[2][2];
    float SPINV_item5 = SP[0][2]*SP[1][1];
    float SPINV_item6 = 1.0F/(SPINV_item0*SP[0][0] - SPINV_item1*SP[0][0] + SPINV_item2*SP[2][0] + SPINV_item3*SP[1][0] - SPINV_item4*SP[1][0] - SPINV_item5*SP[2][0]);
    
    
    float SPINV[3][3];
    SPINV[0][0] = SPINV_item6*(SPINV_item0 - SPINV_item1);
    SPINV[1][0] = -SPINV_item6*(SP[1][0]*SP[2][2] - SP[1][2]*SP[2][0]);
    SPINV[2][0] = SPINV_item6*(SP[1][0]*SP[2][1] - SP[1][1]*SP[2][0]);
    SPINV[0][1] = -SPINV_item6*(-SPINV_item3 + SPINV_item4);
    SPINV[1][1] = SPINV_item6*(SP[0][0]*SP[2][2] - SP[0][2]*SP[2][0]);
    SPINV[2][1] = -SPINV_item6*(SP[0][0]*SP[2][1] - SP[0][1]*SP[2][0]);
    SPINV[0][2] = SPINV_item6*(SPINV_item2 - SPINV_item5);
    SPINV[1][2] = -SPINV_item6*(SP[0][0]*SP[1][2] - SP[0][2]*SP[1][0]);
    SPINV[2][2] = SPINV_item6*(SP[0][0]*SP[1][1] - SP[0][1]*SP[1][0]);

    // kf_state_update.cpp
    float X_new_item0 = velX - gpsVelX;
    float X_new_item1 = velY - gpsVelY;
    float X_new_item2 = velZ - gpsVelZ;
    
    
    float X_new[6];
    X_new[0] = -X_new_item0*(P[0][0]*SPINV[0][0] + P[0][1]*SPINV[1][0] + P[0][2]*SPINV[2][0]) - X_new_item1*(P[0][0]*SPINV[0][1] + P[0][1]*SPINV[1][1] + P[0][2]*SPINV[2][1]) - X_new_item2*(P[0][0]*SPINV[0][2] + P[0][1]*SPINV[1][2] + P[0][2]*SPINV[2][2]) + velX;
    X_new[1] = -X_new_item0*(P[0][1]*SPINV[0][0] + P[1][1]*SPINV[1][0] + P[1][2]*SPINV[2][0]) - X_new_item1*(P[0][1]*SPINV[0][1] + P[1][1]*SPINV[1][1] + P[1][2]*SPINV[2][1]) - X_new_item2*(P[0][1]*SPINV[0][2] + P[1][1]*SPINV[1][2] + P[1][2]*SPINV[2][2]) + velY;
    X_new[2] = -X_new_item0*(P[0][2]*SPINV[0][0] + P[1][2]*SPINV[1][0] + P[2][2]*SPINV[2][0]) - X_new_item1*(P[0][2]*SPINV[0][1] + P[1][2]*SPINV[1][1] + P[2][2]*SPINV[2][1]) - X_new_item2*(P[0][2]*SPINV[0][2] + P[1][2]*SPINV[1][2] + P[2][2]*SPINV[2][2]) + velZ;
    X_new[3] = -X_new_item0*(P[0][3]*SPINV[0][0] + P[1][3]*SPINV[1][0] + P[2][3]*SPINV[2][0]) - X_new_item1*(P[0][3]*SPINV[0][1] + P[1][3]*SPINV[1][1] + P[2][3]*SPINV[2][1]) - X_new_item2*(P[0][3]*SPINV[0][2] + P[1][3]*SPINV[1][2] + P[2][3]*SPINV[2][2]) + a_bx;
    X_new[4] = -X_new_item0*(P[0][4]*SPINV[0][0] + P[1][4]*SPINV[1][0] + P[2][4]*SPINV[2][0]) - X_new_item1*(P[0][4]*SPINV[0][1] + P[1][4]*SPINV[1][1] + P[2][4]*SPINV[2][1]) - X_new_item2*(P[0][4]*SPINV[0][2] + P[1][4]*SPINV[1][2] + P[2][4]*SPINV[2][2]) + a_by;
    X_new[5] = -X_new_item0*(P[0][5]*SPINV[0][0] + P[1][5]*SPINV[1][0] + P[2][5]*SPINV[2][0]) - X_new_item1*(P[0][5]*SPINV[0][1] + P[1][5]*SPINV[1][1] + P[2][5]*SPINV[2][1]) - X_new_item2*(P[0][5]*SPINV[0][2] + P[1][5]*SPINV[1][2] + P[2][5]*SPINV[2][2]) + a_bz;

    // kf_covariance_update.cpp
    float P_new_item0 = P[0][0]*SPINV[0][1] + P[0][1]*SPINV[1][1] + P[0][2]*SPINV[2][1];
    float P_new_item1 = P[0][0]*SPINV[0][2] + P[0][1]*SPINV[1][2] + P[0][2]*SPINV[2][2];
    float P_new_item2 = P[0][0]*SPINV[0][0] + P[0][1]*SPINV[1][0] + P[0][2]*SPINV[2][0] - 1;
    float P_new_item3 = P[0][1]*SPINV[0][0] + P[1][1]*SPINV[1][0] + P[1][2]*SPINV[2][0];
    float P_new_item4 = P[0][1]*SPINV[0][2] + P[1][1]*SPINV[1][2] + P[1][2]*SPINV[2][2];
    float P_new_item5 = P[0][1]*SPINV[0][1] + P[1][1]*SPINV[1][1] + P[1][2]*SPINV[2][1] - 1;
    float P_new_item6 = P[0][2]*SPINV[0][0] + P[1][2]*SPINV[1][0] + P[2][2]*SPINV[2][0];
    float P_new_item7 = P[0][2]*SPINV[0][1] + P[1][2]*SPINV[1][1] + P[2][2]*SPINV[2][1];
    float P_new_item8 = P[0][2]*SPINV[0][2] + P[1][2]*SPINV[1][2] + P[2][2]*SPINV[2][2] - 1;
    float P_new_item9 = P[0][3]*SPINV[0][0] + P[1][3]*SPINV[1][0] + P[2][3]*SPINV[2][0];
    float P_new_item10 = P[0][3]*SPINV[0][1] + P[1][3]*SPINV[1][1] + P[2][3]*SPINV[2][1];
    float P_new_item11 = P[0][3]*SPINV[0][2] + P[1][3]*SPINV[1][2] + P[2][3]*SPINV[2][2];
    float P_new_item12 = P[0][4]*SPINV[0][0] + P[1][4]*SPINV[1][0] + P[2][4]*SPINV[2][0];
    float P_new_item13 = P[0][4]*SPINV[0][1] + P[1][4]*SPINV[1][1] + P[2][4]*SPINV[2][1];
    float P_new_item14 = P[0][4]*SPINV[0][2] + P[1][4]*SPINV[1][2] + P[2][4]*SPINV[2][2];
    float P_new_item15 = P[0][5]*SPINV[0][0] + P[1][5]*SPINV[1][0] + P[2][5]*SPINV[2][0];
    float P_new_item16 = P[0][5]*SPINV[0][1] + P[1][5]*SPINV[1][1] + P[2][5]*SPINV[2][1];
    float P_new_item17 = P[0][5]*SPINV[0][2] + P[1][5]*SPINV[1][2] + P[2][5]*SPINV[2][2];
    
    
    float P_new[6][6];
    P_new[0][0] = -P[0][0]*P_new_item2 - P[0][1]*P_new_item0 - P[0][2]*P_new_item1;
    P_new[1][0] = -P[0][0]*P_new_item3 - P[0][1]*P_new_item5 - P[0][2]*P_new_item4;
    P_new[2][0] = -P[0][0]*P_new_item6 - P[0][1]*P_new_item7 - P[0][2]*P_new_item8;
    P_new[3][0] = -P[0][0]*P_new_item9 - P[0][1]*P_new_item10 - P[0][2]*P_new_item11 + P[0][3];
    P_new[4][0] = -P[0][0]*P_new_item12 - P[0][1]*P_new_item13 - P[0][2]*P_new_item14 + P[0][4];
    P_new[5][0] = -P[0][0]*P_new_item15 - P[0][1]*P_new_item16 - P[0][2]*P_new_item17 + P[0][5];
    P_new[0][1] = -P[0][1]*P_new_item2 - P[1][1]*P_new_item0 - P[1][2]*P_new_item1;
    P_new[1][1] = -P[0][1]*P_new_item3 - P[1][1]*P_new_item5 - P[1][2]*P_new_item4;
    P_new[2][1] = -P[0][1]*P_new_item6 - P[1][1]*P_new_item7 - P[1][2]*P_new_item8;
    P_new[3][1] = -P[0][1]*P_new_item9 - P[1][1]*P_new_item10 - P[1][2]*P_new_item11 + P[1][3];
    P_new[4][1] = -P[0][1]*P_new_item12 - P[1][1]*P_new_item13 - P[1][2]*P_new_item14 + P[1][4];
    P_new[5][1] = -P[0][1]*P_new_item15 - P[1][1]*P_new_item16 - P[1][2]*P_new_item17 + P[1][5];
    P_new[0][2] = -P[0][2]*P_new_item2 - P[1][2]*P_new_item0 - P[2][2]*P_new_item1;
    P_new[1][2] = -P[0][2]*P_new_item3 - P[1][2]*P_new_item5 - P[2][2]*P_new_item4;
    P_new[2][2] = -P[0][2]*P_new_item6 - P[1][2]*P_new_item7 - P[2][2]*P_new_item8;
    P_new[3][2] = -P[0][2]*P_new_item9 - P[1][2]*P_new_item10 - P[2][2]*P_new_item11 + P[2][3];
    P_new[4][2] = -P[0][2]*P_new_item12 - P[1][2]*P_new_item13 - P[2][2]*P_new_item14 + P[2][4];
    P_new[5][2] = -P[0][2]*P_new_item15 - P[1][2]*P_new_item16 - P[2][2]*P_new_item17 + P[2][5];
    P_new[0][3] = -P[0][3]*P_new_item2 - P[1][3]*P_new_item0 - P[2][3]*P_new_item1;
    P_new[1][3] = -P[0][3]*P_new_item3 - P[1][3]*P_new_item5 - P[2][3]*P_new_item4;
    P_new[2][3] = -P[0][3]*P_new_item6 - P[1][3]*P_new_item7 - P[2][3]*P_new_item8;
    P_new[3][3] = -P[0][3]*P_new_item9 - P[1][3]*P_new_item10 - P[2][3]*P_new_item11 + P[3][3];
    P_new[4][3] = -P[0][3]*P_new_item12 - P[1][3]*P_new_item13 - P[2][3]*P_new_item14 + P[3][4];
    P_new[5][3] = -P[0][3]*P_new_item15 - P[1][3]*P_new_item16 - P[2][3]*P_new_item17 + P[3][5];
    P_new[0][4] = -P[0][4]*P_new_item2 - P[1][4]*P_new_item0 - P[2][4]*P_new_item1;
    P_new[1][4] = -P[0][4]*P_new_item3 - P[1][4]*P_new_item5 - P[2][4]*P_new_item4;
    P_new[2][4] = -P[0][4]*P_new_item6 - P[1][4]*P_new_item7 - P[2][4]*P_new_item8;
    P_new[3][4] = -P[0][4]*P_new_item9 - P[1][4]*P_new_item10 - P[2][4]*P_new_item11 + P[3][4];
    P_new[4][4] = -P[0][4]*P_new_item12 - P[1][4]*P_new_item13 - P[2][4]*P_new_item14 + P[4][4];
    P_new[5][4] = -P[0][4]*P_new_item15 - P[1][4]*P_new_item16 - P[2][4]*P_new_item17 + P[4][5];
    P_new[0][5] = -P[0][5]*P_new_item2 - P[1][5]*P_new_item0 - P[2][5]*P_new_item1;
    P_new[1][5] = -P[0][5]*P_new_item3 - P[1][5]*P_new_item5 - P[2][5]*P_new_item4;
    P_new[2][5] = -P[0][5]*P_new_item6 - P[1][5]*P_new_item7 - P[2][5]*P_new_item8;
    P_new[3][5] = -P[0][5]*P_new_item9 - P[1][5]*P_new_item10 - P[2][5]*P_new_item11 + P[3][5];
    P_new[4][5] = -P[0][5]*P_new_item12 - P[1][5]*P_new_item13 - P[2][5]*P_new_item14 + P[4][5];
    P_new[5][5] = -P[0][5]*P_new_item15 - P[1][5]*P_new_item16 - P[2][5]*P_new_item17 + P[5][5];


    Vector3f x = {X_new[3], X_new[4], X_new[5]};
    dcm.transpose();
    Vector3f x_ = dcm*x;
    X_new[3] = x_.x;
    X_new[4] = x_.y;
    X_new[5] = x_.z;
    
    memcpy(X, X_new, sizeof(X));
    memcpy(P, P_new, sizeof(P));
}
