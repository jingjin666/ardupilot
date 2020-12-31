#pragma once

#include <AP_Math/matrixN.h>
#include <AP_Math/AP_Math.h>
#include <AP_AHRS/AP_AHRS.h>

class AP_NavJKF {
public:
    AP_NavJKF(const AP_AHRS *ahrs);
    ~AP_NavJKF();

    void sys_init();
    // calculate the NED earth spin vector in rad/sec
    void calcEarthRateNED(Vector3f &omega, int32_t latitude) const;
    void sys_predict(float delta_t, Vector3f &acc, Matrix3f &dcm);
    void sys_update(Vector3f &gps, Matrix3f &dcm);
    void getVelocity(Vector3f &vel){
        vel.x = X[0];
        vel.y = X[1];
        vel.z = X[2];
    }
    
private:
    float X[6];
    float P[6][6];
    float Q[6][6];
    float R[3][3];

    const AP_AHRS *_ahrs;
};
