#include "velocityplanning.h"

/**
 * @brief VelocityPlanning::TrajectoryTime  计算时间
 * @param time
 * @param Q
 * @param v_0
 * @param v_1
 * @param vmax
 * @param amax
 * @param jmax
 * @param para
 * @return
 */
int VelocityPlanning::TrajectoryTime(double &time, float Q, float v_0, float v_1, float vmax, float amax, float jmax, VectorXd &para)
{
    //judge whether the minimun diaplacement is satisfied
    double Tj[2] = {0.0};
    int index = 0;
    Tj[0] = sqrt(fabs(v_1 - v_0)/jmax);
    Tj[1] = amax/jmax;
    if(Tj[0] > Tj[1]) {
        index = 1;
    }else {
        index = 0;
    }
    if(index == 0) {
        if(Q < Tj[index] * (v_0 + v_1)) {
            return -1;
        }
    }else {
        if(Q < 0.5 * (v_0 + v_1) * (Tj[index] + fabs(v_1 - v_0)/amax)) {
            return -1;
        }
    }
    double Tj1 = 0.0, Tj2 = 0.0, Ta = 0.0, Td = 0.0, Tv = 0.0;
    double a_lima = 0.0, a_limd = 0.0;
    double vlim = 0.0;
    double jmin = -jmax;
    double delta = 1e-6;
    double Sa = 0.0;
    double Sd = 0.0;
    double vl = 0.0;
    double vu = 0.0;
    if(fabs(v_0 - v_1) <= delta && fabs(v_0 - vmax) <= delta) {
        Tv = Q / vmax;
        vlim = vmax;
    } else {
        //amax is not reached
        if ((vmax - v_0) * jmax < pow(amax, 2)) {
            Tj1 = sqrt((vmax - v_0) / jmax);
            Ta = 2 * Tj1;
            a_lima = jmax * Tj1;
        } else {
            Tj1 = amax / jmax;
            Ta = Tj1 + (vmax - v_0) / amax;
            a_lima = amax;
        }
        //amin is not reached
        if ((vmax - v_1) * jmax < pow(amax, 2)) {
            Tj2 = sqrt((vmax - v_1) / jmax);
            Td = 2 * Tj2;
            a_limd = -jmax * Tj2;
        } else {
            Tj2 = amax / jmax;
            Td = Tj2 + (vmax - v_1) / amax;
            a_limd = -amax;
        }
        Sa = 0.5 * (v_0 + vmax) * Ta;
        Sd = 0.5 * (v_1 + vmax) * Td;
        if(Q - Sa - Sd > delta) {
            Tv = (Q - Sa - Sd) / vmax;
            vlim = vmax;
        } else if(fabs(Q - Sa - Sd) <= delta){
            Tv = 0.0;
            vlim = vmax;
        } else {
            Tv = 0.0;
            if(v_0 < v_1) {
                vl = v_1;
            } else {
                vl = v_0;
            }
            vu = vmax;
            vmax = 0.5 * (vl + vu);
            //amax is not reached
            if ((vmax - v_0) * jmax < pow(amax, 2)) {
                Tj1 = sqrt((vmax - v_0) / jmax);
                Ta = 2 * Tj1;
                a_lima = jmax * Tj1;
            } else {
                Tj1 = amax / jmax;
                Ta = Tj1 + (vmax - v_0) / amax;
                a_lima = amax;
            }
            //amin is not reached
            if ((vmax - v_1) * jmax < pow(amax, 2)) {
                Tj2 = sqrt((vmax - v_1) / jmax);
                Td = 2 * Tj2;
                a_limd = -jmax * Tj2;
            } else {
                Tj2 = amax / jmax;
                Td = Tj2 + (vmax - v_1) / amax;
                a_limd = -amax;
            }
            Sa = 0.5 * (v_0 + vmax) * Ta;
            Sd = 0.5 * (v_1 + vmax) * Td;
            while(fabs(Q - Sa - Sd) > delta) {
                if(vu - vl < delta)
                    break;
                if(Q - Sa - Sd > delta) {
                    vl = vmax;
                    vu = vu;
                } else if(Q - Sa - Sd < -delta) {
                    vl = vl;
                    vu = vmax;
                }
                vmax = 0.5 * (vl + vu);
                //amax is not reached
                if ((vmax - v_0) * jmax < pow(amax, 2)) {
                    Tj1 = sqrt((vmax - v_0) / jmax);
                    Ta = 2 * Tj1;
                    a_lima = jmax * Tj1;
                } else {
                    Tj1 = amax / jmax;
                    Ta = Tj1 + (vmax - v_0) / amax;
                    a_lima = amax;
                }
                //amin is not reached
                if ((vmax - v_1) * jmax < pow(amax, 2)) {
                    Tj2 = sqrt((vmax - v_1) / jmax);
                    Td = 2 * Tj2;
                    a_limd = -jmax * Tj2;
                } else {
                    Tj2 = amax / jmax;
                    Td = Tj2 + (vmax - v_1) / amax;
                    a_limd = -amax;
                }
                Sa = 0.5 * (v_0 + vmax) * Ta;
                Sd = 0.5 * (v_1 + vmax) * Td;
            }
            vlim = vmax;
        }
    }
    time = Ta + Tv + Td;
    para << Ta, Tv, Td, Tj1, Tj2, 0, Q, v_0, v_1, vlim, amax, a_lima, a_limd, jmax, jmin;
    //classification
    if(Tv > 0.0) {
        if(Ta > 2*Tj1) {
            if(Td > 2*Tj2) {
                m_accelerationtype = TraZeroTra;
            } else if(Td > 0.0) {
                m_accelerationtype = TraZeroTri;
            } else {
                m_accelerationtype = TraZero;
            }
        } else if(Ta > 0.0) {
            if(Td > 2*Tj2) {
                m_accelerationtype = TriZeroTra;
            } else if(Td > 0.0) {
                m_accelerationtype = TriZeroTri;
            } else {
                m_accelerationtype = TriZero;
            }
        } else {
            if(Td > 2*Tj2) {
                m_accelerationtype = ZeroTra;
            } else if(Td > 0.0) {
                m_accelerationtype = ZeroTri;
            }else {
                m_accelerationtype = Zero;
            }
        }
    } else {
        if(Ta > 2*Tj1) {
            if(Td > 2*Tj2) {
                m_accelerationtype = TraTra;
            } else if(Td > 0.0) {
                m_accelerationtype = TraTri;
            } else {
                m_accelerationtype = TraNone;
            }
        }else if(Ta > 0.0) {
            if(Td > 2*Tj2) {
                m_accelerationtype = TriTra;
            } else if(Td > 0.0) {
                m_accelerationtype = TriTri;
            } else {
                m_accelerationtype = TriNone;
            }
        } else {
            if(Td > 2*Tj2) {
                m_accelerationtype = NoneTra;
            } else if(Td > 0.0) {
                m_accelerationtype = NoneTri;
            } else {
                return -2;
            }
        }
    }
    return 1;
}
