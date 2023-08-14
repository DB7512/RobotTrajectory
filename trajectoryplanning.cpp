#include "trajectoryplanning.h"
#include "mathfunction.h"
#include <QTextStream>
#include <QFile>
#include <QDebug>

TrajectoryPlanning::TrajectoryPlanning(QObject *parent)
    : QObject{parent}
{

}

bool TrajectoryPlanning::LinePlanning(int peroid, PointInformation point_start, PointInformation point_end,
                                      std::vector<std::vector<float> > &trajectory_point, vector<vector<float> > &trajectory_inf)
{
    float pose_start[6] = {0.0};
    float pose_end[6] = {0.0};
    float Q = 0.0; //displacement
    float trajectory_ftime = 0.0; //trajectory time
    int interpolation_peroid_num = 0; //need period number
    float t = 0.0; //time
    for (int i = 0; i < 6; i++) {
        pose_start[i] = point_start.pose[i];
        pose_end[i] = point_end.pose[i];
    }
    //calculate displacement
    float detla_x = pose_end[0] - pose_start[0]; //displacemengt in x-axis
    float detla_y = pose_end[1] - pose_start[1]; //displacemengt in y-axis
    float detla_z = pose_end[2] - pose_start[2]; //displacemengt in z-axis
    Q = sqrt(pow(detla_x,2) + pow(detla_y,2) + pow(detla_z,2))/1000.0;
    //constraint
    float vmax = point_end.vmax;
    float amax = point_end.amax;
    float jmax = point_end.jmax;
    float v_0 = point_start.v;
    float v_1 = point_end.v;
    //calculate time
    int ret = 0;
    ret = TrajectoryTime(trajectory_ftime, Q, v_0, v_1, vmax, amax, jmax, para);
    qDebug()<<"ret"<<ret<<"time"<<trajectory_ftime;
    //corrent parameters
    CorrentionParameters(para, peroid);
    //vmax = para(9);
    //amax = para(10);
    //jmax = para(13);
    //ret = TrajectoryTime(trajectory_ftime, Q, v_0, v_1, vmax, amax, jmax, para);
    float temp = (para(0) + para(1) + para(2)) / (1.0 / peroid);
    interpolation_peroid_num = round(temp);
    if(para(0) + para(1) + para(2) < 1e-6) return false;
    float lambda[interpolation_peroid_num];//归一化参数
    vector<float> interpolation_point(6);//插值结果
    vector<float> infpoint(3);
    float q_dq[3] = {0.0}; //resulets: displacement, velocity, acceleration
    for (int i = 0; i < interpolation_peroid_num; i++) {
        //init vector
        for(int j = 0; j < 6; j++) {
            if(j < 3) {
                interpolation_point[j] = 0.0;
                infpoint[j] = 0.0;
            } else {
                interpolation_point[j] = pose_start[j];
            }
        }
        lambda[i] = 0.0;
        t = i*(1.0/peroid);
        //calculate trajectory
        TrajectoryCalculation(t, para, q_dq);
        for (int var = 0; var < 3; ++var) {
            infpoint[var] = q_dq[var];
        }
        if(Q!=0)
            lambda[i] = q_dq[0] / Q;
        QFile data("data_new.txt");
        if(!data.open(QIODevice::WriteOnly | QIODevice::Text | QIODevice::Append)) return false;
        QTextStream stream(&data);
        stream<<number<<" "<<lambda[i]<<" "<<q_dq[0]<<" "<<q_dq[1]<<"\n";
        data.close();
        //        qDebug()<<lambda[i]<<q_dq[1];
        interpolation_point[0] = pose_start[0] + detla_x * lambda[i];
        interpolation_point[1] = pose_start[1] + detla_y * lambda[i];
        interpolation_point[2] = pose_start[2] + detla_z * lambda[i];
        trajectory_point.push_back(interpolation_point);
        trajectory_inf.push_back(infpoint);
        //        qDebug()<<interpolation_point[0];
    }
    return true;
}

int TrajectoryPlanning::TrajectoryTime(float &time, float Q, float v_0, float v_1, float vmax, float amax, float jmax, VectorXf &para)
{

    //judge whether the minimun diaplacement is satisfied
    float Tj[2] = {0.0};
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
    float Tj1 = 0.0, Tj2 = 0.0, Ta = 0.0, Td = 0.0, Tv = 0.0;
    float a_lima = 0.0, a_limd = 0.0;
    float vlim = 0.0;
    float jmin = -jmax;
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
    Tv = (Q) / vmax - (Ta / 2) * (1 + v_0 / vmax) - (Td / 2) * (1 + v_1 / vmax);
    //vmax is reached
    if (Tv > 0) {
        vlim = vmax;
        time = Ta + Tv + Td;
    } else {
        Tv = 0.0;
        float Tj = 0.0;
        float delta = 0.0;
        Tj = amax / jmax;
        Tj1 = Tj;
        Tj2 = Tj;
        delta = (pow(amax, 4) / pow(jmax, 2)) + 2 * (pow(v_0, 2) + pow(v_1, 2)) + amax * (4 * Q - 2 * (amax / jmax) * (v_0 + v_1));
        Ta = ((pow(amax, 2) / jmax) - 2.0 * v_0 + sqrt(delta)) / (2.0 * amax);
        Td = ((pow(amax, 2) / jmax) - 2.0 * v_1 + sqrt(delta)) / (2.0 * amax);
        if (Ta < 0 || Td < 0) {
            if (Ta < 0) {
                Ta = 0.0;
                Tj1 = 0.0;
                Td = 2.0 * Q / (v_0 + v_1);
                Tj2 = (jmax * Q - sqrt(jmax * (jmax * pow(Q, 2) + pow(v_1 + v_0, 2) * (v_1 - v_0)))) / (jmax * (v_1 + v_0));
                a_lima = 0;
                a_limd = -jmax * Tj2;
                vlim = v_0;
                time = Ta + Tv + Td;
            } else if (Td < 0) {
                Td = 0.0;
                Tj2 = 0.0;
                Ta = 2.0 * Q / (v_0 + v_1);
                Tj1 = (jmax * Q - sqrt(jmax * (jmax * pow(Q, 2)) - pow(v_1 + v_0, 2) * (v_1 - v_0))) / (jmax * (v_1 + v_0));
                a_lima = jmax * Tj1;
                a_limd = 0.0;
                vlim = v_1;
                //vlim = v_0 + a_lima * (Ta - Tj1);
                time = Ta + Tv + Td;
            }
        } else if (Ta >= 2 * Tj1 && Td >= 2 * Tj2) {
            a_lima = amax;
            a_limd = -amax;
            vlim = v_0 + a_lima * (Ta - Tj);
            time = Ta + Tv + Td;
        } else if (Ta < 2 * Tj1 && Td >= 2 * Tj2) {
            float coeff[6] = { 0 };
            float x[5][2] = { {0} };
            coeff[0] = 1;
            coeff[1] = 2 * amax;
            coeff[2] = 2 * jmax * v_0 + pow(amax, 2);
            coeff[3] = 4 * jmax * amax * v_0;
            coeff[4] = jmax * ((pow(v_0, 2) - pow(v_1, 2)) * jmax + pow(amax, 2) * (v_0 + v_1) - 2 * jmax * amax * Q);
            coeff[5] = 0;
            GetMathInstance().Equation_Root5(coeff, x);
            a_lima = 0;
            for (int j = 0; j < 5; j++) {
                if (x[j][0]<0 || x[j][0]>amax) {
                    x[j][0] = 0;
                }
                if (x[j][0] > a_lima) {
                    a_lima = x[j][0];
                }
            }
            Tj1 = a_lima / jmax;
            Tj2 = amax / jmax;
            Ta = 2 * Tj1;
            vlim = Tj1 * a_lima + v_0;
            Td = (vlim - v_1) / amax + Tj2;
            time = Ta + Tv + Td;
        } else if (Ta >= 2 * Tj1 && Td < 2 * Tj2) {
            float coeff[6] = { 0 };
            float x[5][2] = { {0} };
            coeff[0] = 1;
            coeff[1] = 2 * amax;
            coeff[2] = 2 * jmax * v_1 + pow(amax, 2);
            coeff[3] = 4 * jmax * amax * v_1;
            coeff[4] = jmax * ((pow(v_1, 2) - pow(v_0, 2)) * jmax + pow(amax, 2) * (v_0 + v_1) - 2 * jmax * amax * Q);
            coeff[5] = 0;
            GetMathInstance().Equation_Root5(coeff, x);
            a_limd = 0;
            for (int j = 0; j < 5; j++) {
                if (x[j][0]<0 || x[j][0]>amax) {
                    x[j][0] = 0;
                }
                if (x[j][0] > a_limd) {
                    a_limd = x[j][0];
                }
            }
            Tj1 = amax / jmax;
            Tj2 = a_limd / jmax;
            Td = 2 * Tj2;
            vlim = Tj2 * a_limd + v_1;
            Ta = (vlim - v_0) / amax + Tj1;
            time = Ta + Tv + Td;
        } else if (Ta < 2 * Tj1 && Td < 2 * Tj2) {
            if (fabs(v_0 - v_1) > 0) {
                float coeff[5] = { 0 };
                float x[4][2] = { {0} };
                coeff[0] = jmax * (v_0 - v_1);
                coeff[1] = -2 * pow(jmax, 2) * Q;
                coeff[2] = pow(jmax, 2) * pow(v_0 - v_1, 2);
                coeff[3] = -4 * pow(jmax, 3) * v_0 * Q;
                coeff[4] = pow(jmax, 4) * pow(Q, 2) - pow(jmax, 3) * (v_0 - v_1) * pow(v_0 + v_1, 2);
                GetMathInstance().Equation_Root4(coeff, x);
                a_lima = 0;
                a_limd = 0;
                for (int j = 0; j < 4; j++) {
                    if (x[j][0]<0 || x[j][0]>amax) {
                        x[j][0] = 0;
                    }
                    if (x[j][0] > a_lima) {
                        a_lima = x[j][0];
                    }
                }
                a_limd = sqrt(pow(a_lima, 2) + jmax * (v_0 - v_1));
                Tj1 = a_lima / jmax;
                Tj2 = a_limd / jmax;
                Ta = 2 * Tj1;
                Td = 2 * Tj2;
                vlim = Tj2 * a_limd + v_1;
                time = Ta + Tv + Td;
            } else {
                float coeff[4] = { 0 };
                float x[3] = { 0 };
                coeff[0] = 2 * pow(jmax, 2) * Q;
                coeff[1] = 0;
                coeff[2] = 4 * pow(jmax, 3) * v_0 * Q;
                coeff[3] = -pow(jmax, 4) * pow(Q, 2);
                GetMathInstance().Equation_Root3(coeff, x);
                a_lima = 0;
                a_limd = 0;
                for (int j = 0; j < 3; j++) {
                    if (x[j]<0 || x[j]>amax) {
                        x[j] = 0;
                    }
                    if (x[j] > a_lima) {
                        a_lima = x[j];
                    }
                }
                a_limd = sqrt(pow(a_lima, 2) + jmax * (v_0 - v_1));
                Tj1 = a_lima / jmax;
                Tj2 = a_limd / jmax;
                Ta = 2 * Tj1;
                Td = 2 * Tj2;
                vlim = Tj2 * a_limd + v_1;
                time = Ta + Tv + Td;
            }
        }
    }

    para << Ta, Tv, Td, Tj1, Tj2, 0, Q, v_0, v_1, vlim, amax, a_lima, a_limd, jmax, jmin;
    //classification
    if(Tv > 0.0) {
        if(Ta > 2*Tj1) {
            if(Td > 2*Tj2) {
                m_trajectorytype = TraUniTra;
            } else if(Td > 0.0) {
                m_trajectorytype = TraUniTri;
            } else {
                m_trajectorytype = TraUni;
            }
        } else if(Ta > 0.0) {
            if(Td > 2*Tj2) {
                m_trajectorytype = TriUniTra;
            } else if(Td > 0.0) {
                m_trajectorytype = TriUniTri;
            } else {
                m_trajectorytype = TriUni;
            }
        } else {
            if(Td > 2*Tj2) {
                m_trajectorytype = UniTra;
            } else if(Td > 0.0) {
                m_trajectorytype = UniTri;
            }else {
                m_trajectorytype = Uni;
            }
        }
    } else {
        if(Ta > 2*Tj1) {
            if(Td > 2*Tj2) {
                m_trajectorytype = TraTra;
            } else if(Td > 0.0) {
                m_trajectorytype = TraTri;
            } else {
                m_trajectorytype = TraNone;
            }
        }else if(Ta > 0.0) {
            if(Td > 2*Tj2) {
                m_trajectorytype = TriTra;
            } else if(Td > 0.0) {
                m_trajectorytype = TriTri;
            } else {
                m_trajectorytype = TriNone;
            }
        } else {
            if(Td > 2*Tj2) {
                m_trajectorytype = NoneTra;
            } else if(Td > 0.0) {
                m_trajectorytype = NoneTri;
            } else {
                return -2;
            }
        }
    }
    return 1;
}

void TrajectoryPlanning::CorrentionParameters(VectorXf &para, int peroid)
{
    float Ta, Tv, Td, Tj1, Tj2, Q, v_0, v_1, vlim, amax, alima, alimd, jmax, jmin;
    Ta     = para(0);
    Tv     = para(1);
    Td     = para(2);
    Tj1    = para(3);
    Tj2    = para(4);
    Q      = para(6);
    v_0    = para(7);
    v_1    = para(8);
    vlim   = para(9);
    amax   = para(10);
    alima  = para(11);
    alimd  = para(12);
    jmax   = para(13);
    jmin   = para(14);
    float Te = ceil((Ta + Tv + Td) / (1.0 / peroid)) * (1.0 / peroid) - (Ta + Tv + Td);
    float Tj1_new = Tj1, Tj2_new = Tj2, Ta_new = Ta, Td_new = Td, Tv_new = Tv;
    float deltat_T = 0.0;
    if(m_trajectorytype == TraUniTra || m_trajectorytype == TriUniTri) {
        if(v_1 > v_0) { //extend the acceleration time
            deltat_T    = Te * (vlim + v_0) / (vlim - v_0);
            Tv_new      = Tv - deltat_T;
            Ta_new      = Ta + Te + deltat_T;
            Tj1_new     = Tj1 + 0.5 * (Te + deltat_T);
            jmax        = (vlim - v_0) / ((Ta_new - Tj1_new) * Tj1_new);
            alima       = Tj1_new * jmax;
            vlim        = alima * (Ta_new - Tj1_new) + v_0;
        } else if(v_1 < v_0) { //extend the deceleration time
            deltat_T    = Te * (vlim + v_1) / (vlim - v_1);
            Tv_new      = Tv - deltat_T;
            Td_new      = Td + Te + deltat_T;
            Tj2_new     = Tj2 + 0.5 * (Te + deltat_T);
            jmin        = (v_1 - vlim) / ((Td_new - Tj2_new) * Tj2_new);
            alimd       = Tj2_new * jmin;
            vlim        = v_1 - alimd*(Td_new - Tj2_new);
        } else { //extend the acceleration and deceleration time
            deltat_T    = Te * (v_0 + 2 * vlim + v_1) / (2 * vlim - v_0 - v_1);
            Tv_new      = Tv - deltat_T;
            Ta_new      = Ta + 0.5 * (Te + deltat_T);
            Td_new      = Td + 0.5 * (Te + deltat_T);
            Tj1_new     = Tj1 + (Te + deltat_T) * 0.25;
            Tj2_new     = Tj2 + (Te + deltat_T) * 0.25;
            jmax        = (vlim - v_0) / ((Ta_new - Tj1_new) * Tj1_new);
            jmin        = (v_1 - vlim) / ((Td_new - Tj2_new) * Tj2_new);
            alima       = Tj1_new * jmax;
            alimd       = Tj2_new * jmin;
            vlim        = alima * (Ta_new - Tj1_new) + v_0;
        }
    } else if(m_trajectorytype == TriUniTra || m_trajectorytype == UniTra
             || m_trajectorytype == UniTri) { //extend the deceleration time
        deltat_T    = Te * (vlim + v_1) / (vlim - v_1);
        Tv_new      = Tv - deltat_T;
        Td_new      = Td + Te + deltat_T;
        Tj2_new     = Tj2 + 0.5 * (Te + deltat_T);
        jmin        = (v_1 - vlim) / ((Td_new - Tj2_new) * Tj2_new);
        alimd       = Tj2_new * jmin;
        vlim        = v_1 - alimd*(Td_new - Tj2_new);
    } else if(m_trajectorytype == TraUniTri || m_trajectorytype == TraUni
             || m_trajectorytype == TriUni) { //extend the acceleration time
        deltat_T    = Te * (vlim + v_0) / (vlim - v_0);
        Tv_new      = Tv - deltat_T;
        Ta_new      = Ta + Te + deltat_T;
        Tj1_new     = Tj1 + 0.5 * (Te + deltat_T);
        jmax        = (vlim - v_0) / ((Ta_new - Tj1_new) * Tj1_new);
        alima       = Tj1_new * jmax;
        vlim        = alima * (Ta_new - Tj1_new) + v_0;
    } else if(m_trajectorytype == Uni){
        //decrease v_1
        Td_new = cbrt(8*(Tv * vlim - (Tv + Te)*vlim) / jmin);
        Tj2_new = 0.5 * Td_new;
        alimd = 0.5 * jmin * Td_new;
        Tv_new = Tv - Td_new;
        v_1 = vlim + 0.25 * jmin * pow(Td_new,2);
    } else if(m_trajectorytype == TraTra || m_trajectorytype == TriTri) {
        deltat_T    = 0.0;
        Tv_new      = Tv - deltat_T;
        Ta_new      = Ta + 0.5 * (Te + deltat_T);
        Td_new      = Td + 0.5 * (Te + deltat_T);
        Tj1_new     = Tj1 + (Te + deltat_T) * 0.25;
        Tj2_new     = Tj2 + (Te + deltat_T) * 0.25;
        jmax        = (vlim - v_0) / ((Ta_new - Tj1_new) * Tj1_new);
        jmin        = (v_1 - vlim) / ((Td_new - Tj2_new) * Tj2_new);
        alima       = Tj1_new * jmax;
        alimd       = Tj2_new * jmin;
        vlim        = alima * (Ta_new - Tj1_new) + v_0;
    } else if(m_trajectorytype == TraTri || m_trajectorytype == TriNone || m_trajectorytype == TraNone) {
        Tv_new      = Tv - deltat_T;
        Ta_new      = Ta + Te + deltat_T;
        Tj1_new     = Tj1 + 0.5 * (Te + deltat_T);
        jmax        = (vlim - v_0) / ((Ta_new - Tj1_new) * Tj1_new);
        alima       = Tj1_new * jmax;
        vlim        = alima * (Ta_new - Tj1_new) + v_0;
    } else if(m_trajectorytype == TriTra || m_trajectorytype == NoneTra || m_trajectorytype == NoneTri) {
        deltat_T    = 0.0;
        Tv_new      = Tv - deltat_T;
        Td_new      = Td + Te + deltat_T;
        Tj2_new     = Tj2 + 0.5 * (Te + deltat_T);
        jmin        = (v_1 - vlim) / ((Td_new - Tj2_new) * Tj2_new);
        alimd       = Tj2_new * jmin;
        vlim        = v_1 - alimd*(Td_new - Tj2_new);
    }
    para = VectorXf::Zero(15);
    para << Ta_new, Tv_new, Td_new, Tj1_new, Tj2_new, 0, Q, v_0, v_1, vlim, amax, alima, alimd, jmax, jmin;
}

void TrajectoryPlanning::TrajectoryCalculation(float t, VectorXf para, float q_dq[3])
{
    float Ta, Tv, Td, Tj1, Tj2, p1, p2, v_0, v_1, vlim, amax, alima, alimd, jmax, jmin;
    Ta      = para(0);
    Tv      = para(1);
    Td      = para(2);
    Tj1     = para(3);
    Tj2     = para(4);
    p1      = para(5);
    p2      = para(6);
    v_0     = para(7);
    v_1     = para(8);
    vlim    = para(9);
    amax    = para(10);
    alima   = para(11);
    alimd   = para(12);
    jmax    = para(13);
    jmin    = para(14);
    float T = 0.0;
    float q = 0.0, v = 0.0, a = 0.0, j = 0.0;
    T = Ta + Tv + Td;
    if (t >= 0 && t < Tj1) {
        q = v_0 * t + jmax * pow(t, 3) / 6;
        v = v_0 + jmax * (pow(t, 2) * 0.5);
        a = jmax * t;
        j = jmax;
    } else if (t >= Tj1 && t < Ta - Tj1) {
        q = v_0 * t + (alima / 6) * (3 * pow(t, 2) - 3 * Tj1 * t + pow(Tj1, 2));
        v = v_0 + alima * (t - Tj1 * 0.5);
        a = alima;
        j = 0.0;
    } else if (t >= Ta - Tj1 && t < Ta) {
        q = (vlim + v_0) * (Ta / 2) - vlim * (Ta - t) + jmax * (pow(Ta - t, 3) / 6);
        v = vlim - jmax * (pow(Ta - t, 2) * 0.5);
        a = - jmax * (t -Ta);
        j = jmax;
    } else if (t >= Ta && t < Ta + Tv) {
        q = (vlim + v_0) * (Ta * 0.5) + vlim * (t - Ta);
        v = vlim;
        a = 0.0;
        j = 0.0;
    } else if (t >= Ta + Tv && t < T - Td + Tj2) {
        //q = fabs(p2 - p1) - (vlim + v_1) * (Td * 0.5) + vlim * (t - T + Td) - jmax * (pow(t - T + Td, 3) / 6);
        //v = vlim - jmax * (pow(t - T + Td, 2) * 0.5);
        //a = - jmax * (t - T + Td);
        q = fabs(p2 - p1) - (vlim + v_1) * (Td * 0.5) + vlim * (t - T + Td) + jmin * (pow(t - T + Td, 3) / 6);
        v = vlim + jmin * (pow(t - T + Td, 2) * 0.5);
        a = jmin * (t - T + Td);
        j = jmin;
    } else if (t >= T - Td + Tj2 && t < T - Tj2) {
        q = fabs(p2 - p1) - (vlim + v_1) * (Td * 0.5) + vlim * (t - T + Td) + (alimd / 6) * (3 * pow(t - T + Td, 2) - 3 * Tj2 * (t - T + Td) + pow(Tj2, 2));
        v = vlim + alimd * (t - T + Td - Tj2 * 0.5);
        a = alimd;
        j = 0.0;
    } else {
        //q = fabs(p2 - p1) - v_1 * (T - t) - jmax * (pow(T - t, 3) / 6);
        //v = v_1 + jmax * (pow(t - T, 2) * 0.5);
        //a = jmax * (t - T);
        q = fabs(p2 - p1) - v_1 * (T - t) + jmin * (pow(T - t, 3) / 6);
        v = v_1 - jmin * (pow(t - T, 2) * 0.5);
        a = -jmin * (t - T);
        j = jmin;
    }
    q_dq[0] = q;
    q_dq[1] = v;
    q_dq[2] = a;
}
