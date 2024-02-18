#include "trajectoryplanning.h"
#include "mathfunction.h"
#include <QTextStream>
#include <QFile>
#include <QDebug>
#include "cmath"
#include <algorithm>
#include <ctime>


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
    float vmax = point_end.maxVelocity;
    float amax = point_end.maxAcceleration;
    float jmax = point_end.maxJerk;
    float v_0 = point_start.velocity;
    float v_1 = point_end.velocity;
    //calculate time
    int ret = 0;
    ret = TrajectoryTime4(trajectory_ftime, Q, v_0, v_1, vmax, amax, jmax, m_parameters);
    qDebug()<<"ret"<<ret<<"time"<<trajectory_ftime;
    if(ret < 0) return -1;
    //corrent parameters
    CorrentionParameters(m_parameters, peroid);
    float temp = (m_parameters(0) + m_parameters(1) + m_parameters(2)) / (1.0 / peroid);
    interpolation_peroid_num = round(temp);
    if(m_parameters(0) + m_parameters(1) + m_parameters(2) < 1e-6) return false;
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
        TrajectoryCalculation(t, m_parameters, q_dq);
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
        //        interpolation_point[0] = pose_start[0] + detla_x * lambda[i];
        //        interpolation_point[1] = pose_start[1] + detla_y * lambda[i];
        //        interpolation_point[2] = pose_start[2] + detla_z * lambda[i];
        //        trajectory_point.push_back(interpolation_point);
        //        trajectory_inf.push_back(infpoint);
        //        qDebug()<<interpolation_point[0];
    }
    return true;
}

bool TrajectoryPlanning::TimePlanning(int peroid, float Q, float vmax, float amax, float jmax, float v_0, float v_1, std::vector<float> &trajectory_point, vector<vector<float> > &trajectory_inf)
{
    float trajectory_ftime = 0.0; //trajectory time
    //    double trajectory_ftime = 0.0; //trajectory time
    int interpolation_peroid_num = 0; //need period number
    float t = 0.0; //time
    //calculate time
    int ret = 0;
    ret = TrajectoryTime4(trajectory_ftime, Q, v_0, v_1, vmax, amax, jmax, m_parameters);
    qDebug()<<"ret"<<ret<<"time"<<trajectory_ftime;
    if(ret < 0) return -1;
    double temp = (m_parameters(0) + m_parameters(1) + m_parameters(2)) / (1.0 / peroid);
    interpolation_peroid_num = round(temp);
    if(m_parameters(0) + m_parameters(1) + m_parameters(2) < 1e-6) return false;
    double lambda[interpolation_peroid_num];//归一化参数
    vector<float> infpoint(3);
    float q_dq[3] = {0.0}; //resulets: displacement, velocity, acceleration
    for (int i = 0; i < interpolation_peroid_num - 1; i++) {
        //init vector
        lambda[i] = 0.0;
        t = i*(1.0/peroid);
        //calculate trajectory
        TrajectoryCalculation(t, m_parameters, q_dq);
        for (int var = 0; var < 3; ++var) {
            infpoint[var] = q_dq[var];
        }
        if(Q!=0) {
            lambda[i] = q_dq[0] / Q;
            if(lambda[i]<-1e-5 || (lambda[i] - 1.0) > 1e-5) {
                int mamama =0;
            }
        }
        QFile data("data.txt");
        if(!data.open(QIODevice::WriteOnly | QIODevice::Text | QIODevice::Append)) return false;
        QTextStream stream(&data);
        stream<<number<<" "<<lambda[i]<<" "<<q_dq[0]<<" "<<q_dq[1]<<"\n";
        data.close();
        //        qDebug()<<lambda[i]<<q_dq[1];
        trajectory_point.push_back(q_dq[0]);
        trajectory_inf.push_back(infpoint);
    }
    return true;
}

int TrajectoryPlanning::TrajectoryTime1(float &time, float Q, float v_0, float v_1, float vmax, float amax, float jmax, VectorXf &para)
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

int TrajectoryPlanning::TrajectoryTime2(float &time, float Q, float v_0, float v_1, float vmax, float amax, float jmax, VectorXf &para)
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
    float lambda = 0.9;
    float delta;
    //vmax is reached
    if (Tv > -1e-6) {
        vlim = vmax;
        time = Ta + Tv + Td;
    } else {
        Tv = 0.0;
        delta = 0.0;
        Tj1 = amax / jmax;
        Tj2 = amax / jmax;
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
            vlim = v_0 + a_lima * (Ta - Tj1);
            time = Ta + Tv + Td;
        } else {
            while(Ta < 2 * amax / jmax || Td < 2 * amax / jmax) {
                amax = lambda*amax;
                Tv = 0.0;
                delta = 0.0;
                Tj1 = amax / jmax;
                Tj2 = amax / jmax;
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
                        break;
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
                        break;
                    }
                } else if (Ta >= 2 * Tj1 && Td >= 2 * Tj2) {
                    a_lima = amax;
                    a_limd = -amax;
                    vlim = v_0 + a_lima * (Ta - Tj1);
                    time = Ta + Tv + Td;
                    break;
                }
            }
        }
    }
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

int TrajectoryPlanning::TrajectoryTime3(float &time, float Q, float v_0, float v_1, float vmax, float amax, float jmax, VectorXf &para)
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
    float alima = 0.0, alimd = 0.0;
    float vlim = 0.0;
    float jmin = -jmax;
    float delta = 1e-5;
    //amax is not reached
    if ((vmax - v_0) * jmax < pow(amax, 2)) {
        Tj1 = sqrt((vmax - v_0) / jmax);
        Ta = 2 * Tj1;
        alima = jmax * Tj1;
    } else {
        Tj1 = amax / jmax;
        Ta = Tj1 + (vmax - v_0) / amax;
        alima = amax;
    }
    //amin is not reached
    if ((vmax - v_1) * jmax < pow(amax, 2)) {
        Tj2 = sqrt((vmax - v_1) / jmax);
        Td = 2 * Tj2;
        alimd = -jmax * Tj2;
    } else {
        Tj2 = amax / jmax;
        Td = Tj2 + (vmax - v_1) / amax;
        alimd = -amax;
    }
    float Sad = 0.5 * (v_0 + vmax) * Ta + 0.5 * (vmax + v_1) * Td;
    float F = Q - Sad;
    if(F > delta) {
        Tv = F / vmax;
        vlim = vmax;
    } else if(fabs(F) < delta) {
        Tv = 0.0;
        vlim = vmax;
    } else {
        Tv = 0.0;
        if(v_1 - v_0 > delta) {
            //Ta < Td and vlim < vmax
            //if alimd just reach amin and Td = 2 * Tj2, trajectory type is TraTri
            vlim = v_1 + pow(amax,2) / jmax;
            Tj1 = amax / jmax;
            Ta = (vlim - v_0) / amax + Tj1;
            Tj2 = amax / jmax;
            Td = 2 * Tj2;
            Sad = 0.5 * (v_0 + vlim) * Ta + 0.5 * (vlim + v_1) * Td;
            F = Q - Sad;
            if(F > delta) {
                //Td > 2 * Tj2, trajectory type is TraTra, vmax > v > vlim
                float a = 2 * jmax;
                float b = 2 * pow(amax,2);
                float c = pow(amax,2) * (v_0 + v_1) - jmax * (pow(v_0,2) + pow(v_1,2)) - 2 * Q * amax * jmax;
                vlim = (-b + sqrt(pow(b,2) - 4*a*c)) / (2*a);
                Tj1 = amax / jmax;
                Ta = (vlim - v_0) / amax + Tj1;
                Tj2 = amax / jmax;
                Td = (vlim - v_1) / amax + Tj2;
                alima = amax;
                alimd = -amax;
            } else if(fabs(F) < delta) {
                alima = amax;
                alimd = -amax;
            } else {
                //alimd < amin, v_0 + pow(amax,2) / jmax < vmax < v_1 + pow(amax,2) / jmax
                //if alima just reach amax and Ta = 2*Tj2
                vlim = v_0 + pow(amax,2) / jmax;
                float vu = v_1 + pow(amax,2) / jmax;
                float vl = v_0 + pow(amax,2) / jmax;
                if(vlim - v_1 <= delta) {
                    //alima = amax, trajectory type is TraTri
                    vlim = 0.5 * (vu + vl);
                    Tj1 = amax / jmax;
                    Ta = (vlim - v_0) / amax + Tj1;
                    Tj2 = sqrt((vlim - v_1) / jmax);
                    Td = Tj2 * 2;
                    Sad = 0.5 * ((vlim + v_0) * Ta + (vlim + v_1) * Td);
                    F = Q - Sad;
                    while(fabs(F) > delta) {
                        if(F > delta) {
                            vl = vlim;
                            vlim = 0.5 * (vl + vu);
                            Tj1 = amax / jmax;
                            Ta = (vlim - v_0) / amax + Tj1;
                            Tj2 = sqrt((vlim - v_1) / jmax);
                            Td = Tj2 * 2;
                            Sad = 0.5 * ((vlim + v_0) * Ta + (vlim + v_1) * Td);
                            F = Q - Sad;
                        } else if(F < -delta) {
                            vu = vlim;
                            vlim = 0.5 * (vl + vu);
                            Tj1 = amax / jmax;
                            Ta = (vlim - v_0) / amax + Tj1;
                            Tj2 = sqrt((vlim - v_1) / jmax);
                            Td = Tj2 * 2;
                            Sad = 0.5 * ((vlim + v_0) * Ta + (vlim + v_1) * Td);
                            F = Q - Sad;
                        }
                    }
                } else {
                    //trajectory type is TriTri
                    Tj1 = amax / jmax;
                    Ta = 2 * Tj1;
                    Tj2 = sqrt((vlim - v_1) / jmax);
                    Td = 2 * Tj2;
                    Sad = 0.5 * ((vlim + v_0) * Ta + (vlim + v_1) * Td);
                    if(Sad - Q > - delta) {
                        vl = v_1;
                        vu = v_0 + pow(amax,2) / jmax;
                        vlim = 0.5 * (vu + vl);
                        Tj1 = sqrt((vlim - v_0) / jmax);
                        Ta = Tj1 * 2;
                        Tj2 = sqrt((vlim - v_1) / jmax);
                        Td = Tj2 * 2;
                        Sad = 0.5 * ((vlim + v_0) * Ta + (vlim + v_1) * Td);
                        F = Q - Sad;
                        while(fabs(F) > delta) {
                            if(F > delta) {
                                vl = vlim;
                                vlim = 0.5 * (vl + vu);
                                Tj1 = sqrt((vlim - v_0) / jmax);
                                Ta = Tj1 * 2;
                                Tj2 = sqrt((vlim - v_1) / jmax);
                                Td = Tj2 * 2;
                                Sad = 0.5 * ((vlim + v_0) * Ta + (vlim + v_1) * Td);
                                F = Q - Sad;
                            } else if(F < -delta) {
                                vu = vlim;
                                vlim = 0.5 * (vl + vu);
                                Tj1 = sqrt((vlim - v_0) / jmax);
                                Ta = Tj1 * 2;
                                Tj2 = sqrt((vlim - v_1) / jmax);
                                Td = Tj2 * 2;
                                Sad = 0.5 * ((vlim + v_0) * Ta + (vlim + v_1) * Td);
                                F = Q - Sad;
                            }
                        }
                    } else {
                        //alima = amax, trajectory type is TraTri
                        vlim = 0.5 * (vu + vl);
                        Tj1 = amax / jmax;
                        Ta = (vlim - v_0) / amax + Tj1;
                        Tj2 = sqrt((vlim - v_1) / jmax);
                        Td = Tj2 * 2;
                        Sad = 0.5 * ((vlim + v_0) * Ta + (vlim + v_1) * Td);
                        F = Q - Sad;
                        while(fabs(F) > delta) {
                            if(F > delta) {
                                vl = vlim;
                                vlim = 0.5 * (vl + vu);
                                Tj1 = amax / jmax;
                                Ta = (vlim - v_0) / amax + Tj1;
                                Tj2 = sqrt((vlim - v_1) / jmax);
                                Td = Tj2 * 2;
                                Sad = 0.5 * ((vlim + v_0) * Ta + (vlim + v_1) * Td);
                                F = Q - Sad;
                            } else if(F < -delta) {
                                vu = vlim;
                                vlim = 0.5 * (vl + vu);
                                Tj1 = amax / jmax;
                                Ta = (vlim - v_0) / amax + Tj1;
                                Tj2 = sqrt((vlim - v_1) / jmax);
                                Td = Tj2 * 2;
                                Sad = 0.5 * ((vlim + v_0) * Ta + (vlim + v_1) * Td);
                                F = Q - Sad;
                            }
                        }
                    }
                }
            }
        } else {
            GetMathInstance().SwapValue(v_0,v_1);
            //Ta < Td and vlim < vmax
            //if alimd just reach amin and Td = 2 * Tj2, trajectory type is TraTri
            vlim = v_1 + pow(amax,2) / jmax;
            Tj1 = amax / jmax;
            Ta = (vlim - v_0) / amax + Tj1;
            Tj2 = amax / jmax;
            Td = 2 * Tj2;
            Sad = 0.5 * (v_0 + vlim) * Ta + 0.5 * (vlim + v_1) * Td;
            if(Sad - Q < -delta) {
                //Td > 2 * Tj2, trajectory type is TraTra, vmax > v > vlim
                float a = 2 * jmax;
                float b = 2 * pow(amax,2);
                float c = pow(amax,2) * (v_0 + v_1) - jmax * (pow(v_0,2) + pow(v_1,2)) - 2 * Q * amax * jmax;
                vlim = (-b + sqrt(pow(b,2) - 4*a*c)) / (2*a);
                Tj1 = amax / jmax;
                Ta = (vlim - v_0) / amax + Tj1;
                Tj2 = amax / jmax;
                Td = (vlim - v_1) / amax + Tj2;
                alima = amax;
                alimd = -amax;
            } else if(fabs(Sad - Q) < -delta) {
                alima = amax;
                alimd = -amax;
            } else {
                //alimd < amin, v_0 + pow(amax,2) / jmax < vmax < v_1 + pow(amax,2) / jmax
                //if alima just reach amax and Ta = 2*Tj2
                vlim = v_0 + pow(amax,2) / jmax;
                float F = 0.0;
                float vu = v_1 + pow(amax,2) / jmax;
                float vl = v_0 + pow(amax,2) / jmax;
                if(vlim - v_1 <= delta) {
                    //alima = amax, trajectory type is TraTri
                    vlim = 0.5 * (vu + vl);
                    Tj1 = amax / jmax;
                    Ta = (vlim - v_0) / amax + Tj1;
                    Tj2 = sqrt((vlim - v_1) / jmax);
                    Td = Tj2 * 2;
                    Sad = 0.5 * ((vlim + v_0) * Ta + (vlim + v_1) * Td);
                    F = Q - Sad;
                    while(fabs(F) > delta) {
                        if(F > delta) {
                            vl = vlim;
                            vlim = 0.5 * (vl + vu);
                            Tj1 = amax / jmax;
                            Ta = (vlim - v_0) / amax + Tj1;
                            Tj2 = sqrt((vlim - v_1) / jmax);
                            Td = Tj2 * 2;
                            Sad = 0.5 * ((vlim + v_0) * Ta + (vlim + v_1) * Td);
                            F = Q - Sad;
                        } else if(F < -delta) {
                            vu = vlim;
                            vlim = 0.5 * (vl + vu);
                            Tj1 = amax / jmax;
                            Ta = (vlim - v_0) / amax + Tj1;
                            Tj2 = sqrt((vlim - v_1) / jmax);
                            Td = Tj2 * 2;
                            Sad = 0.5 * ((vlim + v_0) * Ta + (vlim + v_1) * Td);
                            F = Q - Sad;
                        }
                    }
                } else {
                    Tj1 = amax / jmax;
                    Ta = 2 * Tj1;
                    Tj2 = sqrt((vlim - v_1) / jmax);
                    Td = 2 * Tj2;
                    Sad = 0.5 * ((vlim + v_0) * Ta + (vlim + v_1) * Td);
                    if(Sad - Q > - delta) {
                        vl = v_1;
                        vu = v_0 + pow(amax,2) / jmax;
                        //trajectory type is TriTri
                        vlim = 0.5 * (vu + vl);
                        Tj1 = sqrt((vlim - v_0) / jmax);
                        Ta = Tj1 * 2;
                        Tj2 = sqrt((vlim - v_1) / jmax);
                        Td = Tj2 * 2;
                        Sad = 0.5 * ((vlim + v_0) * Ta + (vlim + v_1) * Td);
                        F = Q - Sad;
                        while(fabs(F) > delta) {
                            if(F > delta) {
                                vl = vlim;
                                vlim = 0.5 * (vl + vu);
                                Tj1 = sqrt((vlim - v_0) / jmax);
                                Ta = Tj1 * 2;
                                Tj2 = sqrt((vlim - v_1) / jmax);
                                Td = Tj2 * 2;
                                Sad = 0.5 * ((vlim + v_0) * Ta + (vlim + v_1) * Td);
                                F = Q - Sad;
                            } else if(F < -delta) {
                                vu = vlim;
                                vlim = 0.5 * (vl + vu);
                                Tj1 = sqrt((vlim - v_0) / jmax);
                                Ta = Tj1 * 2;
                                Tj2 = sqrt((vlim - v_1) / jmax);
                                Td = Tj2 * 2;
                                Sad = 0.5 * ((vlim + v_0) * Ta + (vlim + v_1) * Td);
                                F = Q - Sad;
                            }
                        }
                    } else {
                        //alima = amax, trajectory type is TraTri
                        vlim = 0.5 * (vu + vl);
                        Tj1 = amax / jmax;
                        Ta = (vlim - v_0) / amax + Tj1;
                        Tj2 = sqrt((vlim - v_1) / jmax);
                        Td = Tj2 * 2;
                        Sad = 0.5 * ((vlim + v_0) * Ta + (vlim + v_1) * Td);
                        F = Q - Sad;
                        while(fabs(F) > delta) {
                            if(F > delta) {
                                vl = vlim;
                                vlim = 0.5 * (vl + vu);
                                Tj1 = amax / jmax;
                                Ta = (vlim - v_0) / amax + Tj1;
                                Tj2 = sqrt((vlim - v_1) / jmax);
                                Td = Tj2 * 2;
                                Sad = 0.5 * ((vlim + v_0) * Ta + (vlim + v_1) * Td);
                                F = Q - Sad;
                            } else if(F < -delta) {
                                vu = vlim;
                                vlim = 0.5 * (vl + vu);
                                Tj1 = amax / jmax;
                                Ta = (vlim - v_0) / amax + Tj1;
                                Tj2 = sqrt((vlim - v_1) / jmax);
                                Td = Tj2 * 2;
                                Sad = 0.5 * ((vlim + v_0) * Ta + (vlim + v_1) * Td);
                                F = Q - Sad;
                            }
                        }
                    }
                }
            }
            GetMathInstance().SwapValue(v_0,v_1);
            GetMathInstance().SwapValue(Ta,Td);
            GetMathInstance().SwapValue(Tj1,Tj2);
            alima = -alima;
            alimd = -alimd;
            GetMathInstance().SwapValue(alima,alimd);
        }
    }
    time = Ta + Tv + Td;
    para << Ta, Tv, Td, Tj1, Tj2, 0, Q, v_0, v_1, vlim, amax, alima, alimd, jmax, jmin;
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

int TrajectoryPlanning::TrajectoryTime4(float &time, float Q, float v_0, float v_1, float vmax, float amax, float jmax, VectorXf &para)
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
    float delta = 1e-5;
    float Sa = 0.0;
    float Sd = 0.0;
    float vl = 0.0;
    float vu = 0.0;
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

int TrajectoryPlanning::TrajectoryTime5(double &time, float Q, float v_0, float v_1, float vmax, float amax, float jmax, VectorXd &para)
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

void TrajectoryPlanning::CorrentionParameters(VectorXf &para, int peroid)
{
    float Ta, Tv, Td, Tj1, Tj2, Q, v_0, v_1, vlim, amax, amin, alima, alimd, jmax, jmin;
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
    float delta_T = 0.0;
    float Tm = Te + Ta + Tv + Td;
    if(m_trajectorytype == TraZeroTra || m_trajectorytype == TriZeroTri) {
        if((v_1 - v_0) > 1e-5) { //extend the acceleration time
            delta_T    = Te * (vlim + v_0) / (vlim - v_0);
            Tv_new      = Tv - delta_T;
            if(Tv_new < -1e-5) {
                //decrease vlim
                Tv_new  = Tv;
                Ta_new  = Ta + Te;
                Tj1_new = Tj1 + 0.5 * Te;
                vlim    = (2 * Q - Ta_new * v_0 - Td_new * v_1) / (2 * Tm - Ta_new - Td_new);
                alima   = (vlim - v_0) / (Ta_new - Tj1_new);
                alimd   = (v_1 - vlim) / (Td_new - Tj2_new);
                jmax    = alima / Tj1_new;
                jmin    = alimd / Tj2_new;
            } else {
                Ta_new  = Ta + Te + delta_T;
                Tj1_new = Tj1 + 0.5 * (Te + delta_T);
                jmax    = (vlim - v_0) / ((Ta_new - Tj1_new) * Tj1_new);
                alima   = Tj1_new * jmax;
                vlim    = alima * (Ta_new - Tj1_new) + v_0;
            }
        } else if((v_1 - v_0) < -1e-5) { //extend the deceleration time
            delta_T    = Te * (vlim + v_1) / (vlim - v_1);
            Tv_new      = Tv - delta_T;
            if(Tv_new < -1e-5) {
                //decrease vlim
                Tv_new  = Tv;
                Td_new  = Td + Te;
                Tj2_new = Tj2 + 0.5 * Te;
                vlim    = (2 * Q - Ta_new * v_0 - Td_new * v_1) / (2 * Tm - Ta_new - Td_new);
                alima   = (vlim - v_0) / (Ta_new - Tj1_new);
                alimd   = (v_1 - vlim) / (Td_new - Tj2_new);
                jmax    = alima / Tj1_new;
                jmin    = alimd / Tj2_new;
            } else {
                Td_new  = Td + Te + delta_T;
                Tj2_new = Tj2 + 0.5 * (Te + delta_T);
                jmin    = (v_1 - vlim) / ((Td_new - Tj2_new) * Tj2_new);
                alimd   = Tj2_new * jmin;
                vlim    = v_1 - alimd*(Td_new - Tj2_new);
            }
        } else { //extend the acceleration and deceleration time
            delta_T    = Te * (v_0 + 2 * vlim + v_1) / (2 * vlim - v_0 - v_1);
            Tv_new      = Tv - delta_T;
            if(Tv_new < -1e-5) {
                //decrease vlim
                Tv_new  = Tv;
                Ta_new  = Ta + 0.5 * Te;
                Tj1_new = Tj1 + 0.25 * Te;
                Td_new  = Td + 0.5 * Te;
                Tj2_new = Tj2 + 0.25 * Te;
                vlim    = (2 * Q - Ta_new * v_0 - Td_new * v_1) / (2 * Tm - Ta_new - Td_new);
                alima   = (vlim - v_0) / (Ta_new - Tj1_new);
                alimd   = (v_1 - vlim) / (Td_new - Tj2_new);
                jmax    = alima / Tj1_new;
                jmin    = alimd / Tj2_new;
            } else {
                Ta_new      = Ta + 0.5 * (Te + delta_T);
                Td_new      = Td + 0.5 * (Te + delta_T);
                Tj1_new     = Tj1 + (Te + delta_T) * 0.25;
                Tj2_new     = Tj2 + (Te + delta_T) * 0.25;
                jmax        = (vlim - v_0) / ((Ta_new - Tj1_new) * Tj1_new);
                jmin        = (v_1 - vlim) / ((Td_new - Tj2_new) * Tj2_new);
                alima       = Tj1_new * jmax;
                alimd       = Tj2_new * jmin;
                vlim        = alima * (Ta_new - Tj1_new) + v_0;
            }
        }
    } else if(m_trajectorytype == TriZeroTra || m_trajectorytype == ZeroTra
               || m_trajectorytype == ZeroTri) { //extend the deceleration time
        delta_T    = Te * (vlim + v_1) / (vlim - v_1);
        Tv_new      = Tv - delta_T;
        if(Tv_new < -1e-5) {
            if(m_trajectorytype == TriZeroTra) {
                //decrease vlim
                Tv_new  = Tv;
                Td_new  = Td + Te;
                Tj2_new = Tj2 + 0.5 * Te;
                vlim    = (2 * Q - Ta_new * v_0 - Td_new * v_1) / (2 * Tm - Ta_new - Td_new);
                alima   = (vlim - v_0) / (Ta_new - Tj1_new);
                alimd   = (v_1 - vlim) / (Td_new - Tj2_new);
                jmax    = alima / Tj1_new;
                jmin    = alimd / Tj2_new;
            } else {
                //decrease vlim and v_0
                Tv_new  = Tv;
                Td_new  = Td + Te;
                Tj2_new = Tj2 + 0.5 * Te;
                vlim    = (2 * Q - Td_new * v_1) / (2 * Tv_new + Td_new);
                alimd   = (v_1 - vlim) / (Td_new - Tj2_new);
                jmin    = alimd / Tj2_new;
                v_0 = vlim;
            }
        } else {
            Td_new      = Td + Te + delta_T;
            Tj2_new     = Tj2 + 0.5 * (Te + delta_T);
            jmin        = (v_1 - vlim) / ((Td_new - Tj2_new) * Tj2_new);
            alimd       = Tj2_new * jmin;
            vlim        = v_1 - alimd*(Td_new - Tj2_new);
        }
    } else if(m_trajectorytype == TraZeroTri || m_trajectorytype == TraZero
               || m_trajectorytype == TriZero) { //extend the acceleration time
        delta_T    = Te * (vlim + v_0) / (vlim - v_0);
        Tv_new      = Tv - delta_T;
        if(Tv_new < -1e-5) {
            if(m_trajectorytype == TraZeroTri) {
                //decrease vlim
                Tv_new  = Tv;
                Ta_new  = Ta + Te;
                Tj1_new = Tj1 + 0.5 * Te;
                vlim    = (2 * Q - Ta_new * v_0 - Td_new * v_1) / (2 * Tm - Ta_new - Td_new);
                alima   = (vlim - v_0) / (Ta_new - Tj1_new);
                alimd   = (v_1 - vlim) / (Td_new - Tj2_new);
                jmax    = alima / Tj1_new;
                jmin    = alimd / Tj2_new;
            } else {
                //decrease vlim and v_1
                Tv_new  = Tv;
                Ta_new  = Ta + Te;
                Tj1_new = Tj1 + 0.5 * Te;
                vlim    = (2 * Q - Ta_new * v_0) / (2 * Tv_new + Ta_new);
                alima   = (vlim - v_0) / (Ta_new - Tj1_new);
                alimd   = (v_1 - vlim) / (Td_new - Tj2_new);
                jmax    = alima / Tj1_new;
                jmin    = alimd / Tj2_new;
                v_1 = vlim;
            }
        } else {
            Ta_new      = Ta + Te + delta_T;
            Tj1_new     = Tj1 + 0.5 * (Te + delta_T);
            jmax        = (vlim - v_0) / ((Ta_new - Tj1_new) * Tj1_new);
            alima       = Tj1_new * jmax;
            vlim        = alima * (Ta_new - Tj1_new) + v_0;
        }
    } else if(m_trajectorytype == Zero){
        //decrease v_0 and v_1
        Td_new = cbrt(8*(Tv * vlim - (Tv + Te)*vlim) / jmin);
        Tv_new = Tv - Td_new;
        if(Tv_new < -1e-5) {
            Tv_new = Tv + Te;
            vlim = Q / Tv_new;
            v_0 = vlim;
            v_1 = vlim;
            Td_new = Td;
        } else {
            Tj2_new = 0.5 * Td_new;
            alimd = 0.5 * jmin * Td_new;
            v_1 = vlim + 0.25 * jmin * pow(Td_new,2);
        }
    } else if(m_trajectorytype == TraTra || m_trajectorytype == TriTri) {
        delta_T    = 0.0;
        Tv_new      = Tv - delta_T;
        Ta_new      = Ta + 0.5 * (Te + delta_T);
        Td_new      = Td + 0.5 * (Te + delta_T);
        Tj1_new     = Tj1 + (Te + delta_T) * 0.25;
        Tj2_new     = Tj2 + (Te + delta_T) * 0.25;
        jmax        = (vlim - v_0) / ((Ta_new - Tj1_new) * Tj1_new);
        jmin        = (v_1 - vlim) / ((Td_new - Tj2_new) * Tj2_new);
        alima       = Tj1_new * jmax;
        alimd       = Tj2_new * jmin;
        vlim        = alima * (Ta_new - Tj1_new) + v_0;
    } else if(m_trajectorytype == TraTri) {
        Tv_new      = Tv - delta_T;
        Ta_new      = Ta + Te + delta_T;
        Tj1_new     = Tj1 + 0.5 * (Te + delta_T);
        jmax        = (vlim - v_0) / ((Ta_new - Tj1_new) * Tj1_new);
        alima       = Tj1_new * jmax;
        vlim        = alima * (Ta_new - Tj1_new) + v_0;
    } else if(m_trajectorytype == TriNone || m_trajectorytype == TraNone) {
        //decrease v_1
        Tv_new  = Tv;
        Ta_new  = Ta + Te;
        Tj1_new = Tj1 + 0.5 * Te;
        vlim    = 2 * Q / Ta_new - v_0;
        alima   = (vlim - v_0) / (Ta_new - Tj1_new);
        jmax    = alima / Tj1_new;
        v_1     = vlim;
    } else if(m_trajectorytype == TriTra) {
        delta_T    = 0.0;
        Tv_new      = Tv - delta_T;
        Td_new      = Td + Te + delta_T;
        Tj2_new     = Tj2 + 0.5 * (Te + delta_T);
        jmin        = (v_1 - vlim) / ((Td_new - Tj2_new) * Tj2_new);
        alimd       = Tj2_new * jmin;
        vlim        = v_1 - alimd * (Td_new - Tj2_new);
    } else if(m_trajectorytype == NoneTra || m_trajectorytype == NoneTri) {
        //decrease v_0
        Tv_new  = Tv;
        Td_new  = Td + Te;
        Tj2_new = Tj2 + 0.5 * Te;
        vlim    = 2 * Q / Td_new - v_1;
        alimd   = (v_1 - vlim) / (Td_new - Tj2_new);
        jmin    = alimd / Tj2_new;
        v_0     = vlim;
    }
    para = VectorXf::Zero(15);
    if(Tv_new < -1e-5) {
        int mannn = 0;
    }
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

void TrajectoryPlanning::TrajectoryCalculationD(double t, VectorXd para, double q_dq[3])
{
    double Ta, Tv, Td, Tj1, Tj2, p1 = 0.0, p2, v_0, v_1, vlim, amax, amin, alima, alimd, jmax, jmin;
    Ta     = para(0);
    Tv     = para(1);
    Td     = para(2);
    Tj1    = para(3);
    Tj2    = para(4);
    p2     = para(6);
    v_0    = para(7);
    v_1    = para(8);
    vlim   = para(9);
    amax   = para(10);
    amin   = para(11);
    alima  = para(12);
    alimd  = para(13);
    jmax   = para(14);
    jmin   = para(15);
    double T = 0.0;
    double q = 0.0, v = 0.0, a = 0.0, j = 0.0;
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

void TrajectoryPlanning::Movep(vector<PointInformation> waypoints, vector<PathInformation> &pathes)
{
    vector<PathInformation> temppath;
    CompoundTrajectory(waypoints, temppath);
//    for (int i = 0; i < temppath.size(); i++) {
//        qDebug()<<temppath[i].startpoint.point[0]<<temppath[i].startpoint.point[1]<<temppath[i].startpoint.point[2];
//        qDebug()<<temppath[i].endpoint.point[0]<<temppath[i].endpoint.point[1]<<temppath[i].endpoint.point[2];
//    }
    //    for (int i = 0; i < waypoints.size(); i++) {
    //        qDebug()<<waypoints[i].point[0]<<waypoints[i].point[1]<<waypoints[i].point[2];
    //    }
    if(temppath.size() < 1) {
        return;
    } else if (temppath.size() == 1) {
        PathInformation path = temppath.front();
    } else {
        vector<PathInformation>::const_iterator path1 = temppath.begin();
        vector<PathInformation>::const_iterator path2 = path1;
        path2++;
        PathInformation finalpath;
        while(path2 != temppath.end() - 1) {
            PathInformation currentpath = *path1;
            PathInformation nextpath = *path2;
            if(nextpath.pathType == Arc) {
                finalpath = nextpath;
                finalpath.pathType = Line2arc;
                finalpath.startpoint = currentpath.startpoint;
                finalpath.endpoint = nextpath.endpoint;
                finalpath.intermediatepoint = currentpath.endpoint;
                finalpath.displacement = currentpath.displacement;
                finalpath.arclength = nextpath.arclength;
                //finalpath.radius = nextpath.radius;
                //finalpath.theta = nextpath.theta;
                finalpath.maxVelocity = max(currentpath.maxVelocity, nextpath.maxVelocity);
                finalpath.maxAcceleration = max(currentpath.maxAcceleration, nextpath.maxAcceleration);
                finalpath.maxJerk = max(currentpath.maxJerk, nextpath.maxJerk);
                CalculatePathParameters(finalpath);
                pathes.push_back(finalpath);
                path2++;
                path1 = path2;
                if(path1 != temppath.end() - 1) {
                    path2++;
                } else {
                    PointInformation pathstart = finalpath.endpoint;
                    finalpath = *path1;
                    finalpath.startpoint = pathstart;
                    finalpath.pathType = Line;
                    CalculatePathParameters(finalpath);
                    pathes.push_back(finalpath);
                    break;
                }
            }
        }
    }
}

bool TrajectoryPlanning::TrajectoryInterpolation(vector<PathInformation> pathes, int peroid)
{
    if(pathes.size() < 1)
        return false;
    double t = 0.0;
    double time = 0.0;
    int interpolation_peroid_num = 0;
    double Q = 0.0;
    for (int i = 0; i < pathes.size(); i++) {
        PathInformation path = pathes[i];
        if(path.pathType == Line) {
            Q = path.displacement;
            time = path.constraints[0] + path.constraints[1] + path.constraints[2];
            if (time < 1e-10) return false;
            interpolation_peroid_num = round(time / (1.0 / peroid));
            double lambda[interpolation_peroid_num];
            vector<Vector3d> interpolation_points;
            Vector3d interpolation_point;
            double q_dq[3] = {0.0};
            for (int j = 0; j < interpolation_peroid_num; j++) {
                t = j * (1.0/peroid);
                //calculate trajectory
                TrajectoryCalculationD(t, path.constraints, q_dq);
                if(Q!=0)
                    lambda[j] = q_dq[0] / Q;
                interpolation_point = path.startpoint.point + lambda[j] * (path.endpoint.point - path.startpoint.point);
                interpolation_points.push_back(interpolation_point);
                QFile data("information.txt");
                if(!data.open(QIODevice::WriteOnly | QIODevice::Text | QIODevice::Append)) return false;
                QTextStream stream(&data);
                stream<<j<<" "<<lambda[j]<<" "<<q_dq[0]<<" "<<q_dq[1]<<" "<<interpolation_point[0]<<" "<<interpolation_point[1]<<" "<<interpolation_point[2]<<"\n";
                data.close();
            }
        }else if(path.pathType == Line2arc) {
            Q = path.displacement + path.arclength;
            time = path.constraints[0] + path.constraints[1] + path.constraints[2];
            if (time < 1e-10) return false;
            interpolation_peroid_num = round(time / (1.0 / peroid));
            double lambda[interpolation_peroid_num];
            vector<Vector3d> interpolation_points;
            Vector3d interpolation_point;
            double q_dq[3] = {0.0};
            for (int j = 0; j < interpolation_peroid_num; j++) {
                t = j*(1.0/peroid);
                //calculate trajectory
                TrajectoryCalculationD(t, path.constraints, q_dq);
                if(Q!=0)
                    lambda[j] = q_dq[0] / Q;
                if(lambda[j] <= (path.displacement/Q)) {
                    interpolation_point = path.startpoint.point + lambda[j] / (path.displacement/Q) * (path.intermediatepoint.point - path.startpoint.point);
                } else {
                    Vector3d yaxis = path.intermediatepoint.point - path.startpoint.point;
                    Vector3d xaxis = path.intermediatepoint.point - path.center;
                    //yaxis.normalized();
                    //xaxis.normalized();
                    double s = lambda[j] * Q - path.displacement;
                    qDebug()<<"s"<<s;
                    interpolation_point = path.center + path.radius * (xaxis/xaxis.norm()*cos(s / path.radius) + yaxis/yaxis.norm()*sin(s / path.radius));
                }
                interpolation_points.push_back(interpolation_point);
                QFile data("information.txt");
                if(!data.open(QIODevice::WriteOnly | QIODevice::Text | QIODevice::Append)) return false;
                QTextStream stream(&data);
                stream<<j<<" "<<lambda[j]<<" "<<q_dq[0]<<" "<<q_dq[1]<<" "<<interpolation_point[0]<<" "<<interpolation_point[1]<<" "<<interpolation_point[2]<<"\n";
                data.close();
            }
        }
    }
}

void TrajectoryPlanning::CompoundTrajectory(vector<PointInformation> &waypoints, vector<PathInformation> &path)
{
    if(waypoints.size() < 2)
        return;
    vector<PointInformation>::const_iterator point1 = waypoints.begin();
    vector<PointInformation>::const_iterator point2 = point1;
    point2++;
    vector<PointInformation>::const_iterator point3;
    PointInformation startpoint = *point1;
    while(point2 != waypoints.end() - 1) { //vector.end()指向最后一个元素的下一个位置，访问最后一个元素应该为vector.end()-1
        point3 = point2;
        point3++;
        PointInformation intermediatepoint; //中间点
        PointInformation endpoint; //终点
        intermediatepoint = *point2;
        endpoint = *point3;
        PathInformation arcpath;
        //计算圆弧段
        ArcSegmentLineToLine(startpoint, intermediatepoint, endpoint, arcpath);
        PathInformation linepath;
        //计算直线段
        SetLinePathSegment(linepath, startpoint, endpoint, arcpath);
        path.push_back(linepath);
        path.push_back(arcpath);
        startpoint = arcpath.endpoint; //圆弧段终点作为下一段规划的起点
        point2 = point3;
    }
    //point2是最后一个路点，此时不用计算圆弧段
    PointInformation endpoint;
    endpoint = *point2;
    PathInformation linepath;
    PathInformation arcpath;
    SetLinePathSegment(linepath, startpoint, endpoint, arcpath);
    path.push_back(linepath);
}

void TrajectoryPlanning::CalculatePathParameters(PathInformation &path)
{
    double time = 0.0;
    double Q = 0.0;
    if(path.pathType == Line2arc) {
        Q = path.displacement + path.arclength;
    } else if(path.pathType == Line) {
        Q = path.displacement;
    } else if(path.pathType == Circular) {
        Q = path.arclength;
    }
    TrajectoryTime(path.accelerationType, time, Q, path.startpoint.velocity, path.endpoint.velocity, path.maxVelocity, path.maxAcceleration, path.maxJerk, path.constraints);
    //修正轨迹约束参数
    CorrentionParameters(path.accelerationType, path.constraints, 100);
}

int TrajectoryPlanning::TrajectoryTime(AccelerationType &accelerationtype, double &time, double Q, double v_0, double v_1, double vmax, double amax, double jmax, VectorXd &para)
{
    //judge whether the minimun diaplacement is satisfied
    double Tj[2] = {0.0};
    int index = 0;
    Tj[0] = sqrt(fabs(v_1 - v_0)/jmax);
    Tj[1] = amax/jmax;
    if(Tj[0] > Tj[1]) {
        index = 1;
    } else {
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
    double delta = 1e-10;
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
    para << Ta, Tv, Td, Tj1, Tj2, 0, Q, v_0, v_1, vlim, amax, -amax, a_lima, a_limd, jmax, jmin;
    //classification
    if(Tv > 0.0) {
        if(Ta - 2*Tj1 > delta) {
            if(Td - 2*Tj2 > delta) {
                accelerationtype = TraZeroTra;
            } else if(Td > 0.0) {
                accelerationtype = TraZeroTri;
            } else {
                accelerationtype = TraZero;
            }
        } else if(Ta > 0.0) {
            if(Td - 2*Tj2 > delta) {
                accelerationtype = TriZeroTra;
            } else if(Td > 0.0) {
                accelerationtype = TriZeroTri;
            } else {
                accelerationtype = TriZero;
            }
        } else {
            if(Td - 2*Tj2 > delta) {
                accelerationtype = ZeroTra;
            } else if(Td > 0.0) {
                accelerationtype = ZeroTri;
            }else {
                accelerationtype = Zero;
            }
        }
    } else {
        if(Ta - 2*Tj1 > delta) {
            if(Td - 2*Tj2 > delta) {
                accelerationtype = TraTra;
            } else if(Td > 0.0) {
                accelerationtype = TraTri;
            } else {
                accelerationtype = TraNone;
            }
        }else if(Ta > 0.0) {
            if(Td - 2*Tj2 > delta) {
                accelerationtype = TriTra;
            } else if(Td > 0.0) {
                accelerationtype = TriTri;
            } else {
                accelerationtype = TriNone;
            }
        } else {
            if(Td - 2*Tj2 > delta) {
                accelerationtype = NoneTra;
            } else if(Td > 0.0) {
                accelerationtype = NoneTri;
            } else {
                return -2;
            }
        }
    }
    return 1;
}

void TrajectoryPlanning::CorrentionParameters(AccelerationType accelerationtype, VectorXd &para, int peroid)
{
    double Ta, Tv, Td, Tj1, Tj2, Q, v_0, v_1, vlim, amax, amin, alima, alimd, jmax, jmin;
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
    amin   = para(11);
    alima  = para(12);
    alimd  = para(13);
    jmax   = para(14);
    jmin   = para(15);
    double Te       = ceil((Ta + Tv + Td) / (1.0 / peroid)) * (1.0 / peroid) - (Ta + Tv + Td);
    double Tj1_new  = Tj1, Tj2_new = Tj2, Ta_new = Ta, Td_new = Td, Tv_new = Tv;
    double delta_T  = 0.0;
    double delta    = 1e-10;
    double Tm       = Te + Ta + Tv + Td;
    if(accelerationtype == TraZeroTra || accelerationtype == TriZeroTri) {
        if((v_1 - v_0) > delta) { //extend the acceleration time
            delta_T = Te * (vlim + v_0) / (vlim - v_0);
            Tv_new  = Tv - delta_T;
            if(Tv_new < -delta) {
                //decrease vlim
                Tv_new  = Tv;
                Ta_new  = Ta + Te;
                Tj1_new = Tj1 + 0.5 * Te;
                vlim    = (2 * Q - Ta_new * v_0 - Td_new * v_1) / (2 * Tm - Ta_new - Td_new);
                alima   = (vlim - v_0) / (Ta_new - Tj1_new);
                alimd   = (v_1 - vlim) / (Td_new - Tj2_new);
                jmax    = alima / Tj1_new;
                jmin    = alimd / Tj2_new;
            } else {
                Ta_new  = Ta + Te + delta_T;
                Tj1_new = Tj1 + 0.5 * (Te + delta_T);
                jmax    = (vlim - v_0) / ((Ta_new - Tj1_new) * Tj1_new);
                alima   = Tj1_new * jmax;
                vlim    = alima * (Ta_new - Tj1_new) + v_0;
            }
        } else if((v_1 - v_0) < -delta) { //extend the deceleration time
            delta_T = Te * (vlim + v_1) / (vlim - v_1);
            Tv_new  = Tv - delta_T;
            if(Tv_new < -delta) {
                //decrease vlim
                Tv_new  = Tv;
                Td_new  = Td + Te;
                Tj2_new = Tj2 + 0.5 * Te;
                vlim    = (2 * Q - Ta_new * v_0 - Td_new * v_1) / (2 * Tm - Ta_new - Td_new);
                alima   = (vlim - v_0) / (Ta_new - Tj1_new);
                alimd   = (v_1 - vlim) / (Td_new - Tj2_new);
                jmax    = alima / Tj1_new;
                jmin    = alimd / Tj2_new;
            } else {
                Td_new  = Td + Te + delta_T;
                Tj2_new = Tj2 + 0.5 * (Te + delta_T);
                jmin    = (v_1 - vlim) / ((Td_new - Tj2_new) * Tj2_new);
                alimd   = Tj2_new * jmin;
                vlim    = v_1 - alimd*(Td_new - Tj2_new);
            }
        } else { //extend the acceleration and deceleration time
            delta_T = Te * (v_0 + 2 * vlim + v_1) / (2 * vlim - v_0 - v_1);
            Tv_new  = Tv - delta_T;
            if(Tv_new < -delta) {
                //decrease vlim
                Tv_new  = Tv;
                Ta_new  = Ta + 0.5 * Te;
                Tj1_new = Tj1 + 0.25 * Te;
                Td_new  = Td + 0.5 * Te;
                Tj2_new = Tj2 + 0.25 * Te;
                vlim    = (2 * Q - Ta_new * v_0 - Td_new * v_1) / (2 * Tm - Ta_new - Td_new);
                alima   = (vlim - v_0) / (Ta_new - Tj1_new);
                alimd   = (v_1 - vlim) / (Td_new - Tj2_new);
                jmax    = alima / Tj1_new;
                jmin    = alimd / Tj2_new;
            } else {
                Ta_new  = Ta + 0.5 * (Te + delta_T);
                Td_new  = Td + 0.5 * (Te + delta_T);
                Tj1_new = Tj1 + (Te + delta_T) * 0.25;
                Tj2_new = Tj2 + (Te + delta_T) * 0.25;
                jmax    = (vlim - v_0) / ((Ta_new - Tj1_new) * Tj1_new);
                jmin    = (v_1 - vlim) / ((Td_new - Tj2_new) * Tj2_new);
                alima   = Tj1_new * jmax;
                alimd   = Tj2_new * jmin;
                vlim    = alima * (Ta_new - Tj1_new) + v_0;
            }
        }
    } else if(accelerationtype == TriZeroTra || accelerationtype == ZeroTra
               || accelerationtype == ZeroTri) { //extend the deceleration time
        delta_T = Te * (vlim + v_1) / (vlim - v_1);
        Tv_new  = Tv - delta_T;
        if(Tv_new < -delta) {
            if(accelerationtype == TriZeroTra) {
                //decrease vlim
                Tv_new  = Tv;
                Td_new  = Td + Te;
                Tj2_new = Tj2 + 0.5 * Te;
                vlim    = (2 * Q - Ta_new * v_0 - Td_new * v_1) / (2 * Tm - Ta_new - Td_new);
                alima   = (vlim - v_0) / (Ta_new - Tj1_new);
                alimd   = (v_1 - vlim) / (Td_new - Tj2_new);
                jmax    = alima / Tj1_new;
                jmin    = alimd / Tj2_new;
            } else {
                //decrease vlim and v_0
                Tv_new  = Tv;
                Td_new  = Td + Te;
                Tj2_new = Tj2 + 0.5 * Te;
                vlim    = (2 * Q - Td_new * v_1) / (2 * Tv_new + Td_new);
                alimd   = (v_1 - vlim) / (Td_new - Tj2_new);
                jmin    = alimd / Tj2_new;
                v_0 = vlim;
            }
        } else {
            Td_new  = Td + Te + delta_T;
            Tj2_new = Tj2 + 0.5 * (Te + delta_T);
            jmin    = (v_1 - vlim) / ((Td_new - Tj2_new) * Tj2_new);
            alimd   = Tj2_new * jmin;
            vlim    = v_1 - alimd*(Td_new - Tj2_new);
        }
    } else if(accelerationtype == TraZeroTri || accelerationtype == TraZero
               || accelerationtype == TriZero) { //extend the acceleration time
        delta_T = Te * (vlim + v_0) / (vlim - v_0);
        Tv_new  = Tv - delta_T;
        if(Tv_new < -delta) {
            if(accelerationtype == TraZeroTri) {
                //decrease vlim
                Tv_new  = Tv;
                Ta_new  = Ta + Te;
                Tj1_new = Tj1 + 0.5 * Te;
                vlim    = (2 * Q - Ta_new * v_0 - Td_new * v_1) / (2 * Tm - Ta_new - Td_new);
                alima   = (vlim - v_0) / (Ta_new - Tj1_new);
                alimd   = (v_1 - vlim) / (Td_new - Tj2_new);
                jmax    = alima / Tj1_new;
                jmin    = alimd / Tj2_new;
            } else {
                //decrease vlim and v_1
                Tv_new  = Tv;
                Ta_new  = Ta + Te;
                Tj1_new = Tj1 + 0.5 * Te;
                vlim    = (2 * Q - Ta_new * v_0) / (2 * Tv_new + Ta_new);
                alima   = (vlim - v_0) / (Ta_new - Tj1_new);
                alimd   = (v_1 - vlim) / (Td_new - Tj2_new);
                jmax    = alima / Tj1_new;
                jmin    = alimd / Tj2_new;
                v_1 = vlim;
            }
        } else {
            Ta_new  = Ta + Te + delta_T;
            Tj1_new = Tj1 + 0.5 * (Te + delta_T);
            jmax    = (vlim - v_0) / ((Ta_new - Tj1_new) * Tj1_new);
            alima   = Tj1_new * jmax;
            vlim    = alima * (Ta_new - Tj1_new) + v_0;
        }
    } else if(accelerationtype == Zero){
        //decrease v_0 and v_1
        Td_new = cbrt(8*(Tv * vlim - (Tv + Te)*vlim) / jmin);
        Tv_new = Tv - Td_new;
        if(Tv_new < -delta) {
            Tv_new  = Tv + Te;
            vlim    = Q / Tv_new;
            v_0     = vlim;
            v_1     = vlim;
            Td_new  = Td;
        } else {
            Tj2_new = 0.5 * Td_new;
            alimd   = 0.5 * jmin * Td_new;
            v_1     = vlim + 0.25 * jmin * pow(Td_new,2);
        }
    } else if(accelerationtype == TraTra || accelerationtype == TriTri) {
        delta_T = 0.0;
        Tv_new  = Tv - delta_T;
        Ta_new  = Ta + 0.5 * (Te + delta_T);
        Td_new  = Td + 0.5 * (Te + delta_T);
        Tj1_new = Tj1 + (Te + delta_T) * 0.25;
        Tj2_new = Tj2 + (Te + delta_T) * 0.25;
        jmax    = (vlim - v_0) / ((Ta_new - Tj1_new) * Tj1_new);
        jmin    = (v_1 - vlim) / ((Td_new - Tj2_new) * Tj2_new);
        alima   = Tj1_new * jmax;
        alimd   = Tj2_new * jmin;
        vlim    = alima * (Ta_new - Tj1_new) + v_0;
    } else if(accelerationtype == TraTri) {
        Tv_new  = Tv - delta_T;
        Ta_new  = Ta + Te + delta_T;
        Tj1_new = Tj1 + 0.5 * (Te + delta_T);
        jmax    = (vlim - v_0) / ((Ta_new - Tj1_new) * Tj1_new);
        alima   = Tj1_new * jmax;
        vlim    = alima * (Ta_new - Tj1_new) + v_0;
    } else if(accelerationtype == TriNone || accelerationtype == TraNone) {
        //decrease v_1
        Tv_new  = Tv;
        Ta_new  = Ta + Te;
        Tj1_new = Tj1 + 0.5 * Te;
        vlim    = 2 * Q / Ta_new - v_0;
        alima   = (vlim - v_0) / (Ta_new - Tj1_new);
        jmax    = alima / Tj1_new;
        v_1     = vlim;
    } else if(accelerationtype == TriTra) {
        delta_T = 0.0;
        Tv_new  = Tv - delta_T;
        Td_new  = Td + Te + delta_T;
        Tj2_new = Tj2 + 0.5 * (Te + delta_T);
        jmin    = (v_1 - vlim) / ((Td_new - Tj2_new) * Tj2_new);
        alimd   = Tj2_new * jmin;
        vlim    = v_1 - alimd * (Td_new - Tj2_new);
    } else if(accelerationtype == NoneTra || accelerationtype == NoneTri) {
        //decrease v_0
        Tv_new  = Tv;
        Td_new  = Td + Te;
        Tj2_new = Tj2 + 0.5 * Te;
        vlim    = 2 * Q / Td_new - v_1;
        alimd   = (v_1 - vlim) / (Td_new - Tj2_new);
        jmin    = alimd / Tj2_new;
        v_0     = vlim;
    }
    para = VectorXd::Zero(16);
    if(Tv_new < -1e-5) {
        int mannn = 0;
    }
    para << Ta_new, Tv_new, Td_new, Tj1_new, Tj2_new, 0, Q, v_0, v_1, vlim, amax, amin, alima, alimd, jmax, jmin;
}

void TrajectoryPlanning::TrajectoryCalculation(double t, VectorXd para, double q_dq[3])
{
    double Ta, Tv, Td, Tj1, Tj2, p1, p2, v_0, v_1, vlim, amax, alima, alimd, jmax, jmin;
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
    double T = 0.0;
    double q = 0.0, v = 0.0, a = 0.0, j = 0.0;
    double delta = 1e-10;
    T = Ta + Tv + Td;
    if (t >= 0.0 &&  t < Tj1) {
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

void TrajectoryPlanning::ArcSegmentLineToLine(PointInformation &startpoint, PointInformation &intermediatepoint, PointInformation &endpoint,
                                              PathInformation &arcpath)
{
    Vector3d spoint = startpoint.point;
    Vector3d epoint = endpoint.point;
    Vector3d ipoint = intermediatepoint.point;
    Vector3d intermediate2start = 0.5 * (spoint - ipoint);
    Vector3d intermediate2end = 0.5 * (epoint - ipoint);
    //两直线夹角
    double theta = acos(intermediate2start.dot(intermediate2end) / (intermediate2start.norm() * intermediate2end.norm()));
    //    double theta = acos((intermediate2start[0]*1000.0*intermediate2end[0]*1000.0 + intermediate2start[1]*1000.0*intermediate2end[1]*1000.0 + intermediate2start[2]*1000.0*intermediate2end[2]*1000.0) / (intermediate2start.norm()*1000.0 * intermediate2end.norm()*1000.0));
    if(theta / M_PI * 180.0 < 1e-10)
        return;
    double l = intermediatepoint.radius / tan(theta * 0.5);
    //判断交融点是否超过直线中点
    if((intermediate2start.norm() - l) < 1e-6 || (intermediate2end.norm() - l)  < 1e-6) {
        if((intermediate2start.norm() - intermediate2end.norm()) < 1e-6) {
            if(intermediate2start.norm() > 1e-6) {
                intermediatepoint.radius = tan(theta * 0.5) * intermediate2start.norm();
            } else {
                intermediatepoint.radius = 0;
                return;
            }
        } else {
            if(intermediate2end.norm() > 1e-6) {
                intermediatepoint.radius = tan(theta * 0.5) * intermediate2end.norm();
            } else {
                intermediatepoint.radius = 0;
                return;
            }
        }
    }
    //更新交融点在直线上的位置
    l = intermediatepoint.radius / tan(theta * 0.5);
    //交融圆弧起点
    Vector3d arcstartpoint = l * intermediate2start.normalized() + intermediatepoint.point;
    //交融圆弧终点
    Vector3d arcendpoint = l * intermediate2end.normalized() + intermediatepoint.point;
    //计算交融圆弧圆心
    Vector3d intermediat2arcstart = arcstartpoint - intermediatepoint.point;
    Vector3d intermediat2arcend = arcendpoint - intermediatepoint.point;
    Vector3d intermediat2middle = 0.5 * (intermediat2arcstart + intermediat2arcend);
    Vector3d intermediat2arccenter = intermediatepoint.radius / sin(theta * 0.5) / intermediat2middle.norm() * intermediat2middle;
    Vector3d arccenter = intermediatepoint.point + intermediat2arccenter;
    //保存圆弧段信息
    SetArcPathSegment(arcpath, intermediatepoint, arcstartpoint, arcendpoint, arccenter, intermediatepoint.radius, M_PI - theta);
}

/**
 * @brief TrajectoryPlanning::SetLinePathSegment
 * @param path
 * @param startpoint
 * @param endpoint      中间点
 * @param point
 */
void TrajectoryPlanning::SetLinePathSegment(PathInformation &path, PointInformation startpoint, PointInformation endpoint,
                                            PathInformation arcpath)
{
    path.startpoint.pathType = Line;
    path.endpoint.pathType = Line;
    path.pathType = Line;
    path.startpoint = startpoint;
    if(arcpath.pathType == Arc) {
        path.endpoint = arcpath.startpoint;
    } else {
        path.endpoint = endpoint;
    }
    Vector3d start2end = path.endpoint.point - path.startpoint.point;
    path.displacement = start2end.norm();
    path.maxVelocity = max(path.startpoint.maxVelocity, path.endpoint.maxVelocity);
    path.maxAcceleration = max(path.startpoint.maxAcceleration, path.endpoint.maxAcceleration);
    path.maxJerk = max(path.startpoint.maxJerk, path.endpoint.maxJerk);
}

void TrajectoryPlanning::SetCricularPathSegment(PathInformation &path, PointInformation startpoint, PointInformation endpoint, PointInformation middlepoint,
                                                Vector3d center, double r, double theta)
{
    path.pathType = Circular;
    path.radius = r;
    path.center = center;
    path.theta = theta;
    path.startpoint.point = startpoint.point;
    path.endpoint.point = endpoint.point;
    path.startpoint.pose = startpoint.pose;
    path.endpoint.pose = endpoint.pose;
    path.maxVelocity = max(startpoint.maxVelocity, middlepoint.maxVelocity);
    path.maxVelocity = max(path.maxVelocity, endpoint.maxVelocity);
    path.maxAcceleration = max(startpoint.maxAcceleration, middlepoint.maxAcceleration);
    path.maxAcceleration = max(path.maxAcceleration, endpoint.maxAcceleration);
    path.maxJerk = max(startpoint.maxJerk, middlepoint.maxJerk);
    path.maxJerk = max(path.maxJerk, endpoint.maxJerk);
    path.startpoint.velocity = startpoint.velocity;
    path.endpoint.velocity = endpoint.velocity;
}

void TrajectoryPlanning::SetArcPathSegment(
    PathInformation &path, PointInformation intermediatepoint,
    Vector3d arcstart, Vector3d arcend,
    Vector3d arccenter, double radius, double theta)
{
    path.startpoint             = intermediatepoint;
    path.endpoint               = intermediatepoint;
    path.waypoint               = intermediatepoint.point;
    path.startpoint.pathType    = Arc;
    path.endpoint.pathType      = Arc;
    path.pathType               = Arc;
    path.radius                 = radius;
    path.center                 = arccenter;
    path.theta                  = theta;
    path.arclength              = path.theta * path.radius;
    path.startpoint.point       = arcstart;
    path.endpoint.point         = arcend;
    path.startpoint.pose        = intermediatepoint.pose;
    path.endpoint.pose          = intermediatepoint.pose;
    path.maxVelocity            = intermediatepoint.maxVelocity;
    path.maxAcceleration        = intermediatepoint.maxAcceleration;
    path.maxJerk                = intermediatepoint.maxJerk;
    path.startpoint.velocity    = intermediatepoint.velocity;
    path.endpoint.velocity      = intermediatepoint.velocity;
}

/**
 * @brief TrajectoryPlanning::ArcParameterCalculate 计算圆弧的参数，包括圆心、半径、旋转轴、旋转角
 * @param p1
 * @param p2
 * @param p3
 * @param center
 * @param radius
 * @param normal
 * @param theta
 */
void TrajectoryPlanning::GetArcParameter(
    Vector3d p1, Vector3d p2, Vector3d p3, Vector3d &center,
    double &radius, Vector3d &normal, double &theta)
{
    // 三个点确定一个平面方程
    double k_11 = (p1(1) - p3(1))*(p2(2) - p3(2)) - (p2(1) - p3(1))*(p1(2) - p3(2));
    double k_12 = (p2(0) - p3(0))*(p1(2) - p3(2)) - (p1(0) - p3(0))*(p2(2) - p3(2));
    double k_13 = (p1(0) - p3(0))*(p2(1) - p3(1)) - (p2(0) - p3(0))*(p1(1) - p3(1));
    double k_14 = -(k_11*p3(0) + k_12*p3(1) + k_13*p3(2));
    // 过p1 p2的中点并且和 p1p2 垂直的平面方程
    double k_21 = p2(0) - p1(0);
    double k_22 = p2(1) - p1(1);
    double k_23 = p2(2) - p1(2);
    double k_24 = -((pow(p2(0),2) - pow(p1(0),2)) + (pow(p2(1),2) - pow(p1(1),2)) + (pow(p2(2), 2) - pow(p1(2),2))) / 2;
    // 过p2 p3的中点并且和 p2p3 垂直的平面方程
    double k_31 = p3(0) - p2(0);
    double k_32 = p3(1) - p2(1);
    double k_33 = p3(2) - p2(2);
    double k_34 = -((pow(p3(0),2) - pow(p2(0),2)) + (pow(p3(1),2) - pow(p2(1),2)) + (pow(p3(2),2) - pow(p2(2),2))) / 2;
    // 圆心肯定在这三个平面上面。利用点法式方程求得圆心
    Matrix3d k_123 = Matrix3d::Zero();
    k_123 << k_11, k_12, k_13,
             k_21, k_22, k_23,
             k_31, k_32, k_33;
    Vector3d k_4 = Vector3d::Zero();
    k_4<< -k_14, -k_24, -k_34;
    // 计算圆心
    center = k_123.inverse() * k_4;
    // 计算半径
    radius = sqrt(pow((center(0) - p1(0)),2) + pow((center(1) - p1(1)),2) + pow((center(2) - p1(2)),2));
    // 计算圆心角及圆弧法向量
    Vector3d center2p1 = p1 - center;
    Vector3d center2p2 = p2 - center;
    Vector3d center2p3 = p3 - center;
    Vector3d normal1 = center2p1.cross(center2p2);
    Vector3d normal2 = center2p1.cross(center2p3);
    double theta1 = acos(center2p1.dot(center2p2) / (center2p1.norm() * center2p2.norm()));
    double theta2 = acos(center2p1.dot(center2p3) / (center2p1.norm() * center2p3.norm()));
    if(normal1.dot(normal2) > 0) {  // 法向量同向
        if(theta1 > theta2) {
            theta = 2*M_PI - theta2;
            normal = - normal1.normalized();
        } else {
            theta = theta2;
            normal = normal1.normalized();
        }
    } else {    // 法向量相反
        theta = 2*M_PI - theta2;
        normal = normal1.normalized();
    }
}

void TrajectoryPlanning::GetArcParameter(PathInformation &arc)
{
    Vector3d p1 = arc.startpoint.point;
    Vector3d p2 = arc.intermediatepoint.point;
    Vector3d p3 = arc.endpoint.point;
    // 三个点确定一个平面方程
    double k_11 = (p1(1) - p3(1))*(p2(2) - p3(2)) - (p2(1) - p3(1))*(p1(2) - p3(2));
    double k_12 = (p2(0) - p3(0))*(p1(2) - p3(2)) - (p1(0) - p3(0))*(p2(2) - p3(2));
    double k_13 = (p1(0) - p3(0))*(p2(1) - p3(1)) - (p2(0) - p3(0))*(p1(1) - p3(1));
    double k_14 = -(k_11*p3(0) + k_12*p3(1) + k_13*p3(2));
    // 过p1 p2的中点并且和 p1p2 垂直的平面方程
    double k_21 = p2(0) - p1(0);
    double k_22 = p2(1) - p1(1);
    double k_23 = p2(2) - p1(2);
    double k_24 = -((pow(p2(0),2) - pow(p1(0),2)) + (pow(p2(1),2) - pow(p1(1),2)) + (pow(p2(2), 2) - pow(p1(2),2))) / 2.0;
    // 过p2 p3的中点并且和 p2p3 垂直的平面方程
    double k_31 = p3(0) - p2(0);
    double k_32 = p3(1) - p2(1);
    double k_33 = p3(2) - p2(2);
    double k_34 = -((pow(p3(0),2) - pow(p2(0),2)) + (pow(p3(1),2) - pow(p2(1),2)) + (pow(p3(2),2) - pow(p2(2),2))) / 2.0;
    // 圆心肯定在这三个平面上面。利用点法式方程求得圆心
    Matrix3d k_123 = Matrix3d::Zero();
    k_123 << k_11, k_12, k_13,
        k_21, k_22, k_23,
        k_31, k_32, k_33;
    Vector3d k_4 = Vector3d::Zero();
    k_4<< -k_14, -k_24, -k_34;
    // 计算圆心
    arc.center = k_123.inverse() * k_4;
    // 计算半径
    arc.radius = sqrt(pow((center(0) - p1(0)),2) + pow((center(1) - p1(1)),2) + pow((center(2) - p1(2)),2));
    // 计算圆心角及圆弧法向量
    Vector3d center2p1 = p1 - arc.center;
    Vector3d center2p2 = p2 - arc.center;
    Vector3d center2p3 = p3 - arc.center;
    Vector3d normal1 = center2p1.cross(center2p2);
    Vector3d normal2 = center2p1.cross(center2p3);
    double theta1 = acos(center2p1.dot(center2p2) / (center2p1.norm() * center2p2.norm()));
    double theta2 = acos(center2p1.dot(center2p3) / (center2p1.norm() * center2p3.norm()));
    if(normal1.dot(normal2) > 0) {  // 法向量同向
        if(theta1 > theta2) {
            arc.theta = 2*M_PI - theta2;
            arc.normal = - normal1.normalized();
        } else {
            arc.theta = theta2;
            arc.normal = normal1.normalized();
        }
    } else {    // 法向量相反
        arc.theta = 2*M_PI - theta2;
        arc.normal = normal1.normalized();
    }
}
/**
 * @brief TrajectoryPlanning::LineIntersectCricular 线段与以端点为圆心，r为半径的圆的交点
 * @param startpoint                                线段起点
 * @param endpoint                                  线段终点
 * @param radius                                    半径
 * @return                                          直线与圆的交点
 */
Vector3d TrajectoryPlanning::LineIntersectCricular(
    Vector3d startpoint, Vector3d endpoint, double radius, ConnectionType type)
{
    Vector3d start2end = endpoint - startpoint;
    if(type == ConnectwithEnd) {
        return endpoint - radius * start2end.normalized();
    } else if(type == ConnectwithStart){
        return startpoint + radius * start2end.normalized();
    }
}

Vector3d TrajectoryPlanning::LineIntersectCricular(PathInformation line, ConnectionType type)
{
    Vector3d start2end = line.endpoint.point - line.startpoint.point;
    Vector3d controlPoint;
    if (type == ConnectwithEnd) {
        controlPoint = line.endpoint.point - line.endpoint.r * start2end.normalized();
        // 若交点越过中点，则设置中点为交点
        if ((line.endpoint.point - controlPoint).norm() > 0.5 * start2end.norm()) {
            controlPoint = line.startpoint.point + 0.5 * start2end;
        }
    } else if (type == ConnectwithStart) {
        controlPoint = line.startpoint.point + line.startpoint.r * start2end.normalized();
        // 若交点越过中点，则设置中点为交点
        if ((line.startpoint.point - controlPoint).norm() > 0.5 * start2end.norm()) {
            controlPoint = line.startpoint.point + 0.5 * start2end;
        }
    }
    return controlPoint;
}

/**
 * @brief TrajectoryPlanning::CricularIntersectCricular 圆弧与以圆弧端点为圆心，r为半径的交融圆的交点
 * @param center                                        圆弧圆心
 * @param circularradius                                圆弧半径
 * @param normal                                        圆弧法线（圆弧起点旋转到圆弧终点）
 * @param intersectpoint                                圆弧端点
 * @param radius                                        交融圆半径
 * @return                                              圆弧与交融圆的交点
 */
Vector3d TrajectoryPlanning::CricularIntersectCricular(
    Vector3d center, double circularradius, Vector3d normal,
    Vector3d intersectpoint, double radius, ConnectionType type)
{
    double theta = 2 * asin(0.5 * radius / circularradius); //计算圆弧与交融圆的交点到圆弧连接点对应的圆心角
    Vector3d center2intersect = intersectpoint - center;    //圆心指向圆弧连接点的矢量
    Eigen::AngleAxisd rotationvectorend(-theta, normal);       //绕圆弧法线旋转-theta角，得到圆弧与圆的交点的矢量
    Eigen::AngleAxisd rotationvectorstart(theta, normal);      //绕圆弧法线旋转theta角，得到圆弧与圆的交点的矢量
    Eigen::Matrix3d rotationmatrix;
    if(type == ConnectwithEnd) {
        rotationmatrix = rotationvectorend.matrix();
    } else if(type == ConnectwithStart){
        rotationmatrix = rotationvectorstart.matrix();
    }
//    rotationmatrix = rotationvector.toRotationMatrix();
    return rotationmatrix * center2intersect + center;      //圆弧与圆的交点
}

Vector3d TrajectoryPlanning::CricularIntersectCricular(PathInformation arc, ConnectionType type)
{
    double theta = 0.0;
    Vector3d center2intersect;
    Eigen::Matrix3d rotationmatrix;
    if(type == ConnectwithEnd) {
        theta = 2 * asin(0.5 * arc.endpoint.r / arc.radius);    //计算圆弧与交融圆的交点到圆弧连接点对应的圆心角
        if (theta > 0.5 * arc.theta) {
            theta = 0.5 * arc.theta;
        }
        Eigen::AngleAxisd rotationvectorend(-theta, arc.normal);    //绕圆弧法线旋转-theta角，得到圆弧与圆的交点的矢量
        rotationmatrix = rotationvectorend.matrix();
        center2intersect = arc.endpoint.point - arc.center; //圆心指向圆弧连接点的矢量
    } else if(type == ConnectwithStart){
        theta = 2 * asin(0.5 * arc.startpoint.r / arc.radius);  //计算圆弧与交融圆的交点到圆弧连接点对应的圆心角
        if (theta > 0.5 * arc.theta) {
            theta = 0.5 * arc.theta;
        }
        Eigen::AngleAxisd rotationvectorstart(theta, arc.normal);   //绕圆弧法线旋转theta角，得到圆弧与圆的交点的矢量
        rotationmatrix = rotationvectorstart.matrix();
        center2intersect = arc.startpoint.point - arc.center;   //圆心指向圆弧连接点的矢量
    }
    return rotationmatrix * center2intersect + arc.center;  //圆弧与圆的交点
}

/**
 * @brief TrajectoryPlanning::LineTangent   线段的单位切向量
 * @param startpoint                        起点
 * @param endpoint                          终点
 * @param type                              连接点类型
 * @return                                  单位切向量（指向远离线段的方向）
 */
Vector3d TrajectoryPlanning::GetLineTangent(
    Vector3d startpoint, Vector3d endpoint, ConnectionType type)
{
    Vector3d start2end = endpoint - startpoint;
    if(type == ConnectwithEnd) {
        return start2end.normalized();
    } else if(type == ConnectwithStart){
        return -start2end.normalized();
    }
}

Vector3d TrajectoryPlanning::GetLineTangent(PathInformation line, ConnectionType type)
{
    Vector3d start2end = line.endpoint.point - line.startpoint.point;
    if(type == ConnectwithEnd) {
        return start2end.normalized();
    } else if(type == ConnectwithStart){
        return - start2end.normalized();
    }
}

/**
 * @brief TrajectoryPlanning::GetArcTangent 圆弧端点的单位切向量（指向远离圆弧的方向）
 * @param startpoint
 * @param endpoint
 * @param center
 * @param radius
 * @param normal
 * @param type
 * @return
 */
Vector3d TrajectoryPlanning::GetArcTangent(
    Vector3d startpoint, Vector3d endpoint, Vector3d center,
    Vector3d normal, ConnectionType type)
{
    Vector3d tangent;
    if(type == ConnectwithEnd) {
        Vector3d center2end = endpoint - center;
        tangent = normal.cross(center2end);
    } else if(type == ConnectwithStart){
        Vector3d center2start = startpoint - center;
        tangent = center2start.cross(normal);
    }
    return tangent.normalized();
}

Vector3d TrajectoryPlanning::GetArcTangent(PathInformation arc, ConnectionType type)
{
    Vector3d tangent;
    if(type == ConnectwithEnd) {
        Vector3d center2end = arc.endpoint.point - arc.center;
        tangent = arc.normal.cross(center2end);
    } else if(type == ConnectwithStart){
        Vector3d center2start = arc.startpoint.point - arc.center;
        tangent = center2start.cross(normal);
    }
    return tangent.normalized();
}

/**
 * @brief TrajectoryPlanning::GetControlPoints  计算直线、圆弧两种基本路径不同组合下的曲线控制点
 * @param path1
 * @param path2
 * @return
 */
vector<PointInformation> TrajectoryPlanning::GetControlPoints(PathInformation path1, PathInformation path2)
{
    vector<PointInformation> controlPoints;
    vector<Vector3d> controlPoint;
    controlPoint.resize(4);
    Vector3d tangent1, tangent2;    // 切向量（远离路径方向）
    // 计算路径1终点与交融圆的交点及交点处的切向量
    if (path1.pathType == Line) {
        controlPoint[0] = LineIntersectCricular(path1, ConnectwithEnd);
        tangent1 = GetLineTangent(path1.startpoint.point, controlPoint[0], ConnectwithEnd);
    } else if (path1.pathType == Arc) {
        controlPoint[0] = CricularIntersectCricular(path1, ConnectwithEnd);
        tangent1 = GetArcTangent(path1.startpoint.point, controlPoint[0], path1.center, path1.normal, ConnectwithEnd);
    }
    // 计算路径2起点与交融圆的交点及交点处的切向量
    if (path2.pathType == Line) {
        controlPoint[3] = LineIntersectCricular(path2, ConnectwithStart);
        tangent2 = GetLineTangent(controlPoint[3], path2.endpoint.point, ConnectwithStart);
    } else if (path2.pathType == Arc) {
        controlPoint[3] = CricularIntersectCricular(path2, ConnectwithStart);
        tangent2 = GetArcTangent(controlPoint[3], path2.endpoint.point, path2.center, path2.normal, ConnectwithStart);
    }
    // 以point1和point4为圆心，以0.5*min{r,|p1p4|}为半径，计算与路径交点处切向量方向的交点point2和point3
    Vector3d tempVector = controlPoint[3] - controlPoint[0];
    double rad = 0.5 * min(tempVector.norm(), path1.endpoint.r);
    controlPoint[1] = rad * tangent1;
    controlPoint[2] = rad * tangent2;
    // p1，p2继承path1.endpoint信息
    PointInformation cpoint;
    cpoint = path1.endpoint;
    cpoint.point = controlPoint[0];
    controlPoints.push_back(cpoint);
    cpoint.point = controlPoint[1];
    controlPoints.push_back(cpoint);
    // p3，p4继承path2.startpoint信息
    cpoint = path2.startpoint;
    cpoint.point = controlPoint[2];
    controlPoints.push_back(cpoint);
    cpoint.point = controlPoint[3];
    controlPoints.push_back(cpoint);
    return controlPoints;
}

/**
 * @brief TrajectoryPlanning::GetCurveParameter 曲线参数计算（弧长等）
 * @param path
 */
void TrajectoryPlanning::GetCurveParameter(PathInformation &path)
{
    vector<Vector3d> ctrlPoint;
}

/**
 * @brief TrajectoryPlanning::GetPointInformation   存放movep路点信息
 * @param point
 * @param pose
 * @param v
 * @param a
 * @param vm
 * @param am
 * @param jm
 * @param r
 * @param type
 * @return
 */
PointInformation TrajectoryPlanning::GetPointInformation(
    Vector3d point, Vector4d pose,
    double v, double a, double vm, double am, double jm,
    double r, PointType type)
{
    PointInformation p;
    p.point = point;
    p.pose = pose;
    p.velocity = v;
    p.acceleration = a;
    p.jerk = jm;
    p.maxVelocity = vm;
    p.maxAcceleration = am;
    p.maxJerk = jm;
    p.r = r;
    p.pointType = type;
    return p;
}

/**
 * @brief TrajectoryPlanning::GetPointInformation 将路点信息转换成没有交融的路径信息
 * @param waypoints
 * @param pathes
 */
int TrajectoryPlanning::GetPathInformation(vector<PointInformation> waypoints, vector<PathInformation> pathes)
{
    if (waypoints.size() < 1) return -1;
    vector<PointInformation>::const_iterator point1 = waypoints.begin();
    vector<PointInformation>::const_iterator point2 = point1;
    point2++;
    vector<PointInformation>::const_iterator point3;
    PointInformation startPoint = *point1;
    PathInformation path;
    while(point2 != waypoints.end() - 1) { //vector.end()指向最后一个元素的下一个位置，访问最后一个元素应该为vector.end()-1
        point3 = point2;
        point3 ++;
        PointInformation nextPoint = *point2;
        if (startPoint.pointType == Line) { // 直线类型点
            path = SetLinePath(startPoint, nextPoint);
            pathes.push_back(path);
            startPoint = *point2;
            point2 = point3;
        } else if (startPoint.pointType == Arc) {   // 圆弧类型点
            if (nextPoint.pointType == Line) {  // 下一个点是直线点
                path = SetLinePath(startPoint, nextPoint);
                pathes.push_back(path);
                startPoint = *point2;
                point2 = point3;
            } else if (nextPoint.pointType == Arc) {    // 下一个点是圆弧点
                PointInformation lastPoint = *point3;
                path = SetArcPath(startPoint, nextPoint, lastPoint);
                pathes.push_back(path);
                startPoint = *point3;
                if (point3 != waypoints.end() - 1) {    // point3不是最后一个点
                    if (point3 != waypoints.end() - 2) {    // point3不是倒数第二个点
                        startPoint = *point3;
                        point3 ++;
                        point2 = point3;
                    } else {    // point3是倒数第二个点
                        startPoint = *point3;
                        point3 ++;
                        nextPoint = *point3;
                        path = SetLinePath(startPoint, nextPoint);
                        pathes.push_back(path);
                        return 0;
                    }
                } else {    // point3是最后一个点
                    return 0;
                }
            }
        }
    }
    // point2是最后一个点
    PointInformation endPoint = *point2;
    path = SetLinePath(startPoint, endPoint);
    pathes.push_back(path);
    return 0;
}

void TrajectoryPlanning::MoveP(vector<PathInformation> path)
{
    // 交融后的路径放入m_movePath
    SetPathSegment(path, m_movePath);
}

/**
 * @brief TrajectoryPlanning::SetPathSegment    将原始路径转换成交融路径
 * @param pathf
 * @param patht
 * @return
 */
int TrajectoryPlanning::SetPathSegment(vector<PathInformation> pathf, vector<PathInformation> patht)
{
    for (int i = 0; i < pathf.size() - 1; ++i) {
        GetPathSegment(pathf[i], pathf[i+1], patht);
    }
}

/**
 * @brief TrajectoryPlanning::SetLinePath 设置原始直线信息
 * @param pf
 * @param pt
 * @return
 */
PathInformation TrajectoryPlanning::SetLinePath(PointInformation pf, PointInformation pt)
{
    PathInformation path;
    path.pathType = Line;
    path.startpoint = pf;
    path.endpoint = pt;
    path.maxVelocity = pt.maxVelocity;
    path.maxAcceleration = pt.maxAcceleration;
    path.maxJerk = pt.maxJerk;
    Vector3d start2end = path.endpoint.point - path.startpoint.point;
    path.displacement = start2end.norm();
    return path;
}

/**
 * @brief TrajectoryPlanning::SetArcPath 设置原始圆弧信息
 * @param pf
 * @param pi
 * @param pt
 * @return
 */
PathInformation TrajectoryPlanning::SetArcPath(PointInformation pf, PointInformation pi, PointInformation pt)
{
    PathInformation path;
    path.pathType = Arc;
    path.startpoint = pf;
    path.intermediatepoint = pi;
    path.endpoint = pt;
    path.maxVelocity = pt.maxVelocity;
    path.maxAcceleration = pt.maxAcceleration;
    path.maxJerk = pt.maxJerk;
    GetArcParameter(path.startpoint.point, path.intermediatepoint.point, path.endpoint.point,
                    path.center, path.radius, path.normal, path.theta);
    path.arclength = path.radius * path.theta;
    return path;
}

/**
 * @brief TrajectoryPlanning::SetCurvePath 设置曲线路径，继承前一路径的约束
 * @param ctrlpoint
 * @param path
 */
void TrajectoryPlanning::SetCurvePath(vector<PointInformation> ctrlpoint, PathInformation &path)
{
    path.pathType = Curve;
    path.controlpoints = ctrlpoint;
    GetCurveParameter(path);
}

/**
 * @brief TrajectoryPlanning::GetPathSegment
 * @param path1
 * @param path2
 * @param path
 */
void TrajectoryPlanning::GetPathSegment(PathInformation path1, PathInformation &path2,
                                        vector<PathInformation> &path)
{
    // 存放计算交融后放入容器的路径
    PathInformation tempPath = path1;
    // 计算曲线4个控制点（包含点信息）
    vector<PointInformation> ctrlPoint = GetControlPoints(path1, path2);
    // p0作为前一路径的终点
    tempPath.endpoint = ctrlPoint[0];
    path.push_back(tempPath);
    // 曲线继承前一路径的约束
    SetCurvePath(ctrlPoint, tempPath);
    path.push_back(tempPath);
    // 将曲线的最后一个控制点信息作为下一路径的起点信息
    path2.startpoint = tempPath.controlpoints[3];
}
