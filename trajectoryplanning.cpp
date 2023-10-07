#include "trajectoryplanning.h"
#include "mathfunction.h"
#include <QTextStream>
#include <QFile>
#include <QDebug>
#include "cmath"

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
    ret = TrajectoryTime4(trajectory_ftime, Q, v_0, v_1, vmax, amax, jmax, para);
    qDebug()<<"ret"<<ret<<"time"<<trajectory_ftime;
    if(ret < 0) return -1;
    //corrent parameters
    CorrentionParameters(para, peroid);
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
//    float trajectory_ftime = 0.0; //trajectory time
    double trajectory_ftime = 0.0; //trajectory time
    int interpolation_peroid_num = 0; //need period number
    float t = 0.0; //time
    //calculate time
    int ret = 0;
    ret = TrajectoryTime5(trajectory_ftime, Q, v_0, v_1, vmax, amax, jmax, parad);
    qDebug()<<"ret"<<ret<<"time"<<trajectory_ftime;
    if(ret < 0) return -1;
    double temp = (parad(0) + parad(1) + parad(2)) / (1.0 / peroid);
    interpolation_peroid_num = round(temp);
    if(parad(0) + parad(1) + parad(2) < 1e-6) return false;
    double lambda[interpolation_peroid_num];//归一化参数
    vector<float> infpoint(3);
    double q_dq[3] = {0.0}; //resulets: displacement, velocity, acceleration
    for (int i = 0; i < interpolation_peroid_num - 1; i++) {
        //init vector
        lambda[i] = 0.0;
        t = i*(1.0/peroid);
        //calculate trajectory
        TrajectoryCalculationD(t, parad, q_dq);
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
    float delta_T = 0.0;
    float Tm = Te + Ta + Tv + Td;
    if(m_trajectorytype == TraUniTra || m_trajectorytype == TriUniTri) {
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
    } else if(m_trajectorytype == TriUniTra || m_trajectorytype == UniTra
               || m_trajectorytype == UniTri) { //extend the deceleration time
        delta_T    = Te * (vlim + v_1) / (vlim - v_1);
        Tv_new      = Tv - delta_T;
        if(Tv_new < -1e-5) {
            if(m_trajectorytype == TriUniTra) {
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
    } else if(m_trajectorytype == TraUniTri || m_trajectorytype == TraUni
               || m_trajectorytype == TriUni) { //extend the acceleration time
        delta_T    = Te * (vlim + v_0) / (vlim - v_0);
        Tv_new      = Tv - delta_T;
        if(Tv_new < -1e-5) {
            if(m_trajectorytype == TraUniTri) {
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
    } else if(m_trajectorytype == Uni){
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
    if(temppath.size() < 1) {
        return;
    } else if (temppath.size() == 1) {
        PathInformation path = temppath.front();
    } else {
        vector<PathInformation>::const_iterator path1 = temppath.begin();
        vector<PathInformation>::const_iterator path2 = path1;
        path2++;
        PathInformation currentpath = *path1;
        PathInformation finalpath;
        while(path2 != temppath.end()) {
            PathInformation nextpath = *path2;
            if(nextpath.pathType == Arc) {
                finalpath = nextpath;
                finalpath.pathType = Line2arc;
                finalpath.startpoint = currentpath.startpoint;
                finalpath.endpoint = nextpath.endpoint;
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
                if(path1 != temppath.end()) {
                    path2++;
                } else {
                    CalculatePathParameters(path1);
                    pathes.push_back(path1);
                    return;
                }
            }
        }
    }
}

void TrajectoryPlanning::CompoundTrajectory(vector<PointInformation> waypoints, vector<PathInformation> &path)
{
    if(waypoints.size() < 2)
        return;
    vector<PointInformation>::const_iterator point1 = waypoints.begin();
    vector<PointInformation>::const_iterator point2 = point1;
    point2++;
    vector<PointInformation>::const_iterator point3;
    PointInformation startpoint = *point1;
    while(point2 != waypoints.end()) {
        point3 = point2;
        point3++;
        PointInformation intermediatepoint;
        PointInformation endpoint;
        intermediatepoint = *point2;
        endpoint = *point3;
        PathInformation arcpath;
        ArcSegmentLineToLine(startpoint, intermediatepoint, endpoint, arcpath);
        PathInformation linepath;
        SetLinePathSegment(linepath, startpoint, endpoint, arcpath);
        path.push_back(linepath);
        path.push_back(arcpath);
        startpoint = arcpath.endpoint;
        point2 = point3;
    }
    PointInformation intermediatepoint;
    PointInformation endpoint;
    intermediatepoint = *point2;
    endpoint = *point3;
    PathInformation linepath;
    PathInformation arcpath;
    SetLinePathSegment(linepath, startpoint, endpoint, arcpath);
    path.push_back(linepath);
}

void TrajectoryPlanning::CalculatePathParameters(PathInformation path)
{
    double time = 0.0;
    double Q = 0.0;
    if(path.pathType = Line2arc) {
        Q = path.displacement + path.arclength;
    } else if(path.pathType == Line) {
        Q = path.displacement;
    } else if(path.pathType == Circular) {
        Q = path.arclength;
    }
    TrajectoryTime(path.accelerationType, time, Q, path.startpoint.velocity, path.endpoint.velocity, path.maxVelocity, path.maxAcceleration, path.maxJerk, path.constraints);
    //添加修正参数函数
    CorrentionParameters(path.constraints, 100);
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

void TrajectoryPlanning::ArcSegmentforLinetoLine(vector<PointInformation> waypoints, vector<PathInformation> pathes)
{
    if(size(waypoints) != 3) return;
    PointInformation waypoint;
    Vector3d startpoint;
    Vector3d endpoint;
    Vector3d middlepoint;
    startpoint = waypoints[0].point;
    endpoint = waypoints[1].point;
    middlepoint = waypoints[2].point;
    double theta = 0.0;
    Vector3d line1;
    Vector3d line2;
    Vector3d startline;
    Vector3d endline;
    startline = startpoint - middlepoint;
    endline = endpoint - middlepoint;
    line1 = 0.5 * startline;
    line2 = 0.5 * endline;
    theta = acos(line1.dot(line2) / (line1.norm() * line2.norm()));
    double l = waypoints[1].radius / tan(theta * 0.5);
    if((line1.norm() - l) < 1e-6 || (line2.norm() - l)  < 1e-6) {
        if((line1.norm() - line2.norm()) < 1e-6) {
            if(line1.norm() > 1e-6) {
                l = line1.norm();
                waypoints[1].radius = tan(theta * 0.5) * l;
            } else {
                waypoints[1].radius = 0;
                return;
            }
        } else {
            if(line2.norm() > 1e-6) {
                l = line2.norm();
                waypoints[1].radius = tan(theta * 0.5) * l;
            } else {
                waypoints[1].radius = 0;
                return;
            }
        }
    }
    Vector3d point1;
    Vector3d point2;
    point1 = startline / startline.norm() * l + startpoint;
    point2 = endline / endline.norm() * l + endpoint;
    Vector3d temppoint;
    double d = l / cos(theta * 0.5);
    temppoint = (0.5 * (point2 - point1) + point1) - middlepoint;
    Vector3d center;
    center = temppoint * d / temppoint.norm() - middlepoint;
    PathInformation path;
    SetLinePathSegment(path,waypoints[0],waypoints[1],point1);
    pathes.push_back(path);
    //添加初始化函数
    SetCricularPathSegment(path,);




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
    path.startPoint = startpoint.point;
    path.endPoint = endpoint.point;
    path.startPose = startpoint.pose;
    path.endPose = endpoint.pose;
    path.maxVelocity = max(startpoint.maxVelocity, middlepoint.maxVelocity, endpoint.maxVelocity);
    path.maxAcceleration = max(startpoint.maxAcceleration, middlepoint.maxAcceleration, endpoint.maxAcceleration);
    path.maxJerk = max(startpoint.maxJerk, middlepoint.maxJerk, endpoint.maxJerk);
    path.startVelocity = startpoint.velocity;
    path.endVelocity = endpoint.velocity;
}

void TrajectoryPlanning::SetArcPathSegment(PathInformation &path, PointInformation intermediatepoint, Vector3d arcstart, Vector3d arcend,
                                           Vector3d arccenter, double radius, double theta)
{
    path.startpoint             = intermediatepoint;
    path.endpoint               = intermediatepoint;
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
