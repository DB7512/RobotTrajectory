#include "velocityplanning.h"

#define EPSILON 1e-6

int VelocityPlanning::velocityPlan(double s, double vs, double ve, double vmax, double amax, double jmax, VectorXd &para, double time)
{
    // 判断是否是匀速
    if (fabs(vs - ve) < EPSILON && fabs(vs - vmax) < EPSILON) {
        para << 0.0, s/vs, 0.0, 0.0, 0.0, s, vs, vs, vs, vs, amax, -amax, amax, -amax, jmax, -jmax;
        return 1;
    }

    double ta, tv, td, tja, tjd, vlim, amin, alima, alimd, jmin;
    double sa, sd;
    DecType decType;
    AccType accType;
    signAmaxJmax(vs, vmax, amax, jmax, amin, jmin); // 判断加速/减速
    // 判断vs-vmax时的变速类型
    if (vmax - vs <= pow(amax, 2) / jmax) {
        ta = 2 * sqrt((vmax - vs) / jmax);  // 欠三角形
        accType = Tri;
    } else {
        ta = (vmax - vs) / amax + amax / jmax;  // 梯形
        accType = Tra;
    }
    sa = 0.5 * ta * (vs + vmax);    // 计算临界三角形位移

    // 判断判断vmax-ve时的变速类型
    signAmaxJmax(vmax, ve, amax, jmax, amin, jmin);
    if (ve - vmax <= pow(amax, 2) / jmax) {
        td = 2 * sqrt((ve - vmax) / jmax);  // 欠三角形
        decType = Tri;
    } else {
        td = (ve - vmax) / amax + amax / jmax;  // 梯形
        decType = Tra;
    }
    sd = 0.5 * td * (vmax + ve);    // 计算临界三角形面积

    if (fabs(s - (sa + sd)) < EPSILON) {
        // vs-vmax-ve都是临界三角形变速
        ta = 2 * amax / jmax;
        tv = 0.0;
        td = 2 * amin / jmin;
        tja = amax / jmax;
        tjd = amin / jmin;
        vlim = vmax;
        alima = amax;
        alimd = amin;
        para << ta,tv,td,tja,tjd,s,vs,ve,vmax,vlim,amax,amin,alima,alimd,jmax,jmin;
        return 1;
    } else if (s - (sa + sd) > EPSILON) {
        // 能达到vmax,存在匀速段
        tv = (s - (sa + sd)) / vmax;
        vlim = vmax;

        if (accType == Tri0) {
            ta = 2 * amax / jmax;
            tja = amax / jmax;
        } else if (accType == Tri) {
            ta = 2 * sqrt((vmax - vs) / jmax);
            tja = 0.5 * ta;
        } else if (accType == Tra) {
            ta = (vmax - vs) / amax + amax / jmax;
            tja = amax / jmax;
        }

        if (decType == Tri0) {
            ta = 2 * amin / jmin;
            tja = amin / jmin;
        } else if (decType == Tri) {
            ta = 2 * sqrt((ve - vmax) / jmin);
            tja = 0.5 * ta;
        } else if (decType == Tra) {
            ta = (ve - vmax) / amin + amin / jmin;
            tja = amin / jmin;
        }
        alima = tja * jmax;
        alimd = tjd * jmin;
        para << ta,tv,td,tja,tjd,s,vs,ve,vmax,vlim,amax,amin,alima,alimd,jmax,jmin;
        return 1;
    } else if ((sa + sd) - s > EPSILON) {
        // 计算vs-ve时s
        double s0;
        signAmaxJmax(vs, ve, amax, jmax, amin, jmin);
        if (ve - vs - pow(amax, 2) / jmax > EPSILON) {
            ta = (ve - vs) / amax + amax / jmax;    // 梯形变速
            accType = Tra;
        } else {
            ta = 2 * sqrt((ve - vs) / jmax);    // 三角形变速
            accType = Tri;
        }
        s0 = 0.5 * ta  * (vs + ve);

        if (s0 - s > EPSILON) {
            // 在[min{vs,ve},max{vs,ve}]中搜索ve'
            double vl = min(vs, ve);
            double vu = max(vs, ve);
            if (accType == Tra) {
                Vector<double> coeff;
                coeff.resize(3);
                vector<double> result;
                result.resize(2);
                coeff[0] = jmax;
                coeff[1] = pow(amax,2) + 2 * jmax * vs;
                coeff[2] = 2 * pow(amax,2) * vs - 2 * amax * jmax * s;
                solve2thEquation(coeff, result);
                for (int i = 0; i < result.size(); ++i) {
                    if ((result[i] + vs >= vl && result[i] + vs <= vu)) {
                        ve = result[i] + vs;
                    }
                }
                if (((result[0] + vs >= vl) && (result[0] + vs <= vu))
                    || ((result[1] + vs >= vl) && (result[1] + vs <= vu))) {
                    if ((result[0] + vs >= vl) && (result[0] + vs <= vu)) {
                        ve = result[0] + vs;
                    } else if ((result[1] + vs >= vl) && (result[1] + vs <= vu)) {
                        ve = result[1] + vs;
                    }
                }




            }
        }
    }

}

int VelocityPlanning::velocityPlan(double s, double vs, double ve, double vmax, double amax, double jmax, VectorXd &para)
{

}

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

void VelocityPlanning::signAmaxJmax(double vf, double vt, double &amax, double &jmax, double &amin, double &jmin)
{
    if (vf < vt) {
        amax = amax;
        jmax = jmax;
    } else {
        amax = -amax;
        jmax = -jmax;
    }
    amin = -amax;
    jmin = -jmax;
}

int VelocityPlanning::solve2thEquation(vector<double> coeff, vector<double> &result)
{
    double a = coeff[0];
    double b = coeff[1];
    double c = coeff[2];
    double delta = pow(b, 2) - 4 * a * c;
    if (delta > 0.0) {
        result.push_back((-b + sqrt(delta) / (2 * a)));
        result.push_back((-b - sqrt(delta) / (2 * a)));
        return 1;
    } else {
        return 0;
    }
}

int VelocityPlanning::solve3thEquation(vector<double> coeff, vector<double> &result)
{
    double a = coeff[0];
    double b = coeff[1];
    double c = coeff[2];
    double d = coeff[3];
    double A = pow(b, 2) - 3 * a * c;
    double B = b * c - 9 * a * d;
    double C = pow(c, 2) - 3 * b * d;
    double delta = pow(B, 2) - 4 * A * C;
    if (fabs(A) < EPSILON && fabs(B) < EPSILON) {
        result[0] = -c / b;
        result[1] = result[0];
        result[2] = result[0];
        return 1;
    }
    if (delta > EPSILON) {
        int zero_1 = 0, zero_2 = 0;
        double Y1 = A * b + 3 * a * (-B + sqrt(delta)) / 2;
        double Y2 = A * b + 3 * a * (-B - sqrt(delta)) / 2;
        if (Y1 < 0) {
            Y1 = -Y1;
            zero_1 = 1;
        }
        if (Y2 < 0) {
            Y2 = -Y2;
            zero_2 = 1;
        }
        double Y_1 = pow(Y1, 1.0 / 3);
        double Y_2 = pow(Y2, 1.0 / 3);
        if (zero_1) {
            Y_1 = -Y_1;
        }
        if (zero_2) {
            Y_2 = -Y_2;
        }
        result[0] = (-b - (Y_1 + Y_2)) / (3 * a);
        double unreal = (sqrt(3) / 2 * (Y_1 - Y_2)) / (3 * a);
        if (fabs(unreal) < EPSILON) {
            result[1] = (-b + 0.5 * (Y_1 + Y_2)) / (3 * a);
            result[2] = result[1];
        } else {
            result[1] = 0;
            result[2] = 0;
        }
        return 1;
    }
    if (fabs(delta) < EPSILON && fabs(A) > EPSILON) {
        double K = B / A;
        result[0] = -b / a + K;
        result[1] = -0.5 * K;
        result[2] = result[1];
        return 1;
    }
    if ((delta < -EPSILON) && (A > EPSILON)) {
        double T = (2 * A * b - 3 * a * B) / (2 * sqrt(pow(A, 3)));
        double theta = acos(T);
        result[0] = (-b - 2 * sqrt(A) * cos(theta / 3)) / (3 * a);
        result[1] = (-b + sqrt(A) * (cos(theta / 3) + sqrt(3) * sin(theta / 3))) / (3 * a);
        result[2] = (-b + sqrt(A) * (cos(theta / 3) - sqrt(3) * sin(theta / 3))) / (3 * a);
        return 1;
    }
}

int VelocityPlanning::solve4thEquation(vector<double> coeff, vector<double> &result)
{
    Matrix<double, 4, 4> matrixCoeff;
    Matrix<complex<double>,Dynamic,Dynamic> matrixSigenvalues;
    VectorXd rl = VectorXd(4);
    VectorXd ig = VectorXd(4);

    matrixCoeff << -coeff[1]/coeff[0], -coeff[2]/coeff[0], -coeff[3]/coeff[0], -coeff[4]/coeff[0],
        1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1, 0;

    matrixSigenvalues = matrixCoeff.eigenvalues();
    rl << matrixSigenvalues.real();
    ig << matrixSigenvalues.imag();

    for (int i = 0; i < 4; ++i) {
        if (fabs(ig(i)) < EPSILON) {
            result[i] = rl(i);
        } else {
            result[i] = 0.0;
        }
    }
}

int VelocityPlanning::solve5thEquation(vector<double> coeff, vector<double> &result)
{
    Matrix<double, 5, 5> matrixCoeff;
    Matrix<complex<double>,Dynamic,Dynamic> matrixSigenvalues;
    VectorXd rl = VectorXd(5);
    VectorXd ig = VectorXd(5);

    matrixCoeff << -coeff[1]/coeff[0], -coeff[2]/coeff[0], -coeff[3]/coeff[0], -coeff[4]/coeff[0], -coeff[5]/coeff[0],
        1, 0, 0, 0, 0,
        0, 1, 0, 0, 0,
        0, 0, 1, 0, 0,
        0, 0, 0, 1, 0;

    matrixSigenvalues = matrixCoeff.eigenvalues();
    rl << matrixSigenvalues.real();
    ig << matrixSigenvalues.imag();

    for (int i = 0; i < 5; ++i) {
        if (fabs(ig(i)) < EPSILON) {
            result[i] = rl(i);
        } else {
            result[i] = 0.0;
        }
    }
}
