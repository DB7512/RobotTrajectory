#include "velocityplanning.h"
#include "trajectoryplanning.h"
#define EPSILON 1e-6

/**
 * @brief VelocityPlanning::PositionPLan
 * @param l
 * @param vs
 * @param ve
 * @param vmax
 * @param amax
 * @param jmax
 * @param para
 * @param time
 * @return 返回值：-1：规划失败；0：规划成功；1：其他
 */
int VelocityPlanning::SetVelocityPlan(double l, double &vs, double &ve,
                                   double &vmax, double amax, double jmax,
                                   VectorXd &para, double &time)
{
    para = VectorXd::Zero(16);
    double ta = 0.0, td = 0.0, tv = 0.0, tja = 0.0, tjd = 0.0;
    double alima = 0.0, alimd = 0.0, vlim = 0.0;
    double ja = 0.0, jd = 0.0;
    double amin = -amax, jmin = -jmax;
    // 处理匀速
    if (fabs(vs - ve) < EPSILON && fabs(vs - vmax) < EPSILON) {
        tv = l / vs;
        ve = vs;
        vlim = vs;
        time = tv;
        para << ta,tv,td,tja,tjd,
            l,vs,ve,vmax,vlim,amax,amin,alima,alimd,jmax,jmin;
        return 0;
    }
    // 计算vs-vmax时的位移s1
    int sign1 = 1;
    if (vmax > vs) {
        sign1 = 1;
    } else {
        sign1 = -1;
    }
    double a = sign1 * amax;
    double j = sign1 * jmax;
    double t1 = 0.0;
    if (fabs(vmax - vs) <= fabs(pow(a, 2) / j)) {
        t1 = 2 * sqrt((vmax - vs) / j);  // (欠)三角形
        tja = 0.5 * t1;
        alima = j * tja;
    } else {
        t1 = (vmax - vs) / a + a / j;  // 梯形
        tja = a / j;
        alima = a;
    }
    ja = j;
    double s1 = 0.5 * t1 * (vs + vmax);    // 计算vs-vmax时的位移
    // 计算vmax-ve时的位移s2
    int sign2 = 1;
    if (ve > vmax) {
        sign2 = 1;
    } else {
        sign2 = -1;
    }
    a = sign2 * amax;
    j = sign2 * jmax;
    double t2 = 0.0;
    if (fabs(ve - vmax) <= fabs(pow(a, 2) / j)) {
        t2 = 2 * sqrt((ve - vmax) / j);  // (欠)三角形
        tjd = 0.5 * t2;
        alimd = j * tjd;
    } else {
        t2 = (ve - vmax) / a + a / j;  // 梯形
        tjd = a / j;
        alimd = a;
    }
    jd = j;
    double s2 = 0.5 * t2 * (vmax + ve);    // 计算vs-vmax时的位移
    double l1 = s1 + s2;
    if (l1 <= l) {  // 存在匀速段
        tv = (l - l1) / vmax;
        time = t1 + tv + t2;
        ta = t1;
        td = t2;
        vlim = vmax;
        jmax = ja;
        jmin = jd;
        para << ta,tv,td,tja,tjd,
            l,vs,ve,vmax,vlim,amax,amin,alima,alimd,jmax,jmin;
        return 0;
    } else {    // l1 > l
        // 计算vs-ve时的位移
        int sign3 = 1;
        if (ve > vs) {
            sign3 = 1;
        } else {
            sign3 = -1;
        }
        a = sign3 * amax;
        j = sign3 * jmax;
        double t = 0.0;
        if (fabs(ve - vs) <= fabs(pow(a, 2) / j)) {
            t = 2 * sqrt((ve - vs) / j);  // (欠)三角形
        } else {
            t = (ve - vs) / a + a / j;  // 梯形
        }
        double l2 = 0.5 * t * (vs + ve);    // 计算vs-ve时的位移

        if (l2 > l) {   // 在vs和ve之间查找ve'
            double left = min(vs, ve);
            double right = max(vs, ve);
            int num = 0;
            while (right - left > EPSILON) {
                if (num > 20) return -1;
                ve = (left + right) / 2.0;
                if (fabs(ve - vs) <= fabs(pow(a, 2) / j)) {
                    t = 2 * sqrt((ve - vs) / j);  // (欠)三角形
                    tja = 0.5 * t;
                    alima = j * tja;
                } else {
                    t = (ve - vs) / a + a / j;  // 梯形
                    tja = a / j;
                    alima = a;
                }
                ja = j;
                double s = 0.5 * t * (vs + ve);    // 计算vs-ve'时的位移
                if (fabs(s - l) < EPSILON) {
                    time = t;
                    ta = t;
                    jmax = ja;
                    para << ta,tv,td,tja,tjd,
                        l,vs,ve,vmax,vlim,amax,amin,alima,alimd,jmax,jmin;
                    return 0;
                } else if (s - l > EPSILON) {
                    right = ve;
                } else if (s - l < -EPSILON) {
                    left = ve;
                }
                num ++;
            }
        } else {    // l2 < l
            if(vmax < vs && vmax < ve) {    // 在(vmax, min{vs,ve})之间查找v'
                double left = vmax;
                double right = min(vs, ve);
                int num = 0;
                while (right - left > EPSILON) {
                    if (num > 20) return -1;
                    vlim = (left + right) / 2.0;
                    // 计算vs-vlim时的位移s1
                    int sign1 = 1;
                    if (vlim > vs) {
                        sign1 = 1;
                    } else {
                        sign1 = -1;
                    }
                    double a = sign1 * amax;
                    double j = sign1 * jmax;
                    double t1 = 0.0;
                    if (fabs(vlim - vs) <= fabs(pow(a, 2) / j)) {
                        t1 = 2 * sqrt((vlim - vs) / j);  // (欠)三角形
                        tja = 0.5 * t1;
                        alima = j * tja;
                    } else {
                        t1 = (vlim - vs) / a + a / j;  // 梯形
                        tja = a / j;
                        alima = a;
                    }
                    ja = j;
                    double s1 = 0.5 * t1 * (vs + vlim);    // 计算vs-vlim时的位移
                    // 计算vlim-ve时的位移s2
                    int sign2 = 1;
                    if (ve > vlim) {
                        sign2 = 1;
                    } else {
                        sign2 = -1;
                    }
                    a = sign2 * amax;
                    j = sign2 * jmax;
                    double t2 = 0.0;
                    if (fabs(ve - vlim) <= fabs(pow(a, 2) / j)) {
                        t2 = 2 * sqrt((ve - vlim) / j);  // (欠)三角形
                        tjd = 0.5 * t2;
                        alimd = j * tjd;
                    } else {
                        t2 = (ve - vlim) / a + a / j;  // 梯形
                        tjd = a / j;
                        alimd = a;
                    }
                    jd = j;
                    double s2 = 0.5 * t2 * (vlim + ve);    // 计算vs-vmax时的位移
                    double s = s1 + s2;
                    if (fabs(s - l) < EPSILON) {
                        time = t1 + t2;
                        ta = t1;
                        td = t2;
                        jmax = ja;
                        jmin = jd;
                        para << ta,tv,td,tja,tjd,
                            l,vs,ve,vmax,vlim,amax,amin,alima,alimd,jmax,jmin;
                        return 0;
                    } else if (s - l > EPSILON) {
                        right = vlim;
                    } else if (s - l < -EPSILON) {
                        left = vlim;
                    }
                    num ++;
                }
            } else if (vmax > vs && vmax > ve) {
                double left = max(vs, ve);
                double right = vmax;
                int num = 0;
                while (right - left > EPSILON) {
                    if (num > 20) return 0;
                    vlim = (left + right) / 2.0;
                    // 计算vs-vlim时的位移s1
                    int sign1 = 1;
                    if (vlim > vs) {
                        sign1 = 1;
                    } else {
                        sign1 = -1;
                    }
                    double a = sign1 * amax;
                    double j = sign1 * jmax;
                    double t1 = 0.0;
                    if (fabs(vlim - vs) <= fabs(pow(a, 2) / j)) {
                        t1 = 2 * sqrt((vlim - vs) / j);  // (欠)三角形
                        tja = 0.5 * t1;
                        alima = j * tja;
                    } else {
                        t1 = (vlim - vs) / a + a / j;  // 梯形
                        tja = a / j;
                        alima = a;
                    }
                    ja = j;
                    double s1 = 0.5 * t1 * (vs + vlim);    // 计算vs-vlim时的位移
                    // 计算vlimx-ve时的位移s2
                    int sign2 = 1;
                    if (ve > vlim) {
                        sign2 = 1;
                    } else {
                        sign2 = -1;
                    }
                    a = sign2 * amax;
                    j = sign2 * jmax;
                    double t2 = 0.0;
                    if (fabs(ve - vlim) <= fabs(pow(a, 2) / j)) {
                        t2 = 2 * sqrt((ve - vlim) / j);  // (欠)三角形
                        tjd = 0.5 * t2;
                        alimd = j * tjd;
                    } else {
                        t2 = (ve - vlim) / a + a / j;  // 梯形
                        tjd = a / j;
                        alimd = a;
                    }
                    jd = j;
                    double s2 = 0.5 * t2 * (vlim + ve);    // 计算vs-vlim时的位移
                    double s = s1 + s2;
                    if (fabs(s - l) < EPSILON) {
                        time = t1 + t2;
                        ta = t1;
                        td = t2;
                        jmax = ja;
                        jmin = jd;
                        para << ta,tv,td,tja,tjd,
                            l,vs,ve,vmax,vlim,
                            amax,amin,alima,alimd,jmax,jmin;
                        return 0;
                    } else if (s - l > EPSILON) {
                        right = vlim;
                    } else if (s - l < -EPSILON) {
                        left = vlim;
                    }
                    num ++;
                }
            } else if (vmax > min(vs, ve) && vmax < max(vs, ve)) {
                if (vs > ve) {
                    vlim = vs;
                    int sign1 = -1;
                    double a = sign1 * amax;
                    double j = sign1 * jmax;
                    double t1 = 0.0;
                    if (fabs(vlim - vs) <= fabs(pow(a, 2) / j)) {
                        t1 = 2 * sqrt((vlim - vs) / j);  // (欠)三角形
                        tjd = 0.5 * t;
                        alimd = j * tjd;
                    } else {
                        t1 = (vlim - vs) / a + a / j;  // 梯形
                        tjd = a / j;
                        alimd = a;
                    }
                    jd = j;
                    double s1 = 0.5 * t1 * (vs + vlim);    // 计算vs-vlim时的位移
                    double tv = (l - s1) / vlim;
                    time = t1 + tv;
                    td = t1;
                    jd = j;
                    jmin = jd;
                    para << ta,tv,td,tja,tjd,
                        l,vs,ve,vmax,vlim,
                        amax,amin,alima,alimd,jmax,jmin;
                    return 0;
                } else {
                    vlim = ve;
                    int sign1 = 1;
                    double a = sign1 * amax;
                    double j = sign1 * jmax;
                    double t1 = 0.0;
                    if (fabs(vlim - vs) <= fabs(pow(a, 2) / j)) {
                        t1 = 2 * sqrt((vlim - vs) / j);  // (欠)三角形
                        tja = 0.5 * t;
                        alima = j * tja;
                    } else {
                        t1 = (vlim - vs) / a + a / j;  // 梯形
                        tja = a / j;
                        alima = a;
                    }
                    ja = j;
                    double s1 = 0.5 * t1 * (vs + vlim);    // 计算vs-vlim时的位移
                    double tv = (l - s1) / vlim;
                    time = t1 + tv;
                    ta = t1;
                    ja = j;
                    jmax = ja;
                    para << ta,tv,td,tja,tjd,
                        l,vs,ve,vmax,vlim,
                        amax,amin,alima,alimd,jmax,jmin;
                    return 0;
                }
            }
        }
    }
}

int VelocityPlanning::TimePLan(double l, double &vs, double &ve, double &vmax, double amax, double jmax, VectorXd &para, double time)
{
    double vavg = l / time;
    if (vavg < vs && vavg < ve) {
        // 计算vs-vavg时的位移
        int sign1 = 1;
        if (vavg > vs) {
            sign1 = 1;
        } else {
            sign1 = -1;
        }
        double a = sign1 * amax;
        double j = sign1 * jmax;
        double t1 = 0.0;
        if (fabs(vavg - vs) <= fabs(pow(a, 2) / j)) {
            t1 = 2 * sqrt((vavg - vs) / j);  // (欠)三角形
        } else {
            t1 = (vavg - vs) / a + a / j;  // 梯形
        }
        double s1 = 0.5 * t1 * (vs + vavg);    // 计算vs-vmax时的位移
        if (s1 > l) {   // 给定参数不合理，给定时间的速度规划失败
            return -1;
        } else {
            // 计算vavg-ve时的位移s2
            int sign2 = 1;
            if (ve > vavg) {
                sign2 = 1;
            } else {
                sign2 = -1;
            }
            a = sign2 * amax;
            j = sign2 * jmax;
            double t2 = 0.0;
            if (fabs(ve - vavg) <= fabs(pow(a, 2) / j)) {
                t2 = 2 * sqrt((ve - vavg) / j);  // (欠)三角形
            } else {
                t2 = (ve - vavg) / a + a / j;  // 梯形
            }
            double s2 = 0.5 * t2 * (vavg + ve);    // 计算vs-vmax时的位移
            // 补充如何同时满足l和t的算法
        }
    }
}

int VelocityPlanning::TrajectoryTime(double &time, float l, float vs, float ve, float vmax, float amax, float jmax, VectorXd &para)
{
    //judge whether the minimun diaplacement is satisfied
    double Tj[2] = {0.0};
    int index = 0;
    Tj[0] = sqrt(fabs(ve - vs)/jmax);
    Tj[1] = amax/jmax;
    if(Tj[0] > Tj[1]) {
        index = 1;
    }else {
        index = 0;
    }
    if(index == 0) {
        if(l < Tj[index] * (vs + ve)) {
            return -1;
        }
    }else {
        if(l < 0.5 * (vs + ve) * (Tj[index] + fabs(ve - vs)/amax)) {
            return -1;
        }
    }
    double tja = 0.0, tjd = 0.0, Ta = 0.0, td = 0.0, tv = 0.0;
    double a_lima = 0.0, a_limd = 0.0;
    double vlim = 0.0;
    double jmin = -jmax;
    double delta = 1e-6;
    double Sa = 0.0;
    double Sd = 0.0;
    double vl = 0.0;
    double vu = 0.0;
    if(fabs(vs - ve) <= delta && fabs(vs - vmax) <= delta) {
        tv = l / vmax;
        vlim = vmax;
    } else {
        //amax is not reached
        if ((vmax - vs) * jmax < pow(amax, 2)) {
            tja = sqrt((vmax - vs) / jmax);
            Ta = 2 * tja;
            a_lima = jmax * tja;
        } else {
            tja = amax / jmax;
            Ta = tja + (vmax - vs) / amax;
            a_lima = amax;
        }
        //amin is not reached
        if ((vmax - ve) * jmax < pow(amax, 2)) {
            tjd = sqrt((vmax - ve) / jmax);
            td = 2 * tjd;
            a_limd = -jmax * tjd;
        } else {
            tjd = amax / jmax;
            td = tjd + (vmax - ve) / amax;
            a_limd = -amax;
        }
        Sa = 0.5 * (vs + vmax) * Ta;
        Sd = 0.5 * (ve + vmax) * td;
        if(l - Sa - Sd > delta) {
            tv = (l - Sa - Sd) / vmax;
            vlim = vmax;
        } else if(fabs(l - Sa - Sd) <= delta){
            tv = 0.0;
            vlim = vmax;
        } else {
            tv = 0.0;
            if(vs < ve) {
                vl = ve;
            } else {
                vl = vs;
            }
            vu = vmax;
            vmax = 0.5 * (vl + vu);
            //amax is not reached
            if ((vmax - vs) * jmax < pow(amax, 2)) {
                tja = sqrt((vmax - vs) / jmax);
                Ta = 2 * tja;
                a_lima = jmax * tja;
            } else {
                tja = amax / jmax;
                Ta = tja + (vmax - vs) / amax;
                a_lima = amax;
            }
            //amin is not reached
            if ((vmax - ve) * jmax < pow(amax, 2)) {
                tjd = sqrt((vmax - ve) / jmax);
                td = 2 * tjd;
                a_limd = -jmax * tjd;
            } else {
                tjd = amax / jmax;
                td = tjd + (vmax - ve) / amax;
                a_limd = -amax;
            }
            Sa = 0.5 * (vs + vmax) * Ta;
            Sd = 0.5 * (ve + vmax) * td;
            while(fabs(l - Sa - Sd) > delta) {
                if(vu - vl < delta)
                    break;
                if(l - Sa - Sd > delta) {
                    vl = vmax;
                    vu = vu;
                } else if(l - Sa - Sd < -delta) {
                    vl = vl;
                    vu = vmax;
                }
                vmax = 0.5 * (vl + vu);
                //amax is not reached
                if ((vmax - vs) * jmax < pow(amax, 2)) {
                    tja = sqrt((vmax - vs) / jmax);
                    Ta = 2 * tja;
                    a_lima = jmax * tja;
                } else {
                    tja = amax / jmax;
                    Ta = tja + (vmax - vs) / amax;
                    a_lima = amax;
                }
                //amin is not reached
                if ((vmax - ve) * jmax < pow(amax, 2)) {
                    tjd = sqrt((vmax - ve) / jmax);
                    td = 2 * tjd;
                    a_limd = -jmax * tjd;
                } else {
                    tjd = amax / jmax;
                    td = tjd + (vmax - ve) / amax;
                    a_limd = -amax;
                }
                Sa = 0.5 * (vs + vmax) * Ta;
                Sd = 0.5 * (ve + vmax) * td;
            }
            vlim = vmax;
        }
    }
    time = Ta + tv + td;
    para << Ta, tv, td, tja, tjd, 0, l, vs, ve, vlim, amax, a_lima, a_limd, jmax, jmin;
    //classification
    if(tv > 0.0) {
        if(Ta > 2*tja) {
            if(td > 2*tjd) {
                m_accelerationtype = TraZeroTra;
            } else if(td > 0.0) {
                m_accelerationtype = TraZeroTri;
            } else {
                m_accelerationtype = TraZero;
            }
        } else if(Ta > 0.0) {
            if(td > 2*tjd) {
                m_accelerationtype = TriZeroTra;
            } else if(td > 0.0) {
                m_accelerationtype = TriZeroTri;
            } else {
                m_accelerationtype = TriZero;
            }
        } else {
            if(td > 2*tjd) {
                m_accelerationtype = ZeroTra;
            } else if(td > 0.0) {
                m_accelerationtype = ZeroTri;
            }else {
                m_accelerationtype = Zero;
            }
        }
    } else {
        if(Ta > 2*tja) {
            if(td > 2*tjd) {
                m_accelerationtype = TraTra;
            } else if(td > 0.0) {
                m_accelerationtype = TraTri;
            } else {
                m_accelerationtype = TraNone;
            }
        }else if(Ta > 0.0) {
            if(td > 2*tjd) {
                m_accelerationtype = TriTra;
            } else if(td > 0.0) {
                m_accelerationtype = TriTri;
            } else {
                m_accelerationtype = TriNone;
            }
        } else {
            if(td > 2*tjd) {
                m_accelerationtype = NoneTra;
            } else if(td > 0.0) {
                m_accelerationtype = NoneTri;
            } else {
                return -2;
            }
        }
    }
    return 1;
}

int VelocityPlanning::GetVelocityType(VectorXd para, int peroid)
{
    double ta, tv, td, tja, tjd;
    double l, vs, ve, vmax, vlim;
    double amax, amin, alima, alimd, jmax, jmin;
    ta = para(0);
    tv = para(1);
    td = para(2);
    tja = para(3);
    tjd = para(4);
    l = para(5);
    vs = para(6);
    ve = para(7);
    vmax = para(8);
    vlim = para(9);
    amax = para(10);
    amin = para(11);
    alima = para(12);
    alimd = para(13);
    jmax = para(14);
    jmin = para(15);
    // m_velocityType = 0;
    if (tv > EPSILON) { // zero
        if (ta > EPSILON) {
            if (jmax > EPSILON) {   // pos
                if (ta - 2 * tja > EPSILON) {   // tra
                    if (td > EPSILON) {
                        if (jmin > EPSILON) {   // pos
                            if (td - 2 * tjd > EPSILON) {   // tra
                                m_velocityType = PosTraZeroPosTra;
                            } else {    // tri
                                m_velocityType = PosTraZeroPosTri;
                            }
                        } else {    // neg
                            if (td - 2 * tjd > EPSILON) {   // tra
                                m_velocityType = PosTraZeroNegTra;
                            } else {    // tri
                                m_velocityType = PosTraZeroNegTri;
                            }
                        }
                    } else {
                        m_velocityType = PosTraZero;
                    }
                } else {    // tri
                    if (td > EPSILON) {
                        if (jmin > EPSILON) {   // pos
                            if (td - 2 * tjd > EPSILON) {   // tra
                                m_velocityType = PosTriZeroPosTra;
                            } else {    // tri
                                m_velocityType = PosTriZeroPosTri;
                            }
                        } else {    // neg
                            if (td - 2 * tjd > EPSILON) {   // tra
                                m_velocityType = PosTriZeroNegTra;
                            } else {    // tri
                                m_velocityType = PosTriZeroNegTri;
                            }
                        }
                    } else {
                        m_velocityType = PosTriZero;
                    }
                }
            } else if (jmax < -EPSILON) {   // neg
                if (ta - 2 * tja > EPSILON) {   // tra
                    if (td > EPSILON) {
                        if (jmin > EPSILON) {   // pos
                            if (td - 2 * tjd > EPSILON) {   // tra
                                m_velocityType = NegTraZeroPosTra;
                            } else {    // tri
                                m_velocityType = NegTraZeroPosTri;
                            }
                        } else {    // neg
                            if (td - 2 * tjd > EPSILON) {   // tra
                                m_velocityType = NegTraZeroNegTra;
                            } else {    // tri
                                m_velocityType = NegTraZeroNegTri;
                            }
                        }
                    } else {
                        m_velocityType = NegTraZero;
                    }
                } else {    // tri
                    if (td > EPSILON) {
                        if (jmin > EPSILON) {   // pos
                            if (td - 2 * tjd > EPSILON) {   // tra
                                m_velocityType = NegTriZeroPosTra;
                            } else {    // tri
                                m_velocityType = NegTriZeroPosTri;
                            }
                        } else {    // neg
                            if (td - 2 * tjd > EPSILON) {   // tra
                                m_velocityType = NegTriZeroNegTra;
                            } else {    // tri
                                m_velocityType = NegTriZeroNegTri;
                            }
                        }
                    } else {
                        m_velocityType = NegTriZero;
                    }
                }
            }
        } else {    // ta = 0.0
            if (td > EPSILON) {
                if (jmin > EPSILON) {   // pos
                    if (td - 2 * tjd > EPSILON) {   // tra
                        m_velocityType = ZeroPosTra;
                    } else {    // tri
                        m_velocityType = ZeroPosTri;
                    }
                } else {    // neg
                    if (td - 2 * tjd > EPSILON) {   // tra
                        m_velocityType = ZeroNegTra;
                    } else {    // tri
                        m_velocityType = ZeroNegTri;
                    }
                }
            } else {    // td = 0.0
                m_velocityType = Uni;
            }
        }
    } else {    // tv = 0.0
        if (ta > EPSILON) {
            if (jmax > EPSILON) {   // pos
                if (ta - 2 * tja > EPSILON) {   // tra
                    if (td > EPSILON) {
                        if (jmin > EPSILON) {   // pos
                            if (td - 2 * tjd > EPSILON) {   // tra
                                m_velocityType = PosTraPosTra;
                            } else {    // tri
                                m_velocityType = PosTraPosTri;
                            }
                        } else {    // neg
                            if (td - 2 * tjd > EPSILON) {   // tra
                                m_velocityType = PosTraNegTra;
                            } else {    // tri
                                m_velocityType = PosTraNegTri;
                            }
                        }
                    } else {
                        m_velocityType = PosTra;
                    }
                } else {    // tri
                    if (td > EPSILON) {
                        if (jmin > EPSILON) {   // pos
                            if (td - 2 * tjd > EPSILON) {   // tra
                                m_velocityType = PosTriPosTra;
                            } else {    // tri
                                m_velocityType = PosTriPosTri;
                            }
                        } else {    // neg
                            if (td - 2 * tjd > EPSILON) {   // tra
                                m_velocityType = PosTriNegTra;
                            } else {    // tri
                                m_velocityType = PosTriNegTri;
                            }
                        }
                    } else {
                        m_velocityType = PosTri;
                    }
                }
            } else if (jmax < -EPSILON) {   // neg
                if (ta - 2 * tja > EPSILON) {   // tra
                    if (td > EPSILON) {
                        if (jmin > EPSILON) {   // pos
                            if (td - 2 * tjd > EPSILON) {   // tra
                                m_velocityType = NegTraPosTra;
                            } else {    // tri
                                m_velocityType = NegTraPosTri;
                            }
                        } else {    // neg
                            if (td - 2 * tjd > EPSILON) {   // tra
                                m_velocityType = NegTraNegTra;
                            } else {    // tri
                                m_velocityType = NegTraNegTri;
                            }
                        }
                    } else {
                        m_velocityType = NegTra;
                    }
                } else {    // tri
                    if (td > EPSILON) {
                        if (jmin > EPSILON) {   // pos
                            if (td - 2 * tjd > EPSILON) {   // tra
                                m_velocityType = NegTriPosTra;
                            } else {    // tri
                                m_velocityType = NegTriPosTri;
                            }
                        } else {    // neg
                            if (td - 2 * tjd > EPSILON) {   // tra
                                m_velocityType = NegTriNegTra;
                            } else {    // tri
                                m_velocityType = NegTriNegTri;
                            }
                        }
                    } else {
                        m_velocityType = NegTri;
                    }
                }
            }
        } else {    // ta = 0.0
            if (td > EPSILON) {
                if (jmin > EPSILON) {   // pos
                    return -1;
                } else {    // neg
                    if (td - 2 * tjd > EPSILON) {   // tra
                        m_velocityType = NegTra;
                    } else {    // tri
                        m_velocityType = NegTri;
                    }
                }
            } else {    // td = 0.0
                return -1;
            }
        }
    }
    return 0;
}

int VelocityPlanning::CorrentParas(VectorXd &para, int peroid)
{
    double ta, tv, td, tja, tjd;
    double l, vs, ve, vmax, vlim;
    double amax, amin, alima, alimd, jmax, jmin;
    ta = para(0);
    tv = para(1);
    td = para(2);
    tja = para(3);
    tjd = para(4);
    l = para(5);
    vs = para(6);
    ve = para(7);
    vmax = para(8);
    vlim = para(9);
    amax = para(10);
    amin = para(11);
    alima = para(12);
    alimd = para(13);
    jmax = para(14);
    jmin = para(15);
    double te       = ceil((ta + tv + td) / (1.0 / peroid)) * (1.0 / peroid) - (ta + tv + td);
    double tja_new  = tja, tjd_new = tjd, ta_new = ta, td_new = td, tv_new = tv;
    double delta_tv = 0.0;
    double delta    = 1e-10;
    double tm       = te + ta + tv + td;
    if(m_velocityType == PosTraZeroPosTra || m_velocityType == PosTraZeroNegTra
        || m_velocityType == NegTraZeroPosTra || m_velocityType == NegTraZeroNegTra
        || m_velocityType == PosTriZeroPosTri || m_velocityType == PosTriZeroNegTri
        || m_velocityType == NegTriZeroPosTri || m_velocityType == NegTriZeroNegTri) {
        if((ve - vs) > delta) { //extend the acceleration time
            delta_tv = te * (vlim + vs) / (vs - vlim);
            tv_new  = tv + delta_tv;
            if(tv_new < -delta) {
                //decrease vlim
                tv_new  = tv;
                ta_new  = ta + te;
                tja_new = tja + 0.5 * te;
                vlim    = (2 * l - ta_new * vs - td_new * ve) / (2 * tm - ta_new - td_new);
                alima   = (vlim - vs) / (ta_new - tja_new);
                alimd   = (ve - vlim) / (td_new - tjd_new);
                jmax    = alima / tja_new;
                jmin    = alimd / tjd_new;
            } else {
                ta_new  = ta + te - delta_tv;
                tja_new = tja + 0.5 * (te - delta_tv);
                jmax    = (vlim - vs) / ((ta_new - tja_new) * tja_new);
                alima   = tja_new * jmax;
                vlim    = alima * (ta_new - tja_new) + vs;
            }
        } else if((ve - vs) < -delta) { //extend the deceleration time
            delta_tv = te * (vlim + ve) / (ve - vlim);
            tv_new  = tv + delta_tv;
            if(tv_new < -delta) {
                //decrease vlim
                tv_new  = tv;
                td_new  = td + te;
                tjd_new = tjd + 0.5 * te;
                vlim    = (2 * l - ta_new * vs - td_new * ve) / (2 * tm - ta_new - td_new);
                alima   = (vlim - vs) / (ta_new - tja_new);
                alimd   = (ve - vlim) / (td_new - tjd_new);
                jmax    = alima / tja_new;
                jmin    = alimd / tjd_new;
            } else {
                td_new  = td + te - delta_tv;
                tjd_new = tjd + 0.5 * (te - delta_tv);
                jmin    = (ve - vlim) / ((td_new - tjd_new) * tjd_new);
                alimd   = tjd_new * jmin;
                vlim    = ve - alimd*(td_new - tjd_new);
            }
        } else { //extend the acceleration and deceleration time
            delta_tv = te * (vs + 2 * vlim + ve) / (vs + ve - 2 * vlim);
            tv_new  = tv + delta_tv;
            if(tv_new < -delta) {
                //decrease vlim
                tv_new  = tv;
                ta_new  = ta + 0.5 * te;
                tja_new = tja + 0.25 * te;
                td_new  = td + 0.5 * te;
                tjd_new = tjd + 0.25 * te;
                vlim    = (2 * l - ta_new * vs - td_new * ve) / (2 * tm - ta_new - td_new);
                alima   = (vlim - vs) / (ta_new - tja_new);
                alimd   = (ve - vlim) / (td_new - tjd_new);
                jmax    = alima / tja_new;
                jmin    = alimd / tjd_new;
            } else {
                ta_new  = ta + 0.5 * (te - delta_tv);
                td_new  = td + 0.5 * (te - delta_tv);
                tja_new = tja + (te - delta_tv) * 0.25;
                tjd_new = tjd + (te - delta_tv) * 0.25;
                jmax    = (vlim - vs) / ((ta_new - tja_new) * tja_new);
                jmin    = (ve - vlim) / ((td_new - tjd_new) * tjd_new);
                alima   = tja_new * jmax;
                alimd   = tjd_new * jmin;
                vlim    = alima * (ta_new - tja_new) + vs;
            }
        }
    } else if(m_velocityType == PosTriZeroPosTra || m_velocityType == PosTriZeroNegTra
               || m_velocityType == NegTriZeroPosTra || m_velocityType == NegTriZeroNegTra
               || m_velocityType == ZeroNegTra || m_velocityType == ZeroNegTri
               || m_velocityType == ZeroPosTra || m_velocityType == ZeroPosTri) { //extend the deceleration time
        delta_tv = te * (vlim + ve) / (ve - vlim);
        tv_new  = tv + delta_tv;
        if(tv_new < -delta) {
            if(m_velocityType == PosTriZeroPosTra || m_velocityType == PosTriZeroNegTra
                || m_velocityType == NegTriZeroPosTra || m_velocityType == NegTriZeroNegTra) {
                //decrease vlim
                tv_new  = tv;
                td_new  = td + te;
                tjd_new = tjd + 0.5 * te;
                vlim    = (2 * l - ta_new * vs - td_new * ve) / (2 * tm - ta_new - td_new);
                alima   = (vlim - vs) / (ta_new - tja_new);
                alimd   = (ve - vlim) / (td_new - tjd_new);
                jmax    = alima / tja_new;
                jmin    = alimd / tjd_new;
            } else {
                //decrease vlim and vs
                tv_new  = tv;
                td_new  = td + te;
                tjd_new = tjd + 0.5 * te;
                vlim    = (2 * l - td_new * ve) / (2 * tv_new + td_new);
                alimd   = (ve - vlim) / (td_new - tjd_new);
                jmin    = alimd / tjd_new;
                vs = vlim;
            }
        } else { //extend the deceleration time
            td_new  = td + te - delta_tv;
            tjd_new = tjd + 0.5 * (te - delta_tv);
            jmin    = (ve - vlim) / ((td_new - tjd_new) * tjd_new);
            alimd   = tjd_new * jmin;
            vlim    = ve - alimd*(td_new - tjd_new);
        }
    } else if(m_velocityType == PosTraZeroPosTri || m_velocityType == PosTraZeroNegTri
               || m_velocityType == NegTraZeroPosTri || m_velocityType == NegTraZeroNegTri
               || m_velocityType == PosTraZero || m_velocityType == PosTriZero
               || m_velocityType == NegTraZero || m_velocityType == NegTriZero) { //extend the acceleration time
        delta_tv = te * (vlim + vs) / (vs - vlim);
        tv_new  = tv + delta_tv;
        if(tv_new < -delta) {
            if(m_velocityType == PosTraZeroPosTri || m_velocityType == PosTraZeroNegTri
                || m_velocityType == NegTraZeroPosTri || m_velocityType == NegTraZeroNegTri) {
                //decrease vlim
                tv_new  = tv;
                ta_new  = ta + te;
                tja_new = tja + 0.5 * te;
                vlim    = (2 * l - ta_new * vs - td_new * ve) / (2 * tm - ta_new - td_new);
                alima   = (vlim - vs) / (ta_new - tja_new);
                alimd   = (ve - vlim) / (td_new - tjd_new);
                jmax    = alima / tja_new;
                jmin    = alimd / tjd_new;
            } else {
                //decrease vlim and ve
                tv_new  = tv;
                ta_new  = ta + te;
                tja_new = tja + 0.5 * te;
                vlim    = (2 * l - ta_new * vs) / (2 * tv_new + ta_new);
                alima   = (vlim - vs) / (ta_new - tja_new);
                alimd   = (ve - vlim) / (td_new - tjd_new);
                jmax    = alima / tja_new;
                jmin    = alimd / tjd_new;
                ve = vlim;
            }
        } else {
            ta_new  = ta + te - delta_tv;
            tja_new = tja + 0.5 * (te - delta_tv);
            jmax    = (vlim - vs) / ((ta_new - tja_new) * tja_new);
            alima   = tja_new * jmax;
            vlim    = alima * (ta_new - tja_new) + vs;
        }
    } else if(m_velocityType == Uni){
        //decrease vs and ve
        td_new = cbrt(8*(tv * vlim - (tv + te)*vlim) / jmin);
        tv_new = tv - td_new;
        if(tv_new < -delta) {
            tv_new  = tv + te;
            vlim    = l / tv_new;
            vs     = vlim;
            ve     = vlim;
            td_new  = td;
        } else {
            tjd_new = 0.5 * td_new;
            alimd   = 0.5 * jmin * td_new;
            ve     = vlim + 0.25 * jmin * pow(td_new,2);
        }
    } else if(m_velocityType == PosTraPosTra || m_velocityType == PosTraNegTra
               || m_velocityType == NegTraPosTra || m_velocityType == NegTraNegTra
               || m_velocityType == PosTriPosTri || m_velocityType == PosTriNegTri
               || m_velocityType == NegTriPosTri || m_velocityType == NegTriNegTri) {
        // extend the acceleration and deceleration time
        delta_tv = 0.0;
        tv_new  = tv + delta_tv;
        ta_new  = ta + 0.5 * (te - delta_tv);
        td_new  = td + 0.5 * (te - delta_tv);
        tja_new = tja + (te - delta_tv) * 0.25;
        tjd_new = tjd + (te - delta_tv) * 0.25;
        jmax    = (vlim - vs) / ((ta_new - tja_new) * tja_new);
        jmin    = (ve - vlim) / ((td_new - tjd_new) * tjd_new);
        alima   = tja_new * jmax;
        alimd   = tjd_new * jmin;
        vlim    = alima * (ta_new - tja_new) + vs;
    } else if(m_velocityType == PosTraPosTri || m_velocityType == PosTraNegTri
               || m_velocityType == NegTraPosTri || m_velocityType == NegTraNegTri) {
        delta_tv = 0.0;
        tv_new  = tv + delta_tv;
        ta_new  = ta + te - delta_tv;
        tja_new = tja + 0.5 * (te - delta_tv);
        jmax    = (vlim - vs) / ((ta_new - tja_new) * tja_new);
        alima   = tja_new * jmax;
        vlim    = alima * (ta_new - tja_new) + vs;
    } else if(m_velocityType == PosTri || m_velocityType == PosTra) {
        //decrease ve
        tv_new  = tv;
        ta_new  = ta + te;
        tja_new = tja + 0.5 * te;
        vlim    = 2 * l / ta_new - vs;
        alima   = (vlim - vs) / (ta_new - tja_new);
        jmax    = alima / tja_new;
        ve     = vlim;
    } else if(m_velocityType == PosTriPosTra || m_velocityType == PosTriNegTra
               || m_velocityType == NegTriPosTra || m_velocityType == NegTriNegTra) {
        delta_tv = 0.0;
        tv_new  = tv + delta_tv;
        td_new  = td + te - delta_tv;
        tjd_new = tjd + 0.5 * (te - delta_tv);
        jmin    = (ve - vlim) / ((td_new - tjd_new) * tjd_new);
        alimd   = tjd_new * jmin;
        vlim    = ve - alimd * (td_new - tjd_new);
    } else if(m_velocityType == NegTra || m_velocityType == NegTri) {
        //decrease vs
        tv_new  = tv;
        td_new  = td + te;
        tjd_new = tjd + 0.5 * te;
        vlim    = 2 * l / td_new - ve;
        alimd   = (ve - vlim) / (td_new - tjd_new);
        jmin    = alimd / tjd_new;
        vs     = vlim;
    }
    if(tv_new < -1e-5) {
        qDebug("tv error!!!");
        return -1;
    }
    // para = VectorXd::Zero(16);
    para << ta_new, tv_new, td_new, tja_new, tjd_new, l, vs, ve, vmax, vlim, amax, amin, alima, alimd, jmax, jmin;
    return 0;
}

double VelocityPlanning::GetPosition(double t, VectorXd para)
{
    double ta, tv, td, tja, tjd;
    double l, vs, ve, vmax, vlim;
    double amax, amin, alima, alimd, jmax, jmin;
    ta = para(0);
    tv = para(1);
    td = para(2);
    tja = para(3);
    tjd = para(4);
    l = para(5);
    vs = para(6);
    ve = para(7);
    vmax = para(8);
    vlim = para(9);
    amax = para(10);
    amin = para(11);
    alima = para(12);
    alimd = para(13);
    jmax = para(14);
    jmin = para(15);
    double tsum = ta + tv + td;
    if (t >= 0.0 && t < tja) {
        return vs * t + jmax * pow(t, 3) / 6.0;
        // v = vs + jmax * pow(t, 2) * 0.5;
        // a = jmax * t;
        // j = jmax;
    } else if (t >= tja && t < ta - tja) {
        return vs * t + alima / 6.0 * (3 * pow(t, 2) - 3 * tja * t + pow(tja, 2));
        // v = vs + alima * (t - tja * 0.5);
        // a = alima;
        // j = 0.0;
    } else if (t >= ta - tja && t < ta) {
        return (vlim + vs) * (ta / 2.0) - vlim * (ta - t) + jmax * pow(ta - t, 3) / 6.0;
        // v = vlim - jmax * pow(ta - t, 2) * 0.5;
        // a = - jmax * (t -ta);
        // j = jmax;
    } else if (t >= ta && t < ta + tv) {
        return (vlim + vs) * ta * 0.5 + vlim * (t - ta);
        // v = vlim;
        // a = 0.0;
        // j = 0.0;
    } else if (t >= ta + tv && t < tsum - td + tjd) {
        return l - (vlim + ve) * td * 0.5 + vlim * (t - tsum + td) + jmin * pow(t - tsum + td, 3) / 6.0;
        // v = vlim + jmin * pow(t - tsum + td, 2) * 0.5;
        // a = jmin * (t - tsum + td);
        // j = jmin;
    } else if (t >= tsum - td + tjd && t < tsum - tjd) {
        return l - (vlim + ve) * td * 0.5 + vlim * (t - tsum + td) + alimd / 6.0 * (3 * pow(t - tsum + td, 2) - 3 * tjd * (t - tsum + td) + pow(tjd, 2));
        // v = vlim + alimd * (t - tsum + td - tjd * 0.5);
        // a = alimd;
        // j = 0.0;
    } else {
        return l - ve * (tsum - t) + jmin * (pow(tsum - t, 3) / 6.0);
        // v = ve - jmin * pow(t - tsum, 2) * 0.5;
        // a = -jmin * (t - tsum);
        // j = jmin;
    }
}

double VelocityPlanning::GetVelocity(double t, VectorXd para)
{
    double ta, tv, td, tja, tjd;
    double l, vs, ve, vmax, vlim;
    double amax, amin, alima, alimd, jmax, jmin;
    ta = para(0);
    tv = para(1);
    td = para(2);
    tja = para(3);
    tjd = para(4);
    l = para(5);
    vs = para(6);
    ve = para(7);
    vmax = para(8);
    vlim = para(9);
    amax = para(10);
    amin = para(11);
    alima = para(12);
    alimd = para(13);
    jmax = para(14);
    jmin = para(15);
    double tsum = ta + tv + td;
    if (t >= 0.0 && t < tja) {
        // return vs * t + jmax * pow(t, 3) / 6.0;
        return vs + jmax * pow(t, 2) * 0.5;
        // a = jmax * t;
        // j = jmax;
    } else if (t >= tja && t < ta - tja) {
        // return vs * t + alima / 6.0 * (3 * pow(t, 2) - 3 * tja * t + pow(tja, 2));
        return vs + alima * (t - tja * 0.5);
        // a = alima;
        // j = 0.0;
    } else if (t >= ta - tja && t < ta) {
        // return (vlim + vs) * (ta / 2.0) - vlim * (ta - t) + jmax * pow(ta - t, 3) / 6.0;
        return vlim - jmax * pow(ta - t, 2) * 0.5;
        // a = - jmax * (t -ta);
        // j = jmax;
    } else if (t >= ta && t < ta + tv) {
        // return (vlim + vs) * ta * 0.5 + vlim * (t - ta);
        return vlim;
        // a = 0.0;
        // j = 0.0;
    } else if (t >= ta + tv && t < tsum - td + tjd) {
        // return l - (vlim + ve) * td * 0.5 + vlim * (t - tsum + td) + jmin * pow(t - tsum + td, 3) / 6.0;
        return vlim + jmin * pow(t - tsum + td, 2) * 0.5;
        // a = jmin * (t - tsum + td);
        // j = jmin;
    } else if (t >= tsum - td + tjd && t < tsum - tjd) {
        // return l - (vlim + ve) * td * 0.5 + vlim * (t - tsum + td) + alimd / 6.0 * (3 * pow(t - tsum + td, 2) - 3 * tjd * (t - tsum + td) + pow(tjd, 2));
        return vlim + alimd * (t - tsum + td - tjd * 0.5);
        // a = alimd;
        // j = 0.0;
    } else {
        // return l - ve * (tsum - t) + jmin * (pow(tsum - t, 3) / 6.0);
        return ve - jmin * pow(t - tsum, 2) * 0.5;
        // a = -jmin * (t - tsum);
        // j = jmin;
    }
}

int VelocityPlanning::Solve2thEquation(vector<double> coeff, vector<double> &result)
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

int VelocityPlanning::Solve3thEquation(vector<double> coeff, vector<double> &result)
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

int VelocityPlanning::Solve4thEquation(vector<double> coeff, vector<double> &result)
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

int VelocityPlanning::Solve5thEquation(vector<double> coeff, vector<double> &result)
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
