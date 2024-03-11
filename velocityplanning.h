#ifndef VELOCITYPLANNING_H
#define VELOCITYPLANNING_H

#include <QObject>
#include "trajectoryplanning.h"

/**
 * @brief The VelocityPlanning class
 * 约定：para有16（0～15）个参数
 * para = {ta,tv,td,tja,tjd,q,vs,ve,vmax,vlim,amax,amin,alima,alimd,jmax,jmin}
 *        { 0, 1, 2,  3,  4,5, 6, 7,   8,   9,  10,  11,  12,   13,   14,  15}
 */

typedef enum {
    // v<vs v<ve
    NegTraZeroPosTra,
    NegTriZeroPosTra,
    NegTraZeroPosTri,
    NegTriZeroPosTri,
    NegTraPosTra,
    NegTraPosTri,
    NegTriPosTra,
    NegTriPosTri,
    // v-{min(vs,ve),max(vs,ve)}
    PosTraZeroPosTra,
    NegTraZeroNegTra,
    PosTriZeroPosTra,
    NegTriZeroNegTra,
    PosTraZeroPosTri,
    NegTraZeroNegTri,
    PosTriZeroPosTri,
    NegTriZeroNegTri,
    PosTraPosTra,
    NegTraNegTra,
    PosTraPosTri,
    NegTriNegTra,
    PosTriPosTri,
    NegTriNegTri,
    NegTraNegTri,
    PosTriPosTra,
    // v>vs,v>ve
    PosTraZeroNegTra,
    PosTriZeroNegTra,
    PosTraZeroNegTri,
    PosTriZeroNegTri,
    PosTraNegTra,
    PosTriNegTri,
    PosTraNegTri,
    PosTriNegTra,
    // vs-ve
    // vs<ve
    PosTraZero,
    PosTriZero,
    ZeroPosTra,
    ZeroPosTri,
    PosTra,
    PosTri,
    // vs>ve
    ZeroNegTra,
    ZeroNegTri,
    NegTraZero,
    NegTriZero,
    NegTra,
    NegTri,
    Uni,
}VelocityType;

class VelocityPlanning : public TrajectoryPlanning
{

public:
    // 速度规划，计算时间，处理三种类型
    int SetVelocityPlan(double l, double &vs, double &ve,
                        double &vmax, double amax, double jmax,
                        VectorXd &para, double &time);
    // xx给定time和l，进行速度规划，算法不完善，目前没有找到可以同时满足l和time的方法
    int TimePLan(double l, double &vs, double &ve, double &vmax, double amax, double jmax, VectorXd &para, double time);
    // xx时间圆整规划，仅可以处理vm > max{vs,ve}情况
    int TrajectoryTime(double &time, float Q, float v_0, float v_1, float vmax, float amax, float jmax, VectorXd &para);
    // 速度规划类型判断
    int GetVelocityType(VectorXd para);
    // 时间周期化之后修正速度规划参数
    int CorrentParas(VectorXd &para, int peroid);
    // 获取t时刻的位置
    double GetPosition(double t, VectorXd para);
    // 获取t时刻的速度
    double GetVelocity(double t, VectorXd para);


    // 求解方程
    int Solve2thEquation(vector<double> coeff, vector<double> &result);
    int Solve3thEquation(vector<double> coeff, vector<double> &result);
    int Solve4thEquation(vector<double> coeff, vector<double> &result);
    int Solve5thEquation(vector<double> coeff, vector<double> &result);

    double m_vmax;
    double m_amax;
    double m_jmax;
    VectorXd m_para;
    VelocityType m_velocityType;

signals:

};

static VelocityPlanning& GetVelPlanInstance()
{
    static VelocityPlanning VelocityPlanningInstance;
    return VelocityPlanningInstance;
}
#endif // VELOCITYPLANNING_H
