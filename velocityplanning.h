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
    int SetVelocityPlan(double l, double &vs, double &ve, double &vmax, double amax, double jmax, VectorXd &para, double &time);
    int TimePLan(double l, double &vs, double &ve, double &vmax, double amax, double jmax, VectorXd &para, double time);
    int TrajectoryTime(double &time, float Q, float v_0, float v_1, float vmax, float amax, float jmax, VectorXd &para);
    int GetVelocityType(VectorXd para, int peroid);
    int CorrentParas(VectorXd &para, int peroid);
    double GetPosition(double t, VectorXd para);
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
