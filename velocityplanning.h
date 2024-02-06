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
    Tri0 = 0,
    Tri,
    Tra,
}AccType;

typedef enum {
    Tri0 = 0,
    Tri,
    Tra,
}DecType;

class VelocityPlanning : public TrajectoryPlanning
{

public:
    int velocityPlan(double s, double vs, double ve, double vmax, double amax, double jmax, VectorXd &para, double time);
    int velocityPlan(double s, double vs, double ve, double vmax, double amax, double jmax, VectorXd &para);
    int TrajectoryTime(double &time, float Q, float v_0, float v_1, float vmax, float amax, float jmax, VectorXd &para);
    void signAmaxJmax(double vf, double vt, double &amax, double &jmax, double &amin, double &jmin);

    // 求解方程
    int solve2thEquation(vector<double> coeff, vector<double> &result);
    int solve3thEquation(vector<double> coeff, vector<double> &result);
    int solve4thEquation(vector<double> coeff, vector<double> &result);
    int solve5thEquation(vector<double> coeff, vector<double> &result);

    double m_vmax;
    double m_amax;
    double m_jmax;

signals:

};

#endif // VELOCITYPLANNING_H
