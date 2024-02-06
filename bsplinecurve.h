#ifndef BSPLINECURVE_H
#define BSPLINECURVE_H

#include <QObject>
#include <iostream>
#include <vector>
#include <cmath>
#include <eigen-3.4.0/Eigen/Dense>

using namespace std;
using namespace Eigen;

class BSplineCurve : public QObject
{
    Q_OBJECT
public:
    BSplineCurve(vector<Vector3d> points, int p);

    double basisFunction(int i, int p, double u, vector<double> U);
    VectorXd controlPointDerivative(int i, int p, int k, vector<double> U);
    VectorXd curveDerivative(int p, int u, int k, vector<double> U);
    double Simpson(int p, double a, double b, vector<double> U);
    double adaptiveSimpson(int p, double a, double b, double eps);
    double calculateCurveLength(double eps);
    double calculateLength(double a, double b, double eps);
    double taylorExpand(double u, int p, double v, double t);
    void calculateBSplinePoint(double u, int p, Vector3d &point);
signals:

public:
    vector<Vector3d> m_controlPoints;    //控制点
    vector<double> m_U;                     //节点向量
    int m_p;                                    //次数（阶数-1）

};

#endif // BSPLINECURVE_H
