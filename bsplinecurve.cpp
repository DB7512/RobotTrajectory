#include "bsplinecurve.h"
#include <iostream>
#include <cmath>

// 参数说明：
// 控制点有{P0,...,Pn}，即有n+1个控制点
// p次B样条曲线，p+1阶
// k阶导数
// 节点矢量 U = {u0,u1,...,um}，则节点有m+1个，m = n+p+1
// i = 0,1,...,n
// 重复度为 p+1

BSplineCurve::BSplineCurve(vector<Vector3d> points, int p)
{
    // 控制点（一共n+1个控制点）
    m_controlPoints = points;
    // 节点数 m = n + p + 1
    int m = m_controlPoints.size() - 1 + p + 1;
    // 次数
    m_p = p;

    // // 相邻控制点距离
    // vector<double> l;
    // // 计算相邻控制点距离
    // for (int i = 0; i < m_controlPoints.size() - 1; ++i) {
    //     l.push_back(sqrt(pow(m_controlPoints[i+1][0] - m_controlPoints[i][0],2)
    //                      + pow(m_controlPoints[i+1][1] - m_controlPoints[i][1],2)
    //                      + pow(m_controlPoints[i+1][2] - m_controlPoints[i][2],2)));
    // }
    // // 相邻控制点距离之和
    // double sum_l = 0.0;
    // for (int i = 0; i < l.size(); ++i) {
    //     sum_l += l[i];
    // }

    // 准均匀有理B样条
    double delta = 1.0 / (double)(m - 2*p);
    // 节点向量，p+1重复度
    m_U.resize(m+1);
    for (int i = 0; i < m + 1; ++i) {
        if(i <= p) {
            m_U[i] = 0.0;
        } else if(i >= m - p) {
            m_U[i] = 1.0;
        } else {
            m_U[i] = m_U[i - 1] + delta;
        }
    }
}

/**
 * @brief BSplineCurve::basisFunction B样条基函数
 * @param i
 * @param p
 * @param u
 * @return
 */
double BSplineCurve::basisFunction(int i, int p, double u, vector<double> U)
{
    if (p == 0)
        return (U[i] <= u && u < U[i+1]) ? 1.0 : 0.0;

    double firstPart = 0.0, secondPart = 0.0;

    if (fabs(U[i+p] - U[i]) > 1e-6)
        firstPart = ((u - U[i]) / (U[i+p] - U[i])) * basisFunction(i, p-1, u, U);
    if (fabs(U[i+p+1] - U[i + 1]) > 1e-6)
        secondPart = ((U[i+p+1] - u) / (U[i+p+1] - U[i+1])) * basisFunction(i+1, p-1, u, U);

    return firstPart + secondPart;
}

/**
 * @brief BSplineCurve::controlPointDerivative  控制点求导Pi(k)
 * @param i
 * @param p     次数
 * @param k     k次导数
 * @param U     节点矢量
 * @return
 */
VectorXd BSplineCurve::controlPointDerivative(int i, int p, int k, vector<double> U)
{
    if (k == 0) {
        return m_controlPoints[i];
    } else {
        double a = (p - k + 1) / (U[i+p+1] - U[i+k]);
        return a * (controlPointDerivative(i+1, p, k-1, U) - controlPointDerivative(i, p, k-1, U));
    }
}

/**
 * @brief BSplineCurve::curveDerivative 曲线在u处k阶的导数
 * @param p     p次曲线
 * @param u
 * @param k     k阶导数
 * @param U     节点矢量
 * @return
 */
VectorXd BSplineCurve::curveDerivative(int p, int u, int k, vector<double> U)
{
    VectorXd sum = Vector3d::Zero(3);
    vector<double> tempU;
    for (int i = 0; i < (int)(U.size() - 2 * k); ++i) {
        tempU.push_back(U[i + k]);
    }
    for (int i = 0; i <= (int)(m_controlPoints.size() - 1 - k); ++i) {
        sum += basisFunction(i, p - k, u, tempU) * controlPointDerivative(i, p, k, tempU);
    }
    return sum;
}

/**
 * @brief BSplineCurve::adaptiveSimpson 辛普森积分计算在[a,b]上的弧长
 * @param p     次数
 * @param a     计算区间
 * @param b     计算区间
 * @return      弧长
 */
double BSplineCurve::Simpson(int p, double a, double b, vector<double> U)
{
    double mid = (a + b) / 2.0;
    double h = (b - a) / 2.0;
    double fa = curveDerivative(p, a, 1, U).norm();
    double fmid = curveDerivative(p, mid, 1, U).norm();
    double fb = curveDerivative(p, b, 1, U).norm();

    return h / 3.0 * (fa + 4 * fmid + fb);
}

/**
 * @brief BSplineCurve::calculateLength 自适应辛普森积分计算[a,b]上的弧长
 * @param p     次数
 * @param a     计算区间
 * @param b     计算区间
 * @param eps   误差
 * @return
 */
double BSplineCurve::adaptiveSimpson(int p, double a, double b, double eps)
{
    double mid = (a + b) / 2.0;
    double l = Simpson(p, a, b, m_U);
    double l1 = Simpson(p, a, mid, m_U);
    double l2 = Simpson(p, mid, b, m_U);
    if(fabs(l1 + l2 -l) < 15 * 1e-6) {
        return l1 + l2 + (l1 + l2 - l) / 15.0;
    } else {
        return adaptiveSimpson(p, a, mid, eps/2.0) + adaptiveSimpson(p, mid, b, eps/2.0);
    }
}

/**
 * @brief BSplineCurve::calculateCurveLength    计算曲线弧长S = S1+S2+...+Sn
 * @param p
 * @return
 */
double BSplineCurve::calculateCurveLength(double eps)
{
    double length = 0.0;
    for (int i = 0; i < (int)(m_U.size() - 1); ++i) {
        length += adaptiveSimpson(m_p, m_U[i], m_U[i+1], eps);
    }
    return length;
}

double BSplineCurve::calculateLength(double a, double b, double eps)
{
    double length = 0.0;
    for (int i = 0; i < (int)(m_U.size() - 1); ++i) {
        length = adaptiveSimpson(m_p, a, b, eps);
    }
    return length;
}

double BSplineCurve::taylorExpand(double u, int p, double v, double t)
{
    return u + v * t / curveDerivative(p, u, 1, m_U).norm();
}

/**
 * @brief BSplineCurve::calculateBSplinePoint   计算B样条上的插值点
 * @param u     自变量
 * @param p     次数
 * @param point 插值点
 */
void BSplineCurve::calculateBSplinePoint(double u, int p, Vector3d &point)
{
    int n = m_controlPoints.size() - 1;           //控制点数量 = n + 1
    point = {0.0, 0.0, 0.0};
    double basis;
    for (int i = 0; i <= n; ++i) {
        basis = basisFunction(i, p, u, m_U);
        for (int j = 0; j < 3; ++j) {
            point[j] += m_controlPoints[i][j] * basis;
        }
    }
}
