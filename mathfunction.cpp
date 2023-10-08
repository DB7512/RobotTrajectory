#include "mathfunction.h"
#include <math.h>
#include <eigen-3.4.0/Eigen/Geometry>
#include <random>
#include <unistd.h>
#include <iostream>
#include <QTime>
#include <QCoreApplication>

MathFunction::MathFunction(QObject *parent)
    : QObject{parent}
{

}

void MathFunction::Equation_Root2(float coffe[3], float Roots_result[2], float Detla)
{
    float D = sqrt(Detla);
    float A = 2 * coffe[0];
    Roots_result[0] = (-coffe[1] + D) / A;
    Roots_result[1] = (-coffe[1] - D) / A;
}

void MathFunction::Equation_Root3(float coffe[4], float Roots_result[3])
{
    float a, b, c, d;
    float A, B, C;
    float delta;

    a = coffe[0];
    b = coffe[1];
    c = coffe[2];
    d = coffe[3];
    A = pow(b, 2) - 3 * a * c;
    B = b * c - 9 * a * d;
    C = pow(c, 2) - 3 * b * d;
    delta = pow(B, 2) - 4 * A * C;
    if (fabs(A) < 1e-6 && fabs(B) < 1e-6) {
        Roots_result[0] = -c / b;
        Roots_result[1] = Roots_result[0];
        Roots_result[2] = Roots_result[0];
        return;
    }
    if (delta > 0) {
        float Y1, Y2;
        float unreal;
        float Y_1, Y_2;
        int zero_1 = 0, zero_2 = 0;
        Y1 = A * b + 3 * a * (-B + sqrt(delta)) / 2;
        Y2 = A * b + 3 * a * (-B - sqrt(delta)) / 2;
        if (Y1 < 0) {
            Y1 = -Y1;
            zero_1 = 1;
        }
        if (Y2 < 0) {
            Y2 = -Y2;
            zero_2 = 1;
        }
        Y_1 = pow(Y1, 1.0 / 3);
        Y_2 = pow(Y2, 1.0 / 3);
        if (zero_1) {
            Y_1 = -Y_1;
        }
        if (zero_2) {
            Y_2 = -Y_2;
        }
        Roots_result[0] = (-b - (Y_1 + Y_2)) / (3 * a);
        unreal = (sqrt(3) / 2 * (Y_1 - Y_2)) / (3 * a);
        if (fabs(unreal) < 1e-6) {
            Roots_result[1] = (-b + 0.5 * (Y_1 + Y_2)) / (3 * a);
            Roots_result[2] = Roots_result[1];
        } else {
            Roots_result[1] = 0;
            Roots_result[2] = 0;
        }
        return;
    }
    if (fabs(delta) < 1e-6 && fabs(A) > 1e-6) {
        float K = B / A;
        Roots_result[0] = -b / a + K;
        Roots_result[1] = -0.5 * K;
        Roots_result[2] = Roots_result[1];
        return;
    }
    if (delta < 0 && A > 0) {
        float T, theta;
        T = (2 * A * b - 3 * a * B) / (2 * sqrt(pow(A, 3)));
        theta = acos(T);
        Roots_result[0] = (-b - 2 * sqrt(A) * cos(theta / 3)) / (3 * a);
        Roots_result[1] = (-b + sqrt(A) * (cos(theta / 3) + sqrt(3) * sin(theta / 3))) / (3 * a);
        Roots_result[2] = (-b + sqrt(A) * (cos(theta / 3) - sqrt(3) * sin(theta / 3))) / (3 * a);
        return;
    }
}

void MathFunction::Equation_Root4(float coffe[5], float Roots_result[4][2])
{
    Matrix<float, 5, 5>matrix_coffe;
    Matrix<complex<float>, Dynamic, Dynamic>matrix_sigenvalues;
    VectorXf rl = VectorXf(5);
    VectorXf ig = VectorXf(5);
    matrix_coffe << -coffe[1] / coffe[0], -coffe[2] / coffe[0], -coffe[3] / coffe[0], -coffe[4] / coffe[0], -coffe[5] / coffe[0],
        1, 0, 0, 0, 0,
        0, 1, 0, 0, 0,
        0, 0, 1, 0, 0,
        0, 0, 0, 1, 0;

    matrix_sigenvalues = matrix_coffe.eigenvalues();
    rl << matrix_sigenvalues.real();
    ig << matrix_sigenvalues.imag();

    for (int i = 0; i < 5; i++) {
        if (ig(i) < 0.000001 && ig(i) > -0.000001) {
            Roots_result[i][0] = rl(i);
            Roots_result[i][1] = 0;
        }
        else {
            Roots_result[i][0] = 0;
            Roots_result[i][1] = 0;
        }
    }
}

void MathFunction::Equation_Root5(float coffe[6], float Roots_result[5][2])
{
    Eigen::Matrix<float, 5, 5>matrix_coffe;
    Eigen::Matrix<complex<float>, Eigen::Dynamic, Eigen::Dynamic>matrix_sigenvalues;
    Eigen::VectorXf rl = VectorXf(5);
    Eigen::VectorXf ig = VectorXf(5);
    matrix_coffe << -coffe[1] / coffe[0], -coffe[2] / coffe[0], -coffe[3] / coffe[0], -coffe[4] / coffe[0], -coffe[5] / coffe[0],
        1, 0, 0, 0, 0,
        0, 1, 0, 0, 0,
        0, 0, 1, 0, 0,
        0, 0, 0, 1, 0;

    matrix_sigenvalues = matrix_coffe.eigenvalues();
    rl << matrix_sigenvalues.real();
    ig << matrix_sigenvalues.imag();

    for (int i = 0; i < 5; i++) {
        if (ig(i) < 0.000001 && ig(i) > -0.000001) {
            Roots_result[i][0] = rl(i);
            Roots_result[i][1] = 0;
        }
        else {
            Roots_result[i][0] = 0;
            Roots_result[i][1] = 0;
        }
    }
}

void MathFunction::SwapValue(float &val1, float &val2)
{
    float temp = 0.0;
    temp = val1;
    val1 = val2;
    val2 = temp;
}

double MathFunction::GetRand(double min, double max)
{
    default_random_engine e;
    //uniform_real_distribution:产生均匀分布的实数
    uniform_real_distribution<double> u(min,max);   //左闭右闭区间
    e.seed(time(0));
    DelayMS(5);
    return u(e);
}

void MathFunction::DelayMS(unsigned int msec)
{
    QTime _Timer = QTime::currentTime().addMSecs(msec);
    while(QTime::currentTime() < _Timer)
        QCoreApplication::processEvents(QEventLoop::AllEvents,100);
}
