#ifndef MATHFUNCTION_H
#define MATHFUNCTION_H

#include <QObject>
#include <eigen-3.4.0/Eigen/Dense>

using namespace std;
using namespace Eigen;
class MathFunction : public QObject
{
    Q_OBJECT
public:
    explicit MathFunction(QObject *parent = nullptr);
    void Equation_Root2(float coffe[3], float Roots_result[2], float Detla);
    void Equation_Root3(float coffe[4], float Roots_result[3]);
    void Equation_Root4(float coffe[5], float Roots_result[4][2]);
    void Equation_Root5(float coffe[6], float Roots_result[5][2]);
    void SwapValue(float &val1, float &val2);

    double GetRand(double min, double max);
    void DelayMS(unsigned int msec);
signals:

};
static MathFunction& GetMathInstance()
{
    static MathFunction MathInstance;
    return MathInstance;
}

#endif // MATHFUNCTION_H
