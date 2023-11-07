#ifndef VELOCITYPLANNING_H
#define VELOCITYPLANNING_H

#include <QObject>
#include "trajectoryplanning.h"

class VelocityPlanning : public TrajectoryPlanning
{

public:
    int TrajectoryTime(double &time, float Q, float v_0, float v_1, float vmax, float amax, float jmax, VectorXd &para);



signals:

};

#endif // VELOCITYPLANNING_H
