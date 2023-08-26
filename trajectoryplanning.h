#ifndef TRAJECTORYPLANNING_H
#define TRAJECTORYPLANNING_H

#include <QObject>
#include <eigen-3.4.0/Eigen/Dense>

using namespace std;
using namespace Eigen;

struct Parameters {
    float Ta;
    float Tv;
    float Td;
    float Tj1;
    float Tj2;
    float q0;
    float q1;
    float v0;
    float v1;
    float vlim;
    float amax;
    float amin;
    float alima;
    float alimd;
    float jmax;
    float jmin;
    Parameters() {
        memset(this,0,sizeof(Parameters));
    }
};
typedef enum {
    TraUniTra = 0,
    TriUniTra,
    TraUniTri,
    TriUniTri,
    TraUni,
    TriUni,
    UniTra,
    UniTri,
    TraTra,
    TraTri,
    TriTra,
    TriTri,
    TraNone,
    TriNone,
    NoneTra,
    NoneTri,
    Uni,
}TrajectoryType;

typedef struct _PointInformation {
    float pose[6];
    float v;
    float a;
    float j;
    float r;
    float vmax;
    float amax;
    float jmax;
    _PointInformation() {
        memset(this,0,sizeof(_PointInformation));
    }
}PointInformation;
typedef struct _TrajectoryInformation {
    float period;
    float vmax;
    float amax;
    float jmax;
    float v_0;
    float v_1;
    VectorXf para;
    TrajectoryType trajectorytype;
    vector<vector<float> > TrajectoryInterpolation;
    vector<vector<float> > TrajectoryVelocity;
    vector<vector<float> > TrajectoryAcceleration;
    _TrajectoryInformation() {
        memset(this,0,sizeof(_TrajectoryInformation));
        para = VectorXf::Zero(15);
    }
}TrajectoryInformation;

class TrajectoryPlanning : public QObject
{
    Q_OBJECT
public:
    explicit TrajectoryPlanning(QObject *parent = nullptr);
    /**
     * @brief LinePlanning      use end point information as vmax, amax and jmax
     * @param peroid
     * @param point_start
     * @param point_end
     * @param trajectory_point
     * @param trajectory_inf
     * @return
     */
    bool LinePlanning(int peroid, PointInformation point_start, PointInformation point_end,
                      std::vector<std::vector<float> > &trajectory_point, vector<vector<float> > &trajectory_inf);
    bool TimePlanning(int peroid, float Q, float vmax, float amax, float jmax, float v_0, float v_1,
                      std::vector<float> &trajectory_point, vector<vector<float> > &trajectory_inf);
    int TrajectoryTime1(float &time, float Q, float v_0, float v_1, float vmax, float amax, float jmax, VectorXf &para);
    int TrajectoryTime2(float &time, float Q, float v_0, float v_1, float vmax, float amax, float jmax, VectorXf &para);
    void CorrentionParameters(VectorXf& para, int peroid);
    void TrajectoryCalculation(float t, VectorXf para, float q_dq[3]);

signals:

public:
    TrajectoryType m_trajectorytype;
    VectorXf para = VectorXf::Zero(15);
    PointInformation m_posestart;
    PointInformation m_poseend;
    int number;

};

static TrajectoryPlanning& GetTrajectoryPlanningInstance()
{
    static TrajectoryPlanning TrajectoryPlanningInstance;
    return TrajectoryPlanningInstance;
}

#endif // TRAJECTORYPLANNING_H
