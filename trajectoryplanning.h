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
    float jmin; //No.16
    Parameters() {
        memset(this,0,sizeof(Parameters));
    }
};
typedef enum {
    TraZeroTra = 0,
    TriZeroTra,
    TraZeroTri,
    TriZeroTri,
    TraZero,
    TriZero,
    ZeroTra,
    ZeroTri,
    TraTra,
    TraTri,
    TriTra,
    TriTri,
    TraNone,
    TriNone,
    NoneTra,
    NoneTri,
    Zero,
}AccelerationType;

typedef enum {
    Position = 0,
    Posture,
}TrajectoryType;

typedef enum {
    None = 0,
    Line,       //直线
    Arc,        //过渡圆弧
    Circular,   //圆弧
    Line2arc,
}PathType;

typedef struct _PointInformation {
    Vector3d point; //位置
    Vector4d pose; //姿态四元数
    double velocity; //速度
    double acceleration; //加速度
    double jerk;
    double maxDeviation; //最大误差
    double maxVelocity;
    double maxAcceleration;
    double maxJerk;
    double radius;
    PathType pathType; //该路点所属路径类型
    _PointInformation() {
        memset(this,0,sizeof(_PointInformation));
        maxJerk = 80;
        pathType = None;
    }
}PointInformation;

typedef struct _PathInformation {
    PointInformation startpoint;
    PointInformation endpoint;
    //Vector3d startPoint;
    //Vector3d endPoint;
    //Vector4d startPose;
    //Vector4d endPose;
    double maxVelocity;
    double maxAcceleration;
    double maxJerk;
    //double startVelocity;
    //double endVelocity;
    PathType pathType; //路径类型
    double radius; //半径
    double theta; //圆心角
    double displacement; //位移
    double arclength; //弧长
    Vector3d center; //圆心
    VectorXd constraints; //16个参数的约束
    AccelerationType accelerationType; //加速度类型
    TrajectoryType trajectoryType; //姿态还是位移
//    vector<vector<double> > TrajectoryInterpolation;
//    vector<vector<double> > TrajectoryVelocity;
//    vector<vector<double> > TrajectoryAcceleration;
//    _TrajectoryInformation() {
//        memset(this,0,sizeof(_TrajectoryInformation));
//        constraints = VectorXd::Zero(16);
//        maxJerk = 80;
//        pathType = None;
//    }
}PathInformation;

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
    int TrajectoryTime3(float &time, float Q, float v_0, float v_1, float vmax, float amax, float jmax, VectorXf &para);
    //dichotomy vmax when Tv < 0.0, low precision
    int TrajectoryTime4(float &time, float Q, float v_0, float v_1, float vmax, float amax, float jmax, VectorXf &para);
    int TrajectoryTime5(double &time, float Q, float v_0, float v_1, float vmax, float amax, float jmax, VectorXd &para);

    void CorrentionParameters(VectorXf& para, int peroid);
    void TrajectoryCalculation(float t, VectorXf para, float q_dq[3]);
    void TrajectoryCalculationD(double t, VectorXd para, double q_dq[3]);

    void Movep(vector<PointInformation> waypoints, vector<PathInformation> &pathes);
    void CompoundTrajectory(vector<PointInformation> &waypoints, vector<PathInformation> &path);
    void CalculatePathParameters(PathInformation path);
    int TrajectoryTime(AccelerationType &accelerationtype, double &time, double Q, double v_0, double v_1, double vmax, double amax, double jmax, VectorXd &para);
    void CorrentionParameters(AccelerationType accelerationtype, VectorXd &para, int peroid);
    void TrajectoryCalculation(double t, VectorXd para, double q_dq[3]);
    void ArcSegmentLineToLine(PointInformation &startpoint, PointInformation &intermediatepoint, PointInformation &endpoint, PathInformation &arcpath);
    void CalculateArcParameterLineToLine(Vector3d intermediate2start, Vector3d intermediate2end, double &radius, Vector3d intermediatepoint, PathInformation ArcSegment);
    void SetLinePathSegment(PathInformation &path, PointInformation startpoint, PointInformation endpoint, PathInformation arcpath);
    void SetCricularPathSegment(PathInformation &path, PointInformation startpoint, PointInformation endpoint, PointInformation middlepoint, Vector3d center, double r, double theta);
    void SetArcPathSegment(PathInformation &path, PointInformation intermediatepoint, Vector3d arcstart, Vector3d arcend, Vector3d arccenter, double radius, double theta);

    double GetRand(double min, double max);

signals:

public:
    TrajectoryType m_trajectorytype;
    AccelerationType m_accelerationtype;
    PointInformation m_posestart;
    PointInformation m_poseend;
    int number;

    vector<VectorXd> m_wayPoints;
    vector<double> m_maxVelocity;
    vector<double> m_maxAcceleration;
    vector<PathInformation> m_pathsegments;
    VectorXf m_parameters = VectorXf::Zero(15);
};

static TrajectoryPlanning& GetTrajectoryPlanningInstance()
{
    static TrajectoryPlanning TrajectoryPlanningInstance;
    return TrajectoryPlanningInstance;
}

#endif // TRAJECTORYPLANNING_H
