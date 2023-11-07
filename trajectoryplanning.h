#ifndef TRAJECTORYPLANNING_H
#define TRAJECTORYPLANNING_H

#include <QObject>
#include <eigen-3.4.0/Eigen/Dense>

using namespace std;
using namespace Eigen;

struct Parameters {
    float Ta;//0
    float Tv;//1
    float Td;//2
    float Tj1;//3
    float Tj2;//4
    float q0;//5
    float q1;//6
    float v0;//7
    float v1;//8
    float vlim;//9
    float amax;//10
    float amin;//11
    float alima;//12
    float alimd;//13
    float jmax;//14
    float jmin; //15
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
    PointInformation intermediatepoint; //中间点
    Vector3d waypoint; //路径点
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
    VectorXd constraints = VectorXd::Zero(16); //16个参数的约束
    AccelerationType accelerationType; //加速度类型
    TrajectoryType trajectoryType; //姿态还是位移
//    vector<vector<double> > TrajectoryInterpolation;
//    vector<vector<double> > TrajectoryVelocity;
//    vector<vector<double> > TrajectoryAcceleration;
//    _TrajectoryInformation() {
////        memset(this,0,sizeof(_TrajectoryInformation));
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
    bool TrajectoryInterpolation(vector<PathInformation> pathes, int peroid);
    void CompoundTrajectory(vector<PointInformation> &waypoints, vector<PathInformation> &path);
    void CalculatePathParameters(PathInformation &path);
    int TrajectoryTime(AccelerationType &accelerationtype, double &time, double Q, double v_0, double v_1, double vmax, double amax, double jmax, VectorXd &para);
    void CorrentionParameters(AccelerationType accelerationtype, VectorXd &para, int peroid);
    void TrajectoryCalculation(double t, VectorXd para, double q_dq[3]);
    void ArcSegmentLineToLine(PointInformation &startpoint, PointInformation &intermediatepoint, PointInformation &endpoint, PathInformation &arcpath);
    void CalculateArcParameterLineToLine(Vector3d intermediate2start, Vector3d intermediate2end, double &radius, Vector3d intermediatepoint, PathInformation ArcSegment);
    void SetLinePathSegment(PathInformation &path, PointInformation startpoint, PointInformation endpoint, PathInformation arcpath);
    void SetCricularPathSegment(PathInformation &path, PointInformation startpoint, PointInformation endpoint, PointInformation middlepoint, Vector3d center, double r, double theta);
    void SetArcPathSegment(PathInformation &path, PointInformation intermediatepoint, Vector3d arcstart, Vector3d arcend, Vector3d arccenter, double radius, double theta);

    void ArcParameterCalculate(Vector3d p1, Vector3d p2, Vector3d p3, Vector3d &center, double radius, Vector3d normal, double theta);
    Vector3d LineIntersectCricular(Vector3d startpoint, Vector3d endpoint, double radius);
    Vector3d CricularIntersectCricular(Vector3d center, double circularradius, Vector3d normal, Vector3d intersectpoint, double radius);

    void LineIntersectSphere(Vector3d startpoint, Vector3d endpoint, Vector3d node, double radius);
    void SolvingQuadratics(double a, double b, double c, vector<double> &t);
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
