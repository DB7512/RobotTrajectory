#ifndef TRAJECTORYPLANNING_H
#define TRAJECTORYPLANNING_H

#include <QObject>
#include <eigen-3.4.0/Eigen/Dense>

using namespace std;
using namespace Eigen;

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
}AccelerationType;  // 由velocityplanning中的velocitytype取代

typedef enum {
    Position = 0,
    Posture,
}TrajectoryType;

typedef enum {
    None = 0,
    Line,       // 重新定义为直线
    Arc,        // 过渡圆弧 // 重新定义为圆弧
    Circular,   // 圆弧   // 重新定义为圆
    LineArc,   // 重新定义为过渡圆弧
    Curve,      // 重新定义为曲线
}PathType;  // 分段路径分类

typedef enum {
    None = 0,
    Line,       //直线
    Arc,        //圆弧
}PointType; // 原始路点分类

typedef enum {
    Fixed,
    Unfixed,
}ArcType; // 原始路点分类

typedef enum {
    ConnectwithStart,
    ConnectwithEnd,
}ConnectionType;    // 与分段路径的起点还是终点连接

typedef struct _PointInformation {
    Vector3d point; // 位置
    Vector4d pose; // 姿态四元数
    double velocity; // 速度
    double acceleration; // 加速度
    double jerk;    // 加加速度
    double maxDeviation; // 最大误差
    double maxVelocity; // 最大速度
    double maxAcceleration; // 最大加速度
    double maxJerk; // 最大加加速度
    // double radius;  // 圆弧半径 !!!可省略
    double r;   // 交融半径
    // PointType pointType; // 该路点所属路径类型
    _PointInformation() {
        memset(this,0,sizeof(_PointInformation));
        maxJerk = 80.0;
        // pointType = None;
    }
}PointInformation;

typedef struct _PathInformation {
    PathType pathType; // 路径类型
    ArcType arcType;    // 圆弧姿态规划类型
    PointInformation startpoint;    // 路径起点
    PointInformation endpoint;  // 路径终点
    PointInformation intermediatepoint; // 圆弧路径中间点
    vector<PointInformation> controlpoints; // 曲线控制点
    Vector3d waypoint; // 预留使用
    double maxVelocity; // 最大速度
    double maxAcceleration; // 最大加速度
    double maxJerk; // 最大加加速度
    double radius; // 半径
    double theta; // 圆心角
    double displacement; // 位移
    double arclength; // 弧长
    Vector3d center; // 圆心
    Vector3d normal; // 圆弧旋转轴
    VectorXd paras = VectorXd::Zero(16); // 16个参数的约束
    AccelerationType accelerationType; // 加速度类型
    TrajectoryType trajectoryType; // 姿态还是位移
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


    // 计算圆弧的参数，包括圆心、半径、旋转轴、旋转角
    void GetArcParameter(Vector3d p1, Vector3d p2, Vector3d p3,
                               Vector3d &center, double &radius, Vector3d &normal, double &theta);
    void GetArcParameter(PathInformation &arc);
    // 线段与以端点为圆心，r为半径的圆的交点
    Vector3d LineIntersectCricular(Vector3d startpoint, Vector3d endpoint, double radius, ConnectionType type);
    Vector3d LineIntersectCricular(PathInformation line, ConnectionType type);
    // 圆弧与以圆弧端点为圆心，r为半径的圆的交点
    Vector3d CricularIntersectCricular(Vector3d center, double circularradius, Vector3d normal,
                                       Vector3d intersectpoint, double radius, ConnectionType type);
    Vector3d CricularIntersectCricular(PathInformation arc, ConnectionType type);
    // 线段的单位切向量
    Vector3d GetLineTangent(Vector3d startpoint, Vector3d endpoint, ConnectionType type);
    Vector3d GetLineTangent(PathInformation line, ConnectionType type);
    // 圆弧端点的单位切向量（指向远离圆弧的方向）
    Vector3d GetArcTangent(Vector3d startpoint, Vector3d endpoint, Vector3d center,
                           Vector3d normal, ConnectionType type);
    Vector3d GetArcTangent(PathInformation arc, ConnectionType type);
    // 计算B样条曲线的控制点
    vector<PointInformation> GetControlPoints(PathInformation path1, PathInformation path2);
    // 计算曲线参数（弧长等）
    void GetCurveParameter(PathInformation &path);


    // 获取movep路点信息
    PointInformation GetPointInformation(Vector3d point, Vector4d pose,
                                         double v, double a, double vm, double am, double jm,
                                         double r, PointType type);
    vector<PointInformation> m_movePoint;   // 用于存放路点信息
    // 将路点信息转换成没有交融的路径信息
    int GetPathInformation(vector<PointInformation> waypoints, vector<PathInformation> pathes);

    void MoveP(vector<PathInformation> path);
    // 1.将原始路径转换成交融路径
    int SetPathSegment(vector<PathInformation> pathf, vector<PathInformation> patht);
    vector<PathInformation> m_movePath; // 2.存放路径信息
    // 设置原始直线信息
    PathInformation SetLinePath(PointInformation pf, PointInformation pt);
    // 设置原始圆弧信息
    PathInformation SetArcPath(PointInformation pf, PointInformation pi, PointInformation pt);
    // 设置曲线信息，约束沿用前一路径的约束信息
    void SetCurvePath(vector<PointInformation> ctrlpoint, PathInformation &path);
    // 路径衔接
    void GetPathSegment(PathInformation path1, PathInformation &path2, vector<PathInformation> &path);


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
