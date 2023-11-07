#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "trajectoryplanning.h"
#include <QTextStream>
#include <QFile>
#include "mathfunction.h"
#include <QDebug>

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::TestTimeCalculation()
{
    float q = 0.05;
    float jmax = 80.0;
    vector<float>interpolation_result;//插值点
    vector<vector<float> >interpolation_inf;
    QFile file("data_inf.txt");
    if(!file.open(QIODevice::WriteOnly | QIODevice::Text | QIODevice::Truncate)) return;
    file.close();
    QFile data("data.txt");
    if(!data.open(QIODevice::WriteOnly | QIODevice::Text | QIODevice::Truncate)) return;
    data.close();
    GetTrajectoryPlanningInstance().number = 0;
    for (int var = 0; var < 20; var++) {
        float v_0 = 0.1;
        for (int m = 0; m < 20; m++) {
            float v_1 = 0.1;
            for(int i = 0; i < 40; i++) {
                float vmax = 0.1;
                for(int j = 0; j < 20; j++) {
                    float amax = 0.2;
                    for(int k = 0; k < 30; k++) {
                        if(vmax < v_0 || vmax < v_1) break;
                        GetTrajectoryPlanningInstance().number += 1;
                        QFile data("data_inf.txt");
                        if(!data.open(QIODevice::WriteOnly | QIODevice::Text | QIODevice::Append)) return;
                        QTextStream stream(&data);
                        stream<<GetTrajectoryPlanningInstance().number<<" pose "<<q<<" v0 "<<v_0<<" v1 "<<v_1<<" vmax "<<vmax<<" a "<<amax<<"\n";
                        data.close();
                        if(GetTrajectoryPlanningInstance().number == 32270) {
                            //                        if(1) {
                            GetTrajectoryPlanningInstance().TimePlanning(100, q, vmax, amax, jmax, v_0, v_1, interpolation_result, interpolation_inf);
                        }
                        amax += 0.2;
                    }
                    vmax += 0.1;
                }
                v_1 += 0.1;
            }
            v_0 += 0.1;
        }
        q += 0.05;
    }
    int a = 0;
    a = 1;
}

void MainWindow::TestInterpolationCalculation()
{
    float start_point[6] = {0,0,0,0,0,0};
    float end_point[6] = {0,0,50.0,0,0,0};
    float jmax = 80.0;
    vector<vector<float> >interpolation_result;//插值点
    vector<vector<float> >interpolation_inf;
    QFile file("data_inf.txt");
    if(!file.open(QIODevice::WriteOnly | QIODevice::Text | QIODevice::Truncate)) return;
    file.close();
    QFile data("data_new.txt");
    if(!data.open(QIODevice::WriteOnly | QIODevice::Text | QIODevice::Truncate)) return;
    data.close();
    GetTrajectoryPlanningInstance().number = 0;
    for (int var = 0; var < 20; var++) {
        for (int i = 0; i < 6; i++) {
            GetTrajectoryPlanningInstance().m_posestart.pose[i] = start_point[i];
            GetTrajectoryPlanningInstance().m_poseend.pose[i] = end_point[i];
        }
        float v_0 = 0.1;
        for (int m = 0; m < 20; m++) {
            float v_1 = 0.1;
            for(int i = 0; i < 40; i++) {
                float vmax = 0.1;
                for(int j = 0; j < 20; j++) {
                    float amax = 0.2;
                    for(int k = 0; k < 30; k++) {
                        if(vmax < v_0 || vmax < v_1) break;
                        GetTrajectoryPlanningInstance().number += 1;
                        GetTrajectoryPlanningInstance().m_posestart.velocity = v_0;
                        GetTrajectoryPlanningInstance().m_posestart.maxVelocity = vmax;
                        GetTrajectoryPlanningInstance().m_posestart.maxAcceleration = amax;
                        GetTrajectoryPlanningInstance().m_posestart.maxJerk = jmax;
                        GetTrajectoryPlanningInstance().m_poseend.maxVelocity = vmax;
                        GetTrajectoryPlanningInstance().m_poseend.maxAcceleration = amax;
                        GetTrajectoryPlanningInstance().m_poseend.maxJerk = jmax;
                        GetTrajectoryPlanningInstance().m_poseend.velocity = v_1;
                        QFile data("data_inf.txt");
                        if(!data.open(QIODevice::WriteOnly | QIODevice::Text | QIODevice::Append)) return;
                        QTextStream stream(&data);
                        stream<<GetTrajectoryPlanningInstance().number<<" pose "<<end_point[2]<<" v0 "<<v_0<<" v1 "<<v_1<<" vmax "<<vmax<<" a "<<amax<<"\n";
                        data.close();
                        //                        if(GetTrajectoryPlanningInstance().number == 39271) {
                        if(1) {
                            GetTrajectoryPlanningInstance().LinePlanning(100, GetTrajectoryPlanningInstance().m_posestart, GetTrajectoryPlanningInstance().m_poseend, interpolation_result, interpolation_inf);
                        }
                        amax += 0.2;
                    }
                    vmax += 0.1;
                }
                v_1 += 0.1;
            }
            v_0 += 0.1;
        }
        end_point[2] += 50.0;
    }
    int a = 0;
    a = 1;
}

void MainWindow::TestMovep()
{
    QFile file("information.txt");
    if(!file.open(QIODevice::WriteOnly | QIODevice::Text | QIODevice::Truncate)) return;
    file.close();
    vector<PointInformation> waypoints;
    vector<PathInformation> pathes;
    PointInformation waypoint;
    Vector3d point; //位置
    //Vector4d pose; //姿态四元数
    double velocity; //速度
    double acceleration; //加速度
    double jerk;
    double maxVelocity;
    double maxAcceleration;
    double radius;
    PathType pathType; //该路点所属路径类型
    point<<0.1,0.2,0.3;
    velocity = 0.1;
    acceleration = 0.6;
    jerk = 80.0;
    maxVelocity = 1.0;
    maxAcceleration = 3.0;
    radius = 0.010;
    pathType = Line;
    double rand;
    waypoint.maxVelocity = maxVelocity;
    waypoint.maxAcceleration = maxAcceleration;
    waypoint.pathType = pathType;
    waypoint.jerk = jerk;
    vector<Vector3d> defwaypoints;
    defwaypoints.push_back(point);
    point<<0.1,0.1,0.2555;
    defwaypoints.push_back(point);
    point<<0.2,0.05,0.2;
    defwaypoints.push_back(point);
    point<<0.3,0.05,0.3;
    defwaypoints.push_back(point);
    for (int i = 0; i < defwaypoints.size(); i++) {
        rand = GetMathInstance().GetRand(0.0,1.0);
        waypoint.point = defwaypoints[i];
        waypoint.velocity = velocity * (1+1*rand);
        waypoint.acceleration = acceleration * (1+2*rand);
        waypoint.radius = radius * (1+0.05*rand);
        waypoints.push_back(waypoint);
    }
    GetTrajectoryPlanningInstance().Movep(waypoints, pathes);
//    for(int i = 0; i < pathes.size(); i++) {
//        qDebug()<<pathes[i].startpoint.point[0]<<pathes[i].startpoint.point[1]<<pathes[i].startpoint.point[2];
//        qDebug()<<pathes[i].endpoint.point[0]<<pathes[i].endpoint.point[1]<<pathes[i].endpoint.point[2];
//        if(pathes[i].pathType != Line)
//            qDebug()<<pathes[i].intermediatepoint.point[0]<<pathes[i].intermediatepoint.point[1]<<pathes[i].intermediatepoint.point[2];
//    }
    GetTrajectoryPlanningInstance().TrajectoryInterpolation(pathes,100);
    int aa = 100;
    int qq = 33;
}


void MainWindow::on_Test_clicked()
{
//    TestTimeCalculation();
//    TestInterpolationCalculation();
    TestMovep();
}


void MainWindow::on_Test_2_clicked()
{

}

