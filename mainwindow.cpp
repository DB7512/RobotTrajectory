#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "trajectoryplanning.h"
#include <QTextStream>
#include <QFile>

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


void MainWindow::on_Test_clicked()
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
                        GetTrajectoryPlanningInstance().m_posestart.v = v_0;
                        GetTrajectoryPlanningInstance().m_posestart.vmax = vmax;
                        GetTrajectoryPlanningInstance().m_posestart.amax = amax;
                        GetTrajectoryPlanningInstance().m_posestart.jmax = jmax;
                        GetTrajectoryPlanningInstance().m_poseend.vmax = vmax;
                        GetTrajectoryPlanningInstance().m_poseend.amax = amax;
                        GetTrajectoryPlanningInstance().m_poseend.jmax = jmax;
                        GetTrajectoryPlanningInstance().m_poseend.v = v_1;
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


void MainWindow::on_Test_2_clicked()
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

