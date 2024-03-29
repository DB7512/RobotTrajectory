#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();
    void TestVelocity();
    void TestTime();
    void TestBSpline();
    void TestTimeCalculation();
    void TestInterpolationCalculation();
    void TestMovep();
private slots:
    void on_Test_clicked();

private:
    Ui::MainWindow *ui;
};
#endif // MAINWINDOW_H
