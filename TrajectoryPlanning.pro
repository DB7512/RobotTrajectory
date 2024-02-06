QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

CONFIG += c++17

# You can make your code fail to compile if it uses deprecated APIs.
# In order to do so, uncomment the following line.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

SOURCES += \
    bsplinecurve.cpp \
    main.cpp \
    mainwindow.cpp \
    mathfunction.cpp \
    pathconstruction.cpp \
    trajectoryplanning.cpp \
    velocityplanning.cpp

HEADERS += \
    bsplinecurve.h \
    mainwindow.h \
    mathfunction.h \
    pathconstruction.h \
    trajectoryplanning.h \
    velocityplanning.h

FORMS += \
    mainwindow.ui

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target

INCLUDEPATH += $$PWD/../
INCLUDEPATH += \home\hua\Trajectory\TrajectoryPlanning\eigen-3.4.0
DEPENDPATH += $$PWD/../
