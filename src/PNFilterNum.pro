#-------------------------------------------------
#
# Project created by QtCreator 2019-01-08T17:39:32
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets printsupport

TARGET = PNFilterNum
TEMPLATE = app


SOURCES += main.cpp\
        mainwindow.cpp \
    complex.cpp \
    ../qcustomplot/qcustomplot.cpp \
    FftComplex.cpp

HEADERS  += mainwindow.h \
    complex.h \
    ../qcustomplot/qcustomplot.h \
    FftComplex.hpp

FORMS    += mainwindow.ui

QMAKE_CXXFLAGS += -Wall
