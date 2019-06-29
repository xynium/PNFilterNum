#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QGraphicsTextItem>
#include <QProgressBar>

#include <qfile.h>
#include <qmath.h>
#include "complex.h"
#include "FftComplex.hpp"

//constante remez
#define BANDPASS  1
#define DIFFERENTIATOR  2
#define HILBERT  3
#define NEGATIVE  0
#define POSITIVE  1
#define GRIDDENSITY  64   //16
#define MAXITERATIONS  160  //40
#define ORDRE  16384      //4096   // ordre FFT
#define Pi             3.1415926535897932
#define Pi2            6.2831853071795865
#define MINPLOTDB -200

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow {
                Q_OBJECT

        public:
                explicit MainWindow(QWidget* parent = 0);
                ~MainWindow();

        private slots:
                void on_pbSyntes_clicked();
                void on_pbSim_clicked();

                void on_pbQuit_clicked();

        private:
                void  remez(long double* h, int numtaps, int numband, long double* bands, long double*  des, long double* weight, int type);
                void CreateDenseGrid(int r, int numtaps, int numband, long double bands[], long double des[],long double weight[], int* gridsize,   long double Grid[], long double D[], long double W[],   int symmetry);
                void InitialGuess(int r, int Ext[], int gridsize);
                void CalcParms(int r, int Ext[], long double Grid[], long double D[], long double W[],   long double ad[],long double x[], long double y[]);
                double ComputeA(long double freq, int r,long  double ad[], long double x[], long double y[]) ;
                void CalcError(int r, long double ad[], long double x[],long  double y[], int gridsize, long double Grid[], long  double D[],long  double W[],long  double E[]);
                short isDone(int r, int Ext[],long  double E[]);
                void FreqSample(int N,long  double A[],long  double h[], int symm);
                void Search(int r, int Ext[], int gridsize, long double E[]);
                Complex* FDFT(Complex* cX) ;
                Complex *FT( Complex* cX);
                void WriteTheFile(long double* h,long double* ht);


        private:
                Ui::MainWindow* ui;
                QProgressBar* progressBar;
                int iM;  // Ordre
                long  double* lovFir;
                long double* hiFir;
                long double dFs, dFc, dFBw;
                QString qsFileName;
                int iHasTmpToDel;

};

#endif // MAINWINDOW_H
