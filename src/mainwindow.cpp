// Remez est desactivé functionnal but noway to select it in UI.
// Remez some time not converge, prefer using OCTAVE or SCILAB

#include "mainwindow.h"
#include "ui_mainwindow.h"

//#define PHHILBERT

MainWindow::MainWindow(QWidget* parent) :
        QMainWindow(parent),
        ui(new Ui::MainWindow) {
        ui->setupUi(this);

        progressBar = new QProgressBar();
        progressBar->setRange(1, 100);
        progressBar->setValue(0);
        progressBar->setTextVisible(true);
        progressBar->setFormat("%p");

        QWidget* w = new QWidget;
        QHBoxLayout* hlayout = new QHBoxLayout();
        hlayout->addSpacerItem(new QSpacerItem(40, 1));
        hlayout->addWidget(progressBar);
        w->setLayout(hlayout);

        ui->statusBar->addPermanentWidget(w, 2);

        iHasTmpToDel = 0;

}

MainWindow::~MainWindow() {
        QFile* mFilterFile;
        QString input;

        if (  iHasTmpToDel == 1) {
                input = qsFileName + "T";
                mFilterFile = new QFile(input);
                mFilterFile->remove();
        }
        delete ui;
}

void MainWindow::on_pbSyntes_clicked() {
        int compt1, iGain, iNbts;
        long double* hb, *hh;
        QString input;
        int i;
        long double FC, SUM;
        QVector<double>  S1;
        QVector<double> vdAx1, vdAx0;
        QString qsTmpa, qsTmpb;
        Complex ctmp;
        long double dhh, dhb;

        dFs = ui->leFSamp->text() .toDouble();
        dFc = ui->leFCut->text().toDouble();
        dFBw = ui->leBand->text().toDouble();
        iGain = ui->leGain->text().toInt();
        iM = ui->leOrder->text().toInt();
        iNbts = ui->leNbts->text().toDouble();
        iM++;
        hb = new long double[iM];
        hh =   new long double[iM];

        if (/*ui->cbRemez->checkState() ==!*/true)    { // debut sincard
                lovFir = new long double[iM];
                hiFir = new long double[iM];

                if ( ui->rbLovPass->isChecked()) {
                        FC = (dFc) / (long double)dFs;    //Calculate the low-pass filter kernel
                        for (i = 1; i < iM; i++) {
                                if ((i == iM / 2) && (iM % 2 == 0)) {
                                        if (ui->cbHilbert->isChecked()) lovFir[i] = 0;
                                        else lovFir[i] = 2.0 * M_PIl * FC;
                                } else if (ui->cbHilbert->isChecked())  lovFir[i] = ( 1.0 - cosl(2.0 * M_PIl * FC * (i - iM / 2.0))) / (i - iM / 2.0);
                                else lovFir[i] = sinl(2.0 * M_PIl * FC * (i - iM / 2.0)) / (i - iM / 2.0);
                                hb[i] = lovFir[i] * (0.42 - 0.5 * cosl(2.0 * M_PIl * i / (iM)) + 0.08 * cosl(4.0 * M_PIl * i / (iM ))); // Blakman window
                                hh[i] = lovFir[i] * (0.54 - 0.46 * cosl(2.0 * M_PIl * i / (iM ))); // hamming window
                        }

                } else {   // Band Pass
                        FC = (dFc - dFBw / 2) / (long double)dFs;    //Calculate the low-pass filter kernel
                        for (i = 1; i < iM; i++) {
                                if ((i == iM / 2) && (iM % 2 == 0)) {
                                        if (ui->cbHilbert->isChecked()) lovFir[i] = 0;
                                        else lovFir[i] = 2.0 * M_PIl * FC;
                                } else if (ui->cbHilbert->isChecked())  lovFir[i] = ( 1.0 - cosl(2.0 * M_PIl * FC * (i - iM / 2.0))) / (i - iM / 2.0);
                                else lovFir[i] = sinl(2.0 * M_PIl * FC * (i - iM / 2.0)) / (i - iM / 2.0);
                                //   lovFir[i] = lovFir[i] * (0.42 - 0.5 * cosl(2.0 * M_PIl * i / iM) + 0.08 * cosl(4.0 * M_PIl * i / iM));
                        }

                        FC = (dFc + dFBw / 2) / (long double)dFs;    //Calculate the second low-pass filter kernel
                        for (i = 1; i < iM; i++) {
                                if ((i == iM / 2) && (iM % 2 == 0)) {
                                        if (ui->cbHilbert->isChecked()) hiFir[i] = 0;
                                        else  hiFir[i] = 2.0 * M_PIl * FC;
                                } else if (ui->cbHilbert->isChecked()) hiFir[i] = (1.0 - cosl(2.0 * M_PIl * FC * (i - iM / 2.0))) / (i - iM / 2.0);
                                else hiFir[i] = sinl(2.0 * M_PIl * FC * (i - iM / 2.0)) / (i - iM / 2.0);
                                //   hiFir[i] = hiFir[i] * (0.42 - 0.5 * cosl(2.0 * M_PIl * i / iM) + 0.08 * cosl(4.0 * M_PIl * i / iM));
                        }

                        for (i = 1; i < iM; i++) {   //Change the low-pass filter kernel in B[ ] into a high-pass
                                hiFir[i] = -hiFir[i];
                        }

                        for (i = 1; i < iM; i++) {      //Add the low-pass filter kernel ], to the high-pass
                                hb[i] = lovFir[i] + hiFir[i];
                        }

                        for (i = 1; i < iM; i++) {    //Change the band-reject filter kernel into a band-pass
                                hb[i] = -hb[i];
                        }

                        for (i = 1; i < iM; i++) {    //Windowing
                                hh[i] = hb[i] * (0.54 - 0.46 * cosl(2.0 * M_PIl * i / (iM ))); // hamming window
                                hb[i] = hb[i] * (0.42 - 0.5 * cosl(2.0 * M_PIl * i / iM) + 0.08 * cosl(4.0 * M_PIl * i / iM));
                        }
                }
        } // fin  sinus cardinal
        else {  //  remez Search   impossible to get here desactivated in UI
                long double*  bands = new long double[6 ];
                long double*   weights = new long double[3 ];
                long double*    desired = new long double[3];
                long double dtrans ;
                int iNb;
                //         dtrans = ui->leRemTrans->text().toDouble();

                weights[0] = 1;
                weights[1] = 1;
                weights[2] = 1;

                iNb = 2;     // cas pass bas
                desired[0] = 1;
                desired[1] = 0;

                bands[0] = 0;
                bands[1] = dFc / dFs;
                bands[2] =  (dFc + dtrans ) / dFs ;
                bands[3] = 0.5;

                /*       if(ui->rbHiPass->isChecked()) {
                               iNb = 2;
                               weights[0] = 1;
                               weights[1] = 1;

                               desired[0] = 0;
                               desired[1] = 1;

                               bands[0] = 0;
                               bands[1] = (dFc - dtrans ) / dFs ;
                               bands[2] = dFc  / dFs;
                               bands[3] = 0.5;
                       }*/
                if (ui->rbBandPass->isChecked()) {
                        iNb = 3;
                        weights[0] = 10;
                        weights[1] = 1;
                        weights[2] = 10;

                        desired[0] = 0;
                        desired[1] = 1;
                        desired[2] = 0;

                        bands[0] = 0;
                        bands[1] = ((dFc - dtrans - dFBw / 2.0) / dFs)  ;
                        bands[2] = (dFc - dFBw / 2.0) / dFs;
                        bands[3] = (dFc + dFBw / 2.0) / dFs;
                        bands[4] = ( (dFc + dtrans + dFBw / 2.0) / dFs)  ;
                        bands[5] = 0.5;
                }
                /*  if (ui->rbNotch->isChecked()) {
                          iNb = 3;
                          weights[0] = 1;
                          weights[1] = 10;
                          weights[2] = 1;

                          desired[0] = 1;
                          desired[1] = 0;
                          desired[2] = 1;

                          bands[0] = 0;
                          bands[1] = ((dFc - dtrans - dFBw / 2.0) / dFs)  ;
                          bands[2] = (dFc - dFBw / 2.0) / dFs;
                          bands[3] = (dFc + dFBw / 2.0) / dFs;
                          bands[4] = ( (dFc + dtrans + dFBw / 2.0) / dFs)  ;
                          bands[5] = 0.5;
                  }*/
                if (ui->cbHilbert->isChecked())  remez(hh, iM - 1, iNb, bands, desired, weights, HILBERT);
                else remez(hh, iM - 1, iNb, bands, desired, weights, BANDPASS);

                for (iNb = iM - 1; iNb > 0; iNb--) hh[iNb] = hh[iNb - 1];

        } //fin remez

        // Normalise le gain
        /*  dhb = 0;
          dhh = 0;
          for (compt1 = 1; compt1 < iM; compt1++) {
                  dhb += hb[compt1] ;
                  dhh += hh[compt1] ;
          }*/

        dhb = iGain / M_PIl;
        dhh = iGain / M_PIl;

        /*   if (ui->cbHilbert->isChecked()) {
                   dhb /= M_PIl;
                   dhh /= M_PIl;
           }*/

        for (compt1 = 1; compt1 < iM; compt1++) {
                hh[compt1] *= dhh;
                hb[compt1] *= dhb;
        }


        if (ui->cbInt->isChecked()) { // formate po int
                dhb = 0; dhh = 0;
                for (int compt1 = 1; compt1 < iM; compt1++) { //cherche max
                        if(fabsl(hb[compt1]) > dhb)  dhb = fabsl(hb[compt1]);
                        if(fabsl(hh[compt1]) > dhh)  dhh = fabsl(hh[compt1]);
                }
                long double dSum = powl(2.0, iNbts);
                dhh = dSum / (2.0 * dhh);
                dhb = dSum / (2.0 * dhb);

                for (compt1 = 1; compt1 < iM; compt1++) {
                        hh[compt1] *= dhh;
                        hb[compt1] *= dhb;
                }
        }

        if (ui->rbBlack->isChecked())
                WriteTheFile(hb, hh);
        else
                WriteTheFile(hh, hb);

        for (compt1 = 0; compt1 < iM; compt1++) {     //trace le resultat
                if (ui->rbBlack->isChecked())
                        S1 << hb[compt1];
                else
                        S1 << hh[compt1];

                vdAx1 << compt1;
        }

        //Affiche la courbe des coeff
        //qsTmpa = "PNFilterNum : Response";
        ui->GraphView->clearGraphs();
        ui->GraphView->addGraph(ui->GraphView->xAxis, ui->GraphView->yAxis);
        ui->GraphView->graph(0)->setPen(QPen(Qt::blue)); // line color blue for first graph
        ui->GraphView->graph(0)->setBrush(QBrush(QColor(0, 0, 255, 20))); // first graph will be filled with translucent blue

        ui->GraphView->xAxis2->setVisible(true);
        ui->GraphView->xAxis2->setTickLabels(true);
        ui->GraphView->yAxis2->setVisible(true);
        ui->GraphView->yAxis2->setTickLabels(true);
        //make left and bottom axes always transfer their ranges to right and top axes:
        connect(ui->GraphView->xAxis, SIGNAL(rangeChanged(QCPRange)), ui->GraphView->xAxis2, SLOT(setRange(QCPRange)));
        connect(ui->GraphView->yAxis, SIGNAL(rangeChanged(QCPRange)), ui->GraphView->yAxis2, SLOT(setRange(QCPRange)));
        //pass data points to graphs:
        ui->GraphView->graph(0)->setData(vdAx1, S1);
        //  ui->GraphView->graph(1)->setData(vdAx1, S0);
        ui->GraphView->graph(0)->rescaleAxes();
        //   ui->GraphView->graph(1)->rescaleAxes(true);
        //ui->GraphView->graph(1)->rescaleAxes(true);
        ui->GraphView->yAxis->setTickLabelColor(Qt::black);
        ui->GraphView->yAxis2->setTickLabelColor(Qt::black);
        ui->GraphView->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
        //ui->GraphView->graph(1)->setName("Phase");
        ui->GraphView->xAxis->setLabel("Number");


        //setWindowTitle(qsTmpa);
        ui->GraphView->replot();
}

/*******************
 * CreateDenseGrid
 *=================
 * Creates the dense grid of frequencies from the specified bands.
 * Also creates the Desired Frequency Response function (D[]) and
 * the Weight function (W[]) on that dense grid
 *
 *
 * INPUT:
 * ------
 * int      r        - 1/2 the number of filter coefficients
 * int      numtaps  - Number of taps in the resulting filter
 * int      numband  - Number of bands in user specification
 * double   bands[]  - User-specified band edges [2*numband]
 * double   des[]    - Desired response per band [numband]
 * double   weight[] - Weight per band [numband]
 * int      symmetry - Symmetry of filter - used for grid check
 *
 * OUTPUT:
 * -------
 * int    gridsize   - Number of elements in the dense frequency grid
 * double Grid[]     - Frequencies (0 to 0.5) on the dense grid [gridsize]
 * double D[]        - Desired response on the dense grid [gridsize]
 * double W[]        - Weight function on the dense grid [gridsize]
 *******************/

void MainWindow::CreateDenseGrid(int r, int numtaps, int numband, long double bands[],  long double des[], long double weight[], int* gridsize,  long double Grid[], long double D[], long double W[],  int symmetry) {
        int i, j, k, band;
        long double delf, lowf, highf;

        delf = 0.5 / (GRIDDENSITY * r);

        /*
         * For differentiator, hilbert,
         *   symmetry is odd and Grid[0] = max(delf, band[0])
         */

        if ((symmetry == NEGATIVE) && (delf > bands[0]))
                bands[0] = delf;

        j = 0;
        for (band = 0; band < numband; band++) {
                Grid[j] = bands[2 * band];
                lowf = bands[2 * band];
                highf = bands[2 * band + 1];
                k = (int)((highf - lowf) / delf + 0.5); /* .5 for rounding */
                for (i = 0; i < k; i++) {
                        D[j] = des[band];
                        W[j] = weight[band];
                        Grid[j] = lowf;
                        lowf += delf;
                        j++;
                }
                Grid[j - 1] = highf;
        }

        /*
         * Similar to above, if odd symmetry, last grid point can't be .5
         *  - but, if there are even taps, leave the last grid point at .5
         */
        if ((symmetry == NEGATIVE) &&
                        (Grid[*gridsize - 1] > (0.5 - delf)) &&
                        (numtaps % 2)) {
                Grid[*gridsize - 1] = 0.5 - delf;
        }
}


/********************
 * InitialGuess
 *==============
 * Places Extremal Frequencies evenly throughout the dense grid.
 *
 *
 * INPUT:
 * ------
 * int r        - 1/2 the number of filter coefficients
 * int gridsize - Number of elements in the dense frequency grid
 *
 * OUTPUT:
 * -------
 * int Ext[]    - Extremal indexes to dense frequency grid [r+1]
 ********************/

void MainWindow::InitialGuess(int r, int Ext[], int gridsize) {
        int i;

        for (i = 0; i <= r; i++)
                Ext[i] = i * (gridsize - 1) / r;
}


/***********************
 * CalcParms
 *===========
 *
 *
 * INPUT:
 * ------
 * int    r      - 1/2 the number of filter coefficients
 * int    Ext[]  - Extremal indexes to dense frequency grid [r+1]
 * double Grid[] - Frequencies (0 to 0.5) on the dense grid [gridsize]
 * double D[]    - Desired response on the dense grid [gridsize]
 * double W[]    - Weight function on the dense grid [gridsize]
 *
 * OUTPUT:
 * -------
 * double ad[]   - 'b' in Oppenheim & Schafer [r+1]
 * double x[]    - [r+1]
 * double y[]    - 'C' in Oppenheim & Schafer [r+1]
 ***********************/

void MainWindow::CalcParms(int r, int Ext[], long double Grid[], long double D[], long double W[],  long double ad[], long double x[], long double y[]) {
        int i, j, k, ld;
        long double sign, xi, delta, denom, numer;

        /*
         * Find x[]
         */
        for (i = 0; i <= r; i++)
                x[i] = cosl(Pi2 * Grid[Ext[i]]);

        /*
         * Calculate ad[]  - Oppenheim & Schafer eq 7.132
         */
        ld = (r - 1) / 15 + 1;     /* Skips around to avoid round errors */
        for (i = 0; i <= r; i++) {
                denom = 1.0;
                xi = x[i];
                for (j = 0; j < ld; j++) {
                        for (k = j; k <= r; k += ld)
                                if (k != i)
                                        denom *= 2.0 * (xi - x[k]);
                }
                if (fabs(denom) < 0.00001)
                        denom = 0.00001;
                ad[i] = 1.0 / denom;
        }

        /*
         * Calculate delta  - Oppenheim & Schafer eq 7.131
         */
        numer = denom = 0;
        sign = 1;
        for (i = 0; i <= r; i++) {
                numer += ad[i] * D[Ext[i]];
                denom += sign * ad[i] / W[Ext[i]];
                sign = -sign;
        }
        delta = numer / denom;
        sign = 1;

        /*
         * Calculate y[]  - Oppenheim & Schafer eq 7.133b
         */
        for (i = 0; i <= r; i++) {
                y[i] = D[Ext[i]] - sign * delta / W[Ext[i]];
                sign = -sign;
        }
}


/*********************
 * ComputeA
 *==========
 * Using values calculated in CalcParms, ComputeA calculates the
 * actual filter response at a given frequency (freq).  Uses
 * eq 7.133a from Oppenheim & Schafer.
 *
 *
 * INPUT:
 * ------
 * double freq - Frequency (0 to 0.5) at which to calculate A
 * int    r    - 1/2 the number of filter coefficients
 * double ad[] - 'b' in Oppenheim & Schafer [r+1]
 * double x[]  - [r+1]
 * double y[]  - 'C' in Oppenheim & Schafer [r+1]
 *
 * OUTPUT:
 * -------
 * Returns double value of A[freq]
 *********************/

double MainWindow::ComputeA(long double freq, int r, long double ad[], long double x[], long double y[]) {
        int i;
        long double xc, c, denom, numer;

        denom = numer = 0;
        xc = cosl(Pi2 * freq);
        for (i = 0; i <= r; i++) {
                c = xc - x[i];
                if (fabsl(c) < 1.0e-7) {
                        numer = y[i];
                        denom = 1;
                        break;
                }
                c = ad[i] / c;
                denom += c;
                numer += c * y[i];
        }
        return numer / denom;
}


/************************
 * CalcError
 *===========
 * Calculates the Error function from the desired frequency response
 * on the dense grid (D[]), the weight function on the dense grid (W[]),
 * and the present response calculation (A[])
 *
 *
 * INPUT:
 * ------
 * int    r      - 1/2 the number of filter coefficients
 * double ad[]   - [r+1]
 * double x[]    - [r+1]
 * double y[]    - [r+1]
 * int gridsize  - Number of elements in the dense frequency grid
 * double Grid[] - Frequencies on the dense grid [gridsize]
 * double D[]    - Desired response on the dense grid [gridsize]
 * double W[]    - Weight function on the desnse grid [gridsize]
 *
 * OUTPUT:
 * -------
 * double E[]    - Error function on dense grid [gridsize]
 ************************/

void MainWindow::CalcError(int r, long double ad[], long  double x[], long double y[],     int gridsize, long double Grid[], long double D[], long double W[], long double E[]) {
        int i;
        long  double A;

        for (i = 0; i < gridsize; i++) {
                A = ComputeA(Grid[i], r, ad, x, y);
                E[i] = W[i] * (D[i] - A);
        }
}

/************************
 * Search
 *========
 * Searches for the maxima/minima of the error curve.  If more than
 * r+1 extrema are found, it uses the following heuristic (thanks
 * Chris Hanson):
 * 1) Adjacent non-alternating extrema deleted first.
 * 2) If there are more than one excess extrema, delete the
 *    one with the smallest error.  This will create a non-alternation
 *    condition that is fixed by 1).
 * 3) If there is exactly one excess extremum, delete the smaller
 *    of the first/last extremum
 *
 *
 * INPUT:
 * ------
 * int    r        - 1/2 the number of filter coefficients
 * int    Ext[]    - Indexes to Grid[] of extremal frequencies [r+1]
 * int    gridsize - Number of elements in the dense frequency grid
 * double E[]      - Array of error values.  [gridsize]
 * OUTPUT:
 * -------
 * int    Ext[]    - New indexes to extremal frequencies [r+1]
 ************************/

void MainWindow::Search(int r, int Ext[],            int gridsize, long double E[]) {
        int i, j, k, l, extra;     /* Counters */
        int up, alt;
        int* foundExt;             /* Array of found extremals */

        /*
         * Allocate enough space for found extremals.
         */
        foundExt = (int*)malloc((2 * r) * sizeof(int));
        k = 0;

        /*
         * Check for extremum at 0.
         */
        if (((E[0] > 0.0) && (E[0] > E[1])) ||
                        ((E[0] < 0.0) && (E[0] < E[1])))
                foundExt[k++] = 0;

        /*
         * Check for extrema inside dense grid
         */
        for (i = 1; i < gridsize - 1; i++) {
                if (((E[i] >= E[i - 1]) && (E[i] > E[i + 1]) && (E[i] > 0.0)) ||
                                ((E[i] <= E[i - 1]) && (E[i] < E[i + 1]) && (E[i] < 0.0)))
                        foundExt[k++] = i;
        }

        /*
         * Check for extremum at 0.5
         */
        j = gridsize - 1;
        if (((E[j] > 0.0) && (E[j] > E[j - 1])) ||
                        ((E[j] < 0.0) && (E[j] < E[j - 1])))
                foundExt[k++] = j;


        /*
         * Remove extra extremals
         */
        extra = k - (r + 1);

        while (extra > 0) {
                if (E[foundExt[0]] > 0.0)
                        up = 1;                /* first one is a maxima */
                else
                        up = 0;                /* first one is a minima */

                l = 0;
                alt = 1;
                for (j = 1; j < k; j++) {
                        if (fabsl(E[foundExt[j]]) < fabsl(E[foundExt[l]]))
                                l = j;               /* new smallest error. */
                        if ((up) && (E[foundExt[j]] < 0.0))
                                up = 0;             /* switch to a minima */
                        else if ((!up) && (E[foundExt[j]] > 0.0))
                                up = 1;             /* switch to a maxima */
                        else {
                                alt = 0;
                                break;              /* Ooops, found two non-alternating */
                        }                      /* extrema.  Delete smallest of them */
                }  /* if the loop finishes, all extrema are alternating */

                /*
                 * If there's only one extremal and all are alternating,
                 * delete the smallest of the first/last extremals.
                 */
                if ((alt) && (extra == 1)) {
                        if (fabsl(E[foundExt[k - 1]]) < fabsl(E[foundExt[0]]))
                                l = foundExt[k - 1]; /* Delete last extremal */
                        else
                                l = foundExt[0];     /* Delete first extremal */
                }

                for (j = l; j < k; j++) {  /* Loop that does the deletion */
                        foundExt[j] = foundExt[j + 1];
                }
                k--;
                extra--;
        }

        for (i = 0; i <= r; i++) {
                Ext[i] = foundExt[i];       /* Copy found extremals to Ext[] */
        }

        free(foundExt);
}


/*********************
 * FreqSample
 *============
 * Simple frequency sampling algorithm to determine the impulse
 * response h[] from A's found in ComputeA
 *
 *
 * INPUT:
 * ------
 * int      N        - Number of filter coefficients
 * double   A[]      - Sample points of desired response [N/2]
 * int      symmetry - Symmetry of desired filter
 *
 * OUTPUT:
 * -------
 * double h[] - Impulse Response of final filter [N]
 *********************/
void MainWindow::FreqSample(int N, long double A[], long double h[], int symm) {
        int n, k;
        long double x, val, M;

        M = (N - 1.0) / 2.0;
        if (symm == POSITIVE) {
                if (N % 2) {
                        for (n = 0; n < N; n++) {
                                val = A[0];
                                x = Pi2 * (n - M) / N;
                                for (k = 1; k <= M; k++)
                                        val += 2.0 * A[k] * cosl(x * k);
                                h[n] = val / N;
                        }
                } else {
                        for (n = 0; n < N; n++) {
                                val = A[0];
                                x = Pi2 * (n - M) / N;
                                for (k = 1; k <= (N / 2 - 1); k++)
                                        val += 2.0 * A[k] * cosl(x * k);
                                h[n] = val / N;
                        }
                }
        } else {
                if (N % 2) {
                        for (n = 0; n < N; n++) {
                                val = 0;
                                x = Pi2 * (n - M) / N;
                                for (k = 1; k <= M; k++)
                                        val += 2.0 * A[k] * sinl(x * k);
                                h[n] = val / N;
                        }
                } else {
                        for (n = 0; n < N; n++) {
                                val = A[N / 2] * sinl(Pi * (n - M));
                                x = Pi2 * (n - M) / N;
                                for (k = 1; k <= (N / 2 - 1); k++)
                                        val += 2.0 * A[k] * sinl(x * k);
                                h[n] = val / N;
                        }
                }
        }
}

/*******************
 * isDone
 *========
 * Checks to see if the error function is small enough to consider
 * the result to have converged.
 *
 * INPUT:
 * ------
 * int    r     - 1/2 the number of filter coeffiecients
 * int    Ext[] - Indexes to extremal frequencies [r+1]
 * double E[]   - Error function on the dense grid [gridsize]
 *
 * OUTPUT:
 * -------
 * Returns 1 if the result converged
 * Returns 0 if the result has not converged
 ********************/

short MainWindow::isDone(int r, int Ext[], long double E[]) {
        int i;
        long double min, max, current;

        min = max = fabs(E[Ext[0]]);
        for (i = 1; i <= r; i++) {
                try {
                        current = fabs(E[Ext[i]]);
                } catch(...) {
                        return 0;
                }
                if (current < min)
                        min = current;
                if (current > max)
                        max = current;
        }
        if (((max - min) / max) < 0.0001)
                return 1;
        return 0;
}

/********************
 * remez
 *=======
 * Calculates the optimal (in the Chebyshev/minimax sense)
 * FIR filter impulse response given a set of band edges,
 * the desired reponse on those bands, and the weight given to
 * the error in those bands.
 *
 * INPUT:
 * ------
 * int     numtaps     - Number of filter coefficients
 * int     numband     - Number of bands in filter specification
 * double  bands[]     - User-specified band edges [2 * numband]
 * double  des[]       - User-specified band responses [numband]
 * double  weight[]    - User-specified error weights [numband]
 * int     type        - Type of filter
 *
 * OUTPUT:
 * -------
 * double h[]      - Impulse response of final filter [numtaps]
 ********************/

void MainWindow::remez(long double h[], int numtaps,  int numband, long double bands[], long double des[], long double weight[],  int type) {
        long double* Grid, *W, *D, *E;
        int    i, iter, gridsize, r, *Ext;
        long  double* taps, c;
        long  double* x, *y, *ad;
        int    symmetry;

        if (type == BANDPASS)
                symmetry = POSITIVE;
        else
                symmetry = NEGATIVE;

        r = numtaps / 2;                /* number of extrema */
        if ((numtaps % 2) && (symmetry == POSITIVE))
                r++;

        /*
         * Predict dense grid size in advance for memory allocation
         *   .5 is so we round up, not truncate
         */
        gridsize = 0;
        for (i = 0; i < numband; i++) {
                gridsize += (int)(2 * r * GRIDDENSITY * (bands[2 * i + 1] - bands[2 * i]) + .5);
        }
        if (symmetry == NEGATIVE) {
                gridsize--;
        }

        /*
         * Dynamically allocate memory for arrays with proper sizes
         */
        Grid = ( long double*)malloc(gridsize * sizeof(long double));
        D = (long double*)malloc(gridsize * sizeof(long double));
        W = (long double*)malloc(gridsize * sizeof(long double));
        E = (long double*)malloc(gridsize * sizeof(long double));
        Ext = (int*)malloc((r + 1) * sizeof(int));
        taps = (long double*)malloc((r + 1) * sizeof(long double));
        x = (long double*)malloc((r + 1) * sizeof(long double));
        y = (long double*)malloc((r + 1) * sizeof(long double));
        ad = (long double*)malloc((r + 1) * sizeof(long double));

        /*
         * Create dense frequency grid
         */
        CreateDenseGrid(r, numtaps, numband, bands, des, weight,
                        &gridsize, Grid, D, W, symmetry);
        InitialGuess(r, Ext, gridsize);

        /*
         * For Differentiator: (fix grid)
         */
        if (type == DIFFERENTIATOR) {
                for (i = 0; i < gridsize; i++) {
                        /* D[i] = D[i]*Grid[i]; */
                        if (D[i] > 0.0001)
                                W[i] = W[i] / Grid[i];
                }
        }

        /*
         * For odd or Negative symmetry filters, alter the
         * D[] and W[] according to Parks McClellan
         */
        if (symmetry == POSITIVE) {
                if (numtaps % 2 == 0) {
                        for (i = 0; i < gridsize; i++) {
                                c = cosl(Pi * Grid[i]);
                                D[i] /= c;
                                W[i] *= c;
                        }
                }
        } else {
                if (numtaps % 2) {
                        for (i = 0; i < gridsize; i++) {
                                c = sinl(Pi2 * Grid[i]);
                                D[i] /= c;
                                W[i] *= c;
                        }
                } else {
                        for (i = 0; i < gridsize; i++) {
                                c = sinl(Pi * Grid[i]);
                                D[i] /= c;
                                W[i] *= c;
                        }
                }
        }

        /*
         * Perform the Remez Exchange algorithm
         */
        for (iter = 0; iter < MAXITERATIONS; iter++) {
                CalcParms(r, Ext, Grid, D, W, ad, x, y);
                CalcError(r, ad, x, y, gridsize, Grid, D, W, E);
                Search(r, Ext, gridsize, E);
                if (isDone(r, Ext, E))
                        break;
        }
        if (iter == MAXITERATIONS) {
                printf("Reached maximum iteration count.\nResults may be bad.\n");
        }

        CalcParms(r, Ext, Grid, D, W, ad, x, y);

        /*
         * Find the 'taps' of the filter for use with Frequency
         * Sampling.  If odd or Negative symmetry, fix the taps
         * according to Parks McClellan
         */
        for (i = 0; i <= numtaps / 2; i++) {
                if (symmetry == POSITIVE) {
                        if (numtaps % 2)
                                c = 1;
                        else
                                c = cosl(Pi * (long double)i / numtaps);
                } else {
                        if (numtaps % 2)
                                c = sinl(Pi2 * (long double)i / numtaps);
                        else
                                c = sinl(Pi * (long double)i / numtaps);
                }
                taps[i] = ComputeA((long double)i / numtaps, r, ad, x, y) * c;
        }

        /*
         * Frequency sampling design with calculated taps
         */
        FreqSample(numtaps, taps, h, symmetry);

        /*         * Delete allocated memory         */
        free(Grid);
        free(W);
        free(D);
        free(E);
        free(Ext);
        free(x);
        free(y);
        free(ad);
}


void MainWindow::on_pbSim_clicked() {
        int n;
        long double dSortie;
        long lFin;
        QString input;
        long double dFs, dFdeb, dFFin, dF, ddf;
        long double dAmpl;
        int iOrdre;
        long double* dFirCoeffP, *dFirCoeffS,*dFirCoeffR;
        long double* dEchantillon;
        long double dmax;
        QVector<double> S0, S1, S2;
        QVector<double> vdAx1, vdAx0;
        QFile*    mFilterFile;
        int iProgVal;
        QString qsTmpa, qsTmpb;
        long double dMult;
        bool bRet;
        int iNPts = 2000, m;
        Complex H,HR;

        iProgVal = 0;
        progressBar->setFormat("%p");
        // Prend les données
        dFs = ui->leFSamp->text() .toDouble();
        dFdeb = ui->leSttF->text().toDouble();
        dFFin = ui->leStpF->text().toDouble();
        if (ui->cbInt->isChecked()) { // Correct gain for int
                n = ui->leNbts->text().toInt();
                dMult = (long double) ((1 << (n - 1)) - 1);
        } else dMult = 1;


        dFirCoeffP = new long double[iM];
        dFirCoeffS = new long double[iM];
         dFirCoeffR = new long double[iM];

        // Charge le tableau de coeff principal
        mFilterFile = new QFile(qsFileName);
        mFilterFile->open(QIODevice::ReadOnly | QIODevice::Text);
        iOrdre = 0;
        while (!mFilterFile->atEnd()) {
                input = mFilterFile->readLine();
                dFirCoeffP[iOrdre] = input.toDouble(&bRet);
                if(bRet)  iOrdre++;
        }
        mFilterFile->close();

        // Charge le tableau de coeff Secondaire
        mFilterFile = new QFile(qsFileName + "T");
        mFilterFile->open(QIODevice::ReadOnly | QIODevice::Text);
        iOrdre = 0;
        while (!mFilterFile->atEnd()) {
                input = mFilterFile->readLine();
                dFirCoeffS[iOrdre] = input.toDouble(&bRet);
                if(bRet)  iOrdre++;
        }
        mFilterFile->close();

#ifdef PHHILBERT
        // Charge le tableau de coeff Hilbert Ref
        mFilterFile = new QFile(qsFileName + "R");  // Faire le meme fir normal et mettre R fin extension
        mFilterFile->open(QIODevice::ReadOnly | QIODevice::Text);
        iOrdre = 0;
        while (!mFilterFile->atEnd()) {
                input = mFilterFile->readLine();
                dFirCoeffR[iOrdre] = input.toDouble(&bRet);
                if(bRet)  iOrdre++;
        }
        mFilterFile->close();
#endif


        dEchantillon = new long  double[iOrdre];
        dF = dFdeb;
        dSortie = 0;
        ddf = (dFFin - dFdeb) / 250;    //250 nombre de pas


        dmax = 1e-29;

#define N2
#ifdef  N0
        //simulation par addition des sinusoides echantillné -> long

        progressBar->setRange(0, 250);
        while (dF < dFFin) {
                for (lFin = 0; lFin < iOrdre; lFin++) { // Charge le tableau d'echantillon
                        dEchantillon[lFin] = sinl(2.0 * M_PI * dF * (long double)lFin / dFs);
                }
                progressBar->setValue( iProgVal++);
                dAmpl = 1e-30;
                while (lFin < 1 * dFs) {
                        dSortie = 0;
                        for (n = 0; n < iOrdre; n++) { //boucle de filtrage
                                dSortie += dEchantillon[n] * dFirCoeff[n];
                        }
                        for (n = 1; n < iOrdre; n++) { //decale les echantillon
                                dEchantillon[n - 1] = dEchantillon[n];
                        }
                        dEchantillon[n - 1] = sinl(2.0 * M_PI * (long double)lFin * dF / dFs);
                        if (dSortie > dAmpl) dAmpl = dSortie;
                        lFin++;
                }
                if (dAmpl > dmax) dmax = dAmpl;
                if (dAmpl > 1e-29) {
                        dAmpl /= dMult;
                        dAmpl /= sqrt(2);
                        dAmpl = 20 * log10l(dAmpl);
                        vdAx0 << dF;
                        S0 << dAmpl;
                }

                dF += ddf;
        }
#endif
#ifdef N2
        // methode module complexe


        dFdeb /= dFs;
        dF = 0;
        progressBar->setRange(0, 2 * iNPts);
        //Trace principal
        m = 0;
        while ((m < iNPts) && ( dF * dFs < dFFin)) {
                dAmpl = 0;
                dF = m  / (2.0 * iNPts);
                dF += dFdeb;
                H = Complex(0, 0);
                for (n = 0; n < iOrdre; n++) {  // fonction transfer ordre pair
                        H +=  Complex (dFirCoeffP[n] * cosl(2.0 * M_PIl * dF * n), dFirCoeffP[n] * sinl( 2.0 * M_PIl * dF * n ));
                }
#ifdef PHHILBERT
                HR = Complex(0, 0);
                for (n = 0; n < iOrdre; n++) {  // fonction transfer ordre pair
                        HR +=  Complex (dFirCoeffR[n] * cosl(2.0 * M_PIl * dF * n), dFirCoeffR[n] * sinl( 2.0 * M_PIl * dF * n ));
                }
#endif
                progressBar->setValue( iProgVal++);
                vdAx0 << dF* dFs;
                dAmpl = H.mag();
                if (dAmpl > dmax) dmax = dAmpl;
                dAmpl = 20 * log10l(dAmpl);
                if (dAmpl<MINPLOTDB) dAmpl=MINPLOTDB;
                S0 << dAmpl ;
#ifdef PHHILBERT
                  dAmpl = H.arg()-HR.arg();
                  if (dAmpl<0) dAmpl+=2*M_PIl;
                 S2 << dAmpl*(180/M_PIl);
#endif
                m++;

        }

        //Trace Secondaire
        m = 0;
        dF = 0;
        while ((m < iNPts) && ( dF * dFs < dFFin)) {
                dAmpl = 0;
                dF = m  / (2.0 * iNPts);
                dF += dFdeb;
                H = Complex(0, 0);
                for (n = 0; n < iOrdre; n++) {  // fonction transfer ordre pair
                        H +=  Complex (dFirCoeffS[n] * cosl(2.0 * M_PIl * dF * n), dFirCoeffS[n] * sinl( 2.0 * M_PIl * dF * n ));
                }
                progressBar->setValue( iProgVal++);
                vdAx1 << dF* dFs;
                dAmpl = H.mag();
                if (dAmpl > dmax) dmax = dAmpl;
                dAmpl = 20 * log10l(dAmpl);
                  if (dAmpl<MINPLOTDB) dAmpl=MINPLOTDB;
                S1 << dAmpl ;
                m++;

        }
        progressBar->setValue( 2 * iNPts);

#endif

        //   ui->GraphView->clearGraphs();
        qsTmpa = "PNFilterNum : Simulation";
        ui->GraphView->clearGraphs();
        ui->GraphView->addGraph(ui->GraphView->xAxis, ui->GraphView->yAxis);
        ui->GraphView->graph(0)->setPen(QPen(Qt::blue)); // line color blue for first graph
        ui->GraphView->graph(0)->setBrush(QBrush(QColor(0, 0, 255, 20))); // first graph will be filled with translucent blue
        ui->GraphView->addGraph(ui->GraphView->xAxis, ui->GraphView->yAxis);
        ui->GraphView->graph(1)->setPen(QPen(Qt::gray)); // line color blue for first graph
        ui->GraphView->graph(1)->setBrush(QBrush(QColor(0, 0, 0, 20))); // first graph will be filled with translucent red

        // set title of plot:
        // ui->GraphView->plotLayout()->insertRow(0);
        // ui->GraphView->plotLayout()->addElement(0, 0, new QCPTextElement(    ui->GraphView, "Simulation Response", QFont("sans", 14, QFont::Bold)));
        ui->GraphView->xAxis2->setVisible(true);
        ui->GraphView->xAxis2->setTickLabels(true);
        ui->GraphView->yAxis2->setVisible(true);
        ui->GraphView->yAxis2->setTickLabels(true);
        // make left and bottom axes always transfer their ranges to right and top axes:
        connect(ui->GraphView->xAxis, SIGNAL(rangeChanged(QCPRange)), ui->GraphView->xAxis2, SLOT(setRange(QCPRange)));
        connect(ui->GraphView->yAxis, SIGNAL(rangeChanged(QCPRange)), ui->GraphView->yAxis2, SLOT(setRange(QCPRange)));
        // pass data points to graphs:
        ui->GraphView->graph(0)->setData(vdAx0, S0);
        ui->GraphView->graph(1)->setData(vdAx1, S1);
        ui->GraphView->graph(0)->rescaleAxes();
        ui->GraphView->graph(1)->rescaleAxes(true);
#ifdef PHHILBERT
        ui->GraphView->addGraph(ui->GraphView->xAxis, ui->GraphView->yAxis);
        ui->GraphView->graph(2)->setPen(QPen(Qt::red)); // line color
         ui->GraphView->graph(2)->setData(vdAx0, S2);
        ui->GraphView->graph(2)->rescaleAxes(true);
#endif
        ui->GraphView->yAxis->setTickLabelColor(Qt::black);
        ui->GraphView->yAxis2->setTickLabelColor(Qt::black);
        ui->GraphView->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
        ui->GraphView->xAxis->setLabel("Frequency [Hz]");
        ui->GraphView->yAxis->setLabel("Magnitude [dB]");

        dmax /= dMult;
        input.sprintf("Max Gain %.2Lf", dmax);
        // setup legend:
        // set graph name, will show up in legend next to icon:
        if (ui->rbBlack->isChecked()) {
                ui->GraphView->graph(0)->setName(QString("Blackman"));
                ui->GraphView->graph(1)->setName(QString("Hamming"));
        } else {
                ui->GraphView->graph(1)->setName(QString("Blackman"));
                ui->GraphView->graph(0)->setName(QString("Hamming"));
        }
        //   ui->GraphView->graph()->setName(qsTmpb);
        ui->GraphView->legend->setFont(QFont(font().family(), 12));
        ui->GraphView->legend->setIconSize(50, 20);
        ui->GraphView->legend->setVisible(true);
        ui->GraphView->axisRect()->insetLayout()->setInsetAlignment(0, Qt::AlignLeft | Qt::AlignBottom);


        progressBar->setFormat(input);

        setWindowTitle(qsTmpa);
        ui->GraphView->replot();
        //   textLabel->disconnect();
}

// Algorithme dspguide chap12
Complex*   MainWindow::FT( Complex* cX) {
        int k, i;
        Complex cS;
        Complex* Sortie = new Complex[ORDRE];

        for (k = 0; k < ORDRE; k++) {
                Sortie[k].im = 0;
                Sortie[k].re = 0;
        }

        for (k = 0; k < ORDRE; k++) {
                for (i = 0; i < ORDRE; i++) {
                        cS.re = cosl(2.0 * M_PIl * k * i / ORDRE);
                        cS.im = -sinl(2.0 * M_PIl * k * i / ORDRE);
                        Sortie[k] = Sortie[k] + cX[i] * cS;
                }
        }
        return Sortie;
}

// Algorithme dspguide chap12
// A verifier
Complex* MainWindow::FDFT(Complex*   cX) {
        int k, i, j, l, m;
        Complex cS, cT;
        //    Complex* Sortie = new Complex[ORDRE];
        Complex cU;
        int iLe, iLe2, jm1, ip;

        m = (int)(log(ORDRE)  / log(2.0));
        cT.re = 0;
        cT.im = 0;
        j = ORDRE / 2;
        for (i = 1; i < ORDRE - 1; i++) {
                if (i >= j) goto lab1;
                cT = cX[j];
                cX[j] = cX[i];
                cX[i] = cT;
lab1:
                k = ORDRE / 2;
lab3:
                if (k > j) goto lab2;
                j -= k;
                k /= 2;
                goto lab3;
lab2:
                j += k;
        }

        for (l = 1; l <= m; l++) {
                iLe = (int)pow(2, l);
                iLe2 = iLe / 2;
                cU.re = 1;
                cU.im = 0;
                cS.re = cosl(M_PIl / iLe2);
                cS.im = -sinl(M_PIl / iLe2);
                for (j = 1; j <= iLe2; j++) {
                        jm1 = j - 1;
                        for (i = jm1; i <= ORDRE - 1; i += iLe) {
                                ip = i + iLe2;
                                cT = cX[ip] * cU;
                                cX[ip] = cX[i] - cT;
                                cX[i] = cX[i] + cT;
                        }
                        cT = cU;
                        cU = cT * cS;
                }
        }

        //   for (l = 0; l < ORDRE; l++) Sortie[l] = cX[l];
        return cX;
}


// Demande le nom du fichier puis ecrit les coeffs
//  Rajoute tmp et ecrit le second set de coeffs ceci pour comparer les 2 fenetrages
// on a 2 fichiers celui avec tmp sera effacé a la fin de la simulation
// Fichier .coe pret pour vivado xilinx

void  MainWindow::WriteTheFile(long double* h, long double* ht) {
        QFile* mFilterFile;
        QFileDialog mFileDlg;
        //    long double dSum;
        QString input;
        long double dh;
        //   long double dGain;
        //   int iNbts;

        // dGain = ui->leGain->text().toDouble();
        //   iNbts = ui->leNbts->text().toDouble();

        if (  iHasTmpToDel == 1) {
                input = qsFileName + "T";
                mFilterFile = new QFile(input);
                mFilterFile->remove();
        }
        qsFileName = mFileDlg.getSaveFileName(this, QString("*.*"), QDir::currentPath() + "/" + "Coeff.coe", QString("Xilinx COE (*.coe);Text files (*.txt)"));
        mFilterFile = new QFile(qsFileName);
        mFilterFile->open(QIODevice::WriteOnly | QIODevice::Text);

        input.sprintf("Radix=10;\r\n");   // po fichier .coe xilinx
        mFilterFile->write(input.toLatin1());
        input.sprintf("CoefData=\r\n");
        mFilterFile->write(input.toLatin1());

        /*   dh=0;
         if (ui->cbInt->isChecked()) {   // formate po int
                    for (int compt1 = 1; compt1 < iM; compt1++) { //cherche max
                            if(fabsl(h[compt1])>dh)  dh=fabsl(h[compt1]);
                    }
                    dSum=powl(2.0,iNbts);
                    dSum/=2.0*dh;
         }*/

        // Fichier principal
        for (int compt1 = 1; compt1 < iM; compt1++) {
                dh = h[compt1] ;
                if (ui->cbInt->isChecked())  input.sprintf("%d", (int)(dh));      //   %04X  po 4 sgn hexa
                else
                        input.sprintf("%Lf", dh  );
                input += "\r\n";
                mFilterFile->write(input.toLatin1());
        }
        mFilterFile->close();

        // Write temp file same name +T at the end
        iHasTmpToDel = 1;
        input = qsFileName + "T";
        mFilterFile = new QFile(input);
        mFilterFile->open(QIODevice::WriteOnly | QIODevice::Text);

        input.sprintf("Radix=10;\r\n");   // po fichier .coe xilinx
        mFilterFile->write(input.toLatin1());
        input.sprintf("CoefData=\r\n");
        mFilterFile->write(input.toLatin1());

        for (int compt1 = 1; compt1 < iM; compt1++) {
                dh = ht[compt1] ;
                if (ui->cbInt->isChecked())  input.sprintf("%d", (int)dh);      //   %04X  po 4 sgn hexa
                else
                        input.sprintf("%Lf", dh  );
                input += "\r\n";
                mFilterFile->write(input.toLatin1());
        }
        mFilterFile->close();

}



void MainWindow::on_pbQuit_clicked()
{
    exit(0);
}
