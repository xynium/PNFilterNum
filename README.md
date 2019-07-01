# PNFilterNum

The purpose of this is to compute coefficients for FIR Filter.
It can also be used for Hilbert filter. 
The output formats for direct use by Xilinx Vivado, it is a .COE file.
Normally Vivado accept real coefficients and arrange for his format, but an integer format is proposed with variable bit number length.

Build: 

Code file are copies of my working directory, they were composed in QT Creator.

Change directory name src for PNFilterNum

Load QCustomPlot library.

Place the file in a directory named qcustomplot at the same level as the PNFilterNum source directory.

Creates a build directory at the same level as PNFilterNum directory.

Copy Makefile in build directory.

From build directory make.
