#!/usr/bin/python2.7

import climaxPyTest as cm
import numpy

import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt

def countCols(filename):

    f = open(filename, 'r');
    nline = 0
    checkCol = False;
    done = False;
    nColLinePrefix = 0;
    iCol = 0;

    colNo   = [];
    colName = [];
    colUnit = [];

    nCol = 0;

    for line in f:

        if not done and "//" in line:
            print 'Found a comment line'
            print line;
            if "Columns are" in line:
                print 'Found a Columns line'
                checkCol = True;

            if checkCol:
                if nColLinePrefix < 3:
                    print 'Incrememnting n = ', nColLinePrefix
                    nColLinePrefix = nColLinePrefix+1;
                else:
                    s = line.split()
                    if "Reduced chi-squared" in line:
                        done = True;
                    else:
                        nCol = nCol+1;
        else:
            continue;

    f.close();
    return nCol

def getCols(filename):
    nCol = countCols(filename);

    f = open(filename, 'r');
    nline = 0
    checkCol = False;
    done = False;
    nColLinePrefix = 0;
    iCol = 0;

    colNo   = [];
    colName = [];
    colUnit = [];

    nCol = 0;

    for line in f:

        if not done and "//" in line:
            print 'Found a comment line'
            print line;
            if "Columns are" in line:
                print 'Found a Columns line'
                checkCol = True;

            if checkCol:
                if nColLinePrefix < 3:
                    print 'Incrememnting n = ', nColLinePrefix
                    nColLinePrefix = nColLinePrefix+1;
                else:
                    s = line.split()
                    if "Reduced chi-squared" in line:
                        done = True;
                    else:
                        colNo.append(s[1]);
                        colName.append(s[2]);
                        colUnit.append(s[3]);
                        nCol = nCol+1;
        else:
            continue;

    f.close();
    return nCol, colNo, colName, colUnit

def countData(filename):

    f = open(filename, 'r');

    nData = 0;
    for line in f:
        if not "//" in line:
            nData = nData+1;
    return nData;

def getVals(filename):

    nCol  = countCols(filename);
    nData = countData(filename);

    f = open(filename, 'r');

    xoffvals = numpy.ndarray((nCol, nData), numpy.double);
    chisq    = numpy.ndarray((nData), numpy.double);
    lnlike   = numpy.ndarray((nData), numpy.double);
    mult     = numpy.ndarray((nData), numpy.double);

    iData = 0;
    print "here 0"
    for line in f:
        if not "//" in line:
            s = line.split();
            for iCol in range(0,nCol):
                xoffvals[iCol, iData] = s[iCol];
            chisq[iData] = s[nCol];
            lnlike[iData] = s[nCol+1];
            mult[iData] = s[nCol+2];
            iData = iData+1;
    return xoffvals, chisq, lnlike, mult


def mantzhistplot(filename, ix, iy, nx,ny, xscale, yscale):
  v = getVals(filename);
  H, xedges, yedges = numpy.histogram2d(v[ix]*xscale, v[iy]*yscale, bins=(nx,ny));
#  extent = [yedges[0], yedges[-1], xedges[-1], xedges[0]]
#  plt.contour(H.transpose(), extent=extent);
  extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
  plt.contour(H.transpose(), extent=extent);
  plt.show()
  return H,xedges,yedges

def histplot2d(filename, ix, iy, nx,ny):
  v = getVals(filename);
  n,no,name,unit = getCols(filename);
  H, xedges, yedges = numpy.histogram2d(v[ix], v[iy], bins=(nx,ny));
#  extent = [yedges[0], yedges[-1], xedges[-1], xedges[0]]
#  plt.contour(H.transpose(), extent=extent);
  extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
  plt.contour(H.transpose(), extent=extent);
  plt.xlabel(name[ix]);
  plt.ylabel(name[iy]);
  plt.show()
  return H,xedges,yedges

def histplot(filename, ix, nx):
  v = getVals(filename);
  n,no,name,unit = getCols(filename);
  plt.hist(v[ix], nx);
  plt.xlabel(name[ix]);
  plt.show()
  return H,xedges,yedges

def compplot(ix, iy, nx,ny):
  plt.clf();
  n,no,name,unit = getCols('arnaud.out');
  histplot2d('arnaud.out', ix, iy, nx, ny);
  plt.hold();
  histplot2d('planck.out', ix, iy, nx, ny);
  histplot2d('sayers.out', ix, iy, nx, ny);
  plt.xlabel(name[ix]);
  plt.ylabel(name[iy]);

def mantzcompplot(ix, iy, nx,ny):
  plt.clf();
  n,no,name,unit = getCols('arnaud.out');
  mantzhistplot('arnaud.out', ix, iy, nx, ny, 1.17, 1.0);
  plt.hold();
  mantzhistplot('planck.out', ix, iy, nx, ny, 1.81, 1.0);
  mantzhistplot('sayers.out', ix, iy, nx, ny, 1.18, 1.0);
  plt.xlabel(name[ix]);
  plt.ylabel(name[iy]);

def plot3dline(file, iX, iY, ax):
    d,chisq,lnlike,mult = getVals(file);
#    fig = plt.figure()
#    ax = fig.gca(projection='3d')
    x = d[iX];
    y = d[iY];
    z = lnlike;
    ax.plot(x,y,z,color='yellow');
    plt.show();
