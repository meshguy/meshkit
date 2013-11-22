#include "Laplace.h"
#include <iostream>
#include <stdlib.h>

#include <assert.h>

int Laplacian::solve2d(int nx, int ny, const vector<double> &x, const vector<double> &y,
                     vector<double> &phi)

{
    int maxIterations = 1000;

    int numSize = nx*ny;

    assert(x.size() == numSize);
    assert(y.size() == numSize);
    assert(phi.size() == numSize);

    vector<double> pold(numSize);
    for (int i = 0; i < numSize; i++)
        pold[i] = phi[i];

    int c, l, r, t, b;
    double dxp, dxm, dyp, dym, dxpm, dypm;
    double cl, cr, ct, cb, cc;

    int iter = 0;

    double maxerror;

    while (1)
    {
        maxerror = 0.0;
        for (int j = 1; j < ny - 1; j++)
        {
            for (int i = 1; i < nx - 1; i++)
            {
                c = j * nx + i;
                l = c - 1;
                r = c + 1;
                t = c + nx;
                b = c - nx;

                dxp = fabs(x[r] - x[c]);
                dxm = fabs(x[l] - x[c]);
                dxpm = fabs(x[r] - x[l]);
                cl = 1.0 / (dxpm * dxm);
                cr = 1.0 / (dxpm * dxp);

                dyp = fabs(y[t] - y[c]);
                dym = fabs(y[b] - y[c]);
                dypm = fabs(y[t] - y[b]);
                ct = 1.0 / (dypm * dyp);
                cb = 1.0 / (dypm * dym);
                cc = 1.0 / (cl + cr + ct + cb);
                phi[c] = cc * (cl * pold[l] + cr * pold[r] + ct * pold[t] + cb * pold[b]);
                maxerror = max(maxerror, fabs(phi[c] - pold[c]));
            }
        }
        if (maxerror < 1.0E-06 || iter == maxIterations) break;

        for (int i = 0; i < numSize; i++)
            pold[i] = phi[i];

        iter++;
    }
//  cout << " Laplacian Solver " << iter << " Max Error : " << maxerror << endl;
}

//////////////////////////////////////////////////////////////////////////////////////

int Laplacian::solve2d(int nx, int ny, vector<double> &phi)
{
    int maxIterations = 100;

    int numSize = nx*ny;

    assert(phi.size() == numSize);

    vector<double> pold(numSize);
    for (int i = 0; i < numSize; i++)
        pold[i] = phi[i];

    int iter = 0;

    double maxerror;

    while (1)
    {
        maxerror = 0.0;
        for (int j = 1; j < ny - 1; j++)
        {
            for (int i = 1; i < nx - 1; i++)
            {
                int idc = j * nx + i;
                int ide = idc + 1;
                int idw = idc - 1;
                int idn = idc + nx;
                int ids = idc - nx;
                phi[idc] = (pold[ide] + pold[idw] +
                        pold[idn] + pold[ids]) / 4.0;
                maxerror = max(maxerror, fabs(phi[idc] - pold[idc]));
            }
        }
        if (maxerror < 1.0E-06 || iter == maxIterations) break;

        for (int i = 0; i < numSize; i++)
            pold[i] = phi[i];

        iter++;
    }
//    cout << " Laplacian Solver " << iter << " Max Error : " << maxerror << endl;
}

////////////////////////////////////////////////////////////////////////////////

int Laplacian::solve3d(int nx, int ny, int nz, vector<double> &phi)

{
    int maxIterations = 10;

    int numSize = nx * ny*nz;

    assert(phi.size() == numSize);

    vector<double> pold(numSize);
    for (int i = 0; i < numSize; i++)
        pold[i] = phi[i];

    int iter = 0;

    double maxerror;

    while (1)
    {
        maxerror = 0.0;
        for (int k = 1; k < nz - 1; k++)
        {
            for (int j = 1; j < ny - 1; j++)
            {
                for (int i = 1; i < nx - 1; i++)
                {
                    int idc = k * nx * ny + j * nx + i;
                    int ide = idc + 1;     // east
                    int idw = idc - 1;     // west
                    int idn = idc + nx;    // north
                    int ids = idc - nx;    // south
                    int idf = idc + nx*ny; // front
                    int idb = idc - nx*ny; // back
                    phi[idc] = (pold[ide] + pold[idw] +
                            pold[idn] + pold[ids] +
                            pold[idf] + pold[idb]) / 6.0;
                    maxerror = max(maxerror, fabs(phi[idc] - pold[idc]));
                }
            }
        }
        if (maxerror < 1.0E-06 || iter == maxIterations) break;

        for (int i = 0; i < numSize; i++)
            pold[i] = phi[i];

        iter++;
    }
//    cout << " Laplacian Solver " << iter << " Max Error : " << maxerror << endl;
}
////////////////////////////////////////////////////////////////////////////////

