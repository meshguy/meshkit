#ifndef SEARCH_UV_H
#define SEARCH_UV_H

#include <math.h>
#include <iostream>
#include <assert.h>

#include <vector>

using namespace std;

#include <iGeom.h>

class SearchUV
{
  public:
    int getUV(const iGeom_Instance &g, const iBase_EntityHandle &face,
              double x, double y, double z, double &u, double &v)
    {
        recursive_depth = 0;
        max_depth = 100;
        geom = g;
        gFace = face;

        xyzSearch[0] = x;
        xyzSearch[1] = y;
        xyzSearch[2] = z;

        int err;
        double umin, vmin, umax, vmax;
        iGeom_getEntUVRange(geom, gFace, &umin, &vmin, &umax, &vmax, &err);
        assert( !err );

        double urange[2], vrange[2];
        urange[0] = umin;
        urange[1] = umax;

        vrange[0] = vmin;
        vrange[1] = vmax;
        int stat = recursive_search(urange, vrange);

        u = uvFound[0];
        v = uvFound[1];

        assert( u >= umin && u <= umax );
        assert( v >= vmin && v <= vmax );

        return stat;
    }
private:
    int recursive_depth, max_depth;
    iGeom_Instance geom;
    iBase_EntityHandle gFace;
    double xyzSearch[3], uvFound[2];

    struct UVPoint
    {
        int index[2];
        double uv[2];
        double dist;

        bool operator<(const UVPoint & rhs) const
        {
            return dist < rhs.dist;
        }
    };

    vector<UVPoint> uvPoints;

    int recursive_search(double *urange, double *vrange)
    {
        double tol = 1.0E-06;
        int err;
        int numU = 10, numV = 10;
        double x, y, z;

        uvPoints.resize(numU * numV);

        double du = (urange[1] - urange[0]) / (double) numU;
        double dv = (vrange[1] - vrange[0]) / (double) numV;

        int index = 0;
        for (int j = 0; j < numV; j++)
        {
            double v = vrange[0] + j*dv;
            for (int i = 0; i < numU; i++)
            {
                double u = urange[0] + i*du;
                iGeom_getEntUVtoXYZ(geom, gFace, u, v, &x, &y, &z, &err);
                assert(!err);
                double dx = fabs(x - xyzSearch[0]);
                double dy = fabs(y - xyzSearch[1]);
                double dz = fabs(z - xyzSearch[2]);
                uvPoints[index].index[0] = i;
                uvPoints[index].index[1] = j;
                uvPoints[index].uv[0] = u;
                uvPoints[index].uv[1] = v;
                uvPoints[index].dist = dx*dx + dy*dy + dz*dz;
                index++;
            }
        }
        sort(uvPoints.begin(), uvPoints.end());
        double mindist = sqrt(uvPoints[0].dist);

        uvFound[0] = uvPoints[0].uv[0];
        uvFound[1] = uvPoints[0].uv[1];

        if (mindist < tol) return 0;

        if (recursive_depth++ == max_depth ) {
            cout << "Warning: Correct {U,V} pair not found : mindist  " << mindist << endl;
            return 1;
        }

        double u0, u1;
        if (uvPoints[0].index[0] == 0)
            u0 = uvPoints[0].uv[0];
        else
            u0 = uvPoints[0].uv[0] - du;

        if (uvPoints[0].index[0] == numU - 1)
            u1 = uvPoints[0].uv[0];
        else
            u1 = uvPoints[0].uv[0] + du;

        double v0, v1;
        if (uvPoints[0].index[1] == 0)
            v0 = uvPoints[0].uv[1];
        else
            v0 = uvPoints[0].uv[1] - dv;

        if (uvPoints[0].index[1] == numV - 1)
            v1 = uvPoints[0].uv[1];
        else
            v1 = uvPoints[0].uv[1] + dv;

        double new_urange[2], new_vrange[2];
        new_urange[0] = u0;
        new_urange[1] = u1;

        new_vrange[0] = v0;
        new_vrange[1] = v1;

        recursive_search(new_urange, new_vrange);
    }
};

#endif
