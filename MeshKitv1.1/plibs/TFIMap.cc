#include "TFIMap.h"

//****************************************************************************//
// Acknowlegement:  The original Code "Blend" from John Burkardt ( Under GNU
//                  lincense) have been significantly modified. Anyhow full
//                  credit is given to the original developers.
//
//                  Document is also copied from the original code.
//
// Modification Date :  7th July 2009;
// Modified By       :  Chaman Singh Verma
//                      Argonne National Lab, Argonne, IL, USA.
// Major Changes:
//                 (1) Parameters range (-1,1) instead of (0,1).
//                 (2) Probably better function names
//                 (3) Better modularity and reuse
//                 (4) Changed ordering (k,j,i) instead of (i,j,k)
//                 (5) Guass point introduction.
//                 (6) Probably much easier to understand than the original code.
//                 (7) Probably better input parameter passing order.
//                 (8) Error checking using assertions introduced.
// Downside:
//                 (1) May be little more expensive ( than the original code) because of
//                     function calling overheads. 
//
//                 (2) Redundant calculations to make code easier to understand.
//     
////////////////////////////////////////////////////////////////////////////////

#include "TFIMap.h"

using namespace std;

class SpectralElements;

//****************************************************************************80

double TFIMap::linear_interpolation(double r, double x0, double x1)
{

    //************************************************************************80
    //
    //  Purpose: extends scalar data at endpoints to a line.
    //
    //  Diagram:
    //
    //    0-----r-----1
    //
    //  Author:
    //
    //    John Burkardt: Modified by Chaman Singh Verma
    //
    //  Parameters:
    //
    //    Input, double R, the coordinate where an interpolated value is desired.
    //
    //    Input, double X0, X1, the data values at the ends of the line.
    //
    //    Output, double *X, the interpolated data value at (R).
    //
    ///////////////////////////////////////////////////////////////////////////

    assert(r >= -1.0 && r <= 1.0);

    double val = (1.0 - r) * x0 + (1.0 + r) * x1;

    val *= 0.5;

    return val;
}

//****************************************************************************80

double TFIMap::bilinear_interpolation(double r, double s, double *valCorners)
{
    assert(r >= -1.0 && r <= 1.0);
    assert(s >= -1.0 && s <= 1.0);

    double val = (1.0 - r) * (1.0 - s) * valCorners[0] +
            (1.0 + r) * (1.0 - s) * valCorners[1] +
            (1.0 + r) * (1.0 + s) * valCorners[2] +
            (1.0 - r) * (1.0 + s) * valCorners[3];

    val *= 0.25;
    return val;
}

//****************************************************************************80

double TFIMap::bilinear_interpolation(double r, double s,
                                      double x00, double x10, double x11, double x01)
{
    double valCorners[4];

    valCorners[0] = x00;
    valCorners[1] = x10;
    valCorners[2] = x11;
    valCorners[3] = x01;

    double val = bilinear_interpolation(r, s, valCorners);

    return val;
}

//****************************************************************************80

double TFIMap::trilinear_interpolation(double r, double s, double t, double *valCorners)
{
    assert(r >= -1.0 && r <= 1.0);
    assert(s >= -1.0 && s <= 1.0);
    assert(t >= -1.0 && t <= 1.0);

    double val = (1.0 - r) * (1.0 - s) * (1.0 - t) * valCorners[0] +
            (1.0 + r) * (1.0 - s) * (1.0 - t) * valCorners[1] +
            (1.0 + r) * (1.0 + s) * (1.0 - t) * valCorners[2] +
            (1.0 - r) * (1.0 + s) * (1.0 - t) * valCorners[3] +
            (1.0 - r) * (1.0 - s) * (1.0 + t) * valCorners[4] +
            (1.0 + r) * (1.0 - s) * (1.0 + t) * valCorners[5] +
            (1.0 + r) * (1.0 + s) * (1.0 + t) * valCorners[6] +
            (1.0 - r) * (1.0 + s) * (1.0 + t) * valCorners[7];

    val *= 0.125;

    return val;
}

//****************************************************************************80

double TFIMap::trilinear_interpolation(double r, double s, double t,
                                       double x000, double x100, double x110, double x010,
                                       double x001, double x101, double x111, double x011)
{
    double valCorners[8];

    valCorners[0] = x000;
    valCorners[1] = x100;
    valCorners[2] = x110;
    valCorners[3] = x010;

    valCorners[4] = x001;
    valCorners[5] = x101;
    valCorners[6] = x111;
    valCorners[7] = x011;

    double val = trilinear_interpolation(r, s, t, valCorners);

    return val;
}

//****************************************************************************80

double TFIMap::transfinite_blend(double r, double s,
                                 double x00, double x10, double x11, double x01,
                                 double xr0, double x1s, double xr1, double x0s)
{

    //****************************************************************************80
    //
    //  Purpose:
    //
    //    BLEND_1D1 extends scalar data along the boundary into a square.
    //
    //  Diagram:
    //
    //    01-----r1-----11
    //     |      .      |
    //     |      .      |
    //    0s.....rs.....1s
    //     |      .      |
    //     |      .      |
    //    00-----r0-----10
    //
    //  Formula:
    //
    //    Written as a polynomial in R and S, the interpolation map has the form
    //
    //      X(R,S) =
    //           1     * ( x0s + xr0 - x00 )
    //         + r     * ( x00 + x1s - x0s - x10 )
    //         + s     * ( x00 + xr0 - x01 - xr1 )
    //         + r * s * ( x01 + x10 - x00 - x11 )
    //
    //    The nonlinear term ( r * s ) has an important role:
    //
    //    If ( x01 + x10 - x00 - x11 ) is zero, then the input data lies in a plane,
    //    and the mapping is affine.  All the interpolated data will lie
    //    on the plane defined by the four corner values.  In particular,
    //    on any line through the square, data values at intermediate points
    //    will lie between the values at the endpoints.
    //
    //    If ( x01 + x10 - x00 - x11 ) is not zero, then the input data does
    //    not lie in a plane, and the interpolation map is nonlinear.  On
    //    any line through the square, data values at intermediate points
    //    may lie above or below the data values at the endpoints.  The
    //    size of the coefficient of r * s will determine how severe this
    //    effect is.
    //
    //  Author:
    //
    //    John Burkardt: Modified by Chaman Singh Verma
    //
    //  Reference:
    //
    //    William Gordon, Charles A Hall,
    //    Construction of Curvilinear Coordinate Systems and Application to
    //    Mesh Generation,
    //    International Journal of Numerical Methods in Engineering,
    //    Volume 7, pages 461-477, 1973.
    //
    //    Joe Thompson, Bharat Soni, Nigel Weatherill,
    //    Handbook of Grid Generation,
    //    CRC Press,
    //    1999.
    //
    //  Parameters:
    //
    //    Input, double R, S, the coordinates where an interpolated value is desired.
    //
    //    Input, double X00, X01, X10, X11, the data values at the corners.
    //
    //    Input, double XR0, XR1, X0S, X1S, the data values at points along the edges.
    //
    //    Output, double *X, the interpolated data value at (R,S).
    //
    ////////////////////////////////////////////////////////////////////////////

    double u = linear_interpolation(r, x0s, x1s);
    double v = linear_interpolation(s, xr0, xr1);
    double uv = bilinear_interpolation(r, s, x00, x10, x11, x01);

    double val = u + v - uv;

    return val;
}
//****************************************************************************80

double TFIMap::transfinite_blend(double r, double s, double t,
                                 double x000, double xr00, double x100,
                                 double x0s0, double xrs0, double x1s0,
                                 double x010, double xr10, double x110,
                                 double x00t, double xr0t, double x10t,
                                 double x0st, double x1st,
                                 double x01t, double xr1t, double x11t,
                                 double x001, double xr01, double x101,
                                 double x0s1, double xrs1, double x1s1,
                                 double x011, double xr11, double x111)
{

    //************************************************************************80
    //
    //  Purpose:  extends scalar data along the surface into a cube.
    //
    //  Diagram:
    //
    //    010-----r10-----110    011-----r11-----111
    //      |       .       |      |       .       |
    //      |       .       |      |       .       |
    //    0s0.....rs0.....1s0    0s1.....rs1.....1s1     S
    //      |       .       |      |       .       |     |
    //      |       .       |      |       .       |     |
    //    000-----r00-----100    001-----r01-----101     +----R
    //       BACK                     FRONT
    //
    //    001-----0s1-----011    101-----1s1-----111
    //      |       .       |      |       .       |
    //      |       .       |      |       .       |
    //    00t.....0st.....01t    10t.....1st.....11t     T
    //      |       .       |      |       .       |     |
    //      |       .       |      |       .       |     |
    //    000-----0s0-----010    100-----1s0-----110     +----S
    //       LEFT                       RIGHT

    //    001-----r01-----101    011-----r11-----111
    //      |       .       |      |       .       |
    //      |       .       |      |       .       |
    //    00t.....r0t.....10t    01t.....r1t.....11t     T
    //      |       .       |      |       .       |     |
    //      |       .       |      |       .       |     |
    //    000-----r0t-----100    010-----r10-----110     +----R
    //       BOTTOM                       TOP
    //
    //  Author:
    //
    //    John Burkardt
    //
    //  Reference:
    //
    //    William Gordon, Charles A Hall,
    //    Construction of Curvilinear Coordinate Systems and Application to
    //    Mesh Generation,
    //    International Journal of Numerical Methods in Engineering,
    //    Volume 7, pages 461-477, 1973.
    //
    //    Joe Thompson, Bharat Soni, Nigel Weatherill,
    //    Handbook of Grid Generation,
    //    CRC Press,
    //    1999.
    //
    //  Parameters:
    //
    //    Input, double R, S, T, the coordinates where an interpolated value is desired.
    //
    //    Input, double X000, X001, X010, X011, X100, X101, X110, X111, the data
    //    values at the corners.
    //
    //    Input, double XR00, XR01, XR10, XR11, X0S0, X0S1, X1S0, X1S1, X00T, X01T,
    //    X10T, X11T, the data values at points along the edges.
    //
    //    Input, double X0ST, X1ST, XR0T, XR1T, XRS0, XRS1, the data values
    //    at points on the faces.
    //
    //    Output, double *X, the interpolated data value at (R,S,T).
    //
    ////////////////////////////////////////////////////////////////////////////

    double u = linear_interpolation(r, x0st, x1st);
    double v = linear_interpolation(s, xr0t, xr1t);
    double w = linear_interpolation(t, xrs0, xrs1);

    double uv = bilinear_interpolation(r, s, x00t, x10t, x11t, x01t);
    double uw = bilinear_interpolation(r, t, x0s0, x1s0, x1s1, x0s1);
    double vw = bilinear_interpolation(s, t, xr00, xr10, xr11, xr01);

    double uvw = trilinear_interpolation(r, s, t,
                                         x000, x100, x110, x010,
                                         x001, x101, x111, x011);

    double val = u + v + w - uw - uv - vw + uvw;

    return val;
}

//****************************************************************************80

void TFIMap::blend_from_corners(vector<double> &x, const vector<double> &glxnodes)
{

    //****************************************************************************80
    //
    //  Purpose: extends indexed scalar data at endpoints along a line.
    //
    //  Diagram:
    //
    //    ( X0, ..., ..., ..., ..., ..., X6 )
    //
    //  Author:
    //
    //    John Burkardt
    //
    //  Reference:
    //
    //    William Gordon, Charles A Hall,
    //    Construction of Curvilinear Coordinate Systems and Application to
    //    Mesh Generation,
    //    International Journal of Numerical Methods in Engineering,
    //    Volume 7, pages 461-477, 1973.
    //
    //    Joe Thompson, Bharat Soni, Nigel Weatherill,
    //    Handbook of Grid Generation,
    //    CRC Press,
    //    1999.
    //
    //  Parameters:
    //
    //    Input/output, double X[M].
    //
    //    On input, X[0] and X[M-1] contain scalar values which are to be
    //    interpolated through the entries X[1] through X[M-2].  It is assumed
    //    that the dependence of the data is linear in the vector index I.
    //
    //    On output, X[1] through X[M-2] have been assigned interpolated values.
    //
    //    Input, int M, the number of entries in X.
    //
    int i;
    double r;

    int n = glxnodes.size();
    assert(x.size() == n);

    for (i = 1; i < n - 1; i++)
    {
        r = glxnodes[i];
        x[i] = linear_interpolation(r, x[0], x[n - 1]);
    }

    return;
}
//****************************************************************************80

void TFIMap::blend_from_edges(vector<double> &x,
                              const vector<double> &glxnodes,
                              const vector<double> &glynodes)
{

    //****************************************************************************80
    //
    //  Purpose: extends indexed scalar data along edges into a table.
    //
    //  Diagram:
    //
    //    ( X00,  X01,  X02,  X03,  X04,  X05,  X06 )
    //    ( X10,  ...,  ...,  ...,  ...,  ...,  X16 )
    //    ( X20,  ...,  ...,  ...,  ...,  ...,  X26 )
    //    ( X30,  ...,  ...,  ...,  ...,  ...,  X36 )
    //    ( X40,  X41,  X42,  X43,  X44,  X45,  X46 )
    //
    //  Licensing:
    //
    //    This code is distributed under the GNU LGPL license.
    //
    //  Modified:
    //
    //    7th July December 2009
    //
    //  Author:
    //
    //    John Burkardt: Modified by Chaman Singh Verma
    //
    //  Reference:
    //
    //    William Gordon, Charles A Hall,
    //    Construction of Curvilinear Coordinate Systems and Application to
    //    Mesh Generation,
    //    International Journal of Numerical Methods in Engineering,
    //    Volume 7, pages 461-477, 1973.
    //
    //    Joe Thompson, Bharat Soni, Nigel Weatherill,
    //    Handbook of Grid Generation,
    //    CRC Press,
    //    1999.
    //
    //  Parameters:
    //
    //    Input/output, double X[NX*NY], a singly dimensioned array that
    //    is "really" doubly dimensioned.  The double dimension index [I][J]
    //    corresponds to the single dimension index J * NX + I.
    //
    //    On input, data is contained in the "edge entries"
    //    X[0][0], X[I][0], X[0][NY-1] and X[NX-1][NY-1],
    //    for I = 0 to NX-1, and J = 0 to NY-1.
    //
    //    On output, all entries in X have been assigned a value.
    //
    //    Input, int NX, NY  the number of rows and columns in X.
    //
    ///////////////////////////////////////////////////////////////////////////////

    int nx = glxnodes.size();
    int ny = glynodes.size();

    assert(x.size() == nx * ny);

    double x0 = x[0];
    double x1 = x[nx - 1];
    double x2 = x[nx * ny - 1];
    double x3 = x[nx * ny - nx];

    for (int j = 1; j < ny - 1; j++)
    {
        double s = glynodes[j];

        for (int i = 1; i < nx - 1; i++)
        {
            double r = glxnodes[i];
            int il = j*nx;
            int ir = j * nx + nx - 1;
            int it = (ny - 1) * nx + i;
            int ib = i;
            int offset = j * nx + i;

            double xr0 = x[ib];
            double x1s = x[ir];
            double xr1 = x[it];
            double x0s = x[il];

            x[offset] = transfinite_blend(r, s, x0, x1, x2, x3, xr0, x1s, xr1, x0s);
        }
    }
    return;
}

////////////////////////////////////////////////////////////////////////////////

void TFIMap::blend_from_corners(vector<double> &x,
                                const vector<double> &glxnodes,
                                const vector<double> &glynodes)
{

    //****************************************************************************80
    //
    //  Purpose: extends indexed scalar data at corners into a table.
    //
    //  Diagram:
    //
    //    ( X00,  ..., ..., ..., ..., ..., X06 )
    //    ( ...,  ..., ..., ..., ..., ..., ... )
    //    ( ...,  ..., ..., ..., ..., ..., ... )
    //    ( ...,  ..., ..., ..., ..., ..., ... )
    //    ( X40,  ..., ..., ..., ..., ..., X46 )
    //
    //  Licensing:
    //
    //    This code is distributed under the GNU LGPL license.
    //
    //  Modified:
    //
    //    19 December 1998
    //
    //  Author:
    //
    //    John Burkardt
    //
    //  Reference:
    //
    //    William Gordon, Charles A Hall,
    //    Construction of Curvilinear Coordinate Systems and Application to
    //    Mesh Generation,
    //    International Journal of Numerical Methods in Engineering,
    //    Volume 7, pages 461-477, 1973.
    //
    //    Joe Thompson, Bharat Soni, Nigel Weatherill,
    //    Handbook of Grid Generation,
    //    CRC Press,
    //    1999.
    //
    //  Parameters:
    //
    //    Input/output, double X[NX*NY], a singly dimensioned array that
    //    is "really" doubly dimensioned.  The double dimension index [I][J]
    //    corresponds to the single dimension index J * NX + I.
    //
    //    On input, data values have been stored in the entries
    //    [0], [NX-1], [NX*NY-NX] and [NX*NY-1], which correspond to the double
    //    dimension entries [0][0], [1][NX-1], [0][NY-1] and [NX-1][NY-1].
    //
    //    On output, all entries in X have been assigned a value.
    //
    //    Input, int NX, NY, the number of rows and columns in the doubly
    //    dimensioned data.
    //
    ///////////////////////////////////////////////////////////////////////////////
    int nx = glxnodes.size();
    int ny = glynodes.size();

    int offset;
    double r, s, x0, x1;

    int nxy = nx*ny;

    for (int i = 1; i < nx - 1; i++)
    {
        r = glxnodes[i];

        x0 = x[0];
        x1 = x[nx - 1];
        offset = i;
        x[offset] = linear_interpolation(r, x0, x1);

        x0 = x[nxy - nx];
        x1 = x[nxy - 1 ];
        offset = (ny - 1) * nx + i;
        x[offset] = linear_interpolation(r, x0, x1);
    }


    for (int j = 1; j < ny - 1; j++)
    {
        s = glynodes[j];

        x0 = x[0];
        x1 = x[nx * ny - nx];
        offset = j*nx;
        x[offset] = linear_interpolation(s, x0, x1);

        x0 = x[nx - 1];
        x1 = x[nxy - 1];
        offset = j * nx + (nx - 1);
        x[offset] = linear_interpolation(s, x0, x1);
    }

    blend_from_edges(x, glxnodes, glynodes);

    return;
}


//****************************************************************************80

void TFIMap::blend_from_faces(vector<double> &x,
                              const vector<double> &glxnodes,
                              const vector<double> &glynodes,
                              const vector<double> &glznodes)
{

    //************************************************************************80
    //
    //  Purpose:
    //
    //    BLEND_FROM_FACES extends indexed scalar data along faces into a cubic table.
    //
    //  Diagram:
    //
    //    ( X000    X010    X020    X030    X040    X050 )
    //    ( X100    X110    X120    X130    X140    X150 )
    //    ( X200    X210    X220    X230    X240    X250 )   Layer 1
    //    ( X300    X310    X320    X330    X340    X350 )
    //    ( X400    X410    X420    X430    X440    X450 )
    //
    //    ( X001    X011    X021    X031    X041    X051 )
    //    ( X101    ...     ....    ....    ....    X151 )
    //    ( X201    ...     ....    ....    ....    X251 )   Layer K
    //    ( X301    ...     ....    ....    ....    X351 )   1 < K < M3
    //    ( X401    X411    X421    X431    X441    X451 )
    //
    //    ( X002    X012    X022    X032    X042    X052 )
    //    ( X102    X112    X122    X132    X142    X152 )
    //    ( X202    X212    X222    X232    X242    X252 )   Layer M3
    //    ( X302    X312    X322    X332    X342    X352 )
    //    ( X402    X412    X422    X432    X442    X452 )
    //
    //  Licensing:
    //
    //    This code is distributed under the GNU LGPL license.
    //
    //  Modified:
    //
    //    22 December 1998
    //
    //  Author:
    //
    //    John Burkardt
    //
    //  Reference:
    //
    //    William Gordon, Charles A Hall,
    //    Construction of Curvilinear Coordinate Systems and Application to
    //    Mesh Generation,
    //    International Journal of Numerical Methods in Engineering,
    //    Volume 7, pages 461-477, 1973.
    //
    //    Joe Thompson, Bharat Soni, Nigel Weatherill,
    //    Handbook of Grid Generation,
    //    CRC Press,
    //    1999.
    //
    //  Parameters:
    //
    //    Input/output, double X[NX*NY*NZ], a singly dimensioned array that
    //    is "really" triply dimensioned.  The triple dimension index
    //    [I][J][K] corresponds to the single dimension index
    //    K * NX*NY + J * NX + I
    //
    //    On input, there is already scalar data in the entries X[I][J][K]
    //    corresponding to "faces" of the table, that is, entries for which
    //    at least one of the three indices I, J and K is equal to their
    //    minimum or maximum possible values.
    //
    //    On output, all entries in X have been assigned a value, using the
    //    table indices as independent variables.
    //
    //    Input, int NX, NY, NY, the number of rows, columns, and layers in X.
    //
    ////////////////////////////////////////////////////////////////////////////
    int nx = glxnodes.size();
    int ny = glynodes.size();
    int nz = glznodes.size();

    double r, s, t;

    int offset, nxy = nx*ny;

    for (int k = 1; k < nz - 1; k++)
    {
        t = glznodes[k];
        for (int j = 1; j < ny - 1; j++)
        {
            s = glynodes[j];
            for (int i = 1; i < nx - 1; i++)
            {
                r = glxnodes[i];

                // Points on Back Plane
                double x000 = x[0];
                double xr00 = x[i];
                double x100 = x[(nx - 1)];

                double x0s0 = x[j * nx];
                double xrs0 = x[i + j * nx];
                double x1s0 = x[(nx - 1) + j * nx];

                double x010 = x[(ny - 1) * nx];
                double xr10 = x[i + (ny - 1) * nx];
                double x110 = x[(ny - 1) * nx + (nx - 1)];

                // Intermediate Plane

                double x00t = x[k * nxy];
                double xr0t = x[i + k * nxy];
                double x10t = x[(nx - 1) + k * nxy];

                double x0st = x[j * nx + k * nxy];
                double x1st = x[(nx - 1) + j * nx + k * nxy];

                double x01t = x[ (ny - 1) * nx + k * nxy];
                double xr1t = x[i + (ny - 1) * nx + k * nxy];
                double x11t = x[(nx - 1) + (ny - 1) * nx + k * nxy];

                // Front Plane
                double x001 = x[(nz - 1) * nxy];
                double xr01 = x[ i + (nz - 1) * nxy];
                double x101 = x[(nz - 1) * nxy + (nx - 1)];

                double x0s1 = x[ j * nx + (nz - 1) * nxy];
                double xrs1 = x[ i + j * nx + (nz - 1) * nxy];
                double x1s1 = x[ (nx - 1) + j * nx + (nz - 1) * nxy];

                double x011 = x[(ny - 1) * nx + (nz - 1) * nxy];
                double xr11 = x[ i + (ny - 1) * nx + (nz - 1) * nxy];
                double x111 = x[(nz - 1) * nxy + (ny - 1) * nx + (nx - 1)];

                offset = k * nxy + j * nx + i;
                x[offset] = transfinite_blend(r, s, t,
                                              x000, xr00, x100,
                                              x0s0, xrs0, x1s0,
                                              x010, xr10, x110,
                                              x00t, xr0t, x10t,
                                              x0st, x1st,
                                              x01t, xr1t, x11t,
                                              x001, xr01, x101,
                                              x0s1, xrs1, x1s1,
                                              x011, xr11, x111);
            }

        }

    }

    return;
}

////////////////////////////////////////////////////////////////////////////////

void TFIMap::blend_from_edges(vector<double> &x,
                              const vector<double> &glxnodes,
                              const vector<double> &glynodes,
                              const vector<double> &glznodes)
{

    //************************************************************************80
    //
    //  Purpose: extends indexed scalar data along "edges" into a cubic table.
    //
    //  Diagram:
    //
    //    ( X000,   X010,   X020,   X030,   X040,   X050 )
    //    ( X100,   ...,    ...,    ...,    ...,    X150 )
    //    ( X200,   ...,    ...,    ...,    ...,    X250 )   Layer 1
    //    ( X300,   ...,    ...,    ...,    ...,    X350 )
    //    ( X400,   X410,   X420,   X430,   X440,   X450 )
    //
    //    ( X001,   ...,    ...,    ...,    ...,    X051 )
    //    ( ....,   ...,    ...,    ...,    ...,    ...  )
    //    ( ....,   ...,    ...,    ...,    ...,    ...  )   Layer K
    //    ( ....,   ...,    ...,    ...,    ...,    ...  )   1 < K < M3
    //    ( X401,   ...,    ...,    ...,    ...,    X451 )
    //
    //    ( X002,   X012,   X022,   X032,   X042,   X052 )
    //    ( X102,   ...,    ...,    ...,    ...,    X152 )
    //    ( X202,   ...,    ...,    ...,    ...,    X252 )   Layer M3
    //    ( X302    ...,    ...,    ...,    ...,    X352 )
    //    ( X402,   X412,   X422,   X432,   X442,   X452 )
    //
    //  Licensing:
    //
    //    This code is distributed under the GNU LGPL license.
    //
    //  Modified:
    //
    //    22 December 1998
    //
    //  Author:
    //
    //    John Burkardt
    //
    //  Reference:
    //
    //    William Gordon, Charles A Hall,
    //    Construction of Curvilinear Coordinate Systems and Application to
    //    Mesh Generation,
    //    International Journal of Numerical Methods in Engineering,
    //    Volume 7, pages 461-477, 1973.
    //
    //    Joe Thompson, Bharat Soni, Nigel Weatherill,
    //    Handbook of Grid Generation,
    //    CRC Press,
    //    1999.
    //
    //  Parameters:
    //
    //    Input/output, double X[NX*NY*NZ], a singly dimensioned array that
    //    is "really" triply dimensioned.  The triple dimension index
    //    [I][J][K] corresponds to the single dimension index
    //    K * NX*NY + J * NX + I
    //
    //    On input, there is already scalar data in the entries X[I][J][K]
    //    corresponding to "edges" of the table, that is, entries for which
    //    at least two of the three indices I, J and K are equal to their
    //    minimum or maximum possible values.
    //
    //    Input, int NX, NY, NZ, the number of rows, columns, and layers in X.
    //
    ////////////////////////////////////////////////////////////////////////////

    int offset;

    int nx = glxnodes.size();
    int ny = glynodes.size();
    int nz = glznodes.size();

    int nxy = nx*ny;

    double r, s, t;
    double x00, x10, x11, x01;
    double xs0, x1t, xs1, x0t;

    for (int k = 1; k < nz - 1; k++)
    {
        t = glznodes[k];
        for (int j = 1; j < ny - 1; j++)
        {
            s = glynodes[j];
            //Left Face ...
            x00 = x[0];
            x10 = x[(ny - 1) * nx];
            x11 = x[(nz - 1) * nxy + (ny - 1) * nx];
            x01 = x[(nz - 1) * nxy];

            xs0 = x[j * nx];
            x1t = x[k * nxy + (ny - 1) * nx];
            xs1 = x[(nz - 1) * nxy + j * nx];
            x0t = x[k * nxy];

            offset = k * nxy + j*nx;

            x[offset] = transfinite_blend(s, t, x00, x10, x11, x01, xs0, x1t, xs1, x0t);

            // Right Face ...
            x00 = x[nx - 1];
            x10 = x[(ny - 1) * nx + (nx - 1)];
            x11 = x[(nz - 1) * nxy + (ny - 1) * nx + (nx - 1)];
            x01 = x[(nz - 1) * nxy + (nx - 1)];

            xs0 = x[j * nx + (nx - 1)];
            x1t = x[k * nxy + (ny - 1) * nx + (nx - 1)];
            xs1 = x[(nz - 1) * nxy + j * nx + (nx - 1)];
            x0t = x[k * nxy + (nx - 1)];

            offset = k * nxy + j * nx + nx - 1;

            x[offset] = transfinite_blend(s, t, x00, x10, x11, x01, xs0, x1t, xs1, x0t);
        }
    }

    double xr0, xr1;

    for (int k = 1; k < nz - 1; k++)
    {
        t = glznodes[k];
        for (int i = 1; i < nx - 1; i++)
        {
            r = glxnodes[i];

            // Bottom Face ...
            x00 = x[0];
            x10 = x[nx - 1];
            x11 = x[(nz - 1) * nxy + (nx - 1)];
            x01 = x[(nz - 1) * nxy];

            xr0 = x[i];
            x1t = x[k * nxy + (nx - 1)];

            xr1 = x[(nz - 1) * nxy + i];
            x0t = x[k * nxy];

            offset = k * nxy + i;

            x[offset] = transfinite_blend(r, t, x00, x10, x11, x01, xr0, x1t, xr1, x0t);

            // Top Face ...
            x00 = x[(ny - 1) * nx];
            x10 = x[(ny - 1) * nx + nx - 1];
            x11 = x[(nz - 1) * nxy + (ny - 1) * nx + nx - 1];
            x01 = x[(nz - 1) * nxy + (ny - 1) * nx];

            xr0 = x[(ny - 1) * nx + i];
            x1t = x[k * nxy + (ny - 1) * nx + nx - 1];
            xr1 = x[(nz - 1) * nxy + (ny - 1) * nx + i];
            x0t = x[k * nxy + (ny - 1) * nx];

            offset = k * nxy + (ny - 1) * nx + i;

            x[offset] = transfinite_blend(r, t, x00, x10, x11, x01, xr0, x1t, xr1, x0t);

        }
    }

    double x0s, x1s;
    for (int j = 1; j < ny - 1; j++)
    {
        s = glynodes[j];
        for (int i = 1; i < nx - 1; i++)
        {
            r = glxnodes[i];

            // Back Face ...
            x00 = x[0];
            x10 = x[nx - 1];
            x11 = x[(ny - 1) * nx + (nx - 1)];
            x01 = x[(ny - 1) * nx];

            xr0 = x[i];
            x1s = x[j * nx + nx - 1];
            xr1 = x[(ny - 1) * nx + i];
            x0s = x[j * nx];

            offset = j * nx + i;

            x[offset] = transfinite_blend(r, s, x00, x10, x11, x01, xr0, x1s, xr1, x0s);

            // Front Face ...
            x00 = x[(nz - 1) * nxy];
            x10 = x[(nz - 1) * nxy + nx - 1];
            x11 = x[(nz - 1) * nxy + (ny - 1) * nx + (nx - 1)];
            x01 = x[(nz - 1) * nxy + (ny - 1) * nx];

            xr0 = x[(nz - 1) * nxy + i];
            x1s = x[(nz - 1) * nxy + j * nx + nx - 1];
            xr1 = x[(nz - 1) * nxy + (ny - 1) * nx + i];
            x0s = x[(nz - 1) * nxy + j * nx];

            offset = (nz - 1) * nx * ny + j * nx + i;

            x[offset] = transfinite_blend(r, s, x00, x10, x11, x01, xr0, x1s, xr1, x0s);
        }
    }

    blend_from_faces(x, glxnodes, glynodes, glznodes);

    return;
}

////////////////////////////////////////////////////////////////////////////////

void TFIMap::blend_from_corners(vector<double> &x,
                                const vector<double> &glxnodes,
                                const vector<double> &glynodes,
                                const vector<double> &glznodes)
{

    //************************************************************************80
    //
    //  Purpose:
    //
    //    BLEND_IJK_0D1 extends indexed scalar data along corners into a cubic table.
    //
    //  Diagram:
    //
    //    ( X000,   ...,  ...,  ...,  ...,  ...,  X060 )
    //    ( ....,   ...,  ...,  ...,  ...,  ...,  ...  )
    //    ( ....,   ...,  ...,  ...,  ...,  ...,  ...  )   First "layer"
    //    ( ....,   ...,  ...,  ...,  ...,  ...,  ...  )
    //    ( X400,   ...,  ...,  ...,  ...,  ...,  X460 )
    //
    //    ( ....,   ...,  ...,  ...,  ...,  ...,  ...  )
    //    ( ....,   ...,  ...,  ...,  ...,  ...,  ...  )
    //    ( ....,   ...,  ...,  ...,  ...,  ...,  ...  )   Middle "layers"
    //    ( ....,   ...,  ...,  ...,  ...,  ...,  ...  )
    //    ( ....,   ...,  ...,  ...,  ...,  ...,  ...  )
    //
    //    ( X003,  ...,  ...,  ...,  ...,  ...,  X063  )
    //    ( ....,   ...,  ...,  ...,  ...,  ...,  ...  )
    //    ( ....,   ...,  ...,  ...,  ...,  ...,  ...  )   Last "layer"
    //    ( ....,   ...,  ...,  ...,  ...,  ...,  ...  )
    //    ( X403,  ...,  ...,  ...,  ...,  ...,  X463  )
    //
    //  Licensing:
    //
    //    This code is distributed under the GNU LGPL license.
    //
    //  Modified:
    //
    //    22 December 1998
    //
    //  Author:
    //
    //    John Burkardt
    //
    //  Reference:
    //
    //    William Gordon, Charles A Hall,
    //    Construction of Curvilinear Coordinate Systems and Application to
    //    Mesh Generation,
    //    International Journal of Numerical Methods in Engineering,
    //    Volume 7, pages 461-477, 1973.
    //
    //    Joe Thompson, Bharat Soni, Nigel Weatherill,
    //    Handbook of Grid Generation,
    //    CRC Press,
    //    1999.
    //
    //  Parameters:
    //
    //    Input/output, double X[NX*NY*NZ], a singly dimensioned array that
    //    is "really" triply dimensioned.  The triple dimension index
    //    [I][J][K] corresponds to the single dimension index
    //    K * NX*NY + J * NX + I
    //
    //    On input, there is already scalar data in the entries X[I][J][K]
    //    corresponding to "cornders" of the table, that is, entries for which
    //    each of the three indices I, J and K is equal to their
    //    minimum or maximum possible values.
    //
    //    Input, int NX, NY, NZ, the number of rows, columns, and layers in X.
    //
    ///////////////////////////////////////////////////////////////////////////



    int offset;
    double r, s, t, x0, x1;

    int nx = glxnodes.size();
    int ny = glynodes.size();
    int nz = glznodes.size();

    int nxy = nx*ny;

    for (int i = 1; i < nx - 1; i++)
    {
        r = glxnodes[i];

        x0 = x[0];
        x1 = x[nx - 1];
        offset = i;
        x[offset] = linear_interpolation(r, x0, x1);

        x0 = x[(ny - 1) * nx];
        x1 = x[(ny - 1) * nx + nx - 1];
        offset = (ny - 1) * nx + + i;
        x[offset] = linear_interpolation(r, x0, x1);

        x0 = x[(nz - 1) * nxy];
        x1 = x[(nz - 1) * nxy + nx - 1];
        offset = (nz - 1) * nxy + i;
        x[offset] = linear_interpolation(r, x0, x1);

        x0 = x[(nz - 1) * nxy + (ny - 1) * nx ];
        x1 = x[(nz - 1) * nxy + (ny - 1) * nx + nx - 1];
        offset = (nz - 1) * nxy + (ny - 1) * nx + i;
        x[offset] = linear_interpolation(r, x0, x1);

    }

    for (int j = 1; j < ny - 1; j++)
    {
        s = glynodes[j];

        x0 = x[0];
        x1 = x[ (ny - 1) * nx ];
        offset = j*nx;
        x[offset] = linear_interpolation(s, x0, x1);

        x0 = x[nx - 1];
        x1 = x[(ny - 1) * nx + (nx - 1)];
        offset = j * nx + nx - 1;
        x[offset] = linear_interpolation(s, x0, x1);

        x0 = x[(nz - 1) * nxy];
        x1 = x[(nz - 1) * nxy + (ny - 1) * nx];
        offset = (nz - 1) * nxy + j*nx;
        x[offset] = linear_interpolation(s, x0, x1);

        x0 = x[ (nz - 1) * nxy + nx - 1];
        x1 = x[ (nz - 1) * nxy + (ny - 1) * nx + nx - 1];
        offset = (nz - 1) * nxy + j * nx + nx - 1;
        x[offset] = linear_interpolation(s, x0, x1);
    }

    for (int k = 1; k < nz - 1; k++)
    {
        t = glznodes[k];

        x0 = x[0];
        x1 = x[(nz - 1) * nxy];
        offset = k*nxy;
        x[offset] = linear_interpolation(t, x0, x1);

        x0 = x[nx - 1];
        x1 = x[(nz - 1) * nxy + nx - 1];
        offset = k * nxy + nx - 1;
        x[offset] = linear_interpolation(t, x0, x1);

        x0 = x[(ny - 1) * nx];
        x1 = x[(nz - 1) * nxy + (ny - 1) * nx ];
        offset = k * nxy + (ny - 1) * nx;
        x[offset] = linear_interpolation(t, x0, x1);

        x0 = x[(ny - 1) * nx + nx - 1];
        x1 = x[(nz - 1) * nxy + (ny - 1) * nx + nx - 1];
        offset = k * nx * ny + (ny - 1) * nx + nx - 1;
        x[offset] = linear_interpolation(t, x0, x1);
    }

    blend_from_edges(x, glxnodes, glynodes, glznodes);

    return;
}

///////////////////////////////////////////////////////////////////////////////

inline Point3D TFIMap::CoEdge::getXYZCoords(double ur)
{
    assert(ur >= 0.0 && ur <= 1.0);

    int err;

    double u, umin, umax;
    iGeom_getEntURange(geometry, edgeHandle, &umin, &umax, &err);
    assert(!err);

    switch (direction)
    {
    case +1 :
        u = umin + ur * (umax - umin);
        break;
    case -1 :
        u = umax + ur * (umin - umax);
        break;
    }

    Point3D p3d;
    iGeom_getEntUtoXYZ(geometry, edgeHandle, u, &p3d[0], &p3d[1], &p3d[2], &err);
    assert(!err);

    return p3d;
}

///////////////////////////////////////////////////////////////////////////////

inline Point2D TFIMap::CoEdge::getUVCoords(double ur)
{
    assert(ur >= 0.0 && ur <= 1.0);

    int err;

    double u, umin, umax;
    iGeom_getEntURange(geometry, edgeHandle, &umin, &umax, &err);
    assert(!err);

    switch (direction)
    {
    case +1 :
        u = umin + ur * (umax - umin);
        break;
    case -1 :
        u = umax + ur * (umin - umax);
        break;
    }

    Point3D p3d;
    iGeom_getEntUtoXYZ(geometry, edgeHandle, u, &p3d[0], &p3d[1], &p3d[2], &err);
    assert(!err);

    Point2D p2d = geomface->getUVCoords(p3d);

    return p2d;
}


///////////////////////////////////////////////////////////////////////////////

int TFIMap::init_coedges2()
{
    return 1;
}

///////////////////////////////////////////////////////////////////////////////

int TFIMap::init_coedges4()
{
    int err;
    SimpleArray<iBase_EntityHandle> faceedges;
    iGeom_getEntAdj(geometry, gFaceHandle, iBase_EDGE, ARRAY_INOUT(faceedges), &err);
    assert(!err);
    assert(faceedges.size() == 4);

    // Need iGeom Extension functions, as we need ordered facenodes
    vector<iBase_EntityHandle> facenodes;
    iGeom_getEntAdj_Ext(geometry, gFaceHandle, iBase_VERTEX, facenodes, &err);
    assert(!err);
    assert(facenodes.size() == 4);

    int side_no, sense;
    SimpleArray<iBase_EntityHandle> edgenodes;
    for (int i = 0; i < 4; i++)
    {
        iGeom_getEntAdj(geometry, faceedges[i], iBase_VERTEX, ARRAY_INOUT(edgenodes), &err);
        assert(!err);
        GEdge::get_quad_edge_number(&facenodes[0], &edgenodes[0], side_no, sense);
        coedge[side_no].edgeHandle = faceedges[i];
        coedge[side_no].direction = sense;
        coedge[side_no].geometry = geometry;
    }

    double dl, eps = 1.0E-06;
    Point3D p1, p2;

    p1 = coedge[0].getXYZCoords(0.0);
    p2 = coedge[3].getXYZCoords(0.0);
    dl = square_length(p1, p2);
    if (dl > eps * eps)
        cout << "Warning: Mismatch at face vertex 0 " << sqrt(dl) << endl;

    p1 = coedge[0].getXYZCoords(1.0);
    p2 = coedge[1].getXYZCoords(0.0);
    dl = square_length(p1, p2);
    if (dl > eps * eps)
        cout << "Warning: Mismatch at face vertex 1 " << sqrt(dl) << endl;

    p1 = coedge[1].getXYZCoords(1.0);
    p2 = coedge[2].getXYZCoords(1.0);
    dl = square_length(p1, p2);
    if (dl > eps * eps)
        cout << "Warning: Mismatch at face vertex 2 " << sqrt(dl) << endl;

    p1 = coedge[2].getXYZCoords(0.0);
    p2 = coedge[3].getXYZCoords(1.0);
    dl = square_length(p1, p2);
    if (dl > eps * eps)
        cout << "Warning: Mismatch at face vertex 3 " << sqrt(dl) << endl;

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

int TFIMap::init_coedges()
{
    int err;
    SimpleArray<iBase_EntityHandle> faceedges;
    iGeom_getEntAdj(geometry, gFaceHandle, iBase_EDGE, ARRAY_INOUT(faceedges), &err);

    // if( faceedges.size() == 2) return init_coedges2();
    if (faceedges.size() == 4) return init_coedges4();

    return 1;
}

///////////////////////////////////////////////////////////////////////////////

int TFIMap::reconstruct_structured_mesh2d()
{
    int err;
    iBase_EntitySetHandle meshSet;

    std::set<iBase_EntityHandle> boundNodes, domainNodes;

    SimpleArray<iBase_EntityHandle> gEdges;
    iGeom_getEntAdj(geometry, gFaceHandle, iBase_EDGE, ARRAY_INOUT(gEdges), &err);

    iBase_TagHandle dim_tag;
    SimpleArray<iBase_EntityHandle> mEdges, edgenodes, facenodes, faceedges;

    const char *tag2 = "GEOM_DIMENSION";
    int namelen = strlen(tag2);
    iMesh_getTagHandle(mesh, tag2, &dim_tag, &err, namelen);
    assert(!err);

    int geom_dim;
    for (int i = 0; i < gEdges.size(); i++)
    {
        iRel_getEntSetAssociation(assoc, rel, gEdges[i], 0, &meshSet, &err);
        assert(!err);

        iMesh_getEntSetIntData(mesh, meshSet, dim_tag, &geom_dim, &err);
        assert(!err);

        if (geom_dim == 1)
        {
            mEdges.clear();

            iMesh_getEntities(mesh, meshSet, iBase_EDGE, iMesh_ALL_TOPOLOGIES, ARRAY_INOUT(mEdges), &err);
            assert(!err);

            for (int j = 0; j < mEdges.size(); j++)
            {

                edgenodes.clear();
                iMesh_getEntAdj(mesh, mEdges[j], iBase_VERTEX, ARRAY_INOUT(edgenodes), &err);
                for (int k = 0; k < edgenodes.size(); k++)
                    boundNodes.insert(edgenodes[k]);
            }
        }
    }

    SimpleArray<iBase_EntityHandle> mFaces;
    iRel_getEntSetAssociation(assoc, rel, gFaceHandle, 0, &meshSet, &err);
    iMesh_getEntities(mesh, meshSet, iBase_FACE, iMesh_ALL_TOPOLOGIES, ARRAY_INOUT(mFaces), &err);

    map<iBase_EntityHandle, vector<iBase_EntityHandle> > relation02Map; // Vertex->Faces Relations
    map<iBase_EntityHandle, set<iBase_EntityHandle> > relation00Map; // Vertex->Vertex Relations

    for (int i = 0; i < mFaces.size(); i++)
    {
        facenodes.clear();
        iMesh_getEntAdj(mesh, mFaces[i], iBase_VERTEX, ARRAY_INOUT(facenodes), &err);
        if (facenodes.size() != 4) return 1;
        for (int j = 0; j < facenodes.size(); j++)
        {
            relation02Map[facenodes[j]].push_back(mFaces[i]);
            domainNodes.insert(facenodes[j]);
        }

        iMesh_getEntAdj(mesh, mFaces[i], iBase_EDGE, ARRAY_INOUT(faceedges), &err);
        for (int j = 0; j < faceedges.size(); j++)
        {
            edgenodes.clear();
            iMesh_getEntAdj(mesh, faceedges[j], iBase_VERTEX, ARRAY_INOUT(edgenodes), &err);
            relation00Map[edgenodes[0]].insert(edgenodes[1]);
            relation00Map[edgenodes[1]].insert(edgenodes[0]);
        }
    }

    std::vector<iBase_EntityHandle> internalNodes;
    std::set_difference(domainNodes.begin(), domainNodes.end(),
                        boundNodes.begin(), boundNodes.end(),
                        back_inserter(internalNodes));

    // All the internal nodes must have four neighboring faces.
    for (int i = 0; i < internalNodes.size(); i++)
        if (relation02Map[internalNodes[i]].size() != 4) return 1;

    int nxedges = 0, nyedges = 0;
    iRel_getEntSetAssociation(assoc, rel, coedge[0].edgeHandle, 0, &meshSet, &err);
    assert(!err);
    iMesh_getEntSetIntData(mesh, meshSet, dim_tag, &geom_dim, &err);
    assert(!err);
    if (geom_dim == 1)
    {
        mEdges.clear();
        iMesh_getEntities(mesh, meshSet, iBase_EDGE, iMesh_ALL_TOPOLOGIES, ARRAY_INOUT(mEdges), &err);
        assert(!err);
        nxedges = mEdges.size();
    }

    ny = 0;
    iRel_getEntSetAssociation(assoc, rel, coedge[3].edgeHandle, 0, &meshSet, &err);
    assert(!err);
    iMesh_getEntSetIntData(mesh, meshSet, dim_tag, &geom_dim, &err);
    assert(!err);
    if (geom_dim == 1)
    {
        mEdges.clear();
        iMesh_getEntities(mesh, meshSet, iBase_EDGE, iMesh_ALL_TOPOLOGIES, ARRAY_INOUT(mEdges), &err);
        assert(!err);
        nyedges = mEdges.size();
    }

    iRel_getEntSetAssociation(assoc, rel, coedge[1].edgeHandle, 0, &meshSet, &err);
    assert(!err);
    iMesh_getEntSetIntData(mesh, meshSet, dim_tag, &geom_dim, &err);
    assert(!err);
    if (geom_dim == 1)
    {
        mEdges.clear();
        iMesh_getEntities(mesh, meshSet, iBase_EDGE, iMesh_ALL_TOPOLOGIES, ARRAY_INOUT(mEdges), &err);
        assert(!err);
        if (mEdges.size() != nyedges)
        {
            cout << "Warning: #Edges on left and right side don't match : TFI skipped " << endl;
            return 1;
        }
    }

    iRel_getEntSetAssociation(assoc, rel, coedge[2].edgeHandle, 0, &meshSet, &err);
    assert(!err);
    iMesh_getEntSetIntData(mesh, meshSet, dim_tag, &geom_dim, &err);
    assert(!err);
    if (geom_dim == 1)
    {
        mEdges.clear();
        iMesh_getEntities(mesh, meshSet, iBase_EDGE, iMesh_ALL_TOPOLOGIES, ARRAY_INOUT(mEdges), &err);
        assert(!err);
        if (mEdges.size() != nxedges)
        {
            cout << "Warning: #Edges on top and bottom side don't match : TFI skipped " << endl;
            return 1;
        }
    }

    nx = nxedges + 1;
    ny = nyedges + 1;

    structured_nodes.resize(nx * ny);
    for (int i = 0; i < nx * ny; i++) structured_nodes[i] = 0;

    //
    // Start with the lowest edge and be sure that it is already structured. Perhaps
    // this cann't be tested (in general ) combinatorically.
    //
    vector<iBase_EntityHandle> gNodes;
    iGeom_getEntAdj_Ext(geometry, gFaceHandle, iBase_VERTEX, gNodes, &err);

    iRel_getEntSetAssociation(assoc, rel, coedge[0].edgeHandle, 0, &meshSet, &err);
    mEdges.clear();
    iMesh_getEntities(mesh, meshSet, iBase_EDGE, iMesh_ALL_TOPOLOGIES, ARRAY_INOUT(mEdges), &err);

    std::list<iBase_EntityHandle> nodechain;
    vector<int> edgeProcessed(mEdges.size());

    for (int i = 0; i < mEdges.size(); i++) edgeProcessed[i] = 0;

    iBase_EntityHandle start_node, next_node;

    iMesh_getEntAdj(mesh, mEdges[0], iBase_VERTEX, ARRAY_INOUT(edgenodes), &err);

    nodechain.push_back(edgenodes[0]);
    nodechain.push_back(edgenodes[1]);

    edgeProcessed[0] = 1;

    // Forward Search;
    start_node = edgenodes[1];
    next_node = edgenodes[1];
    for (int i = 0; i < mEdges.size(); i++)
    {
        if (!edgeProcessed[i])
        {
            iMesh_getEntAdj(mesh, mEdges[i], iBase_VERTEX, ARRAY_INOUT(edgenodes), &err);
            if (edgenodes[0] == next_node)
            {
                next_node = edgenodes[1];
                edgeProcessed[i] = 1;
                nodechain.push_back(next_node);
            }
            else if (edgenodes[1] == next_node)
            {
                next_node = edgenodes[0];
                edgeProcessed[i] = 1;
                nodechain.push_back(next_node);
            }
        }
    }

    iMesh_getEntAdj(mesh, mEdges[0], iBase_VERTEX, ARRAY_INOUT(edgenodes), &err);

    // Back Search;
    next_node = edgenodes[0];
    start_node = next_node;
    for (int i = 0; i < mEdges.size(); i++)
    {
        if (!edgeProcessed[i])
        {
            iMesh_getEntAdj(mesh, mEdges[i], iBase_VERTEX, ARRAY_INOUT(edgenodes), &err);
            if (edgenodes[0] == next_node)
            {
                next_node = edgenodes[1];
                edgeProcessed[i] = 1;
                nodechain.push_front(next_node);
            }
            else if (edgenodes[1] == next_node)
            {
                next_node = edgenodes[0];
                edgeProcessed[i] = 1;
                nodechain.push_front(next_node);
            }
        }
    }

    for (int i = 0; i < mEdges.size(); i++)
    {
        if (!edgeProcessed[i])
        {
            cout << "Warning: Edge is not simple. structured nodes cann't be done" << endl;
            return 1;
        }
    }

    Point3D gxyz, mxyz;
    iGeom_getVtxCoord(geometry, gNodes[0], &gxyz[0], &gxyz[1], &gxyz[2], &err);
    iMesh_getVtxCoord(mesh, nodechain.front(), &mxyz[0], &mxyz[1], &mxyz[2], &err);
    double d0 = square_length(gxyz, mxyz);

    iGeom_getVtxCoord(geometry, gNodes[1], &gxyz[0], &gxyz[1], &gxyz[2], &err);
    iMesh_getVtxCoord(mesh, nodechain.back(), &mxyz[0], &mxyz[1], &mxyz[2], &err);
    double d1 = square_length(gxyz, mxyz);

    if (d1 < d0) reverse(nodechain.begin(), nodechain.end());

    assert(nodechain.size() == nx);

    iBase_EntityHandle thisnode, neigh;

    int index = 0;
    BOOST_FOREACH(thisnode, nodechain)
    structured_nodes[index++] = thisnode;

    for (int j = 0; j < ny - 1; j++)
    {
        for (int i = 0; i < nx; i++)
        {
            int offset = j * nx + i;
            int r = offset + 1;
            int l = offset - 1;
            int t = offset + nx;
            int b = offset - nx;
            thisnode = structured_nodes[offset];
            assert(thisnode != 0);
            if (i > 0)
            {
                neigh = structured_nodes[l]; // erase left node
                relation00Map[thisnode].erase(neigh);
            }
            if (i < nx - 1)
            {
                neigh = structured_nodes[r]; // erase right node
                relation00Map[thisnode].erase(neigh);
            }
            if (j > 0)
            {
                neigh = structured_nodes[b]; // erase bottom node
                relation00Map[thisnode].erase(neigh);
            }

            // Remaining node must be one and that should go to top.
            assert(relation00Map[thisnode].size() == 1);
            neigh = *relation00Map[thisnode].begin();
            structured_nodes[t] = neigh;
        }
    }
    return 0;
}

///////////////////////////////////////////////////////////////////////////////

int TFIMap::search_edge(set<iBase_EntityHandle> &eset,
                        iBase_EntityHandle n0,
                        iBase_EntityHandle n1, iBase_EntityHandle &edgeHandle)
{
    int err;

    SimpleArray<iBase_EntityHandle> edgenodes;

    BOOST_FOREACH(edgeHandle, eset)
    {
        iMesh_getEntAdj(mesh, edgeHandle, iBase_VERTEX, ARRAY_INOUT(edgenodes), &err);
        if (edgenodes[0] == n0 && edgenodes[1] == n1)
        {
            eset.erase(edgeHandle);
            return 1;
        }

        if (edgenodes[1] == n0 && edgenodes[0] == n1)
        {
            eset.erase(edgeHandle);
            return -1;
        }
    }

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

void TFIMap::project_edge_u(iBase_EntityHandle mEdgeHandle, int dir, double u0, double uN, double s)
{
    int err;
    char *tag_val = NULL;
    int tag_val_allocated, tag_val_size;
    iMesh_getData(mesh, mEdgeHandle, horder_tag, &tag_val, &tag_val_allocated, &tag_val_size, &err);
    assert(!err);

    HO_Points *hopoints = (HO_Points *) tag_val;
    iBase_EntityHandle *nodesOnEdge = hopoints->nodeHandles;
    int numHPoints = hopoints->nx;

    Point2D uv0s = coedge[3].getUVCoords(s);
    Point2D uv1s = coedge[1].getUVCoords(s);

    double u00, u10, u11, u01, ur0, u1s, ur1, u0s;
    double v00, v10, v11, v01, vr0, v1s, vr1, v0s;

    u00 = uCorner[0];
    u10 = uCorner[1];
    u11 = uCorner[2];
    u01 = uCorner[3];
    u0s = uv0s[0];
    u1s = uv1s[0];

    v00 = vCorner[0];
    v10 = vCorner[1];
    v11 = vCorner[2];
    v01 = vCorner[3];
    v0s = uv0s[1];
    v1s = uv1s[1];

    double u, r;
    Point3D pon;
    Point2D uv;

    iBase_EntityHandle currvertex;

    for (int i = 1; i < numHPoints - 1; i++)
    {
        u = gllnodes[i];
        r = 0.5 * (1 - u) * u0 + 0.5 * (1 + u) * uN;

        Point2D uvr0 = coedge[0].getUVCoords(r);
        Point2D uvr1 = coedge[2].getUVCoords(r);

        ur0 = uvr0[0];
        ur1 = uvr1[0];
        uv[0] = transfinite_blend(r, s, u00, u10, u11, u01, ur0, u1s, ur1, u0s);

        vr0 = uvr0[1];
        vr1 = uvr1[1];
        uv[1] = transfinite_blend(r, s, v00, v10, v11, v01, vr0, v1s, vr1, v0s);

        if (dir == 1)
            currvertex = nodesOnEdge[i];
        else
            currvertex = nodesOnEdge[numHPoints - 1 - i];

        pon = geomface->getXYZCoords(uv);
        iMesh_setVtxCoord(mesh, currvertex, pon[0], pon[1], pon[2], &err);

    }

}

///////////////////////////////////////////////////////////////////////////////////////////

void TFIMap::project_edge_v(iBase_EntityHandle mEdgeHandle, int dir, double v0, double vN, double r)
{

    int err;
    char *tag_val = NULL;
    int tag_val_allocated, tag_val_size;
    iMesh_getData(mesh, mEdgeHandle, horder_tag, &tag_val, &tag_val_allocated, &tag_val_size, &err);
    assert(!err);

    HO_Points *hopoints = (HO_Points *) tag_val;
    iBase_EntityHandle *nodesOnEdge = hopoints->nodeHandles;
    int numHPoints = hopoints->nx;

    Point2D uvr0 = coedge[0].getUVCoords(r);
    Point2D uvr1 = coedge[2].getUVCoords(r);

    double v, s;
    double u00, u10, u11, u01, ur0, u1s, ur1, u0s;
    double v00, v10, v11, v01, vr0, v1s, vr1, v0s;

    u00 = uCorner[0];
    u10 = uCorner[1];
    u11 = uCorner[2];
    u01 = uCorner[3];
    ur0 = uvr0[0];
    ur1 = uvr1[0];

    v00 = vCorner[0];
    v10 = vCorner[1];
    v11 = vCorner[2];
    v01 = vCorner[3];
    vr0 = uvr0[1];
    vr1 = uvr1[1];

    Point3D pon;
    Point2D uv;
    iBase_EntityHandle currvertex;

    for (int j = 1; j < numHPoints - 1; j++)
    {
        v = gllnodes[j];
        s = 0.5 * (1 - v) * v0 + 0.5 * (1 + v) * vN;
        Point2D uv0s = coedge[3].getUVCoords(s);
        Point2D uv1s = coedge[1].getUVCoords(s);

        u0s = uv0s[0];
        u1s = uv1s[0];
        uv[0] = transfinite_blend(r, s, u00, u10, u11, u01, ur0, u1s, ur1, u0s);

        v0s = uv0s[1];
        v1s = uv1s[1];
        uv[1] = transfinite_blend(r, s, v00, v10, v11, v01, vr0, v1s, vr1, v0s);


        if (dir == 1)
            currvertex = nodesOnEdge[j];
        else
            currvertex = nodesOnEdge[numHPoints - 1 - j];

        pon = geomface->getXYZCoords(uv);
        iMesh_setVtxCoord(mesh, currvertex, pon[0], pon[1], pon[2], &err);

    }
}

////////////////////////////////////////////////////////////////////////////////

void TFIMap::parametricTFI2D()
{
    int err;
    geomface = new GFace(geometry, gFaceHandle);

    cout << " NX NY " << nx << " " << ny << endl;

    vector<double> u, v;
    u.resize(nx * ny);
    v.resize(nx * ny);

    double du = 1.0 / (double) (nx - 1);
    double dv = 1.0 / (double) (ny - 1);

    int offset;
    Point2D uvnear, uv;
    Point3D xyz;

    // In the case of periodic surface, (U,V) coordinates of the corners
    // or on the edge could be ambigous. To remove the ambiguity, we try to
    // select one node on the surface which will be umbiguous and with respect
    // to the unambious point, we find out the (U,V) along an edge.
    //
    // Careful don't change the order of four loops because of dependency of
    // uvnear in 4 loops.

    Point3D p0 = coedge[0].getXYZCoords(du);
    Point3D p1 = coedge[3].getXYZCoords(dv);
    xyz[0] = 0.5 * (p0[0] + p1[0]);
    xyz[1] = 0.5 * (p0[1] + p1[1]);
    xyz[2] = 0.5 * (p0[2] + p1[2]);
    uvnear = geomface->getUVCoords(xyz);

    // Lower Edge ...
    for (int i = 0; i < nx; i++)
    {
        xyz = coedge[0].getXYZCoords(i * du);
        offset = i;
        uv = geomface->getUVCoords(xyz, uvnear);
        u[offset] = uv[0];
        v[offset] = uv[1];
        uvnear = uv;
    }

    // Right Edge ...
    for (int j = 0; j < ny; j++)
    {
        xyz = coedge[1].getXYZCoords(j * dv);
        offset = j * nx + (nx - 1);
        uv = geomface->getUVCoords(xyz, uvnear);
        u[offset] = uv[0];
        v[offset] = uv[1];
        uvnear = uv;
    }

    // Left Edge ...
    uvnear[0] = u[0];
    uvnear[1] = v[0];
    for (int j = 0; j < ny; j++)
    {
        xyz = coedge[3].getXYZCoords(j * dv);
        offset = j*nx;
        uv = geomface->getUVCoords(xyz, uvnear);
        u[offset] = uv[0];
        v[offset] = uv[1];
        uvnear = uv;
    }

    // Top Edge ...
    for (int i = 0; i < nx; i++)
    {
        xyz = coedge[2].getXYZCoords(i * du);
        offset = (ny - 1) * nx + i;
        uv = geomface->getUVCoords(xyz, uvnear);
        u[offset] = uv[0];
        v[offset] = uv[1];
        uvnear = uv;
    }

/*
    vector<double> glxnodes, glynodes;

    SpectralElements::gauss_linear_nodes(nx, glxnodes);
    SpectralElements::gauss_linear_nodes(ny, glynodes);

    blend_from_edges(u, glxnodes, glynodes);
    blend_from_edges(v, glxnodes, glynodes);

    for (int j = 1; j < ny - 1; j++)
    {
        for (int i = 1; i < nx - 1; i++)
        {
            int id = j * nx + i;
            uv[0] = u[id];
            uv[1] = v[id];
            xyz = geomface->getXYZCoords(uv);
            iMesh_setVtxCoord(mesh, structured_nodes[id], xyz[0], xyz[1], xyz[2], &err);
        }
    }

    hoelement_order = gllnodes.size();

    if (hoelement_order < 3) return;

    // Now do the projection of higher order nodes on each using TFI.

    SimpleArray<iBase_EntityHandle> mEdges;
    iBase_EntitySetHandle meshSet;

    iRel_getEntSetAssociation(assoc, rel, gFaceHandle, 0, &meshSet, &err);
    iMesh_getEntities(mesh, meshSet, iBase_EDGE, iMesh_ALL_TOPOLOGIES, ARRAY_INOUT(mEdges), &err);

    set<iBase_EntityHandle> edgeSet;
    for (int i = 0; i < mEdges.size(); i++)
        edgeSet.insert(mEdges[i]);
    mEdges.clear();

    coedge[0].geomface = geomface;
    coedge[1].geomface = geomface;
    coedge[2].geomface = geomface;
    coedge[3].geomface = geomface;

    //
    // Note that Gauss-Lobatto points are symmetric from the ends, therefore, no check is
    // performed on the edge direction in the following code.
    //
    // In the search_edge function, everytime an edge between (n0,n1) is found, it is removed
    // from the set, so that the successive searchs will be faster.

    uCorner.resize(4);
    vCorner.resize(4);

    uCorner[0] = u[0];
    vCorner[0] = v[0];

    uCorner[1] = u[nx - 1];
    vCorner[1] = v[nx - 1];

    uCorner[2] = u[nx * ny - 1];
    vCorner[2] = v[nx * ny - 1];

    uCorner[3] = u[nx * ny - nx];
    vCorner[3] = v[nx * ny - nx];

    iBase_EntityHandle n0, n1, mEdgeHandle;

    double u0, uN, v0, vN;
    for (int j = 1; j < ny - 1; j++)
    {
        v0 = j*dv;
        for (int i = 0; i < nx - 1; i++)
        {
            int id = j * nx + i;
            n0 = structured_nodes[id];
            n1 = structured_nodes[id + 1];
            int dir = search_edge(edgeSet, n0, n1, mEdgeHandle);
            u0 = i*du;
            uN = u0 + du;
            project_edge_u(mEdgeHandle, dir, u0, uN, v0);
        }
    }

    for (int i = 1; i < nx - 1; i++)
    {
        u0 = i*du;
        for (int j = 0; j < ny - 1; j++)
        {
            int id = j * nx + i;
            n0 = structured_nodes[id];
            n1 = structured_nodes[id + nx];
            int dir = search_edge(edgeSet, n0, n1, mEdgeHandle);
            v0 = j*dv;
            vN = v0 + dv;
            project_edge_v(mEdgeHandle, dir, v0, vN, u0);
        }
    }

    SimpleArray<iBase_EntityHandle> mFaces;
    iMesh_getEntities(mesh, meshSet, iBase_FACE, iMesh_ALL_TOPOLOGIES, ARRAY_INOUT(mFaces), &err);

    for (int i = 0; i < mFaces.size(); i++)
        projectHigherOrderNodes2D(mFaces[i], gllnodes);
*/
    cout << " CSV " << endl;

}

////////////////////////////////////////////////////////////////////////////////

void TFIMap::project_face(iBase_EntityHandle mFaceHandle)
{
    int err;

    int offset, nhx, nhy, numHPoints;

    SimpleArray<iBase_EntityHandle> mEdges, faceNodes;
    iMesh_getEntAdj(mesh, mFaceHandle, iBase_EDGE, ARRAY_INOUT(mEdges), &err);
    iMesh_getEntAdj(mesh, mFaceHandle, iBase_VERTEX, ARRAY_INOUT(faceNodes), &err);

    Point3D p3d, pCentroid;
    Point2D uv, uvCentroid;

    pCentroid[0] = 0.0;
    pCentroid[1] = 0.0;
    pCentroid[2] = 0.0;
    for (int i = 0; i < 4; i++)
    {
        iMesh_getVtxCoord(mesh, faceNodes[i], &p3d[0], &p3d[1], &p3d[2], &err);
        pCentroid[0] += p3d[0];
        pCentroid[1] += p3d[1];
        pCentroid[2] += p3d[2];
    }

    pCentroid[0] /= 4.0;
    pCentroid[1] /= 4.0;
    pCentroid[2] /= 4.0;
    uvCentroid = geomface->getUVCoords(pCentroid);

    char *tag_val = NULL;
    int tag_val_allocated, tag_val_size;
    iMesh_getData(mesh, mFaceHandle, horder_tag, &tag_val, &tag_val_allocated, &tag_val_size, &err);
    assert(!err);
    HO_Points *hopoints = (HO_Points *) tag_val;
    nhx = hopoints->nx;
    nhy = hopoints->ny;
    numHPoints = nhx*nhy;

    iBase_EntityHandle *nodeHandles = hopoints->nodeHandles;

    vector<double> u(numHPoints);
    vector<double> v(numHPoints);

    iBase_EntityHandle currvertex;

    offset = 0;
    currvertex = nodeHandles[offset];
    iMesh_getVtxCoord(mesh, currvertex, &p3d[0], &p3d[1], &p3d[2], &err);
    uv = geomface->getUVCoords(p3d, uvCentroid);
    u[offset] = uv[0];
    v[offset] = uv[1];

    offset = nhx - 1;
    currvertex = nodeHandles[offset];
    iMesh_getVtxCoord(mesh, currvertex, &p3d[0], &p3d[1], &p3d[2], &err);
    uv = geomface->getUVCoords(p3d, uvCentroid);
    u[offset] = uv[0];
    v[offset] = uv[1];

    offset = (nhy - 1) * nhx;
    currvertex = nodeHandles[offset];
    iMesh_getVtxCoord(mesh, currvertex, &p3d[0], &p3d[1], &p3d[2], &err);
    uv = geomface->getUVCoords(p3d, uvCentroid);
    u[offset] = uv[0];
    v[offset] = uv[1];

    offset = nhx * nhy - 1;
    currvertex = nodeHandles[offset];
    iMesh_getVtxCoord(mesh, currvertex, &p3d[0], &p3d[1], &p3d[2], &err);
    uv = geomface->getUVCoords(p3d, uvCentroid);
    u[offset] = uv[0];
    v[offset] = uv[1];

    for (int i = 1; i < nhx - 1; i++)
    {
        offset = i;
        currvertex = nodeHandles[offset];
        iMesh_getVtxCoord(mesh, currvertex, &p3d[0], &p3d[1], &p3d[2], &err);
        uv = geomface->getUVCoords(p3d, uvCentroid);
        u[offset] = uv[0];
        v[offset] = uv[1];

        offset = i + (nhy - 1) * nhx;
        currvertex = nodeHandles[offset];
        iMesh_getVtxCoord(mesh, currvertex, &p3d[0], &p3d[1], &p3d[2], &err);
        uv = geomface->getUVCoords(p3d, uvCentroid);
        u[offset] = uv[0];
        v[offset] = uv[1];
    }

    for (int j = 1; j < nhy - 1; j++)
    {
        offset = j*nhx;
        currvertex = nodeHandles[offset];
        iMesh_getVtxCoord(mesh, currvertex, &p3d[0], &p3d[1], &p3d[2], &err);
        uv = geomface->getUVCoords(p3d, uvCentroid);
        u[offset] = uv[0];
        v[offset] = uv[1];

        offset = j * nhx + (nhx - 1);
        currvertex = nodeHandles[offset];
        iMesh_getVtxCoord(mesh, currvertex, &p3d[0], &p3d[1], &p3d[2], &err);
        uv = geomface->getUVCoords(p3d, uvCentroid);
        u[offset] = uv[0];
        v[offset] = uv[1];
    }

    TFIMap::blend_from_edges(u, gllnodes, gllnodes);
    TFIMap::blend_from_edges(v, gllnodes, gllnodes);

    for (int j = 1; j < nhy - 1; j++)
    {
        for (int i = 1; i < nhx - 1; i++)
        {
            offset = j * nhx + i;
            uv[0] = u[offset];
            uv[1] = v[offset];
            Point3D pon = geomface->getXYZCoords(uv);
            iMesh_setVtxCoord(mesh, nodeHandles[offset], pon[0], pon[1], pon[2], &err);
        }
    }
}


///////////////////////////////////////////////////////////////////////////////

void TFIMap::physicalTFI2D()
{
    int err;
    vector<double> x, y, z;

    x.resize(nx * ny);
    y.resize(nx * ny);
    z.resize(nx * ny);
    for (int i = 0; i < nx * ny; i++)
        iMesh_getVtxCoord(mesh, structured_nodes[i], &x[0], &y[0], &z[0], &err);

    vector<double> glxnodes, glynodes;
    SpectralElements::gauss_linear_nodes(nx, glxnodes);
    SpectralElements::gauss_linear_nodes(ny, glynodes);

    blend_from_edges(x, glxnodes, glynodes);
    blend_from_edges(y, glxnodes, glynodes);
    blend_from_edges(z, glxnodes, glynodes);

    for (int j = 1; j < ny - 1; j++)
    {
        for (int i = 1; i < nx - 1; i++)
        {
            int id = j * nx + i;
            iMesh_setVtxCoord(mesh, structured_nodes[id], x[id], y[id], z[id], &err);
        }
    }
}

///////////////////////////////////////////////////////////////////////////////

void TFIMap::projectHigherOrderNodes2D(iBase_EntityHandle mFaceHandle, const vector<double> &glnodes)
{
    int err;
    char *tag_val = NULL;
    int tag_val_allocated, tag_val_size;

    iMesh_getData(mesh, mFaceHandle, horder_tag, &tag_val, &tag_val_allocated, &tag_val_size, &err);
    assert(!err);

    HO_Points *hopoints = (HO_Points *) tag_val;
    int nhx = hopoints->nx;
    int nhy = hopoints->ny;
    int numHPoints = nhx*nhy;
    iBase_EntityHandle *nodeHandles = hopoints->nodeHandles;

    assert(glnodes.size() == nhx);
    assert(glnodes.size() == nhy);

    vector<double> x, y, z;
    x.resize(numHPoints);
    y.resize(numHPoints);
    z.resize(numHPoints);
    for (int i = 0; i < numHPoints; i++)
        iMesh_getVtxCoord(mesh, nodeHandles[i], &x[i], &y[i], &z[i], &err);

    TFIMap::blend_from_edges(x, glnodes, glnodes);
    TFIMap::blend_from_edges(y, glnodes, glnodes);
    TFIMap::blend_from_edges(z, glnodes, glnodes);

    for (int i = 0; i < numHPoints; i++)
        iMesh_setVtxCoord(mesh, nodeHandles[i], x[i], y[i], z[i], &err);
}

///////////////////////////////////////////////////////////////////////////////

void TFIMap::projectHigherOrderNodes3D(iBase_EntityHandle mCellHandle, const vector<double> &glnodes)
{
    int err;
    char *tag_val = NULL;
    int tag_val_allocated, tag_val_size;

    iMesh_getData(mesh, mCellHandle, horder_tag, &tag_val, &tag_val_allocated, &tag_val_size, &err);
    HO_Points *hopoints = (HO_Points *) tag_val;
    int nhx = hopoints->nx;
    int nhy = hopoints->ny;
    int nhz = hopoints->nz;
    int numHPoints = nhx * nhy*nhz;
    iBase_EntityHandle *nodeHandles = hopoints->nodeHandles;

    assert(glnodes.size() == nhx);
    assert(glnodes.size() == nhy);
    assert(glnodes.size() == nhz);

    vector<double> x, y, z;
    x.resize(numHPoints);
    y.resize(numHPoints);
    z.resize(numHPoints);
    for (int i = 0; i < numHPoints; i++)
        iMesh_getVtxCoord(mesh, nodeHandles[i], &x[i], &y[i], &z[i], &err);

    TFIMap::blend_from_faces(x, glnodes, glnodes, glnodes);
    TFIMap::blend_from_faces(y, glnodes, glnodes, glnodes);
    TFIMap::blend_from_faces(z, glnodes, glnodes, glnodes);

    for (int i = 0; i < numHPoints; i++)
        iMesh_setVtxCoord(mesh, nodeHandles[i], x[i], y[i], z[i], &err);
}

///////////////////////////////////////////////////////////////////////////////

int TFIMap::getTFI2D(iBase_EntityHandle gface, const vector<double> &gnodes)
{
    gFaceHandle = gface;

    if (init_coedges()) return 1;

    if (reconstruct_structured_mesh2d()) return 1;

    gllnodes = gnodes;

    switch (tfi_type)
    {
    case PARAMETRIC_TFI:
        parametricTFI2D();
        break;
    case PHYSICAL_TFI:
        physicalTFI2D();
        break;
    }

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

void TFIMap::saveAs(const string &filename)
{
    int err;

    int numHPoints = gllnodes.size();

    int numEdges = ny * (nx - 1) + nx * (ny - 1);
    int numFaces = (nx - 1)*(ny - 1);

    int numNodes = nx*ny + (numHPoints-2)*numEdges + (numHPoints-2)*(numHPoints-2)*numFaces;

    ofstream ofile(filename.c_str(), ios::out);
    if (ofile.fail())
    {
        cout << "Warning: Cann't open file " << filename << endl;
        return;
    }

    ofile << "#Nodes " << numNodes << endl;

    int node_index = 0;
    double x, y, z;

    for (int i = 0; i < nx * ny; i++)
    {
        iMesh_getVtxCoord(mesh, structured_nodes[i], &x, &y, &z, &err);
        assert(!err);
        ofile << node_index++ << " " << x << " " << y << " " << z << endl;
    }


    iBase_EntitySetHandle meshSet;
    SimpleArray<iBase_EntityHandle> mEdges;

    char *tag_val = NULL;
    int tag_val_allocated, tag_val_size;

    SimpleArray<iBase_EntityHandle> gEdges;
    iGeom_getEntAdj(geometry, gFaceHandle, iBase_EDGE, ARRAY_INOUT(gEdges), &err);

    for (int k = 0; k < gEdges.size(); k++)
    {
        mEdges.clear();
        iRel_getEntSetAssociation(assoc, rel, gEdges[k], 0, &meshSet, &err);
        iMesh_getEntities(mesh, meshSet, iBase_EDGE, iMesh_ALL_TOPOLOGIES, ARRAY_INOUT(mEdges), &err);
        for (int j = 0; j < mEdges.size(); j++)
        {
            iMesh_getData(mesh, mEdges[j], horder_tag, &tag_val, &tag_val_allocated, &tag_val_size, &err);
            assert(!err);

            HO_Points *hopoints = (HO_Points *) tag_val;
            iBase_EntityHandle *nodesOnEdge = hopoints->nodeHandles;
            for (int i = 1; i < numHPoints - 1; i++)
            {
                iMesh_getVtxCoord(mesh, nodesOnEdge[i], &x, &y, &z, &err);
                ofile << node_index++ << " " << x << " " << y << " " << z << endl;
            }
        }
    }

    iRel_getEntSetAssociation(assoc, rel, gFaceHandle, 0, &meshSet, &err);

    mEdges.clear();
    iMesh_getEntities(mesh, meshSet, iBase_EDGE, iMesh_ALL_TOPOLOGIES, ARRAY_INOUT(mEdges), &err);

    for (int i = 0; i < mEdges.size(); i++)
    {
        iMesh_getData(mesh, mEdges[i], horder_tag, &tag_val, &tag_val_allocated, &tag_val_size, &err);
        assert(!err);

        HO_Points *hopoints = (HO_Points *) tag_val;
        iBase_EntityHandle *nodesOnEdge = hopoints->nodeHandles;
        for (int j = 1; j < numHPoints - 1; j++)
        {
            iMesh_getVtxCoord(mesh, nodesOnEdge[j], &x, &y, &z, &err);
            ofile << node_index++ << " " << x << " " << y << " " << z << endl;
        }
    }

    SimpleArray<iBase_EntityHandle> mFaces;
    iMesh_getEntities(mesh, meshSet, iBase_FACE, iMesh_ALL_TOPOLOGIES, ARRAY_INOUT(mFaces), &err);

    for (int k = 0; k < mFaces.size(); k++)
    {
        iMesh_getData(mesh, mFaces[k], horder_tag, &tag_val, &tag_val_allocated, &tag_val_size, &err);
        assert(!err);

        HO_Points *hopoints = (HO_Points *) tag_val;
        iBase_EntityHandle *nodesOnEdge = hopoints->nodeHandles;
        for (int j = 1; j < numHPoints - 1; j++)
        {
            for (int i = 1; i < numHPoints - 1; i++)
            {
                iMesh_getVtxCoord(mesh, nodesOnEdge[j * numHPoints + i], &x, &y, &z, &err);
                ofile << node_index++ << " " << x << " " << y << " " << z << endl;
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////
