#include <limits>
#include <iomanip>
#include <string.h>
#include <sstream>
#include <fstream>

#include "GEdge.h"

int GEdge::get_quad_edge_number(const iBase_EntityHandle *facenodes,
                                const iBase_EntityHandle *edgenodes,
                                int &side_no, int &orientation)
{
     // Counter-Clockwise is +ve direction;

     // Edge 0-1 : ID = 0
     if (edgenodes[0] == facenodes[0] && edgenodes[1] == facenodes[1]) {
          side_no = 0;
          orientation = 1;
          return 0;
     }
     if (edgenodes[0] == facenodes[1] && edgenodes[1] == facenodes[0]) {
          side_no = 0;
          orientation = -1;
          return 0;
     }

     // Edge 1-2 : ID = 1
     if (edgenodes[0] == facenodes[1] && edgenodes[1] == facenodes[2]) {
          side_no = 1;
          orientation = 1;
          return 0;
     }
     if (edgenodes[0] == facenodes[2] && edgenodes[1] == facenodes[1]) {
          side_no = 1;
          orientation = -1;
          return 0;
     }

     // Edge 2-3 : ID = 2
     if (edgenodes[0] == facenodes[3] && edgenodes[1] == facenodes[2]) {
          side_no = 2;
          orientation = 1;
          return 0;
     }
     if (edgenodes[0] == facenodes[2] && edgenodes[1] == facenodes[3]) {
          side_no = 2;
          orientation = -1;
          return 0;
     }

     // Edge 0-3 : ID = 3
     if (edgenodes[0] == facenodes[0] && edgenodes[1] == facenodes[3]) {
          side_no = 3;
          orientation = 1;
          return 0;
     }
     if (edgenodes[0] == facenodes[3] && edgenodes[1] == facenodes[0]) {
          side_no = 3;
          orientation = -1;
          return 0;
     }

     cout << " Fatal Error: No Edge found in the face " << endl;
     return 0;
}

////////////////////////////////////////////////////////////////////////////////

int GEdge::get_hex_edge_number(const iBase_EntityHandle *cellnodes,
                               const iBase_EntityHandle *edgenodes,
                               int &side_no, int &orientation)
{

     // Edge 0-1 : ID = 0
     if (edgenodes[0] == cellnodes[0] && edgenodes[1] == cellnodes[1]) {
          side_no = 0;
          orientation = 1;
          return 0;
     }
     if (edgenodes[0] == cellnodes[1] && edgenodes[1] == cellnodes[0]) {
          side_no = 0;
          orientation = -1;
          return 0;
     }

     // Edge 1-2 : ID = 1
     if (edgenodes[0] == cellnodes[1] && edgenodes[1] == cellnodes[2]) {
          side_no = 1;
          orientation = 1;
          return 0;
     }
     if (edgenodes[0] == cellnodes[2] && edgenodes[1] == cellnodes[1]) {
          side_no = 1;
          orientation = -1;
          return 0;
     }

     // Edge 2-3 : ID = 2
     if (edgenodes[0] == cellnodes[3] && edgenodes[1] == cellnodes[2]) {
          side_no = 2;
          orientation = 1;
          return 0;
     }
     if (edgenodes[0] == cellnodes[2] && edgenodes[1] == cellnodes[3]) {
          side_no = 2;
          orientation = -1;
          return 0;
     }

     // Edge 0-3 : ID = 3
     if (edgenodes[0] == cellnodes[0] && edgenodes[1] == cellnodes[3]) {
          side_no = 3;
          orientation = 1;
          return 0;
     }
     if (edgenodes[0] == cellnodes[3] && edgenodes[1] == cellnodes[0]) {
          side_no = 3;
          orientation = -1;
          return 0;
     }

     // Edge 0-4 : ID = 4
     if (edgenodes[0] == cellnodes[0] && edgenodes[1] == cellnodes[4]) {
          side_no = 4;
          orientation = 1;
          return 0;
     }
     if (edgenodes[0] == cellnodes[4] && edgenodes[1] == cellnodes[0]) {
          side_no = 4;
          orientation = -1;
          return 0;
     }

     // Edge 1-5 : ID = 5
     if (edgenodes[0] == cellnodes[1] && edgenodes[1] == cellnodes[5]) {
          side_no = 5;
          orientation = 1;
          return 0;
     }
     if (edgenodes[0] == cellnodes[5] && edgenodes[1] == cellnodes[1]) {
          side_no = 5;
          orientation = -1;
          return 0;
     }

     // Edge 2-6 : ID = 6
     if (edgenodes[0] == cellnodes[2] && edgenodes[1] == cellnodes[6]) {
          side_no = 6;
          orientation = 1;
          return 0;
     }
     if (edgenodes[0] == cellnodes[6] && edgenodes[1] == cellnodes[2]) {
          side_no = 6;
          orientation = -1;
          return 0;
     }

     // Edge 3-7 : ID = 7
     if (edgenodes[0] == cellnodes[3] && edgenodes[1] == cellnodes[7]) {
          side_no = 7;
          orientation = 1;
          return 0;
     }
     if (edgenodes[0] == cellnodes[7] && edgenodes[1] == cellnodes[3]) {
          side_no = 7;
          orientation = -1;
          return 0;
     }

     // Edge 4-5 : ID = 8
     if (edgenodes[0] == cellnodes[4] && edgenodes[1] == cellnodes[5]) {
          side_no = 8;
          orientation = 1;
          return 0;
     }
     if (edgenodes[0] == cellnodes[5] && edgenodes[1] == cellnodes[4]) {
          side_no = 8;
          orientation = -1;
          return 0;
     }

     // Edge 5-6 : ID = 9
     if (edgenodes[0] == cellnodes[5] && edgenodes[1] == cellnodes[6]) {
          side_no = 9;
          orientation = 1;
          return 0;
     }
     if (edgenodes[0] == cellnodes[6] && edgenodes[1] == cellnodes[5]) {
          side_no = 9;
          orientation = -1;
          return 0;
     }

     // Edge 6-7 : ID = 10
     if (edgenodes[0] == cellnodes[7] && edgenodes[1] == cellnodes[6]) {
          side_no = 10;
          orientation = 1;
          return 0;
     }
     if (edgenodes[0] == cellnodes[6] && edgenodes[1] == cellnodes[7]) {
          side_no = 10;
          orientation = -1;
          return 0;
     }

     // Edge 4-7 : ID = 11
     if (edgenodes[0] == cellnodes[4] && edgenodes[1] == cellnodes[7]) {
          side_no = 11;
          orientation = 1;
          return 0;
     }
     if (edgenodes[0] == cellnodes[7] && edgenodes[1] == cellnodes[4]) {
          side_no = 11;
          orientation = -1;
          return 0;
     }

     cout << " Fatal Error: No Edge found in the cell " << endl;
     return 0;
}

//////////////////////////////////////////////////////////////////////////////

GEdge::GEdge(iBase_EntityHandle h, iGeom *g, iMesh *m, iRel *r, iRel::PairHandle *p)
{
     gEdgeHandle = h;
     geom    = g;
     mesh    = m;
     rel     = r;
     relPair = p;
     kdtree = NULL;

     assert( geom );

     int err;
     err = geom->getEntURange(gEdgeHandle, umin, umax);
     assert(!err);

     bool in_u, in_v;
     err = geom->isEntPeriodic(gEdgeHandle, in_u, in_v);
     assert(!err);

     periodic = in_u;

/*
     getXYZCoords(umin, p0);
     getXYZCoords(umax, p1);

     double dx = p0[0] - p1[0];
     double dy = p0[1] - p1[1];
     double dz = p0[2] - p1[2];
     double eps = 1.0E-08;
     if (dx * dx + dy * dy + dz * dz < eps * eps) periodic = 1;

     numNodes = 1000;
     uCoord = NULL;
     numNeighs = 2;
     annIdx.resize(numNeighs);
     anndist.resize(numNeighs);
     build_kdtree();
*/

}


//////////////////////////////////////////////////////////////////////////////

GEdge::~GEdge()
{
     delete_kdtree();
}
//////////////////////////////////////////////////////////////////////////////

/*
int GEdge::getURange( Point2D &p2d) const
{
     p2d[0] = umin;
     p2d[1] = umax;

     return 0;
}

//////////////////////////////////////////////////////////////////////////////

int GEdge::getXYZCoords(double u, Point3D &p3d) const
{
     int err;
     err = geom->getEntUtoXYZ(gEdgeHandle, u, p3d[0], p3d[1], p3d[2]);
     assert(!err);

     return err;
}


//////////////////////////////////////////////////////////////////////////////

double GEdge::getLength() const
{
     double measure;
     int err = geom->measure(&gEdgeHandle, 1, &measure ); assert( !err );
     return measure;
}
*/
//////////////////////////////////////////////////////////////////////////////

double GEdge::getLength(double ustart, double uend) const
{

     int err;
     double u, du;
     double x0, y0, z0;
     double x1, y1, z1;
     double dx, dy, dz, dl, arclen, maxlen;

     // Recursively subdivide the edge and calculate the length.

     int nsub = 2;
     maxlen = 0.0;
     for (int ilevel = 0; ilevel < 10; ilevel++) {
          du = (uend - ustart) / (double) nsub;

          err = geom->getEntUtoXYZ(gEdgeHandle, ustart, x0, y0, z0);
          assert(!err);

          arclen = 0.0;
          for (int i = 0; i < nsub; i++) {
               u = ustart + (i + 1) * du;
               err = geom->getEntUtoXYZ(gEdgeHandle, u, x1, y1, z1);
               assert(!err);
               dx = x1 - x0;
               dy = y1 - y0;
               dz = z1 - z0;
               dl = sqrt(dx * dx + dy * dy + dz * dz);
               arclen += dl;
               x0 = x1;
               y0 = y1;
               z0 = z1;
          }
          if (fabs(arclen - maxlen) < 1.0E-06) break;
          maxlen = arclen;
          nsub *= 2;
     }
     return arclen;

}

//////////////////////////////////////////////////////////////////////////////

double GEdge::getUCoord(const Point3D &qPoint) const
{
     static double p[3];

     assert(kdtree);

     p[0] = qPoint[0];
     p[1] = qPoint[1];
     p[2] = qPoint[2];

     double eps = 0.0;
     kdtree->annkSearch(p, numNeighs, &annIdx[0], &anndist[0], eps);

     int index = annIdx[0];

     return uCoord[index];
}

//////////////////////////////////////////////////////////////////////////////

int GEdge::getUCoord(double ustart, double dist, double &uguess) const
{
     assert( fabs(uguess-ustart) > 1.0E-15 );

     getXYZCoords(ustart, p0);
     getXYZCoords(uguess, p1);

     double dx, dy, dz, dl, u = uguess;

     double tol = 1.0E-12;

     int ntrials = 0;
     while (1) {
          dx = p1[0] - p0[0];
          dy = p1[1] - p0[1];
          dz = p1[2] - p0[2];
          dl = sqrt(dx * dx + dy * dy + dz * dz);
          if ( fabs(dl-dist) < tol) break;
          u = ustart + (u - ustart) * (dist/dl);
          getXYZCoords(u, p1);

          if (ntrials++ == 100) {
               cout << " Warning: Searching for U failed " << endl;
          }
     }

     uguess = u;

     return 0;
}

//////////////////////////////////////////////////////////////////////////////

void GEdge::refine_search(const Point3D &qPoint, double u0, double u1, double &urefine, double &mindist2) const
{
     double um = 0.5 * (u0 + u1);

     double dx, dy, dz, d0, d1, dm;

     getXYZCoords(u0, p0);
     dx = qPoint[0] - p0[0];
     dy = qPoint[1] - p0[1];
     dz = qPoint[2] - p0[2];
     d0 = dx * dx + dy * dy + dz*dz;
     mindist2 = std::min(d0, mindist2);

     getXYZCoords(u1, p1);
     dx = qPoint[0] - p1[0];
     dy = qPoint[1] - p1[1];
     dz = qPoint[2] - p1[2];
     d1 = dx * dx + dy * dy + dz*dz;

     if (d1 < d0)
          cout << " Warning: Nearest neighbour is ambiguous " << endl;

     int N = 1000;
     double du = (u1 - u0) / (double) (N - 1);

     Point3D pm;
     for (int i = 0; i < N - 2; i++) {
          um = u0 + (i + 1) * du;
          getXYZCoords(um, pm);
          dx = qPoint[0] - pm[0];
          dy = qPoint[1] - pm[1];
          dz = qPoint[2] - pm[2];
          dm = dx * dx + dy * dy + dz*dz;
          if (dm > mindist2) break;
          mindist2 = std::min(dm, mindist2);
          urefine = um;
     }
}

//////////////////////////////////////////////////////////////////////////////

double GEdge::getUCoord(const Point3D &qPoint, double unear) const
{
     static double p[3];

     assert(kdtree);

     p[0] = qPoint[0];
     p[1] = qPoint[1];
     p[2] = qPoint[2];

     double eps = 0.0;
     kdtree->annkSearch(p, numNeighs, &annIdx[0], &anndist[0], eps);

     int index = annIdx[0];
     double u = uCoord[index];

     Point3D pnear;
     getXYZCoords(u, pnear);
     double dx = qPoint[0] - pnear[0];
     double dy = qPoint[1] - pnear[1];
     double dz = qPoint[2] - pnear[2];
     double d1 = dx * dx + dy * dy + dz*dz;

     double mindist2 = std::numeric_limits<double>::max();

     double u0, u1;
     if (index < numNodes - 1) {
          u0 = uCoord[index];
          u1 = uCoord[index + 1];
          refine_search(qPoint, u0, u1, u, mindist2);
     }

     if (index > 0) {
          u0 = uCoord[index];
          u1 = uCoord[index - 1];
          refine_search(qPoint, u0, u1, u, mindist2);
     }

     if (periodic) {
          if (index == 0 || index == numNodes - 1)
               u = fabs(unear - umin) < fabs(unear - umax) ? umin : umax;
     }

     return u;
}

//////////////////////////////////////////////////////////////////////////////

int GEdge::getClosestPoint(const Point3D &qPoint, Point3D &psurf) const
{
     int err;
     double xon, yon, zon;

     err = geom->getEntClosestPt(gEdgeHandle, qPoint[0], qPoint[1], qPoint[2],
                                 xon, yon, zon);
     assert(!err);

     if (!err) {
          psurf[0] = xon;
          psurf[1] = yon;
          psurf[2] = zon;
          return 0;
     }

     double u = getUCoord(qPoint);

     if (u == umin || u == umax)
          cout << " Warning: Selecting end points is dangerous " << endl;

     return getXYZCoords(u, psurf);
}
//////////////////////////////////////////////////////////////////////////////

int GEdge::getFirstDer(double u, Vec3D &avec) const
{
     if (u < umin || u > umax)
          cout << "Warning: U value " << u << " Out of Range : " << umin << " " << umax << endl;

     int err;
     double nx, ny, nz;
     double du = 1.0E-06 * (umax - umin);
     double uc, xc, yc, zc; // Center
     double uf, xf, yf, zf; // Forward
     double ub, xb, yb, zb; // backward

     uc = u;
     err = geom->getEntUtoXYZ(gEdgeHandle, uc, xc, yc, zc);

     if (u == umin) {
          uf = u + du;
          err = geom->getEntUtoXYZ(gEdgeHandle, uf, xf, yf, zf);
          nx = (xf - xc) / (uf - uc);
          ny = (yf - yc) / (uf - uc);
          nz = (zf - zc) / (uf - uc);
          avec[0] = nx;
          avec[1] = ny;
          avec[2] = nz;
          return 0;
     }

     if (u == umax) {
          ub = u - du;
          err = geom->getEntUtoXYZ(gEdgeHandle, ub, xb, yb, zb);
          nx = (xb - xc) / (ub - uc);
          ny = (yb - yc) / (ub - uc);
          nz = (zb - zc) / (ub - uc);
          avec[0] = nx;
          avec[1] = ny;
          avec[2] = nz;
          return 0;
     }

     uf = u + du;
     err = geom->getEntUtoXYZ(gEdgeHandle, uf, xf, yf, zf);

     ub = u - du;
     err = geom->getEntUtoXYZ(gEdgeHandle, ub, xb, yb, zb);

     nx = (xf - xb) / (ub - uf);
     ny = (yf - yb) / (ub - uf);
     nz = (zf - zb) / (ub - uf);

     avec[0] = nx;
     avec[1] = ny;
     avec[2] = nz;

     return 0;
}

///////////////////////////////////////////////////////////////////////////////

void GEdge::build_kdtree()
{
     delete_kdtree(); // If Allocated earlier

     int err;

     double x, y, z;

     double du = (umax - umin) / (double) (numNodes - 1);

     kdnodes = annAllocPts(numNodes, 3);

     uCoord = new double[numNodes];

     int index = 0;
     for (int i = 0; i < numNodes; i++) {
          double u = umin + i*du;
          if (u > umax) u = umax;
          err = geom->getEntUtoXYZ(gEdgeHandle, u, x, y, z);
          assert(!err);
          kdnodes[index][0] = x;
          kdnodes[index][1] = y;
          kdnodes[index][2] = z;
          uCoord[index] = u;
          index++;
     }

     kdtree = new ANNkd_tree(kdnodes, numNodes, 3);
     assert(kdtree);
}

///////////////////////////////////////////////////////////////////////////////

void GEdge::delete_kdtree()
{
     if (kdtree == NULL) return;

     annDeallocPts(kdnodes);
     delete [] uCoord;
     delete kdtree;
     kdtree = NULL;
     annClose();
}

///////////////////////////////////////////////////////////////////////////////

bool GEdge::is_u_uniformly_discretized()
{
}

///////////////////////////////////////////////////////////////////////////////

void GEdge::uniform_u_discretization(int n, vector<Point3D> &xyz, vector<double> &u)
{
     assert( n >= 2 );
     double du = (umax - umin) / (double) (n - 1);

     u.resize(n);
     xyz.resize(n);
     for (int i = 0; i < n; i++) {
          u[i] = umin + i*du;
          getXYZCoords(u[i], xyz[i] );
     }
}

///////////////////////////////////////////////////////////////////////////////

void GEdge::uniform_xyz_discretization(int n, vector<Point3D> &xyz, vector<double> &u)
{
     u.resize(n);
     xyz.resize(n);

     double len = getLength();
     double dl = len / (double) (n - 1);
     double du = 0.90 * (umax - umin) / (double) (n - 1);

     u[0] = umin;
     getXYZCoords(umin, xyz[0] );

     for (int i = 1; i < n - 1; i++) {
          double ufrom = u[i - 1];
          double uguess = ufrom + du;
          getUCoord(ufrom, dl, uguess);
          u[i] = uguess;
          getXYZCoords(u[i], xyz[i] );
     }

     u[n - 1] = umax;
     getXYZCoords(umax, xyz[n-1] );
}

///////////////////////////////////////////////////////////////////////////////

void GEdge::projectHigherOrderNodes(const vector<double> &gnodes)
{
     int err;
     gllnodes = gnodes;

     if( mesh == 0 ) {
          cout << " Warning: No mesh is present in geometric edge : Higher order nodes projection not done " << endl;
          return;
     }

     if( relPair == 0 ) {
          cout << " Warning: No assoc is present in geometric edge : Higher order nodes projection not done " << endl;
          return;
     }

     if( rel == 0 ) {
          cout << " Warning: No relation is present in geometric edge : Higher order nodes projection not done " << endl;
          return;
     }

     const char *tag1 = "HO_POINTS";
     err = mesh->getTagHandle(tag1, horder_tag);
     assert(!err);

     int nsize = gllnodes.size();
     arclength_ratio.resize( nsize );
     for (int i = 0; i < nsize; i++)
          arclength_ratio[i] = 0.5 * (1.0 + gllnodes[i]);


     iBase_EntitySetHandle meshSet;
     err = relPair->getEntSetRelation(gEdgeHandle, 0, meshSet);
     assert(!err);

     iBase_TagHandle dim_tag;
     const char *tag2 = "GEOM_DIMENSION";
     err = mesh->getTagHandle(tag2, dim_tag );
     assert(!err);

     int geom_dim;
     err = mesh->getEntSetIntData(meshSet, dim_tag, geom_dim);
     assert(!err);

     if (geom_dim == 1) {
          vector<iBase_EntityHandle> mEdges;
          err = mesh->getEntities(meshSet, iBase_EDGE, iMesh_ALL_TOPOLOGIES, mEdges);
          assert(!err);

          for (int i = 0; i < mEdges.size(); i++) projectHigherOrderNodes(mEdges[i] );
     }
}

///////////////////////////////////////////////////////////////////////////////

void GEdge::projectHigherOrderNodes(iBase_EntityHandle mEdgeHandle)
{
     int err;

     vector<iBase_EntityHandle> edgenodes;
     Point3D pnear, pon;
     double u, u0, uN, unear, arclen, arcdist, dl;

     char *tag_val = NULL;
     err = mesh->getData(mEdgeHandle, horder_tag, &tag_val);
     assert(!err);

     HO_Points *hopoints = (HO_Points *) tag_val;
     iBase_EntityHandle *nodesOnEdge = hopoints->nodeHandles;
     int numHPoints = hopoints->nx;
     assert( gllnodes.size() ==  numHPoints );

     int nhalf = numHPoints / 2;

     edgenodes.clear();
     err = mesh->getEntAdj(mEdgeHandle, iBase_VERTEX, edgenodes);
     assert( !err );
     assert(edgenodes.size() == 2);

     err = mesh->getVtxCoord(edgenodes[0], p0[0], p0[1], p0[2]);
     assert(!err);

     err = mesh->getVtxCoord(edgenodes[1], p1[0], p1[1], p1[2]);
     assert(!err);

     // End points can be ambiguous, but not the internal points
     u = 0.10;
     pnear[0] = (1.0-u)*p0[0] + u*p1[0];
     pnear[1] = (1.0-u)*p0[1] + u*p1[1];
     pnear[2] = (1.0-u)*p0[2] + u*p1[2];
     unear = getUCoord(pnear);
     u0 = getUCoord(p0, unear);

     u = 0.90;
     pnear[0] = (1.0-u)*p0[0] + u*p1[0];
     pnear[1] = (1.0-u)*p0[1] + u*p1[1];
     pnear[2] = (1.0-u)*p0[2] + u*p1[2];
     unear = getUCoord(pnear);
     uN    = getUCoord(p1, unear);

     arclen = getLength(u0, uN);

     arcdist = 0.0;
     unear = u0;
     for (int i = 1; i < nhalf; i++) {
          iBase_EntityHandle currvertex = nodesOnEdge[i];
          dl = arclength_ratio[i] * arclen - arcdist;
          u = gllnodes[i];
          u = 0.5 * (1 - u) * u0 + 0.5 * (1 + u) * uN;
//      getUCoord(unear, dl, u);
          getXYZCoords(u, pon);
          err = mesh->setVtxCoord(currvertex, pon[0], pon[1], pon[2]);
          assert( !err );
          unear = u;
          arcdist += dl;
     }

     // From the end->nhalf: Here is the big assumption that Gauss Node are
     // Symmetric from both the ends.
     arcdist = 0.0;
     unear = uN;
     for (int i = 1; i < nhalf; i++) {

          iBase_EntityHandle currvertex = nodesOnEdge[numHPoints - 1 - i];
          dl = arclength_ratio[i] * arclen - arcdist;
          u = gllnodes[i];
          u = 0.5 * (1 - u) * uN + 0.5 * (1 + u) * u0; // careful, u0, uN swapped
//      getUCoord(unear, dl, u);
          getXYZCoords(u, pon);
          err = mesh->setVtxCoord(currvertex, pon[0], pon[1], pon[2]);
          assert( !err );
          unear = u;
          arcdist += dl;
     }

     if (numHPoints % 2) {
          int midpos = nhalf;
          iBase_EntityHandle currvertex = nodesOnEdge[midpos];
          u = gllnodes[midpos + 1];
          u = 0.5 * (1 - u) * u0 + 0.5 * (1 + u) * uN;
          getXYZCoords(u, pon);
          err = mesh->setVtxCoord(currvertex, pon[0], pon[1], pon[2]);
          assert( !err );
     }

}

///////////////////////////////////////////////////////////////////////////////

