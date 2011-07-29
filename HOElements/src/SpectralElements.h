#ifndef HORDER_H
#define HORDER_H

#include <string>
#include <string.h>
#include <assert.h>
#include <limits>
#include <math.h>
#include <fstream>
#include <sstream>

#include <vector>
#include <map>
#include <set>

#include <boost/foreach.hpp>

#include <meshkit/iGeom.hpp>
#include <meshkit/iMesh.hpp>
#include <meshkit/iRel.hpp>
#include <MBCN.hpp>

#include "GEdge.h"
#include "GFace.h"
#include "tfiblend.h"
#include "TFIMap.h"

////////////////////////////////////////////////////////////////////////////////
// Objective:  Generate Higher order elements ( Gauss-Lobatto ) on an unstructured
//             Quadrilateral and hexahedral mesh.
//
// Notes for Users:
//             Higher order nodes distance is measured in physical space, that
//             means that the spacing is done in physical rather than parametric
//             space.
//
//             If the geometry is available then the nodes are projected on the
//             underlying edges and faces.
//
//             If the geometric face has structured mesh, then TFI is used to
//             modify the mesh and generation of higher order nodes.
//
//             If the geometric region has structured mesh, then TFI is used
//             to modify internal nodes.
//
//             If geometry is not available, then all the calculations are done
//             in physical space and at present, no projection on quadratic
//             surfaces is done.
//
//             Although the mesh is likely to produce correct results, but later
//             some quality certificates modules will be developed to check the
//             final higher order elements.
//
// Notes for Developers:
//
//
// Developed by:  Chaman Singh Verma
//                Argonne National Lab. IL, Chicago.
// Date:          10th August, 2009.
//
// Dependencies:  MOAB, CGM, and Boost Libraries.
//
//
////////////////////////////////////////////////////////////////////////////////

class GFace;

class SpectralElements {
public:

     static void gauss_linear_nodes(int N, vector<double> &x);
     static void gauss_lobatto_nodes(int n, vector<double> &x);
     static void bilinear_weights(double xi, double eta, vector<double> &weight);
     static void trilinear_weights( double x, double y, double z, vector<double> &weight);
     static double linear_interpolation( const vector<double> &x, const vector<double> &w);
     static int canonical_quad_coords( const vector<Point2D> &uvCorners, const Point2D &qPoint,
                                       Point2D &xy);

     SpectralElements();

//  Brick 20 is special subproblem, which can be done much faster than the general
//  solution.
     void brick20(iMesh *m);

     void brick20( iMesh *m,
                   iGeom *g, iRel *r, iRel::PairHandle *a);

//  Generate higher order elements of a given order without geometry.
     void generate( iMesh *n, int order );

//  Generate higher order elements and project on the underlying geometry.
     void generate( iMesh *m, int order,
                    iGeom *g, iRel *r, iRel::PairHandle *a);

//  Should a monolithic hexelements be created at the end ?
     void create_subelements( bool b);

//  Save the mesh. At present, we create a single monolithic hex elements out of
//  higher order elements. The reason is that the most visualization software don't
//  support higher order elements.
     void saveAs(const string &s);

private:

     struct EntityNodes {
          EntityNodes() {
               projected = 0;
          }

          void add( iBase_EntityHandle &h) {
               nodes.push_back(h);
          }
          void clear() {
               nodes.clear();
          }

          int projected;
          vector<iBase_EntityHandle> nodes;
     };

     iMesh  *mesh;
     iGeom  *geom;
     iRel   *rel;
     iRel::PairHandle  *relPair;

     map<iBase_EntityHandle, int>  structured_mesh;

     int  subCellDim[3];
     bool hasGeometry;
     bool subElements;
     bool hasSubStructuredGrid;

     iBase_TagHandle idtag;
     iBase_TagHandle horder_tag;
     iBase_EntitySetHandle meshRootSet;

     vector<double> gllnodes;

     std::map<iBase_EntityHandle, EntityNodes> entityNodesMap;

     void init();

     void project_on_edges( );
     void project_on_faces( );

     void brick20();

     void linear_edge( iBase_EntityHandle h,   int numPoints);
     void bilinear_face( iBase_EntityHandle h, int numPoints);
     void trilinear_cell( iBase_EntityHandle h,int numPoints);

     void project_on_quad_face( GFace &f, iBase_EntityHandle &mHandle);

     void   linear_interpolation( iBase_EntityHandle h, double xi, double eta, double zeta,
                                  double &x, double &y, double &z);
     void arrange_subcell_nodes(iBase_EntityHandle h, vector<iBase_EntityHandle> &n);

     int  getIndex( int *dim, int i, int j, int k);

     int check_edge_discretization();

     void arrange_brick20_nodes(iBase_EntityHandle h, vector<iBase_EntityHandle> &nodes);

     bool hasStructuredMesh2D(iBase_EntityHandle gface) const;
     bool hasStructuredMesh3D(iBase_EntityHandle gcell) const;

     void setTFIMesh();

     int  matching_node( const vector<Point3D> &p, const Point3D &q);
     void get_subcells(iBase_EntityHandle handle, vector<int> &connect);
     void get_subedges(iBase_EntityHandle handle, vector<int> &connect);
};

#endif
