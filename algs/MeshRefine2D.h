#ifndef REFINE_H
#define REFINE_H

#include <assert.h>
#include <math.h>
#include <string.h>
#include <algorithm>

#include <iGeom.h>
#include <iMesh.h>
#include <iRel.h>

#include <vector>
#include <bitset>

using namespace std;

#include "SimpleArray.h"
#include "basic_math.h"

/////////////////////////////////////////////////////////////////////////////////
//   MeshRefine2D:
//   		LongestEdgeRefine2D: ConsistencyRefine2D
//   		ObtuseRefine2D: ConsistencyRefine2D
//   		GradeRefine2D: ConsistencyRefine2D
//   		Refine2D14:  ConsistencyRefine2D
//   		CentroidRefine2D
//   		DelaunayRefine2D
//   		Square3SubDivision
////////////////////////////////////////////////////////////////////////////////
//Class Refine, refines an existing triangulation. This class works for 2D cases
//or the parametric domain of 3D surface triangulation. The user need to provide
//which cells need to be refined. There could be many criteria for deviding the
//cell, so this is not part of the class. 
//
//Input :   Connection :  Connectivity of the triangulation ( 1D array)
//      :   ParamCoords:  parametric coordinates of the vertices.
//      :   markFace   :  Faces which are marked for refinement.
//
//Output:   Connection :
//          ParamCoords:
//          
//If the user is working with 3D triangulated surface, it is assumed that he 
//has some representation ( example NURBS ) from which he can calculate the
//physical Coordinates from the paramCoords. Again this is not the responsiblity
//of the class.
//
//------------------------------------------------------------------------------
//Number of Division     :   New Vertices           :  New Faces
//------------------------------------------------------------------------------
//    2                  :     1                    :  1
//    3                  :     0                    :  2
//    4                  :     3                    :  3
//------------------------------------------------------------------------------
// Default Number of Subdivisions is 2, which is conservative approach, but 
// good at equalizing the aspect ratio with respect to neighbouring triangles.
//
// Programmer : Chaman Singh Verma
// Place      : Argonne National Lab.
//              Argonne, IL, USA
//
// 
////////////////////////////////////////////////////////////////////////////////
//
/**
 * REFINE AREA         : Increase the density where grid cells have high area/volume
 * REFINE ASPECT_RATIO : Increase the aspect Ratio
 * REFINE CURVATURE    : Create high density mesh near high curvature.
 */

enum RefinePolicy { CENTROID_PLACEMENT, CIRCUMCENTER_PLACEMENT, LONGEST_EDGE_BISECTION};

enum RefineObjective {REFINE_AREA, REFINE_ASPECT_RATIO, REFINE_CURVATURE};

///////////////////////////////////////////////////////////////////////////////

//! \brief 2D Mesh Refinement class.
class MeshRefine2D
{
 public:

  MeshRefine2D() { 
     boundary_split_flag = 0; 
  }
  virtual ~MeshRefine2D() {}

  void setMesh(  iMesh_Instance &m ) {  mesh = m; }
  void setGeometry(  const iGeom_Instance &g ) { geom = g; }

  void setBoundarySplitFlag( bool f ) { boundary_split_flag = f; }

  vector<iBase_EntityHandle> getNewNodes() const { return insertedNodes; }
  vector<iBase_EntityHandle> getNewFaces() const { return insertedFaces; }

  size_t  getNumFacesRefined() const { return numfacesRefined; }

  virtual int execute() = 0;

  virtual int initialize();
  virtual int finalize();

 protected:
    typedef iBase_EntityHandle EHandle;
    int   numIterations;

    iBase_TagHandle remove_tag, vertex_on_edge_tag, boundary_tag, globalID_tag;

    iBase_EntitySetHandle rootSet, currSet;
    vector<iBase_EntitySetHandle>  entitySets;

    iMesh_Instance mesh;
    iGeom_Instance geom;

    vector<EHandle>  hangingVertex; 
    vector<EHandle>  insertedNodes;
    vector<EHandle>  insertedFaces;
    
    bool    boundary_split_flag;                
    size_t  numfacesRefined;

    int setVertexOnEdge(EHandle v1, EHandle v2, EHandle &vmid);
    int getVertexOnEdge(EHandle v1, EHandle v2, EHandle &vmid) const;

    bool searchEdge( EHandle v1, EHandle v2, EHandle &edgehandle) const; 
    bool allow_edge_refinement( EHandle &edgehandle ) const;

    double  edge_length( EHandle v1, EHandle v2)  const;
    Point3D edge_centroid( EHandle v1, EHandle v2) const;

    double  face_aspect_ratio( EHandle v1 ) const;
    Point3D face_centroid( EHandle f) const;

    EHandle create_new_edge( EHandle v1, EHandle v2 );

    int prune_mesh();

};

///////////////////////////////////////////////////////////////////////////////

class Sqrt3Refine2D : public MeshRefine2D
{
  public:
   int execute();
  private:
   int create_tags();
};

///////////////////////////////////////////////////////////////////////////////

class CentroidRefine2D : public MeshRefine2D
{
 public:

   int  initialize();
   int  finalize();
   int  execute();

 private:
  int atomicOp( const iBase_EntityHandle &f);
  int refine_tri( const iBase_EntityHandle &f);
  int refine_quad( const iBase_EntityHandle &f);
};

///////////////////////////////////////////////////////////////////////////////

class ConsistencyRefine2D : public MeshRefine2D
{
   public:
     ~ConsistencyRefine2D() {}

     int  initialize();
     int  finalize();
     int  execute();

   private:
     bitset<3>  edge0, edge1, edge2, bitvec;
     typedef iBase_EntityHandle EHandle;

     int create_tags();

     void  atomicOp( const EHandle &f);
     void  refineEdge0(const EHandle &f);
     void  refineEdge1(const EHandle &f);
     void  refineEdge2(const EHandle &f);
     void  subDivideQuad2Tri( const vector<EHandle>  &qnodes);

     void  checkFaceConsistency( const EHandle &f);

     void  makeConsistent();
     void  makeConsistent1( const EHandle &f );
     void  makeConsistent2( const EHandle &f );
     void  makeConsistent3( const EHandle &f );
};


///////////////////////////////////////////////////////////////////////////////

class LongestEdgeRefine2D : public MeshRefine2D 
{
 public:
  LongestEdgeRefine2D()
  { 
    cutOffAspectRatio = 0.50; 
  }

  ~LongestEdgeRefine2D() {}

  void setCutOffAspectRatio(double asp) { cutOffAspectRatio = asp;}

  int  initialize();
  int  finalize();
  int  execute();

 private:
  int create_tags();
  double cutOffAspectRatio;
  void atomicOp( const iBase_EntityHandle &f);
};

///////////////////////////////////////////////////////////////////////////////

class Refine2D14 : public MeshRefine2D 
{
 public:

  ~Refine2D14() {}

  int  initialize();
  int  finalize();
  int  execute();

 private:
  int  atomicOp( const iBase_EntityHandle &f);
  int  refine_tri(const iBase_EntityHandle &f);
  int  refine_quad(const iBase_EntityHandle &f);
};

///////////////////////////////////////////////////////////////////////////////

struct DelaunayRefinement2D : public MeshRefine2D
{
  ~DelaunayRefinement2D() {}

  int  initialize();
  int  finalize();
  int  execute() {}
};

///////////////////////////////////////////////////////////////////////////////

class ObtuseRefine2D : public MeshRefine2D
{
  public:
   ObtuseRefine2D( ) { cutoffAngle = 90.0;}

   void setCutOffAngle( double a ) { cutoffAngle = std::max(90.0, a); }

   int  initialize();
   int  finalize();
   int  execute();

  private:
   double cutoffAngle;
   int  atomicOp(const iBase_EntityHandle &f);
};

///////////////////////////////////////////////////////////////////////////////

class GradeRefine2D : public MeshRefine2D
{
  public:
   int  initialize();
   int  finalize();
   int  execute();

  private:
     int atomicOp( const iBase_EntityHandle &v);
};

///////////////////////////////////////////////////////////////////////////////

#endif

