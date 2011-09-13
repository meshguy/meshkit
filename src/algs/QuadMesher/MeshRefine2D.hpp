#ifndef REFINE_H
#define REFINE_H

#include <bitset>
#include "Mesh.hpp"
#include "SwapTriEdge.hpp"

#include "basic_math.hpp"

using namespace std;

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

namespace Jaal 
{
enum RefinePolicy { CENTROID_PLACEMENT, CIRCUMCENTER_PLACEMENT, LONGEST_EDGE_BISECTION};

enum RefineObjective {REFINE_AREA, REFINE_ASPECT_RATIO, REFINE_CURVATURE};

///////////////////////////////////////////////////////////////////////////////

//! \brief 2D Mesh Refinement class.
class MeshRefine2D
{
 public:

  MeshRefine2D() { 
     boundary_split_flag = 0; 
     numIterations   = 1;
  }
  virtual ~MeshRefine2D() {}

  void setMesh( Mesh *m ) {  mesh = m; }

// void setGeometry(  const iGeom_Instance &g ) { geom = g; }

  void setBoundarySplitFlag( bool f ) { boundary_split_flag = f; }

  const NodeSequence &getNewNodes() const { return insertedNodes; }
  const FaceSequence &getNewFaces() const { return insertedFaces; }

  size_t  getNumFacesRefined() const { return numfacesRefined; }

  virtual int execute() = 0;

  virtual int initialize();

  void setNumOfIterations( int i ) { numIterations = i; }

  // Set Desired Mesh Quality 
  void setAspectRatio( double a )  { desiredAspectRatio = a; }
  void setDesiredArea( double a )  { desiredArea     = a; }
  void setMinimumAngle( double a ) { desiredMinAngle = a; }
  void setMaximumAngle( double a ) { desiredMaxAngle = a; }
  void setFeatureAngle( double a ) { featureAngle    = a; }
  void setMaximumCells( size_t a ) { maxAllowedCells = a; }

 protected:
    Mesh *mesh;

    class RefinedEdgeMap
    {
       public:
        void clear();

        bool hasEdge( Vertex *v1, Vertex *v2) const; 
        bool allow_edge_refinement( const Edge *edge) const;

        Vertex* setVertexOnEdge(Vertex *v1, Vertex *v2);
        Vertex* getVertexOnEdge(Vertex *v1, Vertex *v2) const;
        bool  boundary_split_flag;                

	// Get All the inserted vertices on the edges..
	NodeSequence getInsertedNodes() const
	{
	    NodeSequence result;
            std::map<Vertex*, vector<RefinedEdge> >::const_iterator it;
	    for( it = refined_edges.begin(); it != refined_edges.end(); ++it) 
	    {
	        const vector<RefinedEdge> &refedges = it->second;
		for( size_t i = 0; i < refedges.size(); i++) 
		     result.push_back( refedges[i].midVertex );

            }
	    return result;
	}
       private:
        struct RefinedEdge
        {
            Edge    *edge;
            Vertex  *midVertex;
        };
        Edge* create_new_edge( const Vertex *v1, const Vertex *v2 );
        std::map<Vertex*, vector<RefinedEdge> > refined_edges;
    };

    RefinedEdgeMap *edgemap;

    FaceSequence  insertedFaces;
    NodeSequence  hangingVertex, insertedNodes;
    
    int     numIterations;
    bool    boundary_split_flag;                
    size_t  numfacesRefined;

    double  desiredAspectRatio, desiredArea;
    double  desiredMinAngle, desiredMaxAngle;
    double  featureAngle;
    size_t  maxAllowedCells;

    int finalize();

    void  append_new_node( Vertex *v0 );
    Face* append_new_triangle(Vertex *v0, Vertex *v1, Vertex *v2);
    Face* append_new_quad(Vertex *v0, Vertex *v1, Vertex *v2, Vertex *v3);

    void remove_it(Face *face) {
        face->setStatus( MeshEntity::REMOVE);
    }
};

///////////////////////////////////////////////////////////////////////////////

struct Sqrt3Refine2D : public MeshRefine2D
{
   Sqrt3Refine2D() {}
   Sqrt3Refine2D(Mesh *m) { setMesh(m); }

   int execute();
};

///////////////////////////////////////////////////////////////////////////////

class CentroidRefine2D : public MeshRefine2D
{
 public:
   CentroidRefine2D() {}
   CentroidRefine2D( Mesh *m ) { setMesh(m); }

   void refine( Face *f ) {
        atomicOp( f );
        mesh->prune();
   }

   void refine( Mesh *m ) {
       setMesh( m );

       size_t nSize = mesh->getSize(2);
       for( size_t i = 0; i < nSize; i++) 
            atomicOp( mesh->getFaceAt(i)  );
       mesh->prune();
   }

   int  execute();
 private:
  int atomicOp(  Face *f);
  int refine_tri( Face *f);
  int refine_quad( Face *f);
};

///////////////////////////////////////////////////////////////////////////////

class LongestEdgeRefine2D : public MeshRefine2D 
{
 public:
  LongestEdgeRefine2D()
  { 
    cutOffAspectRatio = 0.50; 
  }

  LongestEdgeRefine2D(Mesh *m) { 
      setMesh(m); 
      cutOffAspectRatio = 0.50; 
  }

  ~LongestEdgeRefine2D() {}

  void setCutOffAspectRatio(double asp) { cutOffAspectRatio = asp;}

  int  execute();

 private:
  double cutOffAspectRatio;
  int  atomicOp( const Face *face);
};

///////////////////////////////////////////////////////////////////////////////

class ConsistencyRefine2D : public MeshRefine2D
{
   public:
     ConsistencyRefine2D() { }
     ConsistencyRefine2D( Mesh *m, RefinedEdgeMap *emap) 
        { setMesh(m); edgemap = emap;}

     ~ConsistencyRefine2D() {}

     int  execute();

   private:
     bitset<3>  edge0, edge1, edge2, bitvec;

     int   atomicOp( Face *f);
     void  refineEdge0(const Face *f);
     void  refineEdge1(const Face *f);
     void  refineEdge2(const Face *f);

     void  subDivideQuad2Tri( const NodeSequence &qnodes);
     void  makeConsistent1( Face *f );
     void  makeConsistent2( Face *f );
     void  makeConsistent3( Face *f );
     void  makeConsistent();
     void  checkFaceConsistency( Face *f);
};

///////////////////////////////////////////////////////////////////////////////

class Refine2D14 : public MeshRefine2D 
{
 public:

  ~Refine2D14() {}

  int  initialize() { return 0; }
  int  execute();

 private:
  int  atomicOp( Face *f);
  int  refine_tri( Face *f);
  int  refine_quad( Face *f);
};

///////////////////////////////////////////////////////////////////////////////

struct DelaunayRefinement2D : public MeshRefine2D
{
  ~DelaunayRefinement2D() {}

  int  initialize();
  int  finalize();
  int  execute() {

       return 0;
   }
};

///////////////////////////////////////////////////////////////////////////////

class ObtuseRefine2D : public MeshRefine2D
{
  public:
   ObtuseRefine2D( ) { cutoffAngle = 90.0;}

   void setCutOffAngle( double a ) { cutoffAngle = std::max(90.0, a); }

   int  initialize();
   int  execute();

  private:
   double cutoffAngle;
   int   atomicOp(const Face *f);
};

///////////////////////////////////////////////////////////////////////////////

class GradeRefine2D : public MeshRefine2D
{
  public:
   int  initialize();
   int  finalize();
   int  execute();

  private:
     int atomicOp( const Vertex *v);
};

///////////////////////////////////////////////////////////////////////////////

}

#endif

