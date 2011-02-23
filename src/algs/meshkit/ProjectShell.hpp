//-----------------------------------C++-------------------------------------//
// File: src/algs/ProjectShell.hpp
// Wednesday February 11 10:50 2011
// Brief: ProjectShell class definition: one scheme is provided to project a 2D
//        surface, the shell, bounding a convex 3D region onto a plane.
//---------------------------------------------------------------------------//

#ifndef MESHKIT_PROJECTSHELL_HPP
#define MESHKIT_PROJECTSHELL_HPP

#include "meshkit/MeshScheme.hpp"
#include "meshkit/iGeom.hpp"
#include "meshkit/Matrix.hpp"
#include "iMesh.h"

#include <vector>
#include <map>
#include <list>
#include <sys/resource.h>
#include <fstream>


  //===========================================================================//
  /*!
   * \class ProjectShell
   * \brief Project the bounding shell of a 3D convex volume onto a plane.
   * 
   * ProjectShell projects a 2D surface (the "shell"), which is assumed to bound
   * a convex 3D domain, orthogonally onto a plane perpendicular to a given direction.
   * The direction is specified by a nonzero vector, with all collinear vectors being
   * equivalent, up to the reflection of the projected mesh.
   */
  //===========================================================================//
namespace MeshKit {
class ProjectShell : public MeshScheme {
public:
  // Constructor
  ProjectShell(MKCore *mk_core, const MEntVector &me_vec = MEntVector());
  
  //return the type of entities in the created mesh
  void mesh_types(std::vector<moab::EntityType> &tps);
  
  //set up the projection parameters: the direction of projection
  virtual void setup_this();
  
  // Generate the projected mesh
  virtual void execute_this();
  
  // Destructor
  virtual ~ProjectShell();
  

  // Parameters controlling the behavior of ProjectShell
  Vector<3>& direction() {return m_direction;};
  bool&      dbg()       {return m_dbg;};

private:
  // These are auxiliary constructs used to check the mesh in 3D
  struct Node {
    int id; // from 0, into vert_adj array
    Vector<3> xyz;
    std::list<int> edges;
  };
  
  struct Edge {
    int id; // from 0, into new created ps_edges array
    int v[2]; // index in adj_vert_coord array // used for construction of the edge array
    // always v[0] < v[1];
    int used;// count how many times is used in the mesh; more than 2 times is troublesome
    int t[2]; // triangle indices ; init with 0
    std::list<int> extraNodes; // they will be intersection nodes, different from original projection nodes
  };
  
  struct Triangle3D
  {
    int v[3];
    int e[3]; // edge indices: can be negative, but not 0 (orientation is important)
    int t[3]; // adjacent triangles
  };
  
  struct Triangle2D 
  {
    int v[3]; // vertices: indices in the m_xy array; start from 0
    int e[3]; // edges that are indices in original ps_edge array; these are used only for blue triangles (or red?)
    // the new extra intersection points will belong to edges of the blue triangles
    // or the alternate is to create points every time, and merge at the end the nodes?
    int t[3]; // neighbors: boundary or other of same color // triangles start from 0 too
    // boundary is flagged with _numPos (neg) + 1;
    // the indices are in this array (from 0 to  
    int oldId; // the index in the ps_triangle array; some triangles are eliminated because they project to 0 area
    double area; // will keep an area of the triangle;
    // it would be used mostly for debugging purposes
  };
  
  struct FinalTriangle
  {
    int v[3];  // these are new vertices that may be formed by intersections
    // they will be indices in the new array of vertices
    // 
    int redTriangle, blueTriangle;// intersection between 2 triangles is a convex polygon with at most 6 edges
    // it can be decomposed in at most 4 triangles; we do not care about quality
  };
  
  
  int checkMeshValidity();
  
  int projectIn2D();
  int computeIntersections();
  
  int getEdge(Node & v1, Node & v2, int & edgeId, Edge * edgesArr, int sizeArr);
  
  // this method computed intersection between 2 triangles: will output n points, area, affected sides
  int computeIntersectionBetweenRedAndBlue(int red, int blue, double * opPoints, int & oNPoints, double & area,
                                           int * opSides); 
  // this method will add extra points (P) for intersection points
  // they will stay on the red edgesi; will create the final triangles (in 2D) and
  //  extra nodes (if needed)
  //
  int findNodes(int red, int blue, double * iP, int nP);

  int edgeIntersections(double * red, double * blue, int mark[3], double * points, int & nPoints);

  double area2D(Vector<2> a, Vector<2> b, Vector<2> c) {
    Matrix<2,2> M;
    M.set_column(1,b-a);
    M.set_column(2,c-a);
    return det(M);
  }

private:
  Vector<3> m_direction;
  Vector<3> m_dirX;
  Vector<3> m_dirY;
  
  int m_numNodes;
  int m_numTriangles;
  double * m_xyz ; // original coordinates
  int * m_triangles; // original triangles
  
  Edge  * ps_edges ;
  int m_numEdges;
  Node * ps_nodes ;
  Triangle3D * ps_triangles;
  
  double * m_xy; // 2d coordinates
  Triangle2D * m_redMesh;
  Triangle2D * m_blueMesh;
  std::vector <FinalTriangle> m_finalMesh;
  int m_numCurrentNodes;
  int m_numFinalTriangles;
  int m_numPos, m_numNeg;  
  int m_num2dPoints;
  int m_2dcapacity;

  bool m_dbg;
  std::ofstream m_dbg_out;
};// class ProjectShell
}// namespace MeshKit

#endif  //  MESHKIT_PROJECTSHELL_HPP
