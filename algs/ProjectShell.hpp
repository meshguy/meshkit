
// Dec 07, 2009
// It reads a shell mesh file 
// makes some tests on the shell mesh
// projects on a plane perpendicular to the direction
// computes the intersected mesh

#ifndef PROJECTSHELL_HPP
#define PROJECTSHELL_HPP

#include <vector>
#include <map>
#include <list>
#include <sys/resource.h>

#include "iMesh.h"

// those are constructs used to check the mesh in 3D
struct PSNode {
  int id; // from 0, into vert_adj array
  double xyz[3];
  std::list<int> edges;
};

struct PSEdge {
  int id; // from 0, into new created ps_edges array
  int v[2]; // index in adj_vert_coord array // used for construction of the edge array
  // always v[0] < v[1];
  int used;// count how many times is used in the mesh; more than 2 times is troublesome
  int t[2]; // triangle indices ; init with 0
  std::list<int> extraNodes; // they will be intersection nodes, different from original projection nodes
};

struct PS3DTriangle
{
  int v[3];
  int e[3]; // edge indices: can be negative, but not 0 (orientation is important)
  int t[3]; // adjacent triangles
};

struct PSTriangle2D 
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

//#include "iBase.h"
//#include "iMesh.h"

class ProjectShell
{
public:

  //explicit ProjectShell(iMesh_Instance mesh,
  ProjectShell(iMesh_Instance mesh,
	      iBase_EntitySetHandle root_set,
	      double direction[3]);
  
  //virtual ~ProjectShell();
  ~ProjectShell();

  iMesh_Instance mesh_impl() const
  {
    return m_mesh;
  }

  int project();

  int writeNewMesh(iMesh_Instance mesh);
  
private:
 
  int  getMeshData(); // read input mesh
  int  checkMeshValidity();

  int  projectIn2D();
  int computeIntersections();

  // this method computed intersection between 2 triangles: will output n points, area, affected sides
  int computeIntersectionBetweenRedAndBlue(int red, int blue, double * opPoints, int & oNPoints, double & area,
					   int * opSides); 
  // this method will add extra points (P) for intersection points
  // they will stay on the red edgesi; will create the final triangles (in 2D) and
  //  extra nodes (if needed)
  //
  int findNodes(int red, int blue, double * iP, int nP);
  iMesh_Instance m_mesh;
  iBase_EntitySetHandle m_hRootSet;

  double m_direction[3];
  double m_dirX[3];
  double m_dirY[3];
  
  int m_numNodes;
  int m_numTriangles;
  double * m_xyz ; // original coordinates
  int * m_triangles; // original triangles

  PSEdge  * ps_edges ;
  int m_numEdges;
  PSNode * ps_nodes ;
  PS3DTriangle * ps_triangles;

  double * m_xy; // 2d coordinates
  PSTriangle2D * m_redMesh;
  PSTriangle2D * m_blueMesh;
  std::vector <FinalTriangle> m_finalMesh;
  int m_numCurrentNodes;
  int m_numFinalTriangles;
  int m_numPos, m_numNeg;  
  int m_num2dPoints;
  int m_2dcapacity;
};

#endif  // PROJECTSHELL_HPP
