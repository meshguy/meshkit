#ifndef MESHKIT_SMOOTH_FACE_EVAL_HPP
#define MESHKIT_SMOOTH_FACE_EVAL_HPP

// do we really need iMesh here; maybe go directly to MOAB
//#include "iMesh.h"
#include "MBInterface.hpp"
#include "MBRange.hpp"
#include "MBCartVect.hpp"
#include "MBTagConventions.hpp"

#define determ3(p1,q1,p2,q2,p3,q3) ((q3)*((p2)-(p1)) + (q2)*((p1)-(p3)) + (q1)*((p3)-(p2)))
#define sqr(a) ((a)*(a))
#define cube(a) (sqr(a) * (a))
#define quart(a) (sqr(a) * sqr(a))
#define blend(x) (-2.0*(x)*(x)*(x) + 3.0*(x)*(x))

#include <vector>
#include <map>
//#include "MBEntityHandle.hpp"

// work only with CAMAL > = 500
// #if CAMAL_VERSION < 500

// #else

class RefFace;
#include "moab/GeomTopoTool.hpp"

#include "CMLSurfEval.hpp"
#include "SmoothBase.hpp"
class SmoothCurveEval;// it is derived from SmoothBase, maybe just need

//! Implement CAMAL geometry callbacks using smooth iMesh 
class SmoothFaceEval : public CMLSurfEval, public SmoothBase
{
public:
  SmoothFaceEval(MBInterface * mb, MBEntityHandle surface_set, MBInterface * mbo,
		  moab::GeomTopoTool * gTool); // entity or entity set

  virtual ~SmoothFaceEval();
  
  virtual double area();

  virtual void bounding_box(double box_min[3], double box_max[3]);

  virtual void move_to_surface(double& x, double& y, double& z);

  virtual void move_to_surface(double& x, double& y, double& z,
			       double& u_guess, double& v_guess);

  virtual bool normal_at(double x, double y, double z, 
                         double& nx, double& ny, double& nz);

  virtual bool normal_at(double x, double y, double z, 
                         double& u_guess, double& v_guess,
                         double& nx, double& ny, double& nz);
  
  virtual bool is_planar();

  virtual bool is_parametric();

  virtual bool is_periodic_in_u(double& u_period);

  virtual bool is_periodic_in_v(double& v_period);

  virtual void get_param_range_u(double& u_low, double& u_high);

  virtual void get_param_range_v(double& v_low, double& v_high);
  
  virtual bool uv_from_position(double x, double y, double z, 
                                double& u, double& v);

  virtual bool uv_from_position(double x, double y, double z, 
                                double& u, double& v,
                                double& cx, double& cy, double& cz);

  virtual void position_from_uv(double u, double v, 
				double& x, double& y, double& z);

  virtual void distortion_at_uv(double u, double v, 
                                double du[3], double dv[3]);

// these are not derived from CMLSurfEval
  void set_mesh_size(double tmp_size);
  double get_mesh_size();
  
  //void set_ref_face(RefFace * refFace);// used to set the ref face from Cholla (or our face)

  //int get_dimension() const { return myDimension; } it should be 2 always

  // this is the setup for all the computation that we need to do
  //
  //int Initialize();

  // initialize normals// they will be stored as tags on nodes
  int init_gradient();

  // some functions for edge evaluations
  MBErrorCode evaluate_smooth_edge(MBEntityHandle eh, double &t,
                              MBCartVect & outv);

  double get_length_loop( int loopIndex);

  MBErrorCode  evaluate_loop_at_u( int loopIndex, double u, MBCartVect & position);

  MBErrorCode  eval_bezier_patch( MBEntityHandle tri,
		  MBCartVect &areacoord,
		  MBCartVect &pt);

  void project_to_facet_plane(
    MBEntityHandle tri,
    MBCartVect &pt,
    MBCartVect &point_on_plane,
    double &dist_to_plane);

  void facet_area_coordinate(
    MBEntityHandle facet,
    MBCartVect & pt_on_plane,
    MBCartVect & areacoord );

  // geom topo tool
  MBErrorCode initialize_geom_topo_tool();

  MBErrorCode project_to_facets_main(MBCartVect &this_point,
                                   bool trim,
                                   bool & outside,
                                   MBCartVect * closest_point_ptr = NULL,  // interested only in normal
                                   MBCartVect * normal_ptr = NULL); // interested only in closest point

  MBErrorCode project_to_facets(
  		std::vector<MBEntityHandle> & facet_list, MBEntityHandle & lastFacet,
  		int interpOrder, double compareTol, MBCartVect &this_point, bool trim,
  		bool & outside, MBCartVect *closest_point_ptr, MBCartVect * normal_ptr);

  MBErrorCode project_to_facet( MBEntityHandle facet,
                                         MBCartVect &pt,
                                         MBCartVect &areacoord,
                                         MBCartVect &close_point,
                                         bool &outside_facet,
                                         double compare_tol);

  bool is_at_vertex(
    MBEntityHandle facet,  // (IN) facet we are evaluating
    MBCartVect &pt,    // (IN) the point
    MBCartVect &ac,    // (IN) the ac of the point on the facet plane
    double compare_tol, // (IN) return TRUE of closer than this
    MBCartVect &eval_pt, // (OUT) location at vertex if TRUE
    MBCartVect *eval_norm_ptr ); // (OUT) normal at vertex if TRUE

  MBErrorCode project_to_patch(
    MBEntityHandle facet,     // (IN) the facet where the patch is defined
    MBCartVect &ac,       // (IN) area coordinate initial guess (from linear facet)
    MBCartVect &pt,       // (IN) point we are projecting to patch
    MBCartVect &eval_pt,  // (OUT) The projected point
    MBCartVect *eval_norm, // (OUT) normal at evaluated point
    bool &outside, // (OUT) the closest point on patch to pt is on an edge
    double compare_tol,    // (IN) comparison tolerance
    int edge_id ) ;         // (IN) only used if this is to be projected to one
                           //      of the edges.  Otherwise, should be -1

  MBErrorCode eval_bezier_patch_normal( MBEntityHandle facet,
  		MBCartVect &areacoord,
  		MBCartVect &normal );

  // this will be called now from driver...
  MBErrorCode compute_tangents_for_each_edge();// they will be used for control points

  MBErrorCode get_normals_for_vertices(const MBEntityHandle * conn2, MBCartVect N[2]);// here we need the gradient tag

  // make this one public, will be called during initialization !!!
  MBErrorCode init_edge_control_points(MBCartVect &P0, MBCartVect &P3,
    		MBCartVect &N0, MBCartVect &N3, MBCartVect &T0, MBCartVect &T3,
    		MBCartVect * Pi);

  // moved from private, because will be called from PaveDriver
  MBErrorCode compute_control_points_on_edges(double min_dot, MBTag edgeCtrlTag, MBTag markTag);

  MBErrorCode compute_internal_control_points_on_facets(
  			double min_dot,  MBTag facetCtrlTag, MBTag facetEdgeCtrlTag);

  // move from private too
  void DumpModelControlPoints();

  //
  MBErrorCode find_loops();

  MBErrorCode mesh_count_total(std::map<MBEntityHandle, SmoothCurveEval*> & mapCurves, int & mesh_count);

  void  evenify(std::map<MBEntityHandle, SmoothCurveEval*> & mapCurves);

  void  mesh(double iMeshSize, std::map<MBEntityHandle, SmoothCurveEval*> & mapCurves);

  int eval_counter() { return _evaluationsCounter;}
private:

  //===========================================================================
  //Function Name: move_ac_inside
  //
  //Member Type:  PRIVATE
  //Description:  find the closest area coordinate to the boundary of the
  //              patch if any of its components are < 0
  //              Return if the ac was modified.
  //===========================================================================
  bool  move_ac_inside( MBCartVect &ac, double tol );

  //===========================================================================
  //Function Name: ac_at_edge
  //
  //Member Type:  PRIVATE
  //Description:  determine the area coordinate of the facet at the edge
  //===========================================================================
  void  ac_at_edge( MBCartVect &fac,  // facet area coordinate
  								MBCartVect &eac,   // edge area coordinate
                                  int edge_id ) ;     // id of edge

  // some local functions that do not need to be public
  MBErrorCode init_bezier_edge(MBEntityHandle edge, double min_dot);
 //

  MBErrorCode find_edges_orientations( MBEntityHandle edges[3],
  		const MBEntityHandle * conn3, int orient[3]); // maybe we will set it?

  MBErrorCode init_facet_control_points(
    MBCartVect N[6],     // vertex normals (per edge)
    MBCartVect  P[3][5],  // edge control points
    MBCartVect  G[6] ) ;   // return internal control points


  //iMesh_Instance _meshIface;
  //iBase_EntitySetHandle _surf_set;

  moab::GeomTopoTool * _my_geomTopoTool;
  MBRange _gsets[4];// redundant here
  MBEntityHandle _obb_root;

  // those are the bounding box limits;
  // they are adjusted for the control points too in each triangle
  void adjust_bounding_box(MBCartVect & vect);
  MBCartVect _minim;
  MBCartVect _maxim;

  MBRange _triangles;
  MBRange _edges;
  MBRange _nodes;

  std::vector<double> _fractions;// they are increasing from 0. to 1., do we need these?
  std::vector<double> _loopLengths;

  // each loop will be actually a vector of MBEntityHandle, paired with a vector of senses
  // number of loops is decided by the size of _loopEnds
  std::vector<MBEntityHandle> _loops; // set1, set 3, set 5, ...
  std::vector<char> _senses; // 0 forward, 1 backward: 0, 0, 1, ...
  std::vector<int> _loopEnds;// the first loop starts at 0 always;

  // number of loops is decided by the size of _loopEnds.size()
  // this ref face will be gone, we will replace it with a new call
  //RefFace * _smooth_face;
  //int myDimension;
  double meshSize;

  // this tag is on edges
  // rval = _mb->tag_create("MARKER", 1, MB_TAG_BIT, _markTag, &value);
  MBTag _markTag; // this is a tag used to mark edges when we look for loops

  // this tag is on nodes
  //MBErrorCode rval = _mb->tag_create("GRADIENT", 3 * sizeof(double),
	// MB_TAG_DENSE, _gradientTag, &defNormal);
  MBTag _gradientTag; // this will be used for normal at nodes

  // this tag is on edges
  //MBErrorCode rval = _mb->tag_create("TANGENTS", 6 * sizeof(double),
  	// MB_TAG_DENSE, _tangentsTag, &defTangents);
  MBTag _tangentsTag; // each edge will have exactly 2 tangents, because there is only
  // one feature edge, and it is periodic
  // the feature edge is exactly on the boundary

  // this tag is on edges
  //MBErrorCode rval = _mb->tag_create("CONTROLEDGE", 9 * sizeof(double),
    	// MB_TAG_DENSE, _edgeCtrlTag, &defControls);
  MBTag _edgeCtrlTag;

  // this tag is on facets (triangles), 6 control points on each facet
  // there are also some 15 points used in evaluation; how to store them?
    //MBErrorCode rval = _mb->tag_create("CONTROLFACE", 18 * sizeof(double),
      	// MB_TAG_DENSE, _facetCtrlTag, &defControls);
  MBTag _facetCtrlTag;

  // these are the 12 points stored for each edge, again
  // it is cheaper this way compared to retrieve the edges every time, determine their orientation, order
  // in triangle, and retrieve the control points from the edge
  // the control points are stored as 12 points, in order : edge 0, 1, and 2, in that order
  //MBErrorCode rval = _mb->tag_create("CONTROLEDGEFACE", 36 * sizeof(double),
       	// MB_TAG_DENSE, _facetEdgeCtrlTag, &defControls);
  MBTag _facetEdgeCtrlTag; //
  // plane of the facet stores as a normal a, b, c and d, distance, for
  // ax+by+cz+d=0
  //MBErrorCode rval = _mb->tag_create("PLANE", 4 * sizeof(double),
        	// MB_TAG_DENSE, _planeTag, &defPlane);
  MBTag _planeTag;

  // counter for calls
  long  _evaluationsCounter;
};
// #endif

#endif /* SMOOTH_FACE_EVAL_HPP*/
