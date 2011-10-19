/*
 * QslimDecimation.cpp
 *
 *  Created on: Mar 10, 2010
 *      Author: iulian
 */

#include <assert.h>
#include "QslimDecimation.hpp"

// for proximity searches
#include "moab/AdaptiveKDTree.hpp"
#include "moab/ReadUtilIface.hpp"

#include "Mat4.h"
#include "defs.h"
#include "quadrics.h"
#include <time.h>
#include <map>

// those are used in model navigation/simplification
#include "primitives.h"

// this is the global thing, that everybody will use
moab::Interface * mb;
moab::Tag uniqIDtag; // this will be used to mark vertices moab::EntityHandles
moab::Tag validTag;

moab::Tag costTag; // simplification induces an error cost at each vertex
// try to keep adding the cost, to see if it is spreading nicely

// this will be used to store plane data for each triangle, as 4 doubles
// will be updated only if needed ( -m option == opts.will_preserve_mesh_quality)
moab::Tag planeDataTag;

moab::Range verts; // original list of vertices, that are part of the original triangles
moab::Range triangles;
moab::Range edgs;
QslimOptions opts; // external
moab::EntityHandle iniSet;

int uniqID(moab::EntityHandle v) {
  int val;
  moab::ErrorCode rval = mb->tag_get_data(uniqIDtag, &v, 1, &val);
  assert(rval==moab::MB_SUCCESS);
  return val;
}
// the vertices are not deleted anymore, just invalidated
// the edges are deleted, though, and triangles
int ehIsValid(moab::EntityHandle v) {
  unsigned char val;
  moab::ErrorCode rval = mb->tag_get_data(validTag, &v, 1, &val);
  assert(rval==moab::MB_SUCCESS);
  return (int) val;
}

// include here the main classes used for decimation

#include "Heap.hpp"
// prox grid is used for proximity grid only
//#include "ProxGrid.h"

class pair_info: public Heapable {
public:
  moab::EntityHandle v0, v1; // Vertex *v0, *v1;

  Vec3 candidate;
  double cost;

  //pair_info ( Vertex *a,Vertex *b ) { v0=a; v1=b; cost=HUGE; }
  pair_info(moab::EntityHandle a, moab::EntityHandle b) {
    v0 = a;
    v1 = b;
    cost = HUGE;
  }

  bool isValid() {
    return ehIsValid(v0) && ehIsValid(v1);
  }//  :v0->isValid() && v1->isValid(); }
};

typedef buffer<pair_info *> pair_buffer;

class vert_info {
public:

  pair_buffer pairs;

  Mat4 Q;
  double norm;

  vert_info() :
    Q(Mat4::zero) {
    pairs.init(2);
    norm = 0.0;
  }
};

// these are
static Heap *heap;
static array<vert_info> vinfo;
static double proximity_limit; // distance threshold squared

int validFaceCount;
int validVertCount;

////////////////////////////////////////////////////////////////////////
//
// Low-level routines for manipulating pairs
//

static inline vert_info& vertex_info(moab::EntityHandle v)//Vertex *v )
{
  //  moab::EntityHandle should return an integer tag with an
  //  index in the big array of vert_info
  //   something like: return tag
  //   for the time being, we can return the simple id...
  return vinfo(uniqID(v));
}

static
bool check_for_pair(moab::EntityHandle v0, moab::EntityHandle v1)//Vertex *v0, Vertex *v1 )
{
  const pair_buffer& pairs = vertex_info(v0).pairs;

  for (int i = 0; i < pairs.length(); i++) {
    if (pairs(i)->v0 == v1 || pairs(i)->v1 == v1)
      return true;
  }

  return false;
}

static pair_info *new_pair(moab::EntityHandle v0, moab::EntityHandle v1)//  Vertex *v0, Vertex *v1 )
{
  vert_info& v0_info = vertex_info(v0);
  vert_info& v1_info = vertex_info(v1);

  pair_info *pair = new pair_info(v0, v1);
  v0_info.pairs.add(pair);
  v1_info.pairs.add(pair);

  return pair;
}

static
void delete_pair(pair_info *pair) {
  vert_info& v0_info = vertex_info(pair->v0);
  vert_info& v1_info = vertex_info(pair->v1);

  v0_info.pairs.remove(v0_info.pairs.find(pair));
  v1_info.pairs.remove(v1_info.pairs.find(pair));

  if (pair->isInHeap())
    heap->kill(pair->getHeapPos());

  delete pair;
}

////////////////////////////////////////////////////////////////////////
//
// The actual guts of the algorithm:
//
//     - pair_is_valid
//     - compute_pair_info
//     - do_contract
//
/*
static
bool pair_is_valid(moab::EntityHandle u, moab::EntityHandle v)// Vertex *u, Vertex *v )
{
  //
  Vec3 vu = getVec3FromMBVertex(mb, u);
  Vec3 vv = getVec3FromMBVertex(mb, v);
  return norm2(vu - vv) < proximity_limit;
  //return  norm2 ( *u - *v ) < proximity_limit;
}*/

static
int predict_face(moab::EntityHandle tria, moab::EntityHandle v1,
    moab::EntityHandle v2, /*Face& F, Vertex *v1, Vertex *v2,*/
    Vec3& vnew, Vec3& f1, Vec3& f2, Vec3& f3) {
  int nmapped = 0;
  const moab::EntityHandle * conn;
  int num_nodes;
  moab::ErrorCode rval = mb->get_connectivity(tria, conn, num_nodes);
  assert(3==num_nodes && rval == moab::MB_SUCCESS);
  if (conn[0] == v1 || conn[0] == v2) {
    f1 = vnew;
    nmapped++;
  } else
    f1 = getVec3FromMBVertex(mb, conn[0]);

  if (conn[1] == v1 || conn[1] == v2) {
    f2 = vnew;
    nmapped++;
  } else
    f2 = getVec3FromMBVertex(mb, conn[1]);

  if (conn[2] == v1 || conn[2] == v2) {
    f3 = vnew;
    nmapped++;
  } else
    f3 = getVec3FromMBVertex(mb, conn[2]);

  // find vertices in face tria
  /*
   if ( F.vertex ( 0 ) == v1 || F.vertex ( 0 ) == v2 )
   { f1 = vnew;  nmapped++; }
   else f1 = *F.vertex ( 0 );

   if ( F.vertex ( 1 ) == v1 || F.vertex ( 1 ) == v2 )
   { f2 = vnew;  nmapped++; }
   else f2 = *F.vertex ( 1 );

   if ( F.vertex ( 2 ) == v1 || F.vertex ( 2 ) == v2 )
   { f3 = vnew;  nmapped++; }
   else f3 = *F.vertex ( 2 );
   */
  return nmapped;
}

#define MESH_INVERSION_PENALTY 1e9

static
double pair_mesh_positivity(/* Model& M,*/moab::EntityHandle v1,
    moab::EntityHandle v2, /*Vertex *v1, Vertex *v2,*/
    Vec3& vnew) {
  std::vector<moab::EntityHandle> changed;

  // :  construct the list of faces influenced by the
  //   moving of vertices v1 and v2 into vnew
  //M.contractionRegion ( v1, v2, changed );
  moab::ErrorCode rval = contractionRegion(mb, v1, v2, changed);
  if (rval != moab::MB_SUCCESS) {
    std::cout << "error in getting adjacency information vs: "
        << mb->id_from_handle(v1) << " " << mb->id_from_handle(v2) << "\n";
  }

  // double Nsum = 0;
  if (opts.logfile && opts.selected_output & OUTPUT_CONTRACTIONS)
    *opts.logfile << " positivity for v1, v2: " << mb->id_from_handle(v1)
        << " " << mb->id_from_handle(v2) << std::endl;

  for (unsigned int i = 0; i < changed.size(); i++) {
    //Face& F = *changed ( i );
    moab::EntityHandle F = changed[i];
    Vec3 f1, f2, f3;

    int nmapped = predict_face(F, v1, v2, vnew, f1, f2, f3);

    //
    // Only consider non-degenerate faces
    if (nmapped < 2) {
      //Plane Pnew ( f1, f2, f3 );
#if 0
      Vec3 normalNew = Pnew.normal();
      if ( normalNew[Z] < positivityMin )
      positivityMin=normalNew[Z]; // Z direction!!!
      if (opts.logfile && opts.selected_output&OUTPUT_CONTRACTIONS )
      *opts.logfile << "Triangle " << mb->id_from_handle(F)
      << " nmapped " << nmapped << std::endl;
      if (opts.logfile && positivityMin<=0 && opts.selected_output&OUTPUT_CONTRACTIONS )
      *opts.logfile << "Triangle " << mb->id_from_handle(F)
      << " normal Z:" << normalNew[Z] << std::endl;
#endif
      double positiv = (f2[0] - f1[0]) * (f3[1] - f1[1]) - (f2[1] - f1[1])
          * (f3[0] - f1[0]);
      if (positiv <= 0) {
        if (opts.logfile && opts.selected_output & OUTPUT_CONTRACTIONS)
          *opts.logfile << "Triangle " << mb->id_from_handle(F) << " nmapped "
              << nmapped << " orient: " << positiv << std::endl;
        return MESH_INVERSION_PENALTY * 10;
      }
    }
  }

  //return (-Nmin) * MESH_INVERSION_PENALTY;

  return 0.0;
}

static
double pair_mesh_penalty( /*Model& M, Vertex *v1, Vertex *v2,*/
moab::EntityHandle v1, moab::EntityHandle v2, Vec3& vnew) {
  std::vector<moab::EntityHandle> changed;

  //   construct the list of faces influenced by the
  //   moving of vertices v1 and v2 into vnew
  //M.contractionRegion ( v1, v2, changed );
  moab::ErrorCode rval = contractionRegion(mb, v1, v2, changed);
  if (rval != moab::MB_SUCCESS) {
    std::cout << "error in getting adjacency information vs: "
        << mb->id_from_handle(v1) << " " << mb->id_from_handle(v2) << "\n";
  }

  // double Nsum = 0;
  double Nmin = 0;

  for (unsigned int i = 0; i < changed.size(); i++) {
    //Face& F = *changed ( i );
    moab::EntityHandle F = changed[i];
    Vec3 f1, f2, f3;

    int nmapped = predict_face(F, v1, v2, vnew, f1, f2, f3);
    //
    // Only consider non-degenerate faces
    if (nmapped < 2) {
      Plane Pnew(f1, f2, f3);
      Plane p = trianglePlane(mb, F);
      double delta = Pnew.normal() * p.normal(); //  Pnew.normal() * F.plane().normal();

      if (Nmin > delta)
        Nmin = delta;
    }
  }

  //return (-Nmin) * MESH_INVERSION_PENALTY;
  if (Nmin < 0.0)
    return MESH_INVERSION_PENALTY;
  else
    return 0.0;
}

static
void compute_pair_info(pair_info *pair /* Model * pM0,*/) {
  moab::EntityHandle v0 = pair->v0;
  moab::EntityHandle v1 = pair->v1;

  // Vertex *v0 = pair->v0;
  // Vertex *v1 = pair->v1;

  vert_info& v0_info = vertex_info(v0);
  vert_info& v1_info = vertex_info(v1);

  Mat4 Q = v0_info.Q + v1_info.Q;
  double norm = v0_info.norm + v1_info.norm;

  pair->cost = quadrix_pair_target(Q, v0, v1, pair->candidate);

  if (opts.will_weight_by_area)
    pair->cost /= norm;

  if (opts.will_preserve_mesh_quality)
    pair->cost += pair_mesh_penalty(/* *pM0,*/v0, v1, pair->candidate);

  if (opts.height_fields)
    pair->cost += pair_mesh_positivity(/* *pM0, */v0, v1, pair->candidate);

  if (opts.logfile && opts.selected_output & OUTPUT_CONTRACTIONS)

    *opts.logfile << "pair ( v" << mb->id_from_handle(v0) << " v"
        << mb->id_from_handle(v1) << " ) cost: " << -pair->cost << std::endl;
  //
  // NOTICE:  In the heap we use the negative cost.  That's because
  //          the heap is implemented as a MAX heap.
  //
  if (pair->isInHeap()) {
    heap->update(pair, -pair->cost);
  } else {
    heap->insert(pair, -pair->cost);
  }
  if (opts.logfile && opts.selected_output & OUTPUT_COST) {
    heap_node *top = heap->top();
    pair_info *pairTop = (pair_info *) top->obj;
    *opts.logfile << " i/u pair (" << uniqID(pair->v0) << "," << uniqID(
        pair->v1) << ")=" << pair->cost << "  min : (" << uniqID(pairTop->v0)
        << "," << uniqID(pairTop->v1) << ") " << pairTop->cost << std::endl;
  }
}
void recomputeChangedPairsCost(std::vector<moab::EntityHandle> & changed,/* Model *pM0, Vertex *v0,*/
moab::EntityHandle v0) {
  //
  for (unsigned int i = 0; i < changed.size(); i++) {

    moab::EntityHandle F = changed[i];
    const moab::EntityHandle * conn;
    int num_nodes;
    mb->get_connectivity(F, conn, num_nodes);
    //if (!F->isValid())
    //   continue;
    // recompute the pair that is not connected to vertex v0
    // loop over all the vertices of F that are not v0, and recompute the costs
    // of all the pairs associated  that do not contain v0
    // we do not have to recreate or delete any pair, we just recompute what we have
    // some will be recomputed 2 times, but it is OK
    for (int k = 0; k < 3; k++) {
      moab::EntityHandle v = conn[k];
      if (v == v0)
        continue;
      vert_info & v_info = vertex_info(v);
      for (int j = 0; j < v_info.pairs.length(); j++) {
        pair_info *p = v_info.pairs(j);
        if (p->v0 == v0 || p->v1 == v0)
          continue; // do not recompute cost of pairs already computed
        if (opts.logfile && (opts.selected_output & OUTPUT_COST))
          *opts.logfile << "recompute cost of pair (v" << uniqID(p->v0) + 1
              << " v" << uniqID(p->v1) + 1 << ")" << std::endl;
        compute_pair_info(p);
      }

    }

  }
}

static
void do_contract(pair_info *pair) {

  moab::EntityHandle v0 = pair->v0;
  moab::EntityHandle v1 = pair->v1;
  // cost of contraction is accumulated at v0
  double costToContract = pair->cost;
  vert_info& v0_info = vertex_info(v0);
  vert_info& v1_info = vertex_info(v1);
  Vec3 vnew = pair->candidate;
  if (opts.logfile && (opts.selected_output & OUTPUT_COST)) {
    *opts.logfile << "---- do contract v0:" << uniqID(v0) + 1 << " v1:"
        << uniqID(v1) + 1 << std::endl;
  }
  int i;

  //
  // Make v0 be the new vertex
  v0_info.Q += v1_info.Q;
  v0_info.norm += v1_info.norm;

  //
  // Perform the actual contraction
  std::vector<moab::EntityHandle> changed;
  moab::ErrorCode rval1 = contract(mb, v0, v1, vnew, changed);
  // this is a list of triangles still connected to v0 (are they valid? probably)
  if (opts.will_preserve_mesh_quality) {
    // recompute normals only in this case, because they are not needed otherwise
    int size = changed.size();
    for (int i = 0; i < size; i++) {
      computeTrianglePlane(mb, changed[i]);
    }
  }
  assert (moab::MB_SUCCESS == rval1);

#ifdef SUPPORT_VCOLOR
  //
  // If the vertices are colored, color the new vertex
  // using the average of the old colors.
  v0->props->color += v1->props->color;
  v0->props->color /= 2;
#endif

  //
  // Remove the pair that we just contracted
  delete_pair(pair);
  //
  // Recalculate pairs associated with v0
  for (i = 0; i < v0_info.pairs.length(); i++) {
    pair_info *p = v0_info.pairs(i);
    compute_pair_info(p);
  }

  //
  // Process pairs associated with now dead vertex

  static pair_buffer condemned(6); // collect condemned pairs for execution
  condemned.reset();

  for (i = 0; i < v1_info.pairs.length(); i++) {
    pair_info *p = v1_info.pairs(i);

    moab::EntityHandle u;
    if (p->v0 == v1)
      u = p->v1;
    else if (p->v1 == v1)
      u = p->v0;
    else
      std::cerr << "YOW!  This is a bogus pair." << std::endl;

    if (!check_for_pair(u, v0)) {
      p->v0 = v0;
      p->v1 = u;
      v0_info.pairs.add(p);
      compute_pair_info(p);
    } else
      condemned.add(p);
  }

  for (i = 0; i < condemned.length(); i++)
    // Do you have any last requests?
    delete_pair(condemned(i));
  // only now you can delete the vertex v1 from database
  // moab::ErrorCode rval = mb->delete_entities(&v1, 1);
  // no, it is better to invalidate the vertex, do not delete it
  // maybe we will delete at the end all that are invalid ??
  int invalid = 0;
  moab::ErrorCode rval = mb->tag_set_data(validTag, &v1, 1, &invalid);

  if (opts.plotCost) {
    double cost_at_v0 = 0; // maybe it is already set before
    rval = mb->tag_get_data(costTag, &v0, 1, &cost_at_v0);
    cost_at_v0 += costToContract;
    rval = mb->tag_set_data(costTag, &v0, 1, &cost_at_v0);
  }

  v1_info.pairs.reset(); // safety precaution
  recomputeChangedPairsCost(changed, v0);

}

////////////////////////////////////////////////////////////////////////
//
// External interface: setup and single step iteration
//

bool decimate_quadric(moab::EntityHandle v, Mat4& Q) {
  if (vinfo.length() > 0) {
    vert_info & vinf = vertex_info(v);
    Q = vinf.Q;
    return true;
  } else
    return false;
}

// it is assumed it is mb, moab interface
void decimate_contract() {
  heap_node *top;
  pair_info *pair;

  for (;;) {
    top = heap->extract();
    if (!top)
      return;
    pair = (pair_info *) top->obj;

    //
    // This may or may not be necessary.  I'm just not quite
    // willing to assume that all the junk has been removed from the
    // heap.
    if (pair->isValid())
      break;

    delete_pair(pair);
  }

  if (opts.logfile && (opts.selected_output & OUTPUT_COST))
    *opts.logfile << "#$cost " << validFaceCount << " before contract: "
        << pair->cost << std::endl;

  do_contract(pair);

  if (opts.logfile && (opts.selected_output & OUTPUT_COST))
    *opts.logfile << "#$cost " << validFaceCount << std::endl;

  validVertCount--; // Attempt to maintain valid vertex information
}

double decimate_error(moab::EntityHandle v) {
  vert_info& info = vertex_info(v);

  double err = quadrix_evaluate_vertex(v, info.Q);

  if (opts.will_weight_by_area)
    err /= info.norm;

  return err;
}

double decimate_min_error() {
  heap_node *top;
  pair_info *pair;

  for (;;) {
    top = heap->top();
    if (!top)
      return -1.0;
    pair = (pair_info *) top->obj;

    if (pair->isValid())
      break;

    top = heap->extract();
    delete_pair(pair);
  }

  return pair->cost;
}
#if 0
double decimate_max_error ( )
{
  real max_err = 0;

  for ( int i=0; i<m.vertCount(); i++ )
  if ( m.vertex ( i )->isValid() )
  {
    max_err = MAX ( max_err, decimate_error ( m.vertex ( i ) ) );
  }

  return max_err;
}
#endif
namespace MeshKit {

QslimDecimation::QslimDecimation(moab::Interface * mbi,
    moab::EntityHandle root_set) {
  _mb = mbi;
  iniSet = root_set;// it is not necessarily the root set; this is external; bad design
}

QslimDecimation::~QslimDecimation() {
}
int QslimDecimation::decimate(QslimOptions & iOpts, moab::Range & oRange) {
  // opts is external
  opts = iOpts;

  mb = _mb; // (reinterpret_cast<MBiMesh *> (m_mesh))->mbImpl;
  // need to get all the triangles from the set
  // also all the edges, and all vertices
  // not  a good design here; mb is extern !
  //
  if (NULL == mb)
    return 1;// error
  //moab::EntityHandle mbSet = reinterpret_cast<moab::EntityHandle>(m_InitialSet);

  clock_t start_time = clock();
  int err = Init();
  if (err) {
    std::cerr << "Error in initialization of decimation. Abort\n";
    return 1;
  }
  clock_t init_time = clock();
  std::cout << "   Initialization  " << (double) (init_time - start_time)
      / CLOCKS_PER_SEC << " s.\n";

  int faces_reduction = validFaceCount - opts.face_target;
  int counter = 0, interval = 0;
  clock_t currentTime = init_time;
  while (validFaceCount > opts.face_target && decimate_min_error()
      < opts.error_tolerance) {
    int initf = validFaceCount;
    // big routine
    decimate_contract();
    counter += (initf - validFaceCount);
    if (counter > faces_reduction / opts.timingIntervals) {
      // print some time stats
      clock_t p10_time = clock();
      std::cout << "     " << ++interval << "/" << opts.timingIntervals
          << " reduce to " << validFaceCount << " faces in "
          << (double) (p10_time - currentTime) / CLOCKS_PER_SEC << " s, total:"
          << (double) (p10_time - init_time) / CLOCKS_PER_SEC << " s.\n";
      counter = 0;
      currentTime = p10_time;
    }
  }

  clock_t finish_time = clock();
  std::cout << "   Decimation: " << (double) (finish_time - init_time)
      / CLOCKS_PER_SEC << " s.\n";

  if (opts.create_range) {
    // the remaining nodes and triangles are copied in the range
    // they are put in the old set, too
    // maybe this has to change
    // count first valid vertices and triangles
    moab::Range::iterator it;
    std::vector<moab::EntityHandle> validVerts;
    std::vector<moab::EntityHandle> validTris;
    for (it = triangles.begin(); it != triangles.end(); it++) {
      if (ehIsValid(*it))
        validTris.push_back(*it);
    }
    for (it = verts.begin(); it != verts.end(); it++) {
      if (ehIsValid(*it))
        validVerts.push_back(*it);
    }
    //
    std::vector<double> coords;
    int numNodes = (int) validVerts.size();
    int numTriangles = (int) validTris.size();
    coords.resize(3 * numNodes);

    moab::ErrorCode rval = mb->get_coords(&(validVerts[0]), numNodes,
        &(coords[0]));
    assert(moab::MB_SUCCESS==rval);

    // create the new vertices, at the same  coordinates
    // put those verts in the range that is output

    // to make sure, the range is cleared
    oRange.clear();
    rval = mb->create_vertices(&coords[0], numNodes, oRange);
    assert(moab::MB_SUCCESS==rval);
    // the first element in the range must be the start of the new vertex sequence
    std::map<moab::EntityHandle, moab::EntityHandle> mapping; // this will be from old
    //  to new vertices, for connectivity
    for (int i = 0; i < numNodes; i++) {
      mapping[validVerts[i]] = oRange[i];
    }

    //get the query interface, which we will use to create the triangles directly
    moab::ReadUtilIface *iface;
    rval = mb -> query_interface(iface);// use the new query interface
    assert(moab::MB_SUCCESS==rval);

    //create the triangles, get a direct ptr to connectivity
    moab::EntityHandle starth, *connect;
    rval = iface -> get_element_connect(numTriangles, 3, moab::MBTRI, 1,
        starth, connect);
    assert(moab::MB_SUCCESS==rval);
    // get first the connectivity of the old triangles
    const moab::EntityHandle * conn3;
    for (int i = 0; i < numTriangles; i++) {
      int num_nodes;
      // get each valid triangle one by one, and set the new connectivity directly
      rval = mb->get_connectivity(validTris[i], conn3, num_nodes);
      assert( (moab::MB_SUCCESS==rval) && (num_nodes==3));
      // directly modify the connect array in database
      for (int j = 0; j < 3; j++)
        connect[j] = mapping[conn3[j]];

      // update adjacencies
      // not very smart...; we would like to update once and for all triangles
      // not in a loop
      rval = iface -> update_adjacencies(starth+i, 1, 3, connect);
      assert(moab::MB_SUCCESS==rval);

      connect += 3; // advance
    }



    // clear completely the initial set, after deleting all elements from it...
    //   ok, we are done, commit to ME ?
    rval = mb->delete_entities(triangles);
    assert(moab::MB_SUCCESS==rval);
    rval = mb->delete_entities(edgs);
    assert(moab::MB_SUCCESS==rval);
    rval = mb->delete_entities(verts);
    assert(moab::MB_SUCCESS==rval);
    // remove everything from the initial set, because we will add the new triangles
    mb->remove_entities(iniSet, triangles);
    mb->remove_entities(iniSet, verts);
    mb->remove_entities(iniSet, edgs);

    //add triangles to output range (for the MESelection)
    oRange.insert(starth, starth + numTriangles - 1);
    // add all entities from the range to the initial set, now
    rval = mb->add_entities(iniSet, oRange);
    assert(moab::MB_SUCCESS==rval);
    //me->commit_mesh(mit->second, COMPLETE_MESH);
    // end copy
  } else {
    moab::Range::const_reverse_iterator rit;
    if (opts.useDelayedDeletion) {

      // put in a range triangles to delete
      moab::Range delRange;
      // delete triangles and edges that are invalid
      for (rit = triangles.rbegin(); rit != triangles.rend(); rit++) {
        moab::EntityHandle tr = *rit;
        // check the validity
        if (ehIsValid(tr))
          continue;
        mb->delete_entities(&tr, 1);
        delRange.insert(tr);
      }
      mb->remove_entities(iniSet, delRange);
      // maybe we should delete all edges, but for now, we will keep them
      for (rit = edgs.rbegin(); rit != edgs.rend(); rit++) {
        moab::EntityHandle v = *rit;
        // check the validity
        if (ehIsValid(v))
          continue;
        mb->delete_entities(&v, 1);
      }

    }
    // delete them one by one
    for (rit = verts.rbegin(); rit != verts.rend(); rit++) {
      moab::EntityHandle v = *rit;
      // check the validity
      if (ehIsValid(v))
        continue;
      mb->delete_entities(&v, 1);
    }
  }
  clock_t delete_vTime = clock();
  std::cout << "   Delete Vertices: " << (double) (delete_vTime - finish_time)
      / CLOCKS_PER_SEC << " s.\n";
  // we need to delete the tags we created; they are artificial
  // list of tags to delete:
  // moab::Tag uniqIDtag; // this will be used to mark vertices moab::EntityHandles
  // moab::Tag validTag;
  // moab::Tag planeDataTag;

  // moab::Tag costTag; // simplification induces an error cost at each vertex
  // try to keep adding the cost, to see if it is spreading nicely

  // keep only the cost, the rest are artificial
  mb->tag_delete(uniqIDtag);
  mb->tag_delete(validTag);
  mb->tag_delete(planeDataTag);
  //

  return 0;
}

int QslimDecimation::Init() {
  int i;
  unsigned int j;

  //moab::EntityHandle * set = reinterpret_cast<moab::EntityHandle *> (&_InitialSet);
  moab::ErrorCode rval = mb->get_entities_by_type(iniSet, moab::MBTRI,
      triangles);
  validFaceCount = triangles.size();// this gets reduced every time we simplify the model

  // store the normals/planes computed at each triangle
  // we may need just the normals, but compute planes, it is about the same job
  double defPlane[] = { 0., 0., 1., 0. };
  rval = mb->tag_get_handle("PlaneTriangleData", 4, moab::MB_TYPE_DOUBLE,
      planeDataTag, moab::MB_TAG_DENSE | moab::MB_TAG_CREAT, &defPlane);

  // compute the triangle plane and store it, for each triangle
  for (moab::Range::iterator itr = triangles.begin(); itr != triangles.end(); itr++) {
    // this index i will be the index in the vinfo array
    moab::EntityHandle tri = *itr;
    computeTrianglePlane(mb, tri);
    // setting the data for the tag/triangle is done in the compute
    //rval = mb->tag_set_data(planeDataTag, &tri, 1, &plane);
  }

  // create all the edges if not existing
  mb->get_adjacencies(triangles, 1, true, edgs, moab::Interface::UNION);

  // moab::Range verts;// the vertices are always there, they do not need to be created
  mb->get_adjacencies(triangles, 0, true, verts, moab::Interface::UNION);
  int numNodes = verts.size();
  validVertCount = numNodes; // this will be kept
  vinfo.init(numNodes);
  // set a unique integer tag with the position in vinfo array
  //  this will be used instead of v->uniqID in the vinfo array
  int def_data = -1;

  rval = mb->tag_get_handle("uniqID", 1, moab::MB_TYPE_INTEGER, uniqIDtag,
      moab::MB_TAG_DENSE | moab::MB_TAG_EXCL, &def_data);
  if (moab::MB_SUCCESS != rval)
    return 1;

  unsigned char def_data_bit = 1;// valid by default

  rval = mb->tag_get_handle("valid", 1, moab::MB_TYPE_BIT, validTag,
      moab::MB_TAG_EXCL | moab::MB_TAG_BIT, &def_data_bit);
  if (moab::MB_SUCCESS != rval)
    return 1;

  if (opts.plotCost) {
    double cost_default = 0.;

    rval = mb->tag_get_handle("costTAG", 1, moab::MB_TYPE_DOUBLE, costTag,
        moab::MB_TAG_DENSE | moab::MB_TAG_CREAT, &cost_default);
    if (moab::MB_SUCCESS != rval)
      return 1;
  }

  // set tag for each vertex; this will not be changed during simplification
  i = 0; // for index
  for (moab::Range::iterator it = verts.begin(); it != verts.end(); it++, i++) {
    // this index i will be the index in the vinfo array
    moab::EntityHandle v = *it;
    rval = mb->tag_set_data(uniqIDtag, &v, 1, &i);
    // the default valid data is 1; set it to 0 only to "mark" the vertex invalid for future deletion
    // is it really necessary if we put the default value as 1 anyway

    //rval = mb->tag_set_data(validTag, &v, 1, &def_data_bit);
  }
  moab::Range::iterator it;
  if (opts.logfile && opts.selected_output & OUTPUT_MODEL_DEFN) {
    for (it = verts.begin(); it != verts.end(); it++) {
      moab::EntityHandle v = *it;
      double coords[3];
      rval = mb->get_coords(&v, 1, coords);
      *opts.logfile << "v: " << uniqID(v) << " " << mb->id_from_handle(v)
          << " " << coords[0] << " " << coords[1] << " " << coords[2]
          << std::endl;
    }
  }
  std::cout << "  Decimate:  Distributing shape constraints." << std::endl;

  if (opts.will_use_vertex_constraint)
    for (it = verts.begin(); it != verts.end(); it++) {
      moab::EntityHandle v = *it;
      vertex_info(v).Q = quadrix_vertex_constraint(v);
    }

  if (opts.will_use_plane_constraint) {
    for (it = triangles.begin(); it != triangles.end(); it++) {

      moab::EntityHandle tr = *it;
      Mat4 Q = quadrix_plane_constraint(tr);
      double norm = 0.0;

      if (opts.will_weight_by_area) {
        norm = 1; // triangle area : m_model->face ( i )->area();
        Q *= norm;
      }
      const moab::EntityHandle * conn;
      int num_nodes;
      rval = mb->get_connectivity(tr, conn, num_nodes);
      for (j = 0; j < 3; j++) {
        vert_info& vj_info = vertex_info(conn[j]);
        vj_info.Q += Q;
        vertex_info(conn[j]).norm += norm;

      }
    }
  }
  // just define (one uniqTag for a triangle, see what is happening)
  moab::EntityHandle tr1 = triangles[0];
  rval = mb->tag_set_data(uniqIDtag, &tr1, 1, &i);// just some value

  if (opts.will_constrain_boundaries) {
    std::cout << "  Decimate:  Accumulating discontinuity constraints."
        << std::endl;
    for (it = edgs.begin(); it != edgs.end(); it++) {
      moab::EntityHandle edg = *it;
      if (is_border(edg)) {
        const moab::EntityHandle * conn;
        int num_nodes;
        rval = mb->get_connectivity(edg, conn, num_nodes);
        if (moab::MB_SUCCESS != rval)
          return 1;// fail
        Mat4 B = quadrix_discontinuity_constraint(edg);
        double norm = 0.0;

        if (opts.will_weight_by_area) {
          Vec3 ve1 = getVec3FromMBVertex(mb, conn[0]);
          Vec3 ve2 = getVec3FromMBVertex(mb, conn[1]);
          norm = norm2(ve1 - ve2);
          B *= norm;
        }

        B *= opts.boundary_constraint_weight;
        vert_info& v0_info = vertex_info(conn[0]);
        vert_info& v1_info = vertex_info(conn[1]);
        v0_info.Q += B;
        v0_info.norm += norm;
        v1_info.Q += B;
        v1_info.norm += norm;
      }
    }
  }

  std::cout << "  Decimate:  Allocating heap." << std::endl;
  heap = new Heap(edgs.size());

  int pair_count = 0;

  std::cout << "  Decimate:  Collecting pairs [edges]." << std::endl;
  for (it = edgs.begin(); it != edgs.end(); it++) {
    moab::EntityHandle edg = *it;
    const moab::EntityHandle * conn;
    int num_nodes;
    rval = mb->get_connectivity(edg, conn, num_nodes);
    if (moab::MB_SUCCESS != rval)
      return 1;// fail
    pair_info *pair = new_pair(conn[0], conn[1]);
    compute_pair_info(pair);
    pair_count++;
  }

  if (opts.pair_selection_tolerance < 0) {
    opts.pair_selection_tolerance = 1;//m_model->bounds.radius * 0.05;
    std::cout << "  Decimate:  Auto-limiting at 5% of model radius."
        << std::endl;
  }
  proximity_limit = opts.pair_selection_tolerance
      * opts.pair_selection_tolerance;
  if (proximity_limit > 0) {
    std::cout << "  Decimate:  Collecting pairs [limit="
        << opts.pair_selection_tolerance << "]." << std::endl;
    // use adaptive kd tree to find proximity vertices
    moab::EntityHandle tree_root = 0;
    moab::AdaptiveKDTree kd(mb, true);
    rval = kd.build_tree(verts, tree_root);
    if (rval != moab::MB_SUCCESS) {
      std::cout << "Can't build tree for vertices" << std::endl;
      return 1;
    }

    for (it = verts.begin(); it != verts.end(); it++) {
      moab::Range closeVertices;
      closeVertices.clear();
      moab::EntityHandle v = *it;
      double coords[3];
      mb->get_coords(&v, 1, coords);
      //moab::CartVect v1(coords);
      std::vector<moab::EntityHandle> leaves; // sets of vertices close by
      kd.leaves_within_distance(tree_root, coords,
          opts.pair_selection_tolerance, leaves);
      // add to the list of close vertices
      for (j = 0; j < leaves.size(); j++) {
        rval = mb->get_entities_by_type(leaves[j], moab::MBVERTEX,
            closeVertices);// add to the list
      }

      for (moab::Range::iterator it2 = closeVertices.begin(); it2
          != closeVertices.end(); it2++) {
        moab::EntityHandle vclose = *it2;
        if (v == vclose)
          continue;
        double coords2[3];
        mb->get_coords(&vclose, 1, coords2);

        //moab::CartVect v2(coords2);
        double dd = (coords[0] - coords2[0]) * (coords[0] - coords2[0])
            + (coords[1] - coords2[1]) * (coords[1] - coords2[1]) + (coords[2]
            - coords2[2]) * (coords[2] - coords2[2]);

        if (dd > proximity_limit)
          continue;

/*
#ifdef SAFETY
        assert ( pair_is_valid ( v1,v2 ) );
#endif
*/
        if (!check_for_pair(v, vclose)) {
          pair_info *pair = new_pair(v, vclose);
          compute_pair_info(pair);
          pair_count++;
        }
      }

    }
  } else
    std::cout << "  Decimate:  Ignoring non-edge pairs [limit=0]." << std::endl;

  std::cout << "  Decimate:  Designated " << pair_count << " pairs."
      << std::endl;

  return 0;// no error
}

} // namespace MeshKit
