#include "meshkit/EdgeMesher.hpp"
#include "meshkit/MKCore.hpp"
#include "meshkit/MeshOpFactory.hpp"
#include "meshkit/SizingFunction.hpp"
#include "moab/ReadUtilIface.hpp"
#include <vector>
#include <math.h>

namespace MeshKit 
{
    
// static registration of this  mesh scheme
moab::EntityType tps[] = {moab::MBVERTEX, moab::MBEDGE};
    
int success = MeshOpFactory::instance()->register_meshop("EdgeMesher", tps, 2, EdgeMesher::factory, NULL);
    
MeshOp *EdgeMesher::factory(MKCore *mkcore, const MEVector &me_vec) 
{
  return new EdgeMesher(mkcore, me_vec);
}

void EdgeMesher::mesh_types(std::vector<moab::EntityType> &tps) 
{
  tps.push_back(moab::MBVERTEX);
  tps.push_back(moab::MBEDGE);
}
    
double EdgeMesher::measure(iGeom::EntityHandle ent, double ustart, double uend) const
{
  double umin, umax;
  iGeom::Error err = mk_core()->igeom_instance()->getEntURange(ent, umin, umax);
  IBERRCHK(err, "Trouble getting parameter range for edge.");
  if (umin == umax) throw Error(MK_BAD_GEOMETRIC_EVALUATION, "Edge evaluated to same parameter umax and umin.");
  
  double measure;
  err = mk_core()->igeom_instance()->measure(&ent, 1, &measure);
  IBERRCHK(err, "Trouble getting edge measure.");
  
  return measure*(uend-ustart)/(umax-umin);
}

void EdgeMesher::setup_this()
{
    // compute the number of intervals for the associated ModelEnts, from the size set on them,
    // the sizing function they point to, or a default sizing function
  for (MESelection::iterator mit = meSelection.begin(); mit != meSelection.end(); mit++) {
    ModelEnt *me = mit->first;

      // first check to see whether entity is meshed
    if (me->get_meshed_state() >= COMPLETE_MESH ||
        me->mesh_intervals() > 0) continue;
    
    SizingFunction *sf = mk_core()->sizing_function(me->sizing_function_index());
    if (!sf && me->mesh_intervals() < 0 && me->interval_firmness() == DEFAULT) {
        // no sizing set, just assume default #intervals as 4
      me->mesh_intervals(4);
      me->interval_firmness(DEFAULT);
    }
    else {
        // check # intervals first, then size, and just those for now
      if (sf->intervals() > 0) {
        me->mesh_intervals(sf->intervals());
        me->interval_firmness(HARD);
      }
      else if (sf->size() > 0) {
        me->mesh_intervals(me->measure() / sf->size());
        me->interval_firmness(SOFT);
      }
      else 
        throw Error(MK_INCOMPLETE_MESH_SPECIFICATION, 
                    "Sizing function for edge had neither positive size nor positive intervals.");
    }
  }
}

void EdgeMesher::execute_this()
{
  std::vector<double> coords;
  std::vector<moab::EntityHandle> nodes, edges;

  for (MESelection::iterator mit = meSelection.begin(); mit != meSelection.end(); mit++) {
    ModelEnt *me = mit->first;

      // resize the coords based on the interval setting
    int num_edges = me->mesh_intervals();
    coords.resize(3*(num_edges+1));
    nodes.resize(num_edges+1);
    edges.resize(num_edges);

      // get bounding mesh entities, use 1st 2 entries of nodes list temporarily
    me->boundary(0, nodes);

      // get coords in list, then move one tuple to last position
    moab::ErrorCode rval = mk_core()->moab_instance()->get_coords(&nodes[0], 2, &coords[0]);
    MBERRCHK(rval, "Trouble getting bounding vertex positions.");
    for (int i = 0; i < 3; i++) coords[3*num_edges+i] = coords[3+i];
    nodes[num_edges] = nodes[1];

    equal_meshing(me, num_edges, coords);
      	
    rval = mk_core()->moab_instance()->create_vertices(&coords[3], num_edges-1, mit->second);
    MBERRCHK(rval, "Couldn't create nodes");
    
      // distribute nodes into vector
    int j = 1;
    for (moab::Range::iterator rit = mit->second.begin(); rit != mit->second.end(); rit++)
      nodes[j++] = *rit;

      // get the query iface, which we'll use to create edges directly
    moab::ReadUtilIface *iface;
    rval = mk_core()->moab_instance()->query_interface("ReadUtilIface", (void**)&iface);
    MBERRCHK(rval, "Couldn't get ReadUtilIface.");
    
      // create the edges, getting a direct ptr to connectivity
    moab::EntityHandle starth, *connect, *tmp_connect;
    rval = iface->get_element_connect(num_edges, 2, moab::MBEDGE, 1, starth, connect);
    MBERRCHK(rval, "Couldn't create edge elements.");
    
      // add edges to range for the MESelection
    mit->second.insert(starth, starth+num_edges-1);
    
      // now set the connectivity array from the nodes
    tmp_connect = &nodes[0];
    for (int i = 0; i < num_edges; i++) {
      connect[0] = tmp_connect[0];
      connect[1] = tmp_connect[1];
        // increment connectivity ptr by 2 to go to next edge
      connect += 2;
        // increment tmp_connect by 1, to go to next node
      tmp_connect++;
    }
    
      // ok, we're done; commit to ME
    me->commit_mesh(mit->second, COMPLETE_MESH);
  }
  
}

void EdgeMesher::equal_meshing(ModelEnt *ent, int num_edges, std::vector<double> &coords)
{
  double umin, umax, measure;
    // get the u range for the edge
  iGeom::Error gerr = mk_core()->igeom_instance()->getEntURange(ent->geom_handle(), umin, umax);
  IBERRCHK(gerr, "Trouble getting parameter range for edge.");
  if (umin == umax) throw Error(MK_BAD_GEOMETRIC_EVALUATION, "Edge evaluated to same parameter umax and umin.");
  
    // get the arc length
  measure = ent->measure();

  int err;
  double x, y, z, u, du;
  du = (umax-umin)/(double)num_edges;

  u=umin;
  for(int i = 1; i < num_edges; i++)
  {
    u = umin + i*du;
    gerr = mk_core()->igeom_instance()->getEntUtoXYZ(ent->geom_handle(), u, coords[3*i], coords[3*i+1], coords[3*i+2]);
    IBERRCHK(gerr, "Trouble getting U from XYZ along an edge.");
  }
}

}

