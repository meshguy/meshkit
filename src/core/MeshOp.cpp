#include <assert.h>

#include "meshkit/MeshOp.hpp"
#include "meshkit/MKCore.hpp"

namespace MeshKit 
{
    
MeshOp::MeshOp(const MeshOp &mesh_op) 
        : mkCore(mesh_op.mk_core())
{
    // create a node for this meshop
  opNode = mkCore->meshop_graph().addNode();

    // set the nodemap value for this meshop
  mkCore->node_map()[opNode] = this;

    // this is a copy; set in/out edges from copied node
  for (lemon::ListDigraph::InArcIt ait(mkCore->meshop_graph(), mesh_op.op_node());
       ait != lemon::INVALID; ++ait) 
    mkCore->meshop_graph().addArc(mkCore->meshop_graph().source(ait), opNode);
  for (lemon::ListDigraph::OutArcIt ait(mkCore->meshop_graph(), mesh_op.op_node());
       ait != lemon::INVALID; ++ait) 
    mkCore->meshop_graph().addArc(opNode, mkCore->meshop_graph().target(ait));
}
  
MeshOp::MeshOp(MKCore *mkcore, const MEVector &me_vec) : mkCore(mkcore) 
{
    // create a node for this meshop
  opNode = mkCore->meshop_graph().addNode();

    // set the nodemap value for this meshop
  mkCore->node_map()[opNode] = this;
  
    // new graph node; need to link to graph root/leaf
  mkCore->meshop_graph().addArc(mkCore->root_node(), opNode);
  mkCore->meshop_graph().addArc(opNode, mkCore->leaf_node());

  if (!me_vec.empty()) {
    for (MEVector::const_iterator vit = me_vec.begin(); vit != me_vec.end(); vit++)
      moab::Range &rit = meSelection[*vit];
  }
}

bool MeshOp::add_modelent(ModelEnt *model_ent) 
{
  MESelection::iterator sit = meSelection.find(model_ent);
  if (sit != meSelection.end()) return false;
  
  moab::Range &rtmp = meSelection[model_ent];
  return true;
}

bool MeshOp::remove_modelent(ModelEnt *model_ent) 
{
  MESelection::iterator sit = meSelection.find(model_ent);
  if (sit != meSelection.end()) {
    meSelection.erase(sit);
    return true;
  }
  else return false;
}

MeshOp::~MeshOp()
{
  mkCore->meshop_graph().erase(opNode);
}

    //! Setup function, called in reverse order before execute
void MeshOp::setup() 
{
  for (lemon::ListDigraph::InArcIt ait = in_arcs(); ait != lemon::INVALID; ++ait) 
    mkCore->source(ait)->setup();
    
  setup_this();
}
  
    //! Execute function, called in forward order after setup
void MeshOp::execute() 
{
  execute_this();

  for (lemon::ListDigraph::OutArcIt ait = out_arcs(); ait != lemon::INVALID; ++ait) 
    mkCore->target(ait)->execute();
}

lemon::ListDigraph::InArcIt MeshOp::in_arcs() 
{
  return lemon::ListDigraph::InArcIt(mkCore->meshop_graph(), opNode);
}

lemon::ListDigraph::OutArcIt MeshOp::out_arcs() 
{
  return lemon::ListDigraph::OutArcIt(mkCore->meshop_graph(), opNode);
}

}

  
