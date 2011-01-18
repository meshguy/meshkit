#ifndef MKCORE_HPP
#define MKCORE_HPP

/** \file MKCore.hpp
 */
#include "meshkit/Types.h"
#include "meshkit/Error.hpp"
#include "iGeom.hh"
#include "moab/Interface.hpp"
#include "iMesh.hh"
#include "iRel.hh"
#include "lemon/list_graph.h"
#include <vector>

// outside the namespace 'cuz it just is
class MBiMesh;

namespace MeshKit {

      //! Forward declare since we store a vector of these
class SizingFunction;

/** \class MKCore MKCore.hpp "meshkit/MKCore.hpp"
 * \brief The core MeshKit instance
 *
 * The MKCore object stores MeshKit-wide data like the geometry/mesh/relations instances,
 * and is the object through which other MeshKit class objects are accessed.  Since it is
 * a child class of MeshOp, the MKCore instance can also be a node in the MeshOp graph.  By
 * convention, the MKCore instance serves as the root of the this (directed) graph.
 *
 * If the MKCore constructor is called with no arguments, then instances of the geometry, mesh,
 * and relations databases are created inside the constructor; these instances are then deleted
 * in the MKCore destructor.  If this is not the desired behavior, either pass in non-NULL
 * instances to the MKCore constructor, or re-set the iCreatedIgeom, iCreatedMoab,
 * and/or iCreatedIrel flags in this class.
 */

class MKCore
{
public:

    /** \brief Constructor
     * \param igeom iGeom instance
     * \param moab MOAB instance
     * \param mbi MBiMesh instance; if non-NULL, should use/point to moab parameter
     * \param irel iRel instance
     * \param construct_missing_ifaces If true, constructs the interfaces whose handles are passed in NULL
     */
  MKCore(iGeom *igeom = NULL, 
         moab::Interface *moab = NULL, 
         MBiMesh *mbi = NULL, 
         iRel *irel = NULL,
         bool construct_missing_ifaces = true);
  
    //! destructor
  ~MKCore();

    //! initialize, creating missing geom/mesh/rel interfaces if requested
  void init(bool construct_missing_ifaces);

    /** \brief Load a geometry model from a file
     * \param filename The file to load
     * \param options File options to be passed to the load function
     */
  void load_geometry(const char *filename, const char *options = NULL);

    /** \brief Load a mesh model from a file
     * \param filename The file to load
     * \param options File options to be passed to the load function
     */
  void load_mesh(const char *filename, const char *options = NULL);

    /** \brief Save a geometry model to a file
     * \param filename The file to save
     * \param options File options to be passed to the save function
     */
  void save_geometry(const char *filename, const char *options = NULL);

    /** \brief Save a mesh model to a file
     * \param filename The file to save
     * \param options File options to be passed to the save function
     */
  void save_mesh(const char *filename, const char *options = NULL);

    /** \brief Populate mesh/relations data for geometric entities
     */
  void populate_mesh();

    /** \brief Get model entities of a given dimension
     * \param dim Dimension of entities to get
     * \param model_ents The list these entities get appended to
     */
  void get_entities_by_dimension(int dim, MEVector &model_ents);

    /** \brief Get all model entities
     * \param model_ents The list these entities get appended to
     */
  void get_entities_by_handle(MEVector &model_ents);

    /** \brief Return the iGeom instance pointer
     */
  iGeom *igeom_instance();
  
    /** \brief Return the MOAB instance pointer
     */
  moab::Interface *moab_instance();
  
    /** \brief Return the iMesh instance pointer
     */
  MBiMesh *mb_imesh();
  
    /** \brief Return the iRel instance pointer
     */
  iRel *irel_instance();

    /** \brief Return the iRel pair handle used to relate geometry/mesh entities
     */
  iRel::PairHandle *irel_pair();

    /** \brief Return the (iGeom) tag used to relate geometry entities to mesh entities
     */
  iGeom::TagHandle igeom_model_tag();
  
    /** \brief Return the (MOAB) tag used to relate mesh entities to model entities
     */
  moab::Tag moab_model_tag();

    //! Get the graph
  lemon::ListDigraph &meshop_graph();

    //! Get root node
  lemon::ListDigraph::Node root_node();
  
    //! Get leaf node
  lemon::ListDigraph::Node leaf_node();

    //! Get access to the node map
  lemon::ListDigraph::NodeMap<MeshOp*> &node_map();

    //! Get the MeshOp corresponding to a graph node
  MeshOp *get_meshop(lemon::ListDigraph::Node node);
  
    //! Get the MeshOp corresponding to the source node of this arc
  MeshOp *source(lemon::ListDigraph::Arc arc);
  
    //! Get the MeshOp corresponding to the target node of this arc
  MeshOp *target(lemon::ListDigraph::Arc arc);

    /** \brief Get sizing function by index
     * If the requested index is outside the range of SizingFunction's currently registered,
     * throws an Error.
     * \param index Index of sizing function requested
     * \return SizingFunction* to requested sizing function, NULL of no SizingFunction with that index
     */
  SizingFunction *sizing_function(int index);
  
    /** \brief Add sizing function to those managed by MeshKit
     *
     * The argument to this function is a SizingFunction*; once added, it is MKCore's
     * responsibility to delete this SizingFunction.  Applications can tell MKCore to delete
     * a given SizingFunction (e.g. if it requires lots of memory) by calling delete_sizing_function.
     * \param sf SizingFunction* to be added
     * \return Index of sizing function in MKCore's list of SizingFunction's
     */
  int add_sizing_function(SizingFunction *sf);

    /** \brief Delete sizing function
     *
     * This function removes the referenced sizing function from MKCore's list (setting the
     * corresponding SizingFunction* to NULL, to keep others at the same index position).  
     * Throws an Error if requested sizing function is NULL.
     * \param index Index of SizingFunction to be removed
     */
  void remove_sizing_function(int index);
  
    //! Run setup on the graph
  void setup();

    //! Run execute on the graph
  void execute();
  
private:
    //! Geometry api instance
  iGeom *iGeomInstance;
  
    //! MOAB instance
  moab::Interface *moabInstance;
  
    //! iMesh api instance, for use in iRel
  MBiMesh *mbImesh;
  
    //! IREL api instance
  iRel *iRelInstance;

    //! iRel pair handle used to relate geometry/mesh entities
  iRel::PairHandle *iRelPair;
  
    //! Tag used to associate geometry entities with model entities
  iGeom::TagHandle iGeomModelTag;
  
    //! Tag used to associate mesh entities with model entities
  moab::Tag moabModelTag;

    //! If true, the corresponding interfaces will get deleted from the destructor
  bool iCreatedIgeom, iCreatedMoab, iCreatedMbimesh, iCreatedIrel;

    //! Model entity sets, in array by topological dimension
  moab::Range modelEnts[4];

    //! The MeshOp graph
  lemon::ListDigraph meshopGraph;

    //! Root and leaf nodes of graph
  lemon::ListDigraph::Node rootNode, leafNode;

    //! Map from graph nodes to MeshOps
  lemon::ListDigraph::NodeMap<MeshOp*> nodeMap;

    //! SizingFunction vector
  std::vector<SizingFunction*> sizingFunctions;
};

inline iGeom *MKCore::igeom_instance() 
{
  return iGeomInstance;
}

inline moab::Interface *MKCore::moab_instance()
{
  return moabInstance;
}

inline MBiMesh *MKCore::mb_imesh()
{
  return mbImesh;
}

inline iRel *MKCore::irel_instance()
{
  return iRelInstance;
}

inline iRel::PairHandle *MKCore::irel_pair()
{
  return iRelPair;
}

inline iGeom::TagHandle MKCore::igeom_model_tag()
{
  return iGeomModelTag;
}

inline moab::Tag MKCore::moab_model_tag()
{
  return moabModelTag;
}

    //! Get root node
inline lemon::ListDigraph &MKCore::meshop_graph() 
{
  return meshopGraph;
}

    //! Get root node
inline lemon::ListDigraph::Node MKCore::root_node() 
{
  return rootNode;
}

    //! Get leaf node
inline lemon::ListDigraph::Node MKCore::leaf_node() 
{
  return leafNode;
}
  
    //! Get access to the node map
inline lemon::ListDigraph::NodeMap<MeshOp*> &MKCore::node_map() 
{
  return nodeMap;
}
  
    //! Get the MeshOp corresponding to a graph node
inline MeshOp *MKCore::get_meshop(lemon::ListDigraph::Node node) 
{
  return nodeMap[node];
}
  
inline MeshOp *MKCore::source(lemon::ListDigraph::Arc arc)
{
  return get_meshop(meshopGraph.source(arc));
}

inline MeshOp *MKCore::target(lemon::ListDigraph::Arc arc)
{
  return get_meshop(meshopGraph.target(arc));
}

inline SizingFunction *MKCore::sizing_function(int index) 
{
    // don't check for NULL here 'cuz sometimes we just want to know there isn't one
    // with that index
  if (index >= (int)sizingFunctions.size())
    throw Error(MK_BAD_INPUT, "Sizing function index outside range of valid indices.");
  return sizingFunctions[index];
}
  
} // namespace MeshKit

#endif

  
