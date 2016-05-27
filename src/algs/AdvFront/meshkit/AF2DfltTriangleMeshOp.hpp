/*
 * AF2DfltTriangleMeshOp.hpp
 *
 * An AF2DfltTriangleMeshOp is a meshing operation that applies a
 * two-dimensional advancing front algorithm.  The algorithm uses a
 * set of rules that add only triangular faces.
 * 
 * Like any advancing front algorithm, it requires at least one
 * mesh edge to start the algorithm.  Thus surfaces like spheres and
 * tori that have no one-dimensional boundary should have an artificial
 * one-dimensional boundary introduced.  If not, the algorithm will fail.
 */
#ifndef AF2DFLTTRIANGLEMESHOP_HPP
#define AF2DFLTTRIANGLEMESHOP_HPP

// C++
#include <list>

// MeshKit
#include "meshkit/MeshScheme.hpp"

namespace MeshKit
{

class AF2DfltTriangleMeshOp : public MeshScheme
{
  private:

    static moab::EntityType meshTypes[];

    /**
     * \brief Get the index of the sizing function that is expected
     *   to be the effective sizing function for the specified ModelEnt
     *   if/when this MeshOp is executed
     *
     * The sizing function set on the ModelEnt will be used only if the
     * size is greater than zero. 
     *
     * \return
     */
    int getSFIndex() const;

  public:

    /**
     * \brief Implement RegisterMeshOp<AF2DfltTriangleMeshOp>::can_mesh(iBase_EntityType)
     *
     * \return true if the dimension is iBase_FACE, false otherwise
     */
    static bool can_mesh(iBase_EntityType dimension);

    /**
     * \brief Implement RegisterMeshOp<AF2DfltTriangleMeshOp>::can_mesh(ModelEnt*)
     *
     * \return true if the ModelEnt is a face, i.e., if MeshOp::canmesh_face returns true, false otherwise
     */
    static bool can_mesh(ModelEnt *me);

    /**
     * \brief Implement RegisterMeshOp<AF2DfltTriangleMeshOp>::name()
     *
     * \return "AF2DfltTriangleMeshOp"
     */
    static const char* name();

    /**
     * \brief Implement RegisterMeshOp<AF2DfltTriangleMeshOp>::output_types()
     *
     * Return an appropriately terminated array indicating that this MeshOp
     * can add vertices and triangles to the mesh.
     *
     * \return a pointer to an array containing the elements
     *   moab::MBVERTEX
     *   moab::MBTRI
     *   moab::MBMAXTYPE
     */
    static const moab::EntityType* output_types();

    /**
     * \brief Constructor
     *
     * Construct an AF2DfltTriangleMeshOp from the MKCore context that
     * it will operate in and the vector of MeshKit entities, which should
     * correspond to geometric surfaces, that it is supposed to mesh.
     *
     * \param meshkitCore a pointer to the MKCore instance
     * \param meshEntVec the vector of entities that the
     *   AF2DfltTriangleMeshOp should attempt to mesh
     */
    AF2DfltTriangleMeshOp(MKCore *meshkitCore, const MEntVector &meshEntVec);

    /**
     * \brief Destructor
     */
    virtual ~AF2DfltTriangleMeshOp();

    /**
     * \brief Copy constructor
     *
     * Initializes a copy of this using the superclass's copy
     * constructor.  In the current implementation (May 2016), this
     * will ultimately create a new graph node and add it to the graph.
     *
     * \param toCopy the AF2DfltTriangleMeshOp that should be
     *   copied as this AF2DfltTriangleMeshOp is constructed
     */
    AF2DfltTriangleMeshOp(const AF2DfltTriangleMeshOp& toCopy);

    /**
     * \brief Assignment operator (throws exception)
     *
     * The superclasses do not implement an assignment operator.
     * Thus this class does not support an assignment operator, and in
     * order to be explicit about that, it throws an exception if
     * an assignment is attempted.
     *
     * \param rhs the AF2DfltTriangleMeshOp that is to be
     *   assigned to this AF2DfltTriangleMeshOp
     * \return a reference to this AF2DfltTriangleMeshOp after
     *   the assignment has occurred
     */
    AF2DfltTriangleMeshOp& operator=(const AF2DfltTriangleMeshOp& rhs);

    /**
     * \brief Execute an advancing front algorithm with the default
     *   triangle rules
     *
     * This method depends on the completed execution of edge meshing
     * for all edges bounding the face.  It extracts the meshes of
     * the bounding geometric edges and uses them to initialize the
     * advancing front.
     *
     * If the algorithm succeeds, the new vertices and triangles are
     * committed to the mesh in MOAB.
     */
    virtual void execute_this();

    /**
     * \brief Get the mesh entity types that this mesh operation might create
     * 
     * \return a pointer to an array containing the elements
     *   moab::MBVERTEX
     *   moab::MBTRI
     *   moab::MBMAXTYPE
     */
    virtual const moab::EntityType* mesh_types_arr() const;

    /**
     * \brief Implement GraphNode::setup_this() to set up edge meshers
     *   on any bounding facets that don't have edge meshing set up
     *
     * This method also ensures that the graph dependencies from this node
     * to whatever edge meshing needs to be done are set up.
     */
    virtual void setup_this();
};

} // end namespace MeshKit

#endif
