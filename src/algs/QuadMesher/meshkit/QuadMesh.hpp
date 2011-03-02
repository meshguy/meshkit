//-----------------------------------C++-------------------------------------//
// File: src/algs/QuadMesher.hpp
// Wednesday February 11 10:50 2011
// Brief: QuadMesher class definition: four schemes are provided: equal meshing,
//        Bias Meshing, Dual Bias Meshing, Curvature-based meshing 
//---------------------------------------------------------------------------//

#ifndef MESHKIT_QUADMESHER_HPP
#define MESHKIT_QUADMESHER_HPP

#include "meshkit/MeshScheme.hpp"

namespace MeshKit
{
//===========================================================================//
  /*!
   * \class QuadMesher
   * \brief Triangle to Quad Transformation.
   * 
   */
//===========================================================================//

using namespace std;

namespace Jaal { class Mesh; }

class QuadMesher : public MeshScheme
{
public:
        enum MeshCleanOps {LOCAL_MESH_CLEANUP, GLOBAL_MESH_CLEANUP };

	//construction function for edge mesher
	QuadMesher(MKCore *mk_core, const MEntVector &me_vec);

	~QuadMesher() {}

	//set up the parameters for edge meshing, e.g. compute the number of intervals
	virtual void setup_this();

	//Generate the edge mesh
	virtual void execute_this();

       /**\brief Get class name */
       static const char* name() { return "QuadMesher"; }

       /**\brief Function returning whether this scheme can mesh entities of t
        *        the specified dimension.
        *\param dim entity dimension
        */
       static bool can_mesh(iBase_EntityType dim)
          { return iBase_FACE == dim; }

       /** \brief Function returning whether this scheme can mesh the specified entity
        * 
        * Used by MeshOpFactory to find scheme for an entity.
        * \param me ModelEnt being queried
        * \return If true, this scheme can mesh the specified ModelEnt
        */
       static bool can_mesh(ModelEnt *me)
          { return canmesh_face(me); }


       /**\brief Get list of mesh entity types that can be generated.
        *\return array terminated with \c moab::MBMAXTYPE
        */
       static const moab::EntityType* output_types();

       /** \brief Return the mesh entity types operated on by this scheme
        * \return array terminated with \c moab::MBMAXTYPE
        */
       virtual const moab::EntityType* mesh_types_arr() const
         { return output_types(); }

       /**  \brief Should we allow boundary modification. By Default: No.
 	*    When boundary is not allowed to modify and if the number of
        *    segments are not even, "All Quads" not possible. 
        */
       void   allow_boundary_steiner_points( bool a = 0 );

       /**   \brief Performs mesh clean to reduce the irreguarity of the
        *    of the quad mesh 
        */
       void   mesh_cleanup( MeshCleanOps mcleanup = GLOBAL_MESH_CLEANUP);

private:
      iMesh_Instance imesh;

      // Base: Everythiing must convert to Jaal format for the time being..
      Jaal::Mesh* tri_quad_conversion (Jaal::Mesh *trimesh);  

      // Interface with iMesh. Internally convert to Jaal format and call the previous function..
      Jaal::Mesh* tri_quad_conversion (iMesh_Instance imesh);
};

}

#endif
