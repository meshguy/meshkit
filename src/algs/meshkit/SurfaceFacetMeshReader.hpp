

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string>
#include <iostream>
#include <fstream>
#include <string.h>
#include <limits.h>

#include <iGeom.h>
#include <iMesh.h>
#include <set>
#include <iRel.h>
#include <vector>
#include "meshkit/MKCore.hpp"
#include "meshkit/MeshScheme.hpp"


namespace MeshKit
{

using namespace std; 

  /** \class SurfaceFacetMeshReader SurfaceFacetMeshReader.hpp "meshkit/SurfaceFacetMeshReader.hpp"
   * \brief A class for generating facet-based mesh of geometric surfaces
   * INPUT: one or more ModelEnts representing geometric surfaces
   * MESH TYPE(S): MBEDGE, MBVERTEX, MBTRI
   * OUTPUT: a set of triangles, edges  and verticesfor each ModelEnt
   * DEPENDENCIES: meshkit build with iGeom interface
   * 
   * This class uses the facets generated for the visualization of solid models in typical CAD 
   * software to represent a geometric surface as mesh. Upon execution, this class will call for
   * the facet data and store it as part of the ModelEnt's mesh.
   */
class SurfaceFacetMeshReader : public MeshScheme
{
public:
       /** \brief Construction function for the Solid CAD Surface Mesher
        */
       SurfaceFacetMeshReader(MKCore *mk, const MEntVector &ments);

      /** \brief Decstructor function for the Solid CAD Surface Mesher
       */
       ~SurfaceFacetMeshReader();

       virtual void setup_this();
  
       virtual void execute_this(); 

       //Functions needed to register the new meshing class in MK
      
       /** \brief Function returning whether this scheme can mesh entities of the 
        *  specified dimension
        *  \param dim entity dimension 
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
 
       /** \brief Get the class name */
       static const char* name()
       {return "SurfaceFacetMeshReader";}

       /** \brief Function that names the output types of this meshing class
        *  \return array terminated with \c moab::MBMAXTYPE
        */
       static const moab::EntityType* output_types();
  
       virtual const moab::EntityType* mesh_types_arr() const
       { return output_types(); }

       void set_mesh_params(double faceting_tolerance=0, double geom_resabs=0);

private:
       double facet_tol; 
       double geom_res; 
       MKCore *mk; 
    
       /** \brief Function to for creating triangles and adding them to the ModelEntity's meshset
        * \param surf ModelEntity for the surface to be meshed.
        */
       void facet(ModelEnt *surf);

       double vtx2vtx_dist( iGeom::EntityHandle vtx1, iMesh::EntityHandle vtx2); 
  };

}
