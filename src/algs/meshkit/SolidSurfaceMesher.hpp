

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

  class SolidSurfaceMesher : public MeshScheme
  {
   public:
    /** \brief Construction function for the Solid CAD Surface Mesher
     */
     SolidSurfaceMesher(MKCore *mk, const MEntVector &ments);

    /** \brief Decstructor function for the Solid CAD Surface Mesher
     */
    ~SolidSurfaceMesher();

    virtual void setup_this();
  
    virtual void execute_this(); 

    //Functions needed to register the new meshing class in MK
    static bool can_mesh(iBase_EntityType dim)
    { return iBase_FACE == dim; }

    static bool can_mesh(ModelEnt *me)
    { return canmesh_face(me); }
 
    static const char* name()
    {return "SolidSurfaceMesher";}

    /** \brief Function that names the output types of this meshing class
     */
    static const moab::EntityType* output_types();
  
    virtual const moab::EntityType* mesh_types_arr() const
    { return output_types(); }

    void set_mesh_params(double faceting_tolerance=0, double geom_resabs=0);

  private:
    double facet_tol; 
    double geom_res; 
    MKCore *mk; 
    MEntVector model_ents; 
    
    /** \brief Function to for creating triangles and adding them to the ModelEntity's meshset
     * \param surf ModelEntity for the surface to be meshed.
     */
    void facet(ModelEnt *surf);

    double vtx2vtx_dist( iGeom::EntityHandle vtx1, iMesh::EntityHandle vtx2); 

  };

}
