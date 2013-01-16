//-----------------------------------C++-------------------------------------//
// File: src/algs/meshkit/Pbl.hpp
//
// Brief: Pbl class definition:
//       
//         class, Pbl
//---------------------------------------------------------------------------// 

#ifndef MESHKIT_Pbl_HPP
#define MESHKIT_Pbl_HPP

#include <cassert>
#include <string>
#include <vector>
#include <set>

#include "meshkit/Types.hpp"
#include "meshkit/Error.hpp"
#include "meshkit/MeshScheme.hpp"
#include "meshkit/ModelEnt.hpp"
#include "meshkit/MKCore.hpp"
#include "meshkit/SizingFunction.hpp"
#include "meshkit/RegisterMeshOp.hpp"

#include "meshkit/LocalSet.hpp"
#include "meshkit/LocalTag.hpp"
#include "meshkit/Matrix.hpp"
#include "../CopyUtils.hpp"

#include "meshkit/iMesh.hpp"
#include "meshkit/iGeom.hpp"
#include "MBCN.h"

#include "SimpleArray.hpp"
#include "../AssyGen/parser.hpp"
#include "../AssyGen/clock.hpp"

#ifdef HAVE_MOAB
#include "MBiMesh.hpp"
#include "MBInterface.hpp"
#include "MBRange.hpp"
#include "MBCartVect.hpp"
#endif



namespace MeshKit {

#define DEFAULT_TEST_PBL  "test_pbl.inp"
  class MKCore;

  class Pbl : public MeshScheme
  {
  public:
    /* \brief Constructor
     *
     * Create a new Pbl instance
     * \param impl the iGeom instance handle for the Geom
     */
    Pbl(MKCore *mk, const MEntVector &me_vec);

    /* \brief Destructor
     */
    virtual ~Pbl();

    /**\brief Get class name */
    static const char* name();

    /**\brief Function returning whether this scheme can mesh entities of t
     *        the specified dimension.
     *\param dim entity dimension
     */
    static bool can_mesh(iBase_EntityType dim);

    /** \brief Function returning whether this scheme can mesh the specified entity
     *
     * Used by MeshOpFactory to find scheme for an entity.
     * \param me ModelEnt being queried
     * \return If true, this scheme can mesh the specified ModelEnt
     */
    static bool can_mesh(ModelEnt *me);

    /**\brief Get list of mesh entity types that can be generated.
     *\return array terminated with \c moab::MBMAXTYPE
     */
    static const moab::EntityType* output_types();

    /** \brief Return the mesh entity types operated on by this scheme
     * \return array terminated with \c moab::MBMAXTYPE
     */
    virtual const moab::EntityType* mesh_types_arr() const;

    /** \brief Re-implemented here so we can check topological dimension of model_ent
     * \param model_ent ModelEnt being added
     */
    virtual bool add_modelent(ModelEnt *model_ent);

    //! Setup is a no-op, but must be provided since it's pure virtual
    virtual void setup_this();

    //! The only setup/execute function we need, since meshing vertices is trivial
    virtual void execute_this();

    /** \brief Prepare input/output files for reading/writing
     *  command line args and testdir for default test case
     */
    void PrepareIO (int argc, char *argv[], std::string TestDir);

    /** \brief get the normals given connectivity of a quad
     *  x, y, z is the normal
     */
    void get_normal_quad (std::vector<EntityHandle>conn, double &x, double &y, double &z);

    /** \brief compute determinant of jacobian of a hex element
     *  input connectivity array
     */
    void get_det_jacobian (std::vector<EntityHandle> conn, int offset, double &detJ);

  private:
    // iGeom Impl for calling geometry creation/manipulation operations
    iGeom *igeom;

    // iMesh Impl for calling mesh creation/manipulation operations
    iMesh *imesh;
    MBInterface *mb;
    // parser related
      // file Input
    std::ifstream m_FileInput; 
    std::string szInputString;
    std::string szComment;
    int MAXCHARS, m_nLineNumber;
    
    // variables to parse
    std::string m_InputFile, m_MeshFile, m_OutFile;
    int m_SurfId, m_NeumannSet;
    double m_Thickness;
    int m_Intervals;
    double m_Bias;

    int m_GD;
    std::string m_Card;
    int err;

    // error handlers
    enum ErrorStates { INVALIDINPUT};
    void IOErrorHandler (ErrorStates) const;
  };

  inline const char* Pbl::name()
  {
    return "Pbl";
  }

  inline bool Pbl::can_mesh(iBase_EntityType)
  {
    return false;
  }

  inline bool Pbl::can_mesh(ModelEnt *)
  {
    return true;
  }

  inline const moab::EntityType* Pbl::mesh_types_arr() const
  {
    return output_types();
  }

} // namespace MeshKit
#endif
