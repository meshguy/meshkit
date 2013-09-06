#ifndef MESHKIT_PostBL_HPP
#define MESHKIT_PostBL_HPP

#include <cassert>
#include <string>
#include <vector>
#include <set>
#include <math.h>
#include <iomanip>

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

#include "meshkit/iMesh.hpp"
#include "meshkit/iGeom.hpp"
#include "MBCN.h"

#include "meshkit/SimpleArray.hpp"
#include "meshkit/parser.hpp"
#include "meshkit/clock.hpp"
#include "mstream.hpp"

#include "MBiMesh.hpp"
#include "moab/Interface.hpp"
#include "moab/Range.hpp"
#include "moab/CartVect.hpp"


#include "moab/Skinner.hpp"
#include "moab/AdaptiveKDTree.hpp"
#include "moab/Matrix3.hpp"

/*!
 * \class PostBL
 * \brief Options and Keywords Used in PostBL Algorithm
 *
 *  RUNNING: Postmesh Boundary Layer Tool can be run using the test_postbl executable in test/algs directory
 example:-  test_postbl <name>.inp, where, <name> is the name of the input file containing the keywords below:
 - bias    \verbatim      <double>    bias b/w different layers of boundary layer is always greater than zero. \endverbatim
 - meshfile \verbatim     <string>    input mesh file name support by moab reader. \endverbatim
 - surfaces  \verbatim    <integer>   id of the surface on which boundary layer needs to be created.  \endverbatim
 - neumannset  \verbatim  <integer>   id of the neumann set on which boundary layer needs to be created.  \endverbatim
 - material    \verbatim  <integer>   material id to which the newly created hexes will be assigned, default is 999. \endverbatim
 - thickness  \verbatim   <double>    boundary layer thickness. \endverbatim
 - debug     \verbatim    <1 or 0>    print all the debug o/p if set to 1. A .log file is always created.  \endverbatim
 - outfile   \verbatim    <string>    name of output mesh file, can be any format that is supported by moab.  \endverbatim
 - end        \verbatim               this marks the end of input file for boundary layer generation.  \endverbatim
 - Sample keyword file can be found here: data/test_postbl.inp
*/

namespace MeshKit {

#define DEFAULT_TEST_POSTBL  "test_postbl.inp"

  class MKCore;

  class PostBL : public MeshScheme
  {
  public:
    /* \brief Constructor
     *
     * Create a new PostBL instance
     * \param impl the iGeom instance handle for the Geom
     */
    PostBL(MKCore *mk, const MEntVector &me_vec);

    /* \brief Destructor
     */
    virtual ~PostBL();

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

    //! execute PostBL graph node
    virtual void execute_this();


    /** \brief Improved algorithm
     */
    void Algo2 ();

    /** \brief Prepare input/output files for reading/writing
     *  command line args and testdir for default test case
     * \param TestDir directory where test will be located and command line arguments.
     */
    void PrepareIO (int argc, char *argv[], std::string TestDir);

    /** \brief get the normals given connectivity of a quad
     * \param conn connectivity array type EntityHandle
     *	\param v return normal vector
     */
    void get_normal_quad (std::vector<moab::EntityHandle>conn, moab::CartVect &v)   ;

    /** \brief get the normals given connectivity of a quad
     * \param conn connectivity array type moab::EntityHandle
     *	\param v return normal vector
     */
    void get_normal_edge (std::vector<moab::EntityHandle>conn, moab::CartVect AB, moab::CartVect &v);

    /** \brief compute determinant of jacobian of a hex element
     *  \param conn connectivity array
     *	\param offset passed when conn array has connectivity of more than one element
     *  \param	detJ is returned
     */
    void get_det_jacobian (std::vector<moab::EntityHandle> conn, int offset, double &detJ);

    /** \brief compute edge length from the vertex specified
     *  \param vector of quads/edges
     *	\param vertex specified
     *  \param	o/p edge length
     */
    void find_min_edge_length (moab::Range, moab::EntityHandle, moab::Range, double &e_len);

  private:
    //! iGeom Impl for calling geometry creation/manipulation operations
    iGeom *igeom;

    //! iMesh Impl for calling mesh creation/manipulation operations
    iMesh *imesh;

    //! MOAB Impl for calling mesh creation/manipulation operations
    moab::Interface *mb;

    // ! parser related
    bool debug, hybrid, check_bl_edge_length;
    // !! file Input
    std::ifstream m_FileInput;
    mstream m_LogFile;
    std::string szInputString;
    std::string szComment;
    int MAXCHARS, m_nLineNumber;
    
    // ! variables to parse
    std::string m_InputFile, m_MeshFile, m_OutFile, m_LogName, m_MeshType;
    int m_SurfId, m_NeumannSet, m_Material, m_HConn;
    double m_Thickness, m_MinEdgeLength;
    int m_Intervals, m_JacCalls, tri_sch, fixmat;
    double m_Bias, m_JLo, m_JHi;
    // ! variable for hex and tet meshes
    int m_Conn, m_BElemNodes;
    int m_GD, m_BLDim;
    std::string m_Card;
    int err;

    // ! error handlers
    enum ErrorStates { INVALIDINPUT};
    void IOErrorHandler (ErrorStates) const;
  };

  inline const char* PostBL::name()
  {
    return "PostBL";
  }

  inline bool PostBL::can_mesh(iBase_EntityType)
  {
    return false;
  }

  inline bool PostBL::can_mesh(ModelEnt *)
  {
    return true;
  }

  inline const moab::EntityType* PostBL::mesh_types_arr() const
  {
    return output_types();
  }

} // namespace MeshKit
#endif
