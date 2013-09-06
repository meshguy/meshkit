//-----------------------------------C++-------------------------------------//
// File: src/algs/meshkit/CoreGen.hpp
//
// Brief: CoreGen class definition:
//        Creates reactor core model from input mesh files
//         class, CoreGen
//---------------------------------------------------------------------------//

#ifndef MESHKIT_COREGEN_HPP
#define MESHKIT_COREGEN_HPP

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

#include "meshkit/iMesh.hpp"
#include "meshkit/iGeom.hpp"
#include "MBCN.h"

#include "meshkit/vectortemplate.hpp"
#include "meshkit/matrixtemplate.hpp"
#include "meshkit/parser.hpp"
#include "meshkit/SimpleArray.hpp"
#include "meshkit/clock.hpp"

#ifdef HAVE_MOAB
#include "iMesh_extensions.h"
#include "MBiMesh.hpp"
#include "MBInterface.hpp"
#include "MBRange.hpp"
#include "MBSkinner.hpp"
#include "MBTagConventions.hpp"
#endif


#ifdef USE_MPI
#include "mpi.h"
#include "iMeshP.h"
#include "moab_mpi.h"
#include "moab/ParallelMergeMesh.hpp"
#include "moab/ParallelComm.hpp"
#include "MBParallelConventions.h"
#endif

#include "../../meshkit/CopyGeom.hpp"
#include "../../meshkit/CopyMesh.hpp"
#include "../../meshkit/MergeMesh.hpp"
#include "../../meshkit/ExtrudeMesh.hpp"
#include "../../meshkit/CESets.hpp"

namespace MeshKit {

#define DEFAULT_TEST_FILE  "CoreGen_default"
#define TEST_FILE_NAME "CoreGen_default"

  class MKCore;

  class CoreGen : public MeshScheme
  {
  public:
    /* \brief Constructor
     *
     * Create a new CoreGen instance
     * \param impl the iGeom instance handle for the Geom
     */
    CoreGen(MKCore *mk, const MEntVector &me_vec);

    /* \brief Destructor
     */
    virtual ~CoreGen();

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

    enum ErrorStates {INVALIDINPUT, ENEGATIVE};
    int prepareIO (int argc, char *argv[], int nrank, int numprocs, std::string  TestDir);
    int load_meshes();
    int load_meshes_parallel(const int, int);
    int distribute_mesh(const int,  int);
    int load_geometries();
    int read_inputs_phase1 ();
    int read_inputs_phase2 ();
    int write_makefile ();
    int write_minfofile ();
    int find_assm(const int i, int &assm_index);
    int banner();
    int copymove(const int nrank, const int numprocs);
    int copymove_all(const int nrank, const int numprocs);
    int set_copymove_coords();
    int merge_nodes();
    int merge_nodes_parallel(const int nrank, const int numprocs);
    int assign_gids();
    int assign_gids_parallel(const int nrank, const int numprocs);
    int save_mesh();
    int save_mesh(int rank);
    int save_mesh_parallel(const int nrank, const int numprocs);
    int save_geometry();
    int close();
    int close_parallel(const int nrank, const int numprocs);
    int extrude();
    int move_verts(iBase_EntitySetHandle set, const double *dx);
    int move_geoms(iBase_EntitySetHandle set, const double *dx);
    int create_neumannset();

    bool extrude_flag;
    bool mem_tflag;
    std::string prob_type, savefiles, info, minfo;
    std::vector<std::string> files, all_meshfiles, mk_files;
    std::vector<int> assm_meshfiles;
    std::vector< std::vector<int> > assm_location;
    std::vector<std::vector<int> > position_core;
    std::vector<int> meshfile_proc;
    std::vector<double> x_coord;
    std::vector<double> y_coord;
    bool nst_flag, nsb_flag, nss_flag;
    std::vector<std::string> core_alias;
    std::vector<double> nsx, nsy, nsc;
    int num_nsside;
    std::string  DIR;
  private:
    //! iGeom Impl for calling geometry creation/manipulation operations
    iGeom *igeom;

    //! iMesh Impl for calling mesh creation/manipulation operations
    iMesh *imesh;

    //! MOAB Impl for calling mesh creation/manipulation operations
    MBInterface *mb;

    std::vector <CopyMesh*> cm;
    std::vector <CopyGeom*> cg;

    iBase_EntitySetHandle root_set;
    std::vector<iBase_EntitySetHandle> assys;
    std::vector<int> assys_index;
    // declare variables read in the inputs
    int rank, procs, err;
    int UNITCELL_DUCT, ASSY_TYPES ;
    int nrings, nringsx, nringsy, pack_type, symm;
    double pitch, pitchx, pitchy;
    bool global_ids, back_mesh;
    std::string outfile, mesh_info;
    int nassys; // the number of mesh files
    int tot_assys; // total no. of assms forming core
    int set_DIM; // default is 3D
    double PII;
    double z_height;    // z_height for extruding surfaces mesh
    int z_divisions; // z_divisions for extruding surface mesh
    int nst_Id, nsb_Id;
    std::vector<int> nss_Id;

    // file related
    std::ifstream file_input;    // File Input
    std::ofstream make_file, info_file, minfo_file;    // File Output
    std::string iname, ifile, mfile, geometry, back_meshfile, geom_engine, nsLoc, infofile, minfofile;
    int linenumber;
    std::string card,geom_type, meshfile, mf_alias, temp_alias;
    std::vector<std::string> assm_alias;

    // parsing related
    std::string input_string;
    std::string comment ;
    int MAXCHARS ;

    // merge related
    double merge_tol;
    int do_merge;
    int update_sets;
    iBase_TagHandle merge_tag;

    // MKUtils obj, assigning gid's etc.
    //  MKUtils *mu;
    // error handler
    void IOErrorHandler (ErrorStates) const;
    CClock Timer;
    std::string szDateTime;
    int run_flag;
    clock_t sTime;
#ifdef USE_MPI
    moab::ParallelComm *pc;
#endif

  };

  inline const char* CoreGen::name()
  {
    return "CoreGen";
  }

  inline bool CoreGen::can_mesh(iBase_EntityType)
  {
    return false;
  }

  inline bool CoreGen::can_mesh(ModelEnt *)
  {
    return true;
  }

  inline const moab::EntityType* CoreGen::mesh_types_arr() const
  {
    return output_types();
  }

} // namespace MeshKit
#endif
