//-----------------------------------C++-------------------------------------//
// File: src/algs/meshkit/CoreGen.hpp
//
// Brief: CoreGen class definition:
//        Creates reactor core model from input mesh files
//         class, CoreGen
//---------------------------------------------------------------------------//

#ifndef MESHKIT_COREGEN_HPP
#define MESHKIT_COREGEN_HPP


#define STRINGIFY_(X) #X
#define STRINGIFY(X) STRINGIFY_(X)
#define COREGEN_DEFAULT_TEST_FILE  "coregen_default"
#define CTEST_FILE_NAME "coregen_default"

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
#include "meshkit/mstream.hpp"

#include "iMesh_extensions.h"
#include "MBiMesh.hpp"
#include "moab/Interface.hpp"
#include "moab/Range.hpp"
#include "moab/Skinner.hpp"
#include "MBTagConventions.hpp"
#include "moab/MergeMesh.hpp"


#ifdef USE_MPI
#include "mpi.h"
#include "iMeshP.h"
#include "moab_mpi.h"
#include "moab/ParallelMergeMesh.hpp"
#include "moab/ParallelComm.hpp"
#include "MBParallelConventions.h"
#endif

#include "meshkit/CopyGeom.hpp"
#include "meshkit/CopyMesh.hpp"
#include "meshkit/ExtrudeMesh.hpp"
#include "meshkit/CESets.hpp"

namespace MeshKit {

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
    int parse_assembly_names(CParser parse, int argc, char *argv[]);
    int load_meshes();
    int load_meshes_more_procs(const int, int);
    int load_meshes_parallel(const int, int);
    int distribute_mesh(const int,  int);
    int load_geometries();
    int read_inputs_phase1 (int argc, char *argv[]);
    int read_inputs_phase2 (int argc, char *argv[]);
    int write_makefile ();
    int write_minfofile ();
    int find_assm(const int i, int &assm_index);
    void banner();
    int copymove(const int nrank, const int numprocs);
    int copymove_all(const int nrank, const int numprocs);
    int set_copymove_coords();
    int save_mesh();
    int save_mesh(int rank);
#ifdef USE_MPI
    int save_mesh_parallel(const int nrank, const int numprocs);
#endif
    int save_geometry();
    int shift_mn_ids(iBase_EntitySetHandle orig_set, int index);
    int extrude();
    int move_verts(iBase_EntitySetHandle set, const double *dx);
    int move_geoms(iBase_EntitySetHandle set, const double *dx);
    int create_neumannset();
    int load_and_compute_meshtogeom(iBase_EntitySetHandle set, std::string filenam);

    bool extrude_flag;
    bool compute_meshtogeom;
    std::vector <int> bsameas;
    bool mem_tflag;
    std::string prob_type, savefiles, info, minfo, same_as, reloading_mf;
    std::vector<std::string> files, all_meshfiles, mk_files;
    std::vector<int> assm_meshfiles,  size_mf, times_loaded;
    std::vector<int> rank_load;
    std::vector<double> load_per_assm;
    std::vector< std::vector<int> > assm_location;
    std::vector<std::vector<int> > position_core;
    std::vector<int> meshfile_proc;
    std::vector<double> x_coord;
    std::vector<double> y_coord;
    bool nst_flag, nsb_flag, nss_flag, nssall_flag;
    std::vector<std::string> core_alias;
    std::vector<double> nsx, nsy, nsc;
    int num_nsside, ms_startid, ns_startid;

  private:
    //! iGeom Impl for calling geometry creation/manipulation operations
    iGeom *igeom;

    //! iMesh Impl for calling mesh creation/manipulation operations
    iMesh *imesh;

    //! MOAB Impl for calling mesh creation/manipulation operations
    moab::Interface *mb;

    std::vector <CopyMesh*> cm;
    std::vector <CopyGeom*> cg;
    ExtrudeMesh *em;

    iBase_EntitySetHandle root_set;
    std::vector<iBase_EntitySetHandle> assys;
    std::vector<int> assys_index;
    // declare variables read in the inputs
    int rank, procs, err;
    int UNITCELL_DUCT, ASSY_TYPES ;
    int nrings, nringsx, nringsy, pack_type, symm;
    double pitch, pitchx, pitchy;
    bool global_ids, back_mesh, have_hex27;
    std::string outfile, mesh_info;
    int nassys; // the number of mesh files
    int tot_assys; // total no. of assms forming core
    int set_DIM; // default is 3D
    double PII;
    double z_height;    // z_height for extruding surfaces mesh
    int z_divisions; // z_divisions for extruding surface mesh
    int nst_Id, nsb_Id, nssall_Id;
    std::vector<int> nss_Id;
    std::string testdir;

    // file related
    std::ifstream file_input;    // File Input
    std::ofstream make_file, info_file, minfo_file;    // File Output
    std::string iname, ifile, mfile, geometry, back_meshfile, geom_engine, nsLoc, meshtogeomfile, infofile, minfofile, logfilename;
    int linenumber;
    std::string card,geom_type, meshfile, mf_alias, temp_alias, etype;
    std::vector<std::string> assm_alias;
    std::vector<int> all_ms_starts, all_ns_starts;

    mstream logfile, meshtogeom_file;

    // parsing related
    std::string input_string;
    std::string comment;
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

    // timing related variables
    double ctload, ctcopymove, ctmerge, ctextrude, ctns, ctgid, ctsave;
    clock_t tload, tcopymove, tmerge, textrude, tns, tgid, tsave;

    // more memory/time related variables
    int ld_t, ld_tload, ld_tcopymove, ld_tsave, ld_tgid, ld_tmerge, ld_tns;
    unsigned long long mem1, mem2, mem3, mem4, mem5, mem6, mem7;
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
