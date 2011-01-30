/*********************************************
June,10
Reactor Assembly Mesh Assembler
Argonne National Laboratory

CCrgen class declaration
*********************************************/
#ifndef __RGG_MESH_H__
#define __RGG_MESH_H__

#include <iostream>
#include <cfloat>
#include <math.h>
#include "MKUtils.hpp"
#include "MergeMesh.hpp"
#include "ExtrudeMesh.hpp"
#include "CopyMesh.hpp"
#include "CopyGeom.hpp"
#include "parser.hpp"
#include "fileio.hpp"
#include "clock.hpp"
#include <sstream>
#include <string>
#include "matrixtemplate.hpp"
#include "utils.hpp"

#ifdef MOAB
#include "MBiMesh.hpp"
#include "MBInterface.hpp"
#include "MBRange.hpp"
#include "MBSkinner.hpp"
#endif

class CCrgen
{
public:
  CCrgen ();    // ctor
  ~CCrgen ();   // dtor
  int prepareIO (int argc, char *argv[]);
  int load_meshes();
  int load_geometries();
  int read_inputs_phase1 ();
  int read_inputs_phase2 ();
  int write_makefile ();
  int find_assm(const int i, int &assm_index);
  int banner();
  int copy_move(); 
  int merge_nodes();
  int assign_gids();
  int save_mesh();
  int save_geometry();
  int close();
  int extrude();
  int move_verts(iBase_EntitySetHandle set, const double *dx);
  int move_geoms(iBase_EntitySetHandle set, const double *dx);
  int create_neumannset();
  int copy_move_hex_flat_assys(CopyMesh **cm,
			       const int nrings, const int pack_type,
			       const double pitch,
			       const int symm,
			       std::vector<std::string> &core_alias,
			       std::vector<iBase_EntitySetHandle> &assys);

  int copy_move_sq_assys(CopyMesh **cm,
			 const int nrings, const int pack_type,
			 const double pitch,
			 const int symm,
			 std::vector<std::string> &core_alias,
			 std::vector<iBase_EntitySetHandle> &assys);

  int copy_move_hex_full_assys(CopyMesh **cm,
			       const int nrings, const int pack_type,
			       const double pitch,
			       const int symm,
			       std::vector<std::string> &core_alias,
			       std::vector<iBase_EntitySetHandle> &assys);

  int copy_move_hex_vertex_assys(CopyMesh **cm,
				 const int nrings, const int pack_type,
				 const double pitch,
				 const int symm,
				 std::vector<std::string> &core_alias,
				 std::vector<iBase_EntitySetHandle> &assys);


  int copy_move_one_twelfth_assys(CopyMesh **cm,
				  const int nrings, const int pack_type,
				  const double pitch,
				  const int symm,
				  std::vector<std::string> &core_alias,
				  std::vector<iBase_EntitySetHandle> &assys);

  // phase 1's
  int copy_move_hex_vertex_assys_p1(CopyMesh **cm,
				    const int nrings, const int pack_type,
				    const double pitch,
				    const int symm,
				    std::vector<std::string> &core_alias,
				    std::vector<iBase_EntitySetHandle> &assys);

  int copy_move_hex_flat_assys_p1(CopyMesh **cm,
				  const int nrings, const int pack_type,
				  const double pitch,
				  const int symm,
				  std::vector<std::string> &core_alias,
				  std::vector<iBase_EntitySetHandle> &assys);

  int copy_move_sq_assys_p1(CopyMesh **cm,
			    const int nrings, const int pack_type,
			    const double pitch,
			    const int symm,
			    std::vector<std::string> &core_alias,
			    std::vector<iBase_EntitySetHandle> &assys);

  int copy_move_hex_full_assys_p1(CopyMesh **cm,
				  const int nrings, const int pack_type,
				  const double pitch,
				  const int symm,
				  std::vector<std::string> &core_alias,
				  std::vector<iBase_EntitySetHandle> &assys);

  int copy_move_one_twelfth_assys_p1(CopyMesh **cm,
				     const int nrings, const int pack_type,
				     const double pitch,
				     const int symm,
				     std::vector<std::string> &core_alias,
				     std::vector<iBase_EntitySetHandle> &assys);

// function for geometries
  int copy_move_hex_flat_assys(CopyGeom **cg,
			       const int nrings, const int pack_type,
			       const double pitch,
			       const int symm,
			       std::vector<std::string> &core_alias,
			       std::vector<iBase_EntitySetHandle> &assys);

  int copy_move_sq_assys(CopyGeom **cg,
			 const int nrings, const int pack_type,
			 const double pitch,
			 const int symm,
			 std::vector<std::string> &core_alias,
			 std::vector<iBase_EntitySetHandle> &assys);

  int copy_move_hex_full_assys(CopyGeom **cg,
			       const int nrings, const int pack_type,
			       const double pitch,
			       const int symm,
			       std::vector<std::string> &core_alias,
			       std::vector<iBase_EntitySetHandle> &assys);

  int copy_move_hex_vertex_assys(CopyGeom **cg,
				 const int nrings, const int pack_type,
				 const double pitch,
				 const int symm,
				 std::vector<std::string> &core_alias,
				 std::vector<iBase_EntitySetHandle> &assys);


  int copy_move_one_twelfth_assys(CopyGeom **cg,
				  const int nrings, const int pack_type,
				  const double pitch,
				  const int symm,
				  std::vector<std::string> &core_alias,
				  std::vector<iBase_EntitySetHandle> &assys);

  // phase 1's
  int copy_move_hex_vertex_assys_p1(CopyGeom **cg,
				    const int nrings, const int pack_type,
				    const double pitch,
				    const int symm,
				    std::vector<std::string> &core_alias,
				    std::vector<iBase_EntitySetHandle> &assys);

  int copy_move_hex_flat_assys_p1(CopyGeom **cg,
				  const int nrings, const int pack_type,
				  const double pitch,
				  const int symm,
				  std::vector<std::string> &core_alias,
				  std::vector<iBase_EntitySetHandle> &assys);

  int copy_move_sq_assys_p1(CopyGeom **cg,
			    const int nrings, const int pack_type,
			    const double pitch,
			    const int symm,
			    std::vector<std::string> &core_alias,
			    std::vector<iBase_EntitySetHandle> &assys);

  int copy_move_hex_full_assys_p1(CopyGeom **cg,
				  const int nrings, const int pack_type,
				  const double pitch,
				  const int symm,
				  std::vector<std::string> &core_alias,
				  std::vector<iBase_EntitySetHandle> &assys);

  int copy_move_one_twelfth_assys_p1(CopyGeom **cg,
				     const int nrings, const int pack_type,
				     const double pitch,
				     const int symm,
				     std::vector<std::string> &core_alias,
				     std::vector<iBase_EntitySetHandle> &assys);


  iMesh_Instance impl;
  iGeom_Instance geom;
  
#ifdef HAVE_MOAB
  MBInterface* mbImpl() {return reinterpret_cast<MBiMesh*> (impl)->mbImpl;};
#endif

  bool extrude_flag;
  bool mem_tflag;
  std::string prob_type;

private:

  CopyMesh **cm;
  MergeMesh *mm;
  CopyGeom **cg;
  iBase_EntitySetHandle root_set;
  std::vector<iBase_EntitySetHandle> assys;

  // declare variables read in the inputs
  int err;
  int UNITCELL_DUCT, ASSY_TYPES ;
  int nrings, nringsx, nringsy, pack_type, symm;
  double pitch, pitchx, pitchy;
  bool global_ids, back_mesh, nst_flag, nsb_flag, nss_flag;
  std::vector<std::string> files;
  std::string outfile;
  int nassys; // the number of mesh files
  int tot_assys; // total no. of assms forming core
  int set_DIM; // default is 3D
  double PI;
  double z_height;    // z_height for extruding surfaces mesh
  int z_divisions; // z_divisions for extruding surface mesh
  int nst_Id, nsb_Id, nss_Id;
  
  // file related
  std::ifstream file_input;    // File Input
  std::ofstream make_file;    // File Output
  std::string iname, ifile, mfile, geometry, back_meshfile, geom_engine, nsLoc;
  int linenumber;
  std::string card,geom_type, meshfile, mf_alias, temp_alias;
  std::vector<std::string> assm_alias;
  std::vector<std::string> core_alias;
  
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
  MKUtils *mu;
};
#endif
