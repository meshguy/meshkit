/** \file SweepWrapper.hpp
 *  \brief 
 */

#ifndef MSQ_SWEEP_WRAPPER_HPP
#define MSQ_SWEEP_WRAPPER_HPP

#include <Wrapper.hpp>
#include <string>

namespace MESQUITE_NS {

class SweepWrapper : public Wrapper {
public:
  /**
   *\param max_termination_vertex_movement Stop optimization when, for a 
   *            single optimization step, no vertex moves more than this
   *            amount.  Select a value approprite for your mesh size.
   *\param src_mesh_coord_tag_name  Name of tag on target mesh containing
   *            coordinates of source mesh vertices.
   */
  SweepWrapper( double max_termination_vertex_movement,
                const char* src_mesh_coord_tag_name )
    : initMeshTag(src_mesh_coord_tag_name),
      maxVtxMovement(max_termination_vertex_movement) {}
protected:
  MESQUITE_EXPORT
  void run_wrapper( Mesh* mesh, ParallelMesh* pmesh,
                    MeshDomain* geom, Settings* settings,
                    QualityAssessor* qa, MsqError& err );
private:
  std::string initMeshTag;
  double maxVtxMovement;
};

} // namespace MESQUITE_NS

#endif
