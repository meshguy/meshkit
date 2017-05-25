
/** \file SweepWrapper.cpp
 *  \brief 
 *  
 */

#include "moab/mesquite/Mesquite.hpp"
#include "SweepWrapper.hpp"

#include "moab/mesquite/InstructionQueue.hpp"
#include "moab/mesquite/PMeanPTemplate.hpp"
#include "moab/mesquite/TrustRegion.hpp"
#include "moab/mesquite/QualityAssessor.hpp"
#include "moab/mesquite/RefMeshTargetCalculator.hpp"
#include "moab/mesquite/ReferenceMesh.hpp"
#include "moab/mesquite/TagVertexMesh.hpp"
#include "moab/mesquite/TQualityMetric.hpp"
#include "moab/mesquite/ConjugateGradient.hpp"
#include "moab/mesquite/SteepestDescent.hpp"
#include "moab/mesquite/QuasiNewton.hpp"
#include "moab/mesquite/NonSmoothDescent.hpp"


//#include "TRel2DShape.hpp"
//#include "TRel3DShape.hpp"

#include "moab/mesquite/TShapeNB1.hpp"

#include "moab/mesquite/MsqError.hpp"

#include "moab/mesquite/MsqMatrix.hpp"

namespace MESQUITE_NS {

class TargetFlipper : public TargetCalculator
{
  public:
    TargetFlipper( TargetCalculator* unflipped ) : realTargets(unflipped) {}
    
    bool get_3D_target( PatchData& pd, 
                        size_t element,
                        Sample sample,
                        MsqMatrix<3,3>& W_out,
                        MsqError& err );

    bool get_2D_target( PatchData& pd, 
                        size_t element,
                        Sample sample,
                        MsqMatrix<2,2>& W_out,
                        MsqError& err );
    bool have_surface_orient() const;

    bool get_surface_target( PatchData& pd, 
                                   size_t element,
                                   Sample sample,
                                   MsqMatrix<3,2>& W_out,
                                   MsqError& err );
  private:
    TargetCalculator* realTargets;
};
bool TargetFlipper::get_3D_target( PatchData& pd, 
                        size_t element,
                        Sample sample,
                        MsqMatrix<3,3>& W_out,
                        MsqError& err )
{
	return realTargets->get_3D_target(pd, element, sample, W_out, err);

}
    
bool TargetFlipper::get_2D_target( PatchData& pd, 
                                   size_t element,
                                   Sample sample,
                                   MsqMatrix<2,2>& W_out,
                                   MsqError& err )
{
  return realTargets->get_2D_target( pd, element, sample, W_out, err );
}

bool TargetFlipper::have_surface_orient() const
{ return false; }

bool TargetFlipper::get_surface_target( PatchData& pd, 
                                   size_t element,
                                   Sample sample,
                                   MsqMatrix<3,2>& W_out,
                                   MsqError& err )
{
  bool valid = realTargets->get_surface_target( pd, element, sample, W_out, err );
  if (MSQ_CHKERR(err) || !valid) return false;
  
  if (pd.domain_set()) {
    Vector3D norm;
    pd.get_domain_normal_at_sample( element, sample, norm, err );
    MSQ_ERRZERO(err);
    MsqVector<3> c0(W_out.column(0)), c1(W_out.column(1));
    if ((c0 * c1) % MsqVector<3>(&norm[0]) < 0.0) {
      W_out.set_column(0, c1);
      W_out.set_column(1, c0);
    }
  }
  return true;
}

void SweepWrapper::run_wrapper( Mesh* mesh,
                                ParallelMesh* pmesh,
                                MeshDomain* domain,
                                Settings* settings,
                                QualityAssessor* qa,
                                MsqError& err )
{
    // construct target calculator
  TagVertexMesh source_mesh( err, mesh, false, initMeshTag.c_str() ); MSQ_ERRRTN(err);
  ReferenceMesh ref_mesh( &source_mesh );
  RefMeshTargetCalculator ref_tc( &ref_mesh );
  TargetFlipper flipped_targets( &ref_tc );
  
    // construct objective function
  //TRel2DShape shape2d;
  //TRel3DShape shape3d;
  TShapeNB1 shape;
  //TShapeNB1 shape3d;
  TQualityMetric metric(&flipped_targets, &shape);
  PMeanPTemplate objfunc( 2.0, &metric );
  
    // construct smoother
  //TrustRegion solver(&objfunc);
  //ConjugateGradient solver(&objfunc);
  SteepestDescent solver(&objfunc);
  //QuasiNewton solver(&objfunc);  


  TerminationCriterion terminate;
  terminate.add_absolute_vertex_movement( maxVtxMovement );
  terminate.write_iterations( "mesquite.gpt", err );
  solver.set_inner_termination_criterion( &terminate );
  
    // construct instruction queue
  qa->add_quality_assessment( &metric );
  InstructionQueue q;
  q.add_quality_assessor( qa, err ); 
  MSQ_ERRRTN(err);
  q.set_master_quality_improver( &solver, err ); 
  MSQ_ERRRTN(err);
  q.add_quality_assessor( qa, err ); 
  MSQ_ERRRTN(err);
  
    // smooth mesh
  MeshDomainAssoc mesh_and_domain(mesh, domain);
  q.run_common( &mesh_and_domain, pmesh, settings, err );
  MSQ_ERRRTN(err);
}

} // namespace MESQUITE_NS
