#include "meshkit/ParPostRecv.hpp"
#include "moab/Types.hpp"
#include "RefEntity.hpp"
#include "TDParallel.hpp"

namespace MeshKit 
{
// static registration of this mesh scheme
  moab::EntityType ParPostRecv_tps[] = {moab::MBVERTEX, moab::MBEDGE, moab::MBTRI, moab::MBTET, moab::MBMAXTYPE};
const moab::EntityType* ParPostRecv::output_types()
{ return ParPostRecv_tps; }
  
ParPostRecv::ParPostRecv(MKCore *mkcore, const MEntVector &me_vec)
  : MeshScheme(mkcore, me_vec)
{
  // get information related to MOAB parallel communication
  m_mpcomm = moab::ParallelComm::get_pcomm(mk_core()->moab_instance(), 0);
  if (NULL == m_mpcomm) {
    throw Error(iBase_FAILURE, "Parallel communication should be already created.");
  }
  m_rank = m_mpcomm->proc_config().proc_rank();
}

ParPostRecv::~ParPostRecv()
{
}

void ParPostRecv::setup_this()
{
}

void ParPostRecv::execute_this()
{
  iMesh::Error err;
  for (MEntSelection::iterator mit = mentSelection.begin(); mit != mentSelection.end(); mit++) {
    ModelEnt *me = mit->first;
    RefEntity* entity = reinterpret_cast<RefEntity*> (me->geom_handle());
    TDParallel *td_par = (TDParallel *) entity->get_TD(&TDParallel::is_parallel);
    if (td_par == NULL) ECERRCHK(MK_FAILURE, "Exchange entity should have partitioned information.");
    int charge_p = td_par->get_charge_proc();

    if (m_rank == charge_p) { // send proc
      // put processors sharing this interface vertex to buffer
      DLIList<int>* proc_list = td_par->get_shared_proc_list();
      int n_proc_list = proc_list->size();
      proc_list->reset();
      
      for (int i = 0; i < n_proc_list; i++) {
        int proc = proc_list->get_and_step();
        if (proc != m_rank) {
          int ind = get_shared_list(proc);
        }
      }
    }
    else {
      get_shared_list(charge_p); // receive proc: just create shared entity/proc list
    }
  }

  // post irecv
  moab::ErrorCode rval = m_mpcomm->post_irecv(m_shared_procs);
  MBERRCHK(rval, mk_core()->moab_instance());
}

int ParPostRecv::get_shared_list(const int proc)
{
  int ind = -1;
  std::vector<unsigned int>::iterator vit = 
    std::find(m_shared_procs.begin(), m_shared_procs.end(), proc);
  if (vit == m_shared_procs.end()) {
    ind = m_shared_procs.size();
    m_shared_procs.push_back(proc);
  }
  else ind = vit - m_shared_procs.begin();

  return ind;
}
} // namespace MeshKit
