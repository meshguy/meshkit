/** \file test_graph.cpp
 *
 * Test the MeshOp graph code in MeshKit.  Do this by building a graph
 * which looks like:
 *
 *            R
 *          /   \
 *         A     B
 *        / \    |
 *       C   D  /
 *         \ | /
 *          \|/
 *           L
 *  and traversing it using BFS and RBFS (reverse BFS).  Do this by registering
 *  a MeshOp that does nothing but print its name during the setup/execute
 *  functions.
 *
 */

#include "meshkit/MKCore.hpp"
#include "meshkit/MeshOp.hpp"
#include "meshkit/RegisterMeshOp.hpp"

using namespace MeshKit;

class MyScheme;

class MyScheme : public MeshKit::MeshOp 
{
public:
  MyScheme(MKCore*, const MEntVector &);
  
  inline void setup_this() 
      {std::cout << "myScheme (setup), node " << get_name() << std::endl;
      }
  
  inline void execute_this() 
      {std::cout << "myScheme (execute), node " << get_name() << std::endl;
      }

  inline void mesh_types(std::vector<moab::EntityType> &tps) 
      {
      }
      
  static const char* name() { return "MyScheme"; }
  static bool can_mesh(iBase_EntityType dim) { return dim == iBase_REGION; }
  static bool can_mesh(ModelEnt* ent) { return canmesh_region(ent); }
  static const moab::EntityType* output_types() 
    { 
      static moab::EntityType end = moab::MBMAXTYPE;
      return &end;
    }
  const moab::EntityType* mesh_types_arr() const
    { return output_types(); }
  
};

inline MyScheme::MyScheme(MKCore *mk_core, const MEntVector & me_vec) 
        : MeshOp(mk_core, me_vec)
{}

//---------------------------------------------------------------------------//
RegisterMeshOp<MyScheme> INIT;
//---------------------------------------------------------------------------//

int main(int argc, char **argv) 
{
    // start up MK and register my scheme with it
  MKCore mk;

    // create the scheme objects
  MeshOp *A = mk.construct_meshop("MyScheme"),
      *B = mk.construct_meshop("MyScheme"),
      *C = mk.construct_meshop("MyScheme"),
      *D = mk.construct_meshop("MyScheme");

  A->set_name("A");
  B->set_name("B");
  C->set_name("C");
  D->set_name("D");
  
    // put them in the graph
  mk.get_graph().addArc(mk.root_node()->get_node(), A->get_node());
  mk.get_graph().addArc(mk.root_node()->get_node(), B->get_node());
  mk.get_graph().addArc(A->get_node(), C->get_node());
  mk.get_graph().addArc(A->get_node(), D->get_node());
  mk.get_graph().addArc(B->get_node(), mk.leaf_node()->get_node());
  mk.get_graph().addArc(C->get_node(), mk.leaf_node()->get_node());
  mk.get_graph().addArc(D->get_node(), mk.leaf_node()->get_node());
  
    // now traverse
  mk.setup();
  
  mk.execute();
}

  
