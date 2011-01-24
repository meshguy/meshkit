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
#include "meshkit/MeshOpFactory.hpp"

using namespace MeshKit;

class MyScheme;

class MyScheme : public MeshKit::MeshOp 
{
public:
  static MeshOp *factory(MKCore *, const MEVector &);
  
  MyScheme(MKCore*, const MEVector &);
  
  inline void setup_this() 
      {std::cout << "myScheme (setup), node " << get_name() << std::endl;
      }
  
  inline void execute_this() 
      {std::cout << "myScheme (execute), node " << get_name() << std::endl;
      }

  inline void mesh_types(std::vector<moab::EntityType> &tps) 
      {
      }
  
};

inline MyScheme::MyScheme(MKCore *mk_core, const MEVector & me_vec) 
        : MeshOp(mk_core, me_vec)
{}
  
inline MeshOp *MyScheme::factory(MKCore *mk_core, const MEVector &me_vec)
{
  return new MyScheme(mk_core, me_vec);
}

int main(int argc, char **argv) 
{
    // start up MK and register my scheme with it
  MeshOpFactory::instance()->register_meshop("MyScheme", moab::MBMAXTYPE, MyScheme::factory);
  MKCore *mk = MeshOpFactory::instance()->mk_core();

    // create the scheme objects
  MeshOp *A = MeshOpFactory::instance()->construct_meshop("MyScheme"),
      *B = MeshOpFactory::instance()->construct_meshop("MyScheme"),
      *C = MeshOpFactory::instance()->construct_meshop("MyScheme"),
      *D = MeshOpFactory::instance()->construct_meshop("MyScheme");

  A->set_name("A");
  B->set_name("B");
  C->set_name("C");
  D->set_name("D");
  
    // put them in the graph
  mk->get_graph().addArc(mk->root_node()->get_node(), A->get_node());
  mk->get_graph().addArc(mk->root_node()->get_node(), B->get_node());
  mk->get_graph().addArc(A->get_node(), C->get_node());
  mk->get_graph().addArc(A->get_node(), D->get_node());
  mk->get_graph().addArc(B->get_node(), mk->leaf_node()->get_node());
  mk->get_graph().addArc(C->get_node(), mk->leaf_node()->get_node());
  mk->get_graph().addArc(D->get_node(), mk->leaf_node()->get_node());
  
    // now traverse
  mk->setup();
  
  mk->execute();

  MeshOpFactory::instance()->destroy_instance();
}

  
