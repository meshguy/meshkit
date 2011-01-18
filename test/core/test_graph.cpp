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
      {std::cout << "myScheme (setup), node " << opName << std::endl;
      }
  
  inline void execute_this() 
      {std::cout << "myScheme (execute), node " << opName << std::endl;
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
  mk->meshop_graph().addArc(mk->root_node(), A->op_node());
  mk->meshop_graph().addArc(mk->root_node(), B->op_node());
  mk->meshop_graph().addArc(A->op_node(), C->op_node());
  mk->meshop_graph().addArc(A->op_node(), D->op_node());
  mk->meshop_graph().addArc(B->op_node(), mk->leaf_node());
  mk->meshop_graph().addArc(C->op_node(), mk->leaf_node());
  mk->meshop_graph().addArc(D->op_node(), mk->leaf_node());
  
    // now traverse
  mk->setup();
  
  mk->execute();

  MeshOpFactory::instance()->destroy_instance();
}

  
