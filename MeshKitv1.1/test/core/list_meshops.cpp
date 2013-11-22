#include "meshkit/MKCore.hpp"
#include "meshkit/MeshOpProxy.hpp"
#include "moab/CN.hpp"
#include <iostream>
#include <iomanip>

using namespace MeshKit;
using namespace std;

int main()
{
  const unsigned n = MKCore::num_meshops();
  if (n == 0)
    return 1; // fail
  
  for (unsigned i = 0; i < n; ++i) {
    MeshOpProxy* p = MKCore::meshop_proxy(i);
    
    cout << i << " " << p->name() << " { ";
    bool first = true;
    for (unsigned j = 0; j < 4; ++j) {
      if (p->can_mesh((iBase_EntityType)j)) {
        if (first)
          first = false;
        else
          cout << ", ";
        cout << j;
      }
    }
    
    cout << " } { ";

    const moab::EntityType* types = p->output_types();
    first = true;
    for (unsigned j = 0; types[j] != moab::MBMAXTYPE; ++j) {
      if (first)
        first = false;
      else 
        cout << ", ";
      cout << moab::CN::EntityTypeName(types[j]);
    }
    
    cout << " }" << endl;
    
  }
  
  return 0;
}

  
