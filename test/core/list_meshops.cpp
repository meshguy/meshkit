#include "meshkit/MKCore.hpp"
#include "moab/CN.hpp"
#include <iostream>
#include <iomanip>

using namespace MeshKit;
using namespace std;

const char* str( iBase_EntityType t ) 
{
  switch (t) {
    case iBase_REGION: return "Region";
    case iBase_FACE:   return "Face";
    case iBase_EDGE:   return "Edge";
    case iBase_VERTEX: return "Vertex";
    default:           return "(unknown)";
  }
}

const char* str( moab::EntityType t ) 
{
  return moab::CN::EntityTypeName(t);
}

int main()
{
  MKCore::MeshOpFactory* f = MKCore::op_factory();
  std::vector<MKCore::OpInfo>::const_iterator i;
  std::vector<iBase_EntityType>::const_iterator b;
  std::vector<moab::EntityType>::const_iterator m;
  int have_none = 1;
  for (i = f->registeredOps.begin(); i != f->registeredOps.end(); ++i) {
    have_none = 0;
    
    cout << i->opIndex << " " << i->opName << " " << "{ ";
    if (!i->modelEntTypes.empty()) {
      b = i->modelEntTypes.begin();
      cout << str(*b); 
      for (++b; b != i->modelEntTypes.end(); ++b)
        cout << ", " << str(*b);
    }
    cout << " } { ";
    if (!i->meshEntTypes.empty()) {
      m = i->meshEntTypes.begin();
      cout << str(*m); 
      for (++m; m != i->meshEntTypes.end(); ++m)
        cout << ", " << str(*m);
    }
    
    cout << " }" << endl;
    
  }
  
  return have_none;
}

  
