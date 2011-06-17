#include "meshkit/MeshScheme.hpp"
#include "meshkit/ModelEnt.hpp"
#include "meshkit/Types.hpp"

namespace MeshKit 
{

void MeshScheme::constrain_even() 
{
    // constrain all edges to be even
  MEntVector edges;
  for (MEntSelection::iterator sit = mentSelection.begin(); sit != mentSelection.end(); sit++)
  {
    (*sit).first->get_adjacencies(1, edges);

    for (MEntVector::iterator vit = edges.begin(); vit != edges.end(); vit++)
      (*vit)->constrain_even(true);
  }
}

} // namespace MeshKit
