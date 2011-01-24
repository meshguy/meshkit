#ifndef __EDGEMESHER_HPP
#define __EDGEMESHER_HPP

#include "meshkit/MeshScheme.hpp"
#include "iGeom.hh"
#include <vector>

namespace MeshKit 
{
    
struct Point3D
{
  double px;
  double py;
  double pz;	
};

class EdgeMesher : public MeshScheme
{
public:
	
  enum EdgeSchemeType {equalMesh=0, biasMesh, dualMesh, curvatureMesh};
	
public:
  EdgeMesher(MKCore *mk_core, const MEVector &me_vec);

  virtual ~EdgeMesher();

  static MeshOp *factory(MKCore *mkcore, const MEVector &me_vec);
  
  void mesh_types(std::vector<moab::EntityType> &tps);

  double measure(iGeom::EntityHandle ent, double ustart, double uend) const;
  
  void set_edge_scheme(EdgeSchemeType scheme);

  EdgeSchemeType get_edge_scheme() const;

    // no setup or execute functions, can use parent class versions without modification

    //! Adds bounding vertices to (single) vertexmesher instance
  virtual void setup_this();

    //! Generates the edge mesh
  virtual void execute_this();
	
private:
	
  void equal_meshing(ModelEnt *ent, int num_edges, std::vector<double> &coords);

  EdgeSchemeType schemeType;
};

inline EdgeMesher::EdgeMesher(MKCore *mk_core, const MEVector &me_vec) 
        : MeshScheme(mk_core, me_vec)
{
}

inline EdgeMesher::~EdgeMesher()
{
}

inline void EdgeMesher::set_edge_scheme(EdgeMesher::EdgeSchemeType scheme)
{
  schemeType = scheme;
}

inline EdgeMesher::EdgeSchemeType EdgeMesher::get_edge_scheme() const
{
  return schemeType;
}

} // namespace MeshKit

#endif
