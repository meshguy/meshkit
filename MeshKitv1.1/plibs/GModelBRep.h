#ifndef _GMODEL_OCC_H_
#define _GMODEL_OCC_H_

#include "gmsh/GmshConfig.h"
#include "gmsh/GModel.h"
#include "gmsh/OCCIncludes.h"

class OCC_Options {
 private:
  int _dummy;
 public:
  OCC_Options(int dummy) : _dummy(dummy){}
};

class OCC_Reader {
 protected :
  TopoDS_Shape shape;
  TopTools_IndexedMapOfShape fmap, emap, vmap, somap, shmap, wmap;
 public:
  enum BooleanOperator { Add, Cut }; 
  OCC_Reader()
  {
    somap.Clear();
    shmap.Clear();
    fmap.Clear();
    wmap.Clear();
    emap.Clear();
    vmap.Clear();
  }
  void healGeometry(double tolerance, bool fixsmalledges, 
                    bool fixspotstripfaces, bool sewfaces, 
                    bool makesolids=false);
  void loadBREP(const char *);  
  void loadShape(const TopoDS_Shape *);
  void buildGModel(GModel *gm);
  void buildLists();
  void removeAllDuplicates(const double &tolerance);

  void Sphere(const SPoint3 &center, const double &radius, const BooleanOperator &op);
  void Cylinder(const SPoint3 &bottom_center, const SVector3 &dir, const BooleanOperator &op);
  void applyBooleanOperator(TopoDS_Shape tool, const BooleanOperator &op);
};

GModel* read_BREP_Model(const std::string &fn);

#endif
