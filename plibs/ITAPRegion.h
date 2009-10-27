#ifndef ITAP_REGION_H
#define ITAP_REGION_H

#include <gmsh/GRegion.h>

#include <iGeom.h>

#include "ITAPFace.h"

class ITAPRegion : public GRegion 
{
 public:
  ITAPRegion(GModel *m, iGeom_Instance &g, iBase_EntityHandle *eHandle);

  virtual ~ITAPRegion() {}

  void * getNativePtr() const { return (void*)&cellHandle; }

 private:
    static iBase_TagHandle geom_id_tag;
    static bool tag_available;

    void setGeomTag();

    iBase_EntityHandle *cellHandle;
    iGeom_Instance   geometry;
};

#endif
