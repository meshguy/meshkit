#include "ITAPRegion.h"

iBase_TagHandle ITAPRegion:: geom_id_tag = 0; 
bool ITAPRegion:: tag_available = 0; 

void ITAPRegion :: setGeomTag()
{
   if( !tag_available )  {
      int err;
      const char *tag = "GLOBAL_ID";
      int namelen = strlen(tag);
      iGeom_getTagHandle(geometry, tag, &geom_id_tag, &err, namelen);
      tag_available = 1;
    }
}

ITAPRegion::ITAPRegion( GModel *model, iGeom_Instance &g, iBase_EntityHandle *eHandle)
  : GRegion(model, 0), geometry(g), cellHandle(eHandle)
{
  int err;
  int cellID;

  setGeomTag();
  iGeom_getIntData(geometry, *cellHandle, geom_id_tag, &cellID, &err);
  setTag( cellID );

  SimpleArray<iBase_EntityHandle>  faceHandles;
  iGeom_getEntAdj(geometry, *cellHandle, iBase_FACE, ARRAY_INOUT(faceHandles), &err);

  int faceID;
  for( int i = 0; i < faceHandles.size(); i++) 
  {
     iGeom_getIntData(geometry, faceHandles[i], geom_id_tag, &faceID, &err);
     GFace *f = model->getFaceByTag(faceID);
     l_faces.push_back(f);
     f->addRegion(this);
   }
}

