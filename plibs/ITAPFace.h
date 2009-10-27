#ifndef ITAP_FACE_H
#define ITAP_FACE_H

#include <boost/foreach.hpp>

#include <gmsh/GModel.h>
#include <gmsh/GFace.h>
#include "gmsh/ANN/ANN.h"

#include "ITAPVertex.h"
#include "ITAPEdge.h"


#include "SimpleArray.h"

class ITAPFace : public GFace 
{

 public:
  ITAPFace(GModel *m, iGeom_Instance &g, iBase_EntityHandle *f);
  virtual ~ITAPFace(){}

  void setModelEdges( std::list<GEdge*> &k);

  Range<double> parBounds(int i) const; 
  virtual GPoint point(double par1, double par2) const; 
  virtual GPoint closestPoint(const SPoint3 & queryPoint, const double initialGuess[2]) const; 
  virtual SVector3 normal(const SPoint2 &param) const; 
  virtual Pair<SVector3,SVector3> firstDer(const SPoint2 &param) const; 

  void * getNativePtr() const { return faceHandle; }

  virtual SPoint2 parFromPoint(const SPoint3 &) const;
  virtual bool periodic( int dim) const { return _periodic[dim]; }

  bool  hasSeam() const;

  void  create_kdtree();
  void  delete_kdtree();

private:
    static iBase_TagHandle geom_id_tag;
    static bool tag_available;

    double xmin, xmax, xlength;
    double ymin, ymax, ylength;
    double zmin, zmax, zlength;
    double maxlength;
    bool   _periodic[2];

    double *uvCoords;
    int     numNeighs;
    ANNkd_tree  *kdtree;
    ANNpointArray kdnodes;
    mutable vector<ANNidx>   annIdx;
    mutable vector<ANNdist>  anndist;
   
    double umin, umax, vmin, vmax;

    void setGeomTag();

    vector<std::list<GEdge*> > get_icc_list( std::list<GEdge*> &l);
    void addEdgeLoop( std::list<GEdge*> &l);

    iBase_EntityHandle *faceHandle;
    iGeom_Instance geometry;
};
void save( GFace *f);

#endif
