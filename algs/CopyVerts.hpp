#ifndef COPYVERTS_HPP
#define COPYVERTS_HPP

#include "iBase.h"
#include "iMesh.h"

class CopyVerts
{
public:
  void operator ()(int n,iBase_EntityHandle *src, int src_size,
                   iBase_EntityHandle **dest, int *dest_alloc, int *dest_size)
    const;

  void operator ()(iBase_EntityHandle *src, int src_size,
                   iBase_EntityHandle **dest, int *dest_alloc, int *dest_size)
    const
  {
    operator ()(max_,src,src_size,dest,dest_alloc,dest_size);
  }
protected:
  CopyVerts(iMesh_Instance impl,int max);
  virtual void transform(int n,int i,double *coords) const = 0;
private:
  iMesh_Instance impl_;
  int max_;
};

class CopyMoveVerts : public CopyVerts
{
public:
  CopyMoveVerts(iMesh_Instance impl, const double *dv, int max=1);
protected:
  virtual void transform(int n, int i, double *coords) const;
private:
  double dv_[3];
};

class CopyRotateVerts : public CopyVerts
{
public:
  CopyRotateVerts(iMesh_Instance impl, const double *origin, const double *z,
                  double theta, int max=1);
protected:
  virtual void transform(int n, int i, double *coords) const;
private:
  double origin_[3];
  double z_[3];
  double theta_;
};

#endif
