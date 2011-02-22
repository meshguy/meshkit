#ifndef MESHKIT_TRANSFORM_HPP
#define MESHKIT_TRANSFORM_HPP

#include "iBase.h"
#include "iMesh.h"

namespace copy {
  class Transform
  {
  public:
    void operator ()(
      iMesh_Instance impl, iBase_EntityHandle *src, int src_size,
      iBase_EntityHandle **dest, int *dest_alloc, int *dest_size) const;
  protected:
    virtual void transform(double *coords) const = 0;
  };

  class Identity : public Transform
  {
  public:
    Identity();
  protected:
    virtual void transform(double *coords) const;
  };

  class Translate : public Transform
  {
  public:
    Translate(const double *dv);
  protected:
    virtual void transform(double *coords) const;
  private:
    double dv_[3];
  };

  class Rotate : public Transform
  {
  public:
    Rotate(const double *origin, const double *z, double dtheta);
  protected:
    virtual void transform(double *coords) const;
  private:
    double origin_[3];
    double z_[3];
    double dtheta_;
  };
}

namespace extrude {
  class Transform
  {
  public:
    void operator ()(
      int step, iMesh_Instance impl, iBase_EntityHandle *src, int src_size,
      iBase_EntityHandle **dest, int *dest_alloc, int *dest_size) const;

    int steps() const { return steps_; }
  protected:
    Transform(int steps);
    virtual void transform(int step, double *coords) const = 0;
  private:
    int steps_;
  };

  class Translate : public Transform
  {
  public:
    Translate(const double *dv, int steps);
  protected:
    virtual void transform(int step, double *coords) const;
  private:
    double dv_[3];
  };

  class Rotate : public Transform
  {
  public:
    Rotate(const double *origin, const double *z, double dtheta, int steps);
  protected:
    virtual void transform(int step, double *coords) const;
  private:
    double origin_[3];
    double z_[3];
    double dtheta_;
  };
}

#endif
