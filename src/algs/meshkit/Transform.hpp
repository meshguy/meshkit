#ifndef TRANSFORM_HPP
#define TRANSFORM_HPP

#include "meshkit/Error.hpp"
#include "meshkit/iMesh.hpp"
#include "src/core/meshkit/Matrix.hpp"
#include <vector>

namespace MeshKit {
  namespace Copy {
    class Transform
    {
    public:
      // TODO: change this to use vectors
      void operator ()(iMesh &impl, iMesh::EntityHandle *src, int src_size,
                       iMesh::EntityHandle **dest, int *dest_alloc,
                       int *dest_size) const;
    protected:
      virtual void transform(Vector<3> &coords) const = 0;
    };

    class Identity : public Transform
    {
    public:
      Identity();
    protected:
      virtual void transform(Vector<3> &coords) const;
    };

    class Translate : public Transform
    {
    public:
      Translate(const Vector<3> &dv);
    protected:
      virtual void transform(Vector<3> &coords) const;
    private:
      Vector<3> dv_;
    };

    class Rotate : public Transform
    {
    public:
      Rotate(const Vector<3> &origin, const Vector<3> &z, double dtheta);
    protected:
      virtual void transform(Vector<3> &coords) const;
    private:
      Vector<3> origin_;
      Vector<3> z_;
      double dtheta_;
    };
  }

  namespace Extrude {
    class Transform
    {
    public:
      // TODO: change this to use vectors
      void operator ()(int step, iMesh &impl, iMesh::EntityHandle *src,
                       int src_size, iMesh::EntityHandle **dest,
                       int *dest_alloc, int *dest_size) const;

      int steps() const { return steps_; }
    protected:
      Transform(int steps);
      virtual void transform(int step, Vector<3> &coords) const = 0;
    private:
      int steps_;
    };

    class Translate : public Transform
    {
    public:
      Translate(const Vector<3> &dv, int steps);
    protected:
      virtual void transform(int step, Vector<3> &coords) const;
    private:
      Vector<3> dv_;
    };

    class Rotate : public Transform
    {
    public:
      Rotate(const Vector<3> &origin, const Vector<3> &z, double dtheta,
             int steps);
    protected:
      virtual void transform(int step, Vector<3> &coords) const;
    private:
      Vector<3> origin_;
      Vector<3> z_;
      double dtheta_;
    };
  }
}

#endif
