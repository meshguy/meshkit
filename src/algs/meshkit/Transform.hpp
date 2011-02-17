#ifndef TRANSFORM_HPP
#define TRANSFORM_HPP

#include "meshkit/Error.hpp"
#include "meshkit/iMesh.hpp"
#include "meshkit/Matrix.hpp"
#include <vector>
#include <tr1/functional>

namespace MeshKit {
  namespace Copy {
    typedef std::tr1::function<void(iMesh &, iMesh::EntityHandle*, int,
                                    iMesh::EntityHandle**, int*, int*)>
            AnyTransform;

    /** \class Transform Transform.hpp "meshkit/Transform.hpp"
     * \brief A base class for transforming copied meshes
     *
     * This is the common base class used to transform vertices in copied
     * meshes. Subclasses of this type implement particular transformation
     * functions, e.g. translation or rotation.
     */
    class Transform
    {
    public:
      // TODO: change this to use vectors

      /** \brief Transform the selected vertices
       * \param impl the iMesh implementation holding the vertices
       * \param src a pointer to an array of the source vertices
       * \param src_size the number of source vertices
       * \param dest a pointer to a pointer-to-array of the destination
       *        vertices
       * \param dest_alloc the amount of memory allocated for dest
       * \param dest_size, the number of destination vertices
       */
      void operator ()(iMesh &impl, iMesh::EntityHandle *src, int src_size,
                       iMesh::EntityHandle **dest, int *dest_alloc,
                       int *dest_size) const;
    protected:
      virtual void transform(Vector<3> &coords) const = 0;
    };

    /** \class Identity Transform.hpp "meshkit/Transform.hpp"
     * \brief The identity transformation
     */
    class Identity : public Transform
    {
    public:
      Identity();
    protected:
      virtual void transform(Vector<3> &coords) const;
    };

    /** \class Translate Transform.hpp "meshkit/Transform.hpp"
     * \brief A translation function
     */
    class Translate : public Transform
    {
    public:
      Translate(const Vector<3> &dv);
    protected:
      virtual void transform(Vector<3> &coords) const;
    private:
      Vector<3> dv_;
    };

    /** \class Rotate Transform.hpp "meshkit/Transform.hpp"
     * \brief A rotation function, using Rodrigues' formula
     */
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
    typedef std::tr1::function<void(int, iMesh &, iMesh::EntityHandle*, int,
                                    iMesh::EntityHandle**, int*, int*)>
            AnyTransform;

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
