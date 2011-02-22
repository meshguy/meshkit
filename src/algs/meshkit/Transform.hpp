#ifndef MESHKIT_TRANSFORM_HPP
#define MESHKIT_TRANSFORM_HPP

#include "meshkit/Error.hpp"
#include "meshkit/iMesh.hpp"
#include "meshkit/Matrix.hpp"
#include <vector>

namespace MeshKit {
  inline double * vec2ptr(std::vector< Vector<3> > &v) {
    return reinterpret_cast<double *>(&v[0]);
  }

  Vector<3> rodrigues(const Vector<3> &pt, const Vector<3> &z, double dtheta);

  namespace Copy {
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
      /** \brief Transform the selected vertices
       * \param mesh the iMesh implementation holding the vertices
       * \param src an array of the source vertices
       * \param dest an array of the destination vertices
       */
      virtual void transform(iMesh *mesh,
                             const std::vector<iMesh::EntityHandle> &src,
                             std::vector<iMesh::EntityHandle> &dest) const = 0;

      /** \brief Clone this transform object
       */
      virtual Transform * clone() const = 0;
    };

    /** \class BasicTransform Transform.hpp "meshkit/Transform.hpp"
     * \brief A utility class for transforming copied meshes
     *
     * This class template simplifies the creation of new transformation
     * functions. To use this, create a new class and inherit from this one,
     * passing the new class name as the template parameter (the curiously-
     * recurring template pattern).
     *
     * Example:
     *    class Example : public BasicTransform<Example>
     *    {
     *      friend class BasicTransform<Example>;
     *    public:
     *      Example() { ... }
     *    protected:
     *      void transform_one(Vector<3> &coords) const { ... }
     *    };
     */
    template<typename T>
    class BasicTransform : public Transform
    {
    public:
      /** \brief Transform the selected vertices
       * \param mesh the iMesh implementation holding the vertices
       * \param src an array of the source vertices
       * \param dest an array of the destination vertices
       */
      virtual void transform(iMesh *mesh,
                             const std::vector<iMesh::EntityHandle> &src,
                             std::vector<iMesh::EntityHandle> &dest) const
      {
        std::vector< Vector<3> > coords(src.size());
        IBERRCHK(mesh->getVtxArrCoords(&src[0], src.size(), iBase_INTERLEAVED,
                                       vec2ptr(coords)), *mesh);

        for (size_t i=0; i<coords.size(); i++)
          static_cast<const T*>(this)->transform_one(coords[i]);

        if (&src == &dest)
        {
          IBERRCHK(mesh->setVtxArrCoords(&src[0], src.size(), iBase_INTERLEAVED,
                                         vec2ptr(coords)), *mesh);
        }
        else {
          dest.resize(src.size());
          IBERRCHK(mesh->createVtxArr(src.size(), iBase_INTERLEAVED,
                                      vec2ptr(coords), &dest[0]), *mesh);
        }
      }

      /** \brief Clone this transform object
       */
      virtual Transform * clone() const
      {
        return new T(*static_cast<const T*>(this));
      }
    protected:
      BasicTransform() {}
    };

    /** \class Identity Transform.hpp "meshkit/Transform.hpp"
     * \brief The identity transformation
     */
    class Identity : public BasicTransform<Identity>
    {
      friend class BasicTransform<Identity>;
    public:
      Identity() {}
    protected:
      void transform_one(Vector<3> &coords) const {}
    };

    /** \class Translate Transform.hpp "meshkit/Transform.hpp"
     * \brief A translation function
     */
    class Translate : public BasicTransform<Translate>
    {
      friend class BasicTransform<Translate>;
    public:
      Translate(const Vector<3> &dv) : dv_(dv) {}
    protected:
      void transform_one(Vector<3> &coords) const
      {
        coords += dv_;
      }
    private:
      Vector<3> dv_;
    };

    /** \class Rotate Transform.hpp "meshkit/Transform.hpp"
     * \brief A rotation function, using Rodrigues' formula
     */
    class Rotate : public BasicTransform<Rotate>
    {
      friend class BasicTransform<Rotate>;
    public:
      Rotate(const Vector<3> &origin, const Vector<3> &z, double dtheta)
        : origin_(origin), z_(z/length(z)), dtheta_(dtheta)
      {}
    protected:
      void transform_one(Vector<3> &coords) const
      {
        coords = rodrigues(coords-origin_, z_, dtheta_) + origin_;
      }
    private:
      Vector<3> origin_;
      Vector<3> z_;
      double dtheta_;
    };
  }

  namespace Extrude {
    /** \class Transform Transform.hpp "meshkit/Transform.hpp"
     * \brief A base class for transforming extruded meshes
     *
     * This is the common base class used to transform vertices in extruded
     * meshes. Subclasses of this type implement particular transformation
     * functions, e.g. translation or rotation.
     */
    class Transform
    {
    public:
      /** \brief Transform the selected vertices
       * \param step the step number for the extrusion, with 0 being the
       *        already-existing mesh data
       * \param mesh the iMesh implementation holding the vertices
       * \param src an array of the source vertices
       * \param dest an array of the destination vertices
       */
      virtual void transform(int step, iMesh *mesh,
                             const std::vector<iMesh::EntityHandle> &src,
                             std::vector<iMesh::EntityHandle> &dest) const = 0;
      
      /** \brief The number of steps in this extrusion
       */
      virtual int steps() const = 0;

      /** \brief Clone this transform object
       */
      virtual Transform * clone() const = 0;
    };

    /** \class BasicTransform Transform.hpp "meshkit/Transform.hpp"
     * \brief A utility class for transforming extruded meshes
     *
     * This class template simplifies the creation of new transformation
     * functions. To use this, create a new class and inherit from this one,
     * passing the new class name as the template parameter (the curiously-
     * recurring template pattern).
     *
     * Example:
     *    class Example : public BasicTransform<Example>
     *    {
     *      friend class BasicTransform<Example>;
     *    public:
     *      Example() { ... }
     *    protected:
     *      void transform_one(int step, Vector<3> &coords) const { ... }
     *    };
     */
    template<typename T>
    class BasicTransform : public Transform
    {
    public:
      /** \brief Transform the selected vertices
       * \param step the step number for the extrusion, with 0 being the
       *        already-existing mesh data
       * \param mesh the iMesh implementation holding the vertices
       * \param src an array of the source vertices
       * \param dest an array of the destination vertices
       */
      virtual void transform(int step, iMesh *mesh,
                             const std::vector<iMesh::EntityHandle> &src,
                             std::vector<iMesh::EntityHandle> &dest) const
      {
        std::vector< Vector<3> > coords(src.size());
        IBERRCHK(mesh->getVtxArrCoords(&src[0], src.size(), iBase_INTERLEAVED,
                                      vec2ptr(coords)), *mesh);

        for (size_t i=0; i<coords.size(); i++)
          static_cast<const T*>(this)->transform_one(step, coords[i]);

        dest.resize(src.size());
        IBERRCHK(mesh->createVtxArr(src.size(), iBase_INTERLEAVED,
                                    vec2ptr(coords), &dest[0]), *mesh);
      }

      /** \brief The number of steps in this extrusion
       */
      virtual int steps() const
      {
        return steps_;
      }

      /** \brief Clone this transform object
       */
      virtual Transform * clone() const
      {
        return new T(*static_cast<const T*>(this));
      }
    protected:
      BasicTransform(int steps) : steps_(steps)
      {}

      int steps_;
    };

    /** \class Translate Transform.hpp "meshkit/Transform.hpp"
     * \brief A translation function
     */
    class Translate : public BasicTransform<Translate>
    {
      friend class BasicTransform<Translate>;
    public:
      Translate(const Vector<3> &dv, int steps)
        : BasicTransform<Translate>(steps), dv_(dv/steps)
      {}
    protected:
      virtual void transform_one(int step, Vector<3> &coords) const
      {
        coords += dv_*step;
      }
    private:
      Vector<3> dv_;
    };

    /** \class Rotate Transform.hpp "meshkit/Transform.hpp"
     * \brief A rotation function, using Rodrigues' formula
     */
    class Rotate : public BasicTransform<Rotate>
    {
      friend class BasicTransform<Rotate>;
    public:
      Rotate(const Vector<3> &origin, const Vector<3> &z, double dtheta,
             int steps)
        : BasicTransform<Rotate>(steps), origin_(origin), z_(z/length(z)),
        dtheta_(dtheta/steps)
      {}
    protected:
      virtual void transform_one(int step, Vector<3> &coords) const
      {
        coords = rodrigues(coords-origin_, z_, dtheta_*step) + origin_;
      }
    private:
      Vector<3> origin_;
      Vector<3> z_;
      double dtheta_;
    };
  }
}

#endif
