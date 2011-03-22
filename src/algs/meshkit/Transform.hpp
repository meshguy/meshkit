#ifndef MESHKIT_TRANSFORM_HPP
#define MESHKIT_TRANSFORM_HPP

#include "meshkit/Error.hpp"
#include "meshkit/iMesh.hpp"
#include "meshkit/Matrix.hpp"
#include "meshkit/TransformBase.hpp"
#include <vector>

namespace MeshKit {
  Vector<3> rodrigues(const Vector<3> &pt, const Vector<3> &z, double dtheta);

  namespace Copy {
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
