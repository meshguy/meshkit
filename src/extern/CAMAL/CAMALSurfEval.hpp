#ifndef CAMAL_GEOM_EVAL_HPP
#define CAMAL_GEOM_EVAL_HPP

/** \file CAMALSurfEval.hpp
 */

#include "CMLSurfEval.hpp"

namespace MeshKit 
{

class ModelEnt;
    
/** \class CAMALSurfEval CAMALSurfEval.hpp "meshkit/CAMALSurfEval.hpp"
 * \brief The MeshKit-based surface evaluator for CAMAL meshing algorithms
 *
 * This class implements surface evaluation functions required by CAMAL algorithms in terms of functions
 * on ModelEnt objects (and, in some cases, functions in iGeom).
 */
class CAMALSurfEval : public CMLSurfEval
{
public:
    /** \brief Constructor
     * \param ment The ModelEnt to which this evaluator applies
     */
  CAMALSurfEval(ModelEnt *ment);

    /** \brief Destructor
     */
  virtual ~CAMALSurfEval();

    /** \brief Area of the surface
     */
  virtual double area();

    /** \brief Surface bounding box
     * \param box_min Lower corner of the Cartesian-aligned box
     * \param box_max Upper corner of the Cartesian-aligned box
     */
  virtual void bounding_box(double box_min[3], double box_max[3]);

    /** \brief Move the point to the closest point on the surface
     * \param x X point
     * \param y Y point
     * \param z Z point
     */
  virtual void move_to_surface(double& x, double& y, double& z);

    /** \brief Move the point to the closest point on the surface, based on a guess for the UV parameters
     * \param x X point
     * \param y Y point
     * \param z Z point
     * \param u_guess U hint
     * \param v_guess V hint
     */
  virtual void move_to_surface(double& x, double& y, double& z,
			       double& u_guess, double& v_guess);

    /** \brief Compute the normal to the surface at a point
     * \param x X point
     * \param y Y point
     * \param z Z point
     * \param nx X component of the normal
     * \param ny Y component of the normal
     * \param nz Z component of the normal
     * \return Not sure, maybe whether or not the point is on the surface?
     */
  virtual bool normal_at(double x, double y, double z, 
                         double& nx, double& ny, double& nz);

    /** \brief Compute the normal to the surface at a point, starting from a hint on the UV parameters
     * \param x X point
     * \param y Y point
     * \param z Z point
     * \param u_guess U hint
     * \param v_guess V hint
     * \param nx X component of the normal
     * \param ny Y component of the normal
     * \param nz Z component of the normal
     * \return \a true if normal has unit length (normalized),
     * \a false otherwise
     */
  virtual bool normal_at(double x, double y, double z, 
                         double& u_guess, double& v_guess,
                         double& nx, double& ny, double& nz);
  
    /** \brief Return whether the underlying surface is planar
     * \return \a true if the surface is planar
     */
  virtual bool is_planar();

    /** \brief Return whether the surface is parametric
     * \return \a true if the surface is parametric
     */
  virtual bool is_parametric();

    /** \brief Return whether the surface is periodic in the U direction
     * \param u_period Period in the U direction if it's periodic in that direction
     * \return \a true if the surface is periodic in the U direction
     */
  virtual bool is_periodic_in_u(double& u_period);

    /** \brief Return whether the surface is periodic in the V direction
     * \param v_period Period in the V direction if it's periodic in that direction
     * \return \a true if the surface is periodic in the V direction
     */
  virtual bool is_periodic_in_v(double& v_period);

    /** \brief Return the parameter range in U
     * \param u_low Lower bound on the parameter range
     * \param u_high Upper bound on the parameter range
     */
  virtual void get_param_range_u(double& u_low, double& u_high);

    /** \brief Return the parameter range in V
     * \param v_low Lower bound on the parameter range
     * \param v_high Upper bound on the parameter range
     */
  virtual void get_param_range_v(double& v_low, double& v_high);
  
    /** \brief Get the parametric coordinates from the position
     * \param x X point
     * \param y Y point
     * \param z Z point
     * \param u U coordinates
     * \param v V coordinates
     * \return \a true if successful, \a false otherwise
     */
  virtual bool uv_from_position(double x, double y, double z, 
                                double& u, double& v);

    /** \brief Get the parametric coordinates and closest point on a surface from a position
     * \param x X point
     * \param y Y point
     * \param z Z point
     * \param u U coordinates
     * \param v V coordinates
     * \param cx Closest point on the surface to the X point
     * \param cy Closest point on the surface to the Y point
     * \param cz Closest point on the surface to the Z point
     * \return \a true if successful, \a false otherwise
     */
  virtual bool uv_from_position(double x, double y, double z, 
                                double& u, double& v,
                                double& cx, double& cy, double& cz);

    /** \brief Get the spatial coordinates given parametric coordinates
     * \param u U coordinates
     * \param v V coordinates
     * \param x X point
     * \param y Y point
     * \param z Z point
     */
  virtual void position_from_uv(double u, double v, 
				double& x, double& y, double& z);

    /** \brief Get the parametric derivative vectors at a point
     * \param x X point
     * \param y Y point
     * \param z Z point
     * \param du[3] du/d(xyz)
     * \param dv[3] dv/d(xyz)
     */
  virtual void distortion_at_uv(double u, double v, 
                                double du[3], double dv[3]);

private:

    //! ModelEnt associated with this evaluator
  ModelEnt *modelEnt;
};

} // namespace MeshKit

#endif
