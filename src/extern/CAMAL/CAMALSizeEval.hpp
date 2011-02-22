#ifndef MESHKIT_CAMAL_SIZE_EVAL_HPP
#define MESHKIT_CAMAL_SIZE_EVAL_HPP

#include "CMLSizeEval.hpp"

namespace MeshKit 
{

/** \class CAMALSizeEval CAMALSizeEval.hpp "CAMALSizeEval.hpp"
 * \brief The MeshKit-based size evaluator for CAMAL meshing algorithms
 */
class CAMALSizeEval : public ::CMLSizeEval
{
public:
    /** \brief Constructor
     * \param size Size setting for this evaluator
     */
  CAMALSizeEval(double size);
  
    /** \brief Destructor
     */
  virtual ~CAMALSizeEval();
  
    /** \brief Get the desired mesh size at a point in space
     *
     * \param x The x coordinate of the point
     * \param y The y coordinate of the point
     * \param z The z coordinate of the point
     * \param size The desired mesh size at the point.
     * \param level The number of elements between the boundary and
     * this point.
     *
     * \return \a true if a valid size is returned, \a false otherwise.
     */
  virtual bool size_at_point(double x, double y, double z, double& size, 
                             int level = -1);

    /** \brief Get the desired anisotropic mesh size at a point in space in 
     * the given direction
     *
     * \param loc_x The x coordinate of the location
     * \param loc_y The y coordinate of the location
     * \param loc_z The z coordinate of the location
     * \param vec_x The x component of the vector to be stretched, 
     * returned modified
     * \param vec_y The y component of the vector to be stretched, 
     * returned modified
     * \param vec_z The z component of the vector to be stretched, 
     * returned modified
     *
     * \return \a true if a valid size is returned, \a false otherwise.
     */
  virtual bool stretch_vector(double loc_x, double loc_y, double loc_z,
                              double& vec_x, double& vec_y, double& vec_z);

    /** \brief Get the size tensor at a point in space
     *
     * \param x The x coordinate of the point
     * \param y The y coordinate of the point
     * \param z The z coordinate of the point
     * \param size The size tensor at the point, xx, yy, zz, xy, yz, xz. 
     * Pointer only changed if a NULL is passed in.
     *
     * \return \a true if a valid size is returned, \a false otherwise.
     */
  virtual bool tensor_at_point(double x, double y, double z, 
                               double*& size);

    /** \brief Get the desired mesh size at a point on the surface
     *
     * \param u The u parametric coordinate of the point
     * \param v The v parametric coordinate of the point
     * \param size The desired mesh size at the point.
     * \param level The number of elements between the boundary and
     * this point.
     *
     * \return \a true if a valid size is returned, \a false otherwise.
     */
  virtual bool size_at_point(double u, double v, double& size, 
                             int level = -1);

    /** \brief Get the desired anisotropic mesh size at a point in space in 
     * the given direction
     *
     * \param loc_u The u parametric coordinate of the point in the 
     * sizing field
     * \param loc_v The v parametric coordinate of the point in the 
     * sizing field
     * \param vec_u The u component of the direction to be stretched, 
     * returned modified
     * \param vec_v The v component of the direction to be stretched, 
     * returned modified
     *
     * \return \a true if a valid size is returned, \a false otherwise.
     */
  virtual bool stretch_vector(double loc_u, double loc_v,
			      double& vec_u, double& vec_v);

    /** \brief Get the size tensor at a point in space
     *
     * \param u The u parametric coordinate of the point
     * \param v The v parametric coordinate of the point
     * \param size The size tensor at the point.  Pointer only changed
     * if a NULL is passed in.
     *
     * \return \a true if a valid size is returned, \a false otherwise.
     */
  virtual bool tensor_at_point(double u, double v, 
                               double*& size);

  /** \brief Get the desired mesh size at a point on a curve
     *
     * \param u The u parametric coordinate of the point
     * \param size The desired mesh size at the point.
     * \param level The number of elements between the boundary and
     * this point.
     *
     * \return \a true if a valid size is returned, \a false otherwise.
     */
  virtual bool size_at_point(double u, double& size, 
                             int level = -1);

    /** \brief Get whether the sizing function is anisotropic or not.
     * \return \a true if sizing function is anisotropic
     */
  virtual bool is_anisotropic();

    /** \brief Get the mesh size set on this sizeeval
     * \return Mesh size
     */
  double get_size();
  
    /** \brief Set the mesh size set on this sizeeval
     * \param mesh_size Mesh size being set
     */
  void set_size(double mesh_size);
  
private:
  double meshSize;
};


inline double CAMALSizeEval::get_size() 
{
  return meshSize;
}

inline void CAMALSizeEval::set_size(double mesh_size)
{
  meshSize = mesh_size;
}
  
} // namespace MeshKit

#endif
