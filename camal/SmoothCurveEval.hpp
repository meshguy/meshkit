/*
 * SmoothCurveEval.hpp
 *
 *  Created on: Jun 9, 2010
 *      Author: iulian
 */

#ifndef SMOOTHCURVEEVAL_HPP_
#define SMOOTHCURVEEVAL_HPP_

#include "CMLCurveEval.hpp"
#include "RefEdge.hpp"

class SmoothMeshEval;
// evaluator for Camal Edge Mesher
// this is really copying what Cubit is doing

class SmoothCurveEval : public CMLCurveEval
{
public:
	SmoothCurveEval(RefEdge * edge, SmoothMeshEval * smoothFaceEval, int loopIndex);
	virtual ~SmoothCurveEval();

	virtual double arc_length() ;

	    //! \brief Get the parametric status of the curve.
	    //!
	    //! \return \a true if curve is parametric, \a false otherwise.
	  virtual bool is_parametric() ;

	    //! \brief Get the periodic status of the curve.
	    //!
	    //! \param period The period of the curve if periodic.
	    //!
	    //! \return \a true if curve is periodic, \a false otherwise.
	  virtual bool is_periodic(double& period) ;

	    //! \brief Get the parameter range of the curve.
	    //!
	    //! \param u_start The beginning curve parameter
	    //! \param u_end The ending curve parameter
	    //!
	    //! \note The numerical value of \a u_start may be greater
	    //! than the numerical value of \a u_end.
	  virtual void get_param_range(double& u_start, double& u_end) ;

	    //! Compute the parameter value at a specified distance along the curve.
	    //!
	    //! \param u_root The start parameter from which to compute the distance
	    //! along the curve.
	    //! \param arc_length The distance to move along the curve.
	    //!
	    //! \note For positive values of \a arc_length the distance will be
	    //! computed in the direction of increasing parameter value along the
	    //! curve.  For negative values of \a arc_length the distance will be
	    //! computed in the direction of decreasing parameter value along the
	    //! curve.
	    //!
	    //! \return The parametric coordinate u along the curve
	  virtual double u_from_arc_length(double u_root, double arc_length);


	    //! \brief Evaluate the curve at a specified parameter value.
	    //!
	    //! \param u The parameter at which to evaluate the curve
	    //! \param x The x coordinate of the evaluated point
	    //! \param y The y coordinate of the evaluated point
	    //! \param z The z coordinate of the evaluated point
	  virtual bool position_from_u(double u,
	                               double& x, double& y, double& z ) ;

	    //! \brief Move a point near the curve to the closest point on the curve.
	    //!
	    //! \param x The x coordinate of the point
	    //! \param y The y coordinate of the point
	    //! \param z The z coordinate of the point
	  virtual void move_to_curve(double& x, double& y, double& z);

	    //! Get the u parameter value on the curve closest to x,y,z
	    //! and the point on the curve.
	    //!
	    //! \param x The x coordinate of the point
	    //! \param y The y coordinate of the point
	    //! \param z The z coordinate of the point
	    //!
	    //! \return The parametric coordinate u on the curve
	  virtual double u_from_position(double x, double y, double z) ;

	    //! \brief Get the starting point of the curve.
	    //!
	    //! \param x The x coordinate of the start point
	    //! \param y The y coordinate of the start point
	    //! \param z The z coordinate of the start point
	  virtual void start_coordinates(double& x, double& y, double& z) ;

	    //! \brief Get the ending point of the curve.
	    //!
	    //! \param x The x coordinate of the start point
	    //! \param y The y coordinate of the start point
	    //! \param z The z coordinate of the start point
	  virtual void end_coordinates(double& x, double& y, double& z) ;

private:
	  RefEdge * _ref_edge;
	  SmoothMeshEval * _smoothFaceEval; // just store the face pointer here,
	  int _loopIndex;
};

#endif /* SMOOTHCURVEEVAL_HPP_ */
