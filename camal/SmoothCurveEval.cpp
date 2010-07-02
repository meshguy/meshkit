/*
 * SmoothCurveEval.cpp
 *
 *  Created on: Jun 9, 2010
 *      Author: iulian
 */

#include "SmoothCurveEval.hpp"
#include "CubitVector.hpp"
#include "SmoothMeshEval.hpp"

// most of these functions will simply pass the call to the ref edge
// (delegate)
// it is most likely what CUBIT does with CAMAL
SmoothCurveEval::SmoothCurveEval(RefEdge * edge, SmoothMeshEval * smoothFaceEval, int loop)
{
	_ref_edge = edge;
	_smoothFaceEval = smoothFaceEval;
	_loopIndex = loop;
}

SmoothCurveEval::~SmoothCurveEval() {
	// TODO Auto-generated destructor stub
}

double SmoothCurveEval::arc_length()
{
	assert(_ref_edge);
	// return the length of the loop i
	double len = _smoothFaceEval->get_length_loop( _loopIndex);
	double leng = _ref_edge->get_arc_length();
	return leng;
}

    //! \brief Get the parametric status of the curve.
    //!
    //! \return \a true if curve is parametric, \a false otherwise.
bool SmoothCurveEval::is_parametric()
{
	return true;
}

    //! \brief Get the periodic status of the curve.
    //!
    //! \param period The period of the curve if periodic.
    //!
    //! \return \a true if curve is periodic, \a false otherwise.
bool SmoothCurveEval::is_periodic(double& period)
{
	assert(_ref_edge);
	return _ref_edge->is_periodic(   period);
}

    //! \brief Get the parameter range of the curve.
    //!
    //! \param u_start The beginning curve parameter
    //! \param u_end The ending curve parameter
    //!
    //! \note The numerical value of \a u_start may be greater
    //! than the numerical value of \a u_end.
void SmoothCurveEval::get_param_range(double& u_start, double& u_end)
{
    assert(_ref_edge);
    u_start = 0;
    u_end = 1.;
    u_start = _ref_edge->start_param();
    u_end = _ref_edge->end_param();
    return;
}

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
double SmoothCurveEval::u_from_arc_length(double u_root, double arc_length)
{
    assert(_ref_edge);
    return _ref_edge->u_from_arc_length( u_root, arc_length);
}


    //! \brief Evaluate the curve at a specified parameter value.
    //!
    //! \param u The parameter at which to evaluate the curve
    //! \param x The x coordinate of the evaluated point
    //! \param y The y coordinate of the evaluated point
    //! \param z The z coordinate of the evaluated point
bool SmoothCurveEval::position_from_u(double u,
                               double& x, double& y, double& z )
{
	// we will call the evaluator from smooth face evaluator
	assert (this->_smoothFaceEval);
	MBCartVect position;
	_smoothFaceEval->evaluate_loop_at_u(_loopIndex, u, position);
	assert(_ref_edge);
	CubitVector output_position;
	CubitStatus status = _ref_edge->position_from_u ( u, output_position);
	if (CUBIT_SUCCESS==status)
	{
		output_position.get_xyz( x, y, z ) ;
		return true;
	}
	return false;

}

    //! \brief Move a point near the curve to the closest point on the curve.
    //!
    //! \param x The x coordinate of the point
    //! \param y The y coordinate of the point
    //! \param z The z coordinate of the point
void SmoothCurveEval::move_to_curve(double& x, double& y, double& z)
{
	assert(_ref_edge);
	CubitVector a(x, y, z);
	_ref_edge->move_to_curve ( a );
	a.get_xyz( x, y, z ) ;
	return;
}

    //! Get the u parameter value on the curve closest to x,y,z
    //! and the point on the curve.
    //!
    //! \param x The x coordinate of the point
    //! \param y The y coordinate of the point
    //! \param z The z coordinate of the point
    //!
    //! \return The parametric coordinate u on the curve
double SmoothCurveEval::u_from_position(double x, double y, double z)
{
	assert(_ref_edge);
	CubitVector inp(x, y, z);
	return  _ref_edge->u_from_position (inp);
}

    //! \brief Get the starting point of the curve.
    //!
    //! \param x The x coordinate of the start point
    //! \param y The y coordinate of the start point
    //! \param z The z coordinate of the start point
void SmoothCurveEval::start_coordinates(double& x, double& y, double& z)
{
	assert(_ref_edge);
	CubitVector st = _ref_edge->start_coordinates();
	st.get_xyz( x, y, z ) ;
	return;
}

    //! \brief Get the ending point of the curve.
    //!
    //! \param x The x coordinate of the start point
    //! \param y The y coordinate of the start point
    //! \param z The z coordinate of the start point
void SmoothCurveEval::end_coordinates(double& x, double& y, double& z)
{
	assert(_ref_edge);
	CubitVector en = _ref_edge->end_coordinates();
	en.get_xyz( x, y, z ) ;
	return;
}
