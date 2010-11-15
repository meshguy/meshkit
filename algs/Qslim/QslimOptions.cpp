/*
 * QslimOptions.cpp
 *
 *  Created on: Mar 10, 2010
 *      Author: iulian
 */

#include "QslimOptions.h"
// for HUGES
#include "math.h"

QslimOptions::QslimOptions() {
	// TODO Auto-generated constructor stub
	// default options:
	error_tolerance=0;
	face_target = 0;
	error_tolerance = HUGE;
	will_use_plane_constraint = true;
	will_use_vertex_constraint = false;

	will_preserve_boundaries = false;
	will_preserve_mesh_quality = false;
	will_constrain_boundaries = false;
	boundary_constraint_weight = 1.0;

	will_weight_by_area = false;
	height_fields = false;
	timingIntervals = 10;

	useDelayedDeletion = false;

	placement_policy = PLACE_OPTIMAL;

	pair_selection_tolerance = 0.0;
	logfile = 0;
	selected_output = 0;
	plotCost = 0;


}

QslimOptions::~QslimOptions() {
	// TODO Auto-generated destructor stub
}
