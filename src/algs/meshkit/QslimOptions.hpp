/*
 * QslimOptions.hpp
 *
 */

#ifndef QSLIMOPTIONS_H_
#define QSLIMOPTIONS_H_

#define PLACE_ENDPOINTS 0
#define PLACE_ENDORMID  1
#define PLACE_LINE      2
#define PLACE_OPTIMAL   3

#define OUTPUT_NONE         0x00
#define OUTPUT_CONTRACTIONS 0x01
#define OUTPUT_QUADRICS     0x02
#define OUTPUT_COST         0x04
#define OUTPUT_VERT_NOTES   0x08
#define OUTPUT_FACE_NOTES   0x10
#define OUTPUT_MODEL_DEFN   0x20
#define OUTPUT_ALL          0xFFFFFFFF

#include <fstream>

//===========================================================================//
  /*!
   * \class QslimOptions
   * \brief Place holder for decimation options

   *  Options for decimation can be passed with the test driver.
   *
   *   - -o \<file\>      Output final model to given file, using moab output conventions
   *   - -s \<count\>      Set the target number of faces.
   *   - -i \<timings\>   Intervals for timing reports.
   *   - -e \<thresh\> Set the maximum error tolerance.
   *   - -t \<t\>        Set pair selection tolerance.
   *   - -Q[pv]        Select what constraint quadrics to use [default=p].
   *   - -On       Optimal placement policy.
   *   0=endpoints, 1=endormid, 2=line, 3=optimal [default]
   *   - -B \<weight\> Use boundary preservation planes with given weight.
   *   - -b       preserve boundary (do not use with -B option)
   *   - -m            Preserve mesh quality.
   *   - -a            Enable area weighting.
   *   - -p          Height fields positivity. Used for height fields, assume
   *   triangles are originally positively oriented.
   *   - -d          Use delayed deletion, as opposed to merging
   *   - -c           keep costs in a (sparse!!!!) tag
   *   - -l           log file
   *   - -L           output level in the log file

 */
class QslimOptions {
public:
	QslimOptions();
	virtual ~QslimOptions();
	//! decimation stops when number of triangle is less than face_target; default 0
	int face_target;

	//! decimation can stop if the error is bigger than the given tolerance; default: HUGE (very large, so it is not stopping because of this)
	double error_tolerance;

	//! will use plane constraint;  default 1 (true)
	int     will_use_plane_constraint;

	//! will use vertex constraint; default : 0 (false)
	int     will_use_vertex_constraint;

	//! to preserve boundary edges; default 0 (false)
	int     will_preserve_boundaries;

	//! cost of creating inverted elements is artificially increased, so it should not happen;  default: false
	int     will_preserve_mesh_quality;

	//! triggered by -B option; it will add an additional cost when boundary edges are collapsed used in conjunction with the boundary constraint weight: default false
	int     will_constrain_boundaries;

	//! see the will_constrain_boundaries option
	double boundary_constraint_weight;

	//! used if areas are planar ; default 0 (no)
	int     will_weight_by_area;

	//! for terrain-type data, do not allow reverse of the normal in z direction; default 0 (no check); triggered by -p option
	int     height_fields;

	//! performance decimation intervals; default 10 intervals
	int     timingIntervals;

	//! create a tag with the error cost at each vertex in the final mesh (default 0, no)
	int     plotCost;

	//! default false. If used, it will delay deletion of edges and triangles until the end  of decimation
	int useDelayedDeletion;

	//! default 3, optimal placement policy, for the vertex when an edge is collapsed; other options are 0=end points, 1=end points or middle point, 2 = on the line
	int placement_policy;

	//! create the pairs based on this minimum distance (if the nodes are not connected); by default, the pairs are created only between connected vertices
	double pair_selection_tolerance;

	//! output a debug file
	std::ostream * logfile;

	//! debug level in output
	int selected_output;
};

#endif /* QSLIMOPTIONS_H_ */
