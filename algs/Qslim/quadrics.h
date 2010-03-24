// -*- C++ -*-

//#include "AdjModel.h"
#include "std.h"
#include "Mat4.h"
#include "QslimOptions.h"
#include "MBEntityHandle.h"
#include "MBInterface.hpp"
extern MBInterface * mb;
extern QslimOptions opts;

extern Mat4 quadrix_vertex_constraint(const Vec3&);
extern Mat4 quadrix_vertex_constraint(MBEntityHandle vert);
extern Mat4 quadrix_plane_constraint(double a, double b, double c, double d);
extern Mat4 quadrix_plane_constraint(MBEntityHandle triangle);// (Face& T);
extern Mat4 quadrix_plane_constraint(const Vec3& n, double);
extern Mat4 quadrix_plane_constraint(const Vec3&, const Vec3&, const Vec3&);
extern double quadrix_evaluate_vertex(const Vec3& v, const Mat4& K);


extern bool is_border(MBEntityHandle mbedge);//(Edge *);
extern bool check_for_discontinuity(MBEntityHandle mbedge);//(Edge *);
extern Mat4 quadrix_discontinuity_constraint(MBEntityHandle mbedge, const Vec3&);//(Edge *, const Vec3&);
extern Mat4 quadrix_discontinuity_constraint(MBEntityHandle mbedge);//(Edge *);


extern bool quadrix_find_local_fit(const Mat4& K,
	    const Vec3& v1, const Vec3& v2,
	    Vec3& candidate);
extern bool quadrix_find_line_fit(const Mat4& Q,
				  const Vec3& v1, const Vec3& v2,
				  Vec3& candidate);
extern bool quadrix_find_best_fit(const Mat4& Q, Vec3& candidate);
extern double quadrix_pair_target(const Mat4& Q,
				MBEntityHandle v1, //Vertex *v1,
				MBEntityHandle v2, //Vertex *v2,
				Vec3& candidate);
