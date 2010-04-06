// $Id: quadrics.cxx,v 1.5 1997/10/01 14:07:28 garland Exp $

//#include "qslim.h"
#include "quadrics.h"
#include "3D.h"
#include "primitives.h"

////////////////////////////////////////////////////////////////////////
//
// Primitive quadric construction and evaluation routines
//

//
// Construct a quadric to evaluate the squared distance of any point
// to the given point v.  Naturally, the iso-surfaces are just spheres
// centered at v.
//
Mat4 quadrix_vertex_constraint(const Vec3& v)
{
    Mat4 L(Mat4::identity);

    L(0,3) = -v[0];
    L(1,3) = -v[1];
    L(2,3) = -v[2];
    L(3,3) = v*v;

    L(3,0) = L(0,3);
    L(3,1) = L(1,3);
    L(3,2) = L(2,3);

    return L;
}
Mat4 quadrix_vertex_constraint(MBEntityHandle vert)
{
    Mat4 L(Mat4::identity);
    double v[3];
    mb->get_coords(&vert, 1 , v);
    L(0,3) = -v[0];
    L(1,3) = -v[1];
    L(2,3) = -v[2];
    L(3,3) = v[0]*v[0] + v[1]*v[1] + v[2]*v[2] ; //  v*v;

    L(3,0) = L(0,3);
    L(3,1) = L(1,3);
    L(3,2) = L(2,3);

    return L;
}
//
// Construct a quadric to evaluate the squared distance of any point
// to the given plane [ax+by+cz+d = 0].  This is the "fundamental error
// quadric" discussed in the paper.
//
Mat4 quadrix_plane_constraint(double a, double b, double c, double d)
{
    Mat4 K(Mat4::zero);

    K(0,0) = a*a;   K(0,1) = a*b;   K(0,2) = a*c;  K(0,3) = a*d;
    K(1,0) =K(0,1); K(1,1) = b*b;   K(1,2) = b*c;  K(1,3) = b*d;
    K(2,0) =K(0,2); K(2,1) =K(1,2); K(2,2) = c*c;  K(2,3) = c*d;
    K(3,0) =K(0,3); K(3,1) =K(1,3); K(3,2) =K(2,3);K(3,3) = d*d;

    return K;
}

//
// Define some other convenient ways for constructing these plane quadrics.
//
Mat4 quadrix_plane_constraint(const Vec3& n, double d)
{
    return quadrix_plane_constraint(n[X], n[Y], n[Z], d);
}

//Mat4 quadrix_plane_constraint(Face& T)
Mat4 quadrix_plane_constraint(MBEntityHandle triangle)
{
    // const Plane& p = T.plane();
    Plane p=trianglePlane(mb, triangle);
    double a,b,c,d;
    p.coeffs(&a, &b, &c, &d);

    return quadrix_plane_constraint(a, b, c, d);
}

Mat4 quadrix_plane_constraint(const Vec3& v1, const Vec3& v2, const Vec3& v3)
{
    Plane P(v1,v2,v3);
    double a,b,c,d;
    P.coeffs(&a, &b, &c, &d);

    return quadrix_plane_constraint(a, b, c, d);
}

double quadrix_evaluate_vertex(const Vec3& v, const Mat4& K)
{
    double x=v[X], y=v[Y], z=v[Z];

#ifndef VECTOR_COST_EVALUATION
    //
    // This is the fast way of computing (v^T Q v).
    // 
    return x*x*K(0,0) + 2*x*y*K(0,1) + 2*x*z*K(0,2) + 2*x*K(0,3)
	              + y*y*K(1,1)   + 2*y*z*K(1,2) + 2*y*K(1,3)
	                             + z*z*K(2,2)   + 2*z*K(2,3)
	                                            + K(3,3);
#else
    //
    // The equivalent thing using matrix/vector operations.
    // It's a lot clearer, but it's also slower.
    //
    Vec4 v2(x,y,z,1);
    return v2*(K*v2);
#endif
}



////////////////////////////////////////////////////////////////////////
//
// Routines for computing discontinuity constraints
//

//bool is_border(Edge *e )
bool is_border(MBEntityHandle eh )
{
   //  use nb of triangles connected to check if number of tri is 1
    std::vector<MBEntityHandle> adjTri;
    MBErrorCode rval = mb->get_adjacencies(&eh, 1, 2, false, adjTri,
                                      MBInterface::UNION);
    if (MB_SUCCESS==rval  && adjTri.size()==1)
         return true;
    return false;
}

bool check_for_discontinuity(MBEntityHandle eh) //Edge *e)
{
    return is_border(eh);
}

Mat4 quadrix_discontinuity_constraint( MBEntityHandle mbe
    /*Edge *edge*/, const Vec3& n)
{
    //Vec3& org = *edge->org();
    //Vec3& dest = *edge->dest();
    // to do: Vec3 origin of edge
	// is the orientation important for an edge?
	const MBEntityHandle * conn;
	int num_nodes;
	mb->get_connectivity(mbe, conn, num_nodes);
    Vec3 dest = getVec3FromMBVertex(mb, conn[1]); //
    Vec3 org = getVec3FromMBVertex(mb, conn[0]);
    Vec3 e = dest - org;

    Vec3 n2 = e ^ n;
    unitize(n2);

    double d = -n2 * org;
    return quadrix_plane_constraint(n2, d);
}


Mat4 quadrix_discontinuity_constraint(MBEntityHandle mbedge)// Edge *edge)
{
    Mat4 D(Mat4::zero);

/*
    face_buffer& faces = edge->faceUses();

    for(int i=0; i<faces.length(); i++)
	D += quadrix_discontinuity_constraint(edge,faces(i)->plane().normal());
*/
    // to do : to get the connected faces to the edge

	std::vector<MBEntityHandle> adjFaces;
	MBErrorCode rval = mb->get_adjacencies(&mbedge, 1, 2, false,  adjFaces,
								  MBInterface::UNION);

	for (int i=0; i<adjFaces.size(); i++)
	{
		MBEntityHandle F = adjFaces[i];
		Plane p = trianglePlane( mb, F);
		D += quadrix_discontinuity_constraint(mbedge, p.normal());
	}
    return D;
}



////////////////////////////////////////////////////////////////////////
//
// Routines for computing contraction target
//

bool quadrix_find_local_fit(const Mat4& K,
			    const Vec3& v1, const Vec3& v2,
			    Vec3& candidate)
{

    Vec3 v3 = (v1 + v2) / 2;

    bool try_midpoint = opts.placement_policy > PLACE_ENDPOINTS;

    double c1 = quadrix_evaluate_vertex(v1, K);
    double c2 = quadrix_evaluate_vertex(v2, K);
    double c3;
    if( try_midpoint ) c3 = quadrix_evaluate_vertex(v3, K);

    if( c1<c2 )
    {
	if( try_midpoint && c3<c1 )
	    candidate=v3;
	else
	    candidate=v1;
    }
    else
    {
	if( try_midpoint && c3<c2 )
	    candidate=v3;
	else
	    candidate=v2;
    }

    return true;
}

bool quadrix_find_line_fit(const Mat4& Q,
			   const Vec3& v1, const Vec3& v2,
			   Vec3& candidate)
{
    Vec3 d = v1-v2;

    Vec3 Qv2 = Q*v2;
    Vec3 Qd  = Q*d;

    double denom = 2*d*Qd;

    if( denom == 0.0 )
	return false;

    double a = (d*Qv2 + v2*Qd) / denom;

    if( a<0.0 ) a=0.0;
    if( a>1.0 ) a=1.0;


    candidate = a*d + v2;
    return true;
}

bool quadrix_find_best_fit(const Mat4& Q, Vec3& candidate)
{
    Mat4 K = Q;
    K(3,0) = K(3,1) = K(3,2) = 0.0;  K(3,3) = 1;


    Mat4 M;
    double det = K.inverse(M);
    if( FEQ(det, 0.0, 1e-12) )
	return false;


#ifdef SAFETY
    //
    // The homogeneous division SHOULDN'T be necessary.
    // But, when we're being SAFE, we do it anyway just in case.
    //
    candidate[X] = M(0,3)/M(3,3);
    candidate[Y] = M(1,3)/M(3,3);
    candidate[Z] = M(2,3)/M(3,3);
#else
    candidate[X] = M(0,3);
    candidate[Y] = M(1,3);
    candidate[Z] = M(2,3);
#endif

    return true;
}


double quadrix_pair_target(const Mat4& Q,
			 MBEntityHandle v1, //Vertex *v1,
			 MBEntityHandle v2, // Vertex *v2,
			 Vec3& candidate)
{
    int policy = opts.placement_policy;

    //
    // This analytic boundary preservation isn't really necessary.  The
    // boundary constraint quadrics are quite effective.  But, I've left it
    // in anyway.
    //
    Vec3 vec1 = getVec3FromMBVertex(mb, v1);
    Vec3 vec2 = getVec3FromMBVertex(mb, v2);
    if( opts.will_preserve_boundaries )
    {
    	int c1 = classifyVertex(mb, v1);
    	int c2 = classifyVertex(mb, v2);

    	// if both are on boundary, put a high penalty cost
    	if (c1>0 && c2>0)
    		return 1.e11;// greater than quality error
    	if( c1 > c2 )
    	{
    		candidate = vec1;
    		return quadrix_evaluate_vertex(candidate, Q);
    	}
    	else if( c2 > c1 )
		{
			//candidate = *v2;
			candidate = vec2;
			return quadrix_evaluate_vertex(candidate, Q);
		}
		else if( c1>0 && policy>PLACE_LINE )
			policy = PLACE_LINE;

    	if( policy == PLACE_OPTIMAL ) assert(c1==0 && c2==0);
    }

    switch( policy )
    {
    case PLACE_OPTIMAL:
    	if( quadrix_find_best_fit(Q, candidate) )
    		break;

    case PLACE_LINE:
    	if( quadrix_find_line_fit(Q, vec1, vec2, candidate) )
    		break;

    default:
    	quadrix_find_local_fit(Q, vec1, vec2, candidate);
		break;
    }

    return quadrix_evaluate_vertex(candidate, Q);
}
