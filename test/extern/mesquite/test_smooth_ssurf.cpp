/*! \file test_mesquiteopt.cpp \test
 *
 * test_mesquite.cpp
 *
 */

#include "meshkit/MKCore.hpp"
#include "meshkit/ModelEnt.hpp"
#include "meshkit/MesquiteOpt.hpp"
#include "meshkit/iGeom.hpp"
#include "meshkit/iMesh.hpp"
#include "meshkit/Matrix.hpp"
#include "TestUtil.hpp"

using namespace MeshKit;

void test_smooth_spherical_surf();

ModelEnt* create_surface_mesh( moab::EntityHandle vertices[4][4] );
//ModelEnt* create_volume_mesh( moab::EntityHandle vertices[4][4][4] );
ModelEnt* get_model_ent( iBase_EntityHandle geom );

bool verbose = false;
bool write_vtk = false;
MKCore* core;
int main(int argc, char* argv[])
{
  for (int i = 1; i < argc; ++i) {
    if (!strcmp(argv[i],"-v"))
      verbose = true;
    else if (!strcmp(argv[i],"-w"))
      write_vtk = true;
    else if (!strcmp(argv[i],"-h"))
      std::cout << argv[0] << " [-v] [-w]" << std::endl;
    else
      std::cerr << "Invalid option: \"" << argv[i] << '"' << std::endl;
  }
  
  MKCore my_core;
  core = &my_core;

  int result = 0;
  result += RUN_TEST( test_smooth_spherical_surf );
 
  return result;
}

struct PointSort
{
  bool operator()(ModelEnt* pt1, ModelEnt* pt2)
  {
    double c1[3], c2[3];
    pt1->evaluate( 0, 0, 0, c1 );
    pt2->evaluate( 0, 0, 0, c2 );
    if (fabs(c1[2] - c2[2]) > 1e-6)
      return c1[2] < c2[2];
    else if (fabs(c1[1] - c2[1]) > 1e-6)
      return c1[1] < c2[1];
    else if (fabs(c1[0] - c2[0]) > 1e-6)
      return c1[0] < c2[0];
    else
      return false;
  }
};
    
// Generate meshed square face such that positions of the vertices
// passed back at each index are (values below are indices into 
// passed back 'vertices' array, not coordinates):
//
// (0,3), (1,3), (2,3), (3,3)  Y
// (0,2), (1,2), (2,2), (3,2)  ^
// (0,1), (1,1), (2,1), (3,1)  |
// (0,0), (1,0), (2,0), (3,0)  |
//                             +---->X
ModelEnt* create_surface_mesh( moab::EntityHandle vertices[4][4] )
{
    // create a 3x3x3 brick
  iGeom* geom = core->igeom_instance();
  iBase_EntityHandle brick;
  iGeom::Error err = geom->createBrick( 3, 3, 3, brick );
  CHECK_EQUAL( iBase_SUCCESS, err );
  
    // find ModelEnt for a surface parallel to the Z-plane and normal
    // in positive-Z direction
  ModelEnt* ent = get_model_ent( brick );
  MEntVector surfaces;
  ent->get_adjacencies( 2, surfaces );
  double normal[3];
  double point[3];
  size_t i;
  for (i = 0; i < surfaces.size(); ++i)
  {
    surfaces[i]->evaluate( 0, 0, 0, point, normal );
    if (fabs(normal[0]) < 1e-6 &&
        fabs(normal[1]) < 1e-6 &&
             normal[2]  > 1e-6)
      break;
  }
  CHECK( i != surfaces.size() );
  ModelEnt* surf = surfaces[i];
  
    // find sub-entities such that indices are:
    // 2-----3-----3  Y
    // |           |  ^
    // |           |  |
    // 1           2  |
    // |           |  +---->X
    // |           |
    // 0-----0-----1
    
  MEntVector curves, points;
  surf->get_adjacencies( 1, curves );
  surf->get_adjacencies( 0, points );
  CHECK_EQUAL( (size_t)4, curves.size() );
  CHECK_EQUAL( (size_t)4, points.size() );
  std::sort( points.begin(), points.end(), PointSort() );
  
  int ends[4][2] = { { 0, 1 }, {0, 2}, {1, 3}, {2, 3} };
  bool reversed[4];
  MEntVector endpts;
  for (i = 0; i < curves.size(); ++i) {
    ModelEnt* exp1 = points[ends[i][0]];
    ModelEnt* exp2 = points[ends[i][1]];
    size_t j;
    for (j = i; j < curves.size(); ++j) {
      endpts.clear();
      curves[j]->get_adjacencies( 0, endpts );
      CHECK_EQUAL( (size_t)2, endpts.size() );
      if (endpts[0] == exp1 && endpts[1] == exp2) {
        reversed[i] = false;
        break;
      }
      else if (endpts[0] == exp2 && endpts[1] == exp1) {
        reversed[i] = true;
        break;
      }
    }
    CHECK(j != 4);
    if (j != i)
      std::swap( curves[i], curves[j] );
  }
  
    // create mesh vertices
    
    // first get corner coords
  Vector<3> corners[4];
  for (i = 0; i < 4; ++i) 
    points[i]->evaluate( 0, 0, 0, corners[i].data() );
  
    // now create all vertices
  moab::ErrorCode rval;
  for (int u = 0; u < 4; ++u) {
    for (int v = 0; v < 4; ++v) {
      const double u3 = u/3.0;
      const double v3 = v/3.0;
      Vector<3> coords = (1-u3)*(1-v3)*corners[0] +
                         (  u3)*(1-v3)*corners[1] +
                         (1-u3)*(  v3)*corners[2] +
                         (  u3)*(  v3)*corners[3];
      rval = core->moab_instance()->create_vertex( coords.data(), vertices[u][v] );
      CHECK_EQUAL( moab::MB_SUCCESS, rval );
    }
  }
  
    // assign mesh vertices to geometric vertices
  moab::EntityHandle verts[] = { vertices[0][0], vertices[3][0],
                                 vertices[0][3], vertices[3][3] };
  for (i = 0; i < 4; ++i) {
      // check that we created vertices with correct coords
    Vector<3> coords, vcoords;
    points[i]->evaluate( 0, 0, 0, coords.data() );
    rval = core->moab_instance()->get_coords( &verts[i], 1, vcoords.data() );
    CHECK_EQUAL( moab::MB_SUCCESS, rval );
    CHECK( length(coords - vcoords) < 1e-6 );
      // assign to vertex
    points[i]->commit_mesh( &verts[i], 1, COMPLETE_MESH );
  }
  
    // mesh edges
  for (int i = 0; i < 4; ++i) {
      // get list of four vertices
    int u, v, u_stride, v_stride;
    if (i < 2) {
      u = v = 0;
      u_stride = 1-i;
      v_stride = i;
    }
    else {
      u_stride = i-2;
      v_stride = 3-i;
      v = u_stride;
      u = v_stride;
    }
    moab::EntityHandle verts[4];
    for (int j = 0; j < 4; ++j)
      verts[j] = vertices[3*u + j*u_stride][3*v + j*v_stride];
    if (reversed[i]) {
      std::swap( verts[0], verts[3] );
      std::swap( verts[1], verts[2] );
    }
    
    moab::EntityHandle curve_mesh[5]; // three edges and two interior vertices
    curve_mesh[1] = verts[1];
    curve_mesh[3] = verts[2];
    rval = core->moab_instance()->create_element( moab::MBEDGE, &verts[0], 2, curve_mesh[0] );
    CHECK_EQUAL( moab::MB_SUCCESS, rval );
    rval = core->moab_instance()->create_element( moab::MBEDGE, &verts[1], 2, curve_mesh[2] );
    CHECK_EQUAL( moab::MB_SUCCESS, rval );
    rval = core->moab_instance()->create_element( moab::MBEDGE, &verts[2], 2, curve_mesh[4] );
    CHECK_EQUAL( moab::MB_SUCCESS, rval );
    
    curves[i]->commit_mesh( curve_mesh, 5, COMPLETE_MESH );
  }
      
    // mesh surface
  moab::Range surf_mesh;
  for (int u = 0; u < 3; ++u) {
    for (int v = 0; v < 3; ++v) {
      moab::EntityHandle conn[4] = { vertices[u  ][v  ],
                                     vertices[u+1][v  ],
                                     vertices[u+1][v+1],
                                     vertices[u  ][v+1] };
      moab::EntityHandle quad;
      rval = core->moab_instance()->create_element( moab::MBQUAD, conn, 4, quad );
      CHECK_EQUAL( moab::MB_SUCCESS, rval );
      surf_mesh.insert(quad);
    }
  }
    // add interior vertices to range
  for (int u = 1; u < 3; ++u)
    for (int v = 1; v < 3; ++v)
      surf_mesh.insert( vertices[u][v] );
  surf->commit_mesh( surf_mesh, COMPLETE_MESH );
  
  return surf;
}

ModelEnt* get_model_ent( iBase_EntityHandle geom )
{
  core->populate_model_ents();
  iBase_EntityType type;
  iGeom::Error err = core->igeom_instance()->getEntType( geom, type );
  CHECK_EQUAL(iBase_SUCCESS,err);
  CHECK(type < 4);
  MEntVector ents;
  core->get_entities_by_dimension(type, ents);
  for (MEntVector::iterator i = ents.begin(); i != ents.end(); ++i)
    if ((*i)->geom_handle() == geom)
      return *i;
  return 0;
}

void check_surf_coords( Vector<3> coords[4][4], moab::EntityHandle verts[4][4], double epsilon )
{
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      Vector<3> new_coords;
      moab::ErrorCode rval = core->moab_instance()->get_coords( &verts[i][j], 1, new_coords.data() );
      CHECK_EQUAL( moab::MB_SUCCESS, rval );
      CHECK_REAL_EQUAL( 0.0, length(new_coords - coords[i][j]), epsilon );
    }
  }
}

void test_smooth_spherical_surf()
{
    // create sphere
  const double rad = std::sqrt(3.0);
  iGeom* geom = core->igeom_instance();
  iBase_EntityHandle vol_handle;
  iGeom::Error err = geom->createSphere( rad, vol_handle );
  CHECK_EQUAL( iBase_SUCCESS, err );
  ModelEnt* vol = get_model_ent( vol_handle );
  
    // check expected topology and get surface pointer
  MEntVector list;
  vol->get_adjacencies( 2, list );
  CHECK_EQUAL( (size_t)1, list.size() );
  ModelEnt* surf = list.front();
  //list.clear();
  //surf->get_adjacencies( 1, list );
  //CHECK( list.empty() );
  
    // define quad faces of a single hex circumscribed by the sphere
  const double xp = 0.5;
  const double yp = std::sqrt( 2 - xp * xp );
  const double coords[][3] = { { -1, -1, -1 },
                               {  1, -1, -1 },
                               { xp, yp, -1 }, // perturb one vertex
                               { -1,  1, -1 },
                               { -1, -1,  1 },
                               {  1, -1,  1 },
                               {  1,  1,  1 },
                               { -1,  1,  1 } };
  const int conn[][4] = { { 0, 1, 5, 4 },
                          { 1, 2, 6, 5 },
                          { 2, 3, 7, 6 },
                          { 3, 0, 4, 7 },
                          { 3, 2, 1, 0 },
                          { 4, 5, 6, 7 } };
  
    // create mesh entities
  moab::ErrorCode rval;
  moab::EntityHandle mesh[8+6], *verts=mesh, *quads=verts+8;
  for (int i = 0; i < 8; ++i) {
    rval = core->moab_instance()->create_vertex( coords[i], verts[i] );
    CHECK_EQUAL( moab::MB_SUCCESS, rval );
  }
  for (int i = 0; i < 6; ++i) {
    moab::EntityHandle mconn[4];
    for (int j = 0; j < 4; ++j)
      mconn[j] = verts[conn[i][j]];
    rval = core->moab_instance()->create_element( moab::MBQUAD, mconn, 4, quads[i] );
    CHECK_EQUAL( moab::MB_SUCCESS, rval );
  }
  surf->commit_mesh( mesh, 8+6, COMPLETE_MESH );
  moab::EntityHandle set = surf->mesh_handle();
  if (write_vtk)
    core->moab_instance()->write_file( "test_smooth_spherical_surf-init.vtk", "VTK", 0, &set, 1 );

  list.clear();
  list.push_back(surf);
  MesquiteOpt op( core, list );
  op.print_quality(verbose);
  op.smooth_with_fixed_boundary();
  if (verbose)
    std::cout << std::endl << "Fixed smooth of mesh of spherical surface" << std::endl << std::endl;
  op.setup_this();
  op.execute_this();
  if (write_vtk)
    core->moab_instance()->write_file( "test_smooth_spherical_surf-opt.vtk", "VTK", 0, &set, 1 );

    // check that all vertices lie on sphere
  for (int i = 0; i < 8; ++i) {
    Vector<3> c;
    rval = core->moab_instance()->get_coords( &verts[i], 1, c.data() );
    CHECK_REAL_EQUAL( rad, length(c), 1e-3 );
  }

    // check that all quad faces are squares
  for (int i = 0; i < 6; ++i) {
    Vector<3> c[4];
    moab::EntityHandle mconn[4];
    for (int j = 0; j < 4; ++j)
      mconn[j] = verts[conn[i][j]];
    rval = core->moab_instance()->get_coords( mconn, 4, c[0].data() );
    CHECK_EQUAL( moab::MB_SUCCESS, rval );
    for (int j = 0; j < 4; ++j) {
      Vector<3> e1 = c[(j+1)%4] - c[j];
      Vector<3> e2 = c[(j+2)%4] - c[(j+1)%4];
      CHECK_REAL_EQUAL( 2.0, length(e1),    5e-2 );
      CHECK_REAL_EQUAL( 4.0, length(e1*e2), 5e-2 );
    }
  }
}


