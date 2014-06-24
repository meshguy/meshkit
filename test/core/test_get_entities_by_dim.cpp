/** \file test_copymesh.cpp \test
 *
 * Test Loading of a file w/ Geom Group Sets
 *
 */

#include "meshkit/MKCore.hpp"
#include <stdlib.h>

using namespace MeshKit;

#include "TestUtil.hpp"

#define DEFAULT_TEST_FILE "cyl_grps.sat"

MKCore *mk;

void test_get_entities_by_dim();

int main(int argc, char **argv)
{
    mk = new MKCore();
    int num_fail = 0; 
    #ifdef HAVE_ACIS
    num_fail += RUN_TEST(test_get_entities_by_dim);
    #endif
    return num_fail;
}

void test_get_entities_by_dim()
{
  std::cout << TestDir << std::endl;
  std::string filename = TestDir + "/" + DEFAULT_TEST_FILE;
  // load the test file
  mk->load_geometry(filename.c_str()); 

  // check that we find our group
  MEntVector groups;
  mk->get_entities_by_dimension(4, groups);
  int num_of_grps = groups.size();
  CHECK_EQUAL(1, num_of_grps);

  // check that we get the right number of volumes
  MEntVector vols;
  mk->get_entities_by_dimension(3, vols);
  int num_of_vols = vols.size();
  CHECK_EQUAL(1, num_of_vols);

  // check surfaces
  MEntVector surfs;
  mk->get_entities_by_dimension(2, surfs);
  int num_of_surfs = surfs.size();
  CHECK_EQUAL(3, num_of_surfs);

  MEntVector curves;
  mk->get_entities_by_dimension(1, curves);
  int num_of_curves = curves.size();
  CHECK_EQUAL(2, num_of_curves);

  MEntVector verts;
  mk->get_entities_by_dimension(0, verts);
  int num_of_verts = verts.size();
  CHECK_EQUAL(2, num_of_verts);

}
