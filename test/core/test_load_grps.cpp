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

void test_load_groups();

int main(int argc, char **argv)
{
    mk = new MKCore();
    int num_fail = 0; 
    num_fail += RUN_TEST(test_load_groups);
    return num_fail;
}

void test_load_groups()
{
  std::cout << TestDir << std::endl;
  std::string filename = TestDir + "/" + DEFAULT_TEST_FILE;
  // load the test file
  mk->load_geometry(filename.c_str()); 
}
