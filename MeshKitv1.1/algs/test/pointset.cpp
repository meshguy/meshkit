#include <iostream>

#include <fstream>
#include <vector>
#define STRINGIFY_(X) #X
#define STRINGIFY(X) STRINGIFY_(X)

#define DEFAULT_TEST_FILE STRINGIFY(SRCDIR) "/gridData.txt"
#define DEFAULT_OUT_FILE "output.h5m"

#include "moab/Core.hpp"

int main(int argc, char* argv[])
{
  // check command line arg
  std::string filename_data;

  std::string output_h5m_file;

  bool save = false;
  if (argc == 3) {
    filename_data = argv[1];
    output_h5m_file = argv[2];
    save = true;

  } else {
    printf("Usage: %s <vertex_filename> <output> \n", argv[0]);
    if (argc != 1)
      return 1;
    printf("No files specified.  Defaulting to: %s %s\n", DEFAULT_TEST_FILE,
        DEFAULT_OUT_FILE);
    filename_data = DEFAULT_TEST_FILE;
    output_h5m_file = DEFAULT_OUT_FILE;
  }

  // read the file with the vertex data
  std::ifstream datafile(filename_data.c_str(), std::ifstream::in);
  if (!datafile) {
    std::cout << "can't read file\n";
    return 1;
  }

  char temp[100];
  datafile.getline(temp, 100);// first line has no data in it

  std::vector<double> coords;
  while (!datafile.eof()) {
    datafile.getline(temp, 100);
    double x, y, z;
    int nr = sscanf(temp, "%lf %lf %lf", &x, &y, &z);

    if (3 == nr) {
      coords.push_back(x);
      coords.push_back(y);
      coords.push_back(z);
    }
  }
  datafile.close();
  int numNodes = coords.size() / 3;
  std::cout << "num nodes " << numNodes << "\n";

  moab::Core mb;
  moab::Range verts;
  moab::ErrorCode rval = mb.create_vertices(&(coords[0]), numNodes, verts);
  if (rval != moab::MB_SUCCESS)
    return 1;

  // put all verts in an entity set of dim 2
  moab::EntityHandle newSet;
  rval = mb.create_meshset(moab::MESHSET_SET, newSet);
  if (rval != moab::MB_SUCCESS)
    return 1;

  rval = mb.add_entities(newSet, verts);
  if (rval != moab::MB_SUCCESS)
    return 1;

  //give geom dimension 2
  moab::Tag gtag;
  rval = mb.tag_get_handle("GEOM_DIMENSION", 1, moab::MB_TYPE_INTEGER, gtag,
      moab::MB_TAG_SPARSE|moab::MB_TAG_CREAT);
  if (rval != moab::MB_SUCCESS)
    return 1;
  int dimension =2;
  rval = mb.tag_set_data(gtag, &newSet, 1, &dimension);
  if (rval != moab::MB_SUCCESS)
    return 1;

  rval = mb.write_mesh(output_h5m_file.c_str());
  if (rval != moab::MB_SUCCESS)
      return 1;

  if (!save)
    remove(output_h5m_file.c_str());
  return 0;
}
