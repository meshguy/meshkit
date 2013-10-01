#ifndef READPOLYLINE_H
#define READPOLYLINE_H

// this should not be compiled separately, included only after some stuff
/*
#include <vector>
#include <fstream>
#include <iostream>
*/
// this is just a utility function used too many times in tests for
// mesh-based geometry, so I put it separately

int ReadPolyLineFromFile(const char * filename, double direction[3],
    std::vector<double> & points)
{
  // get the direction, and the polygon/ polyline points
  std::ifstream datafile(filename, std::ifstream::in);
  if (!datafile) {
    std::cout << "can't read polyline file\n";
    return 1;
  }
  //
  char temp[100];
  //double direction[3];// normalized
  double gridSize;
  datafile.getline(temp, 100);// first line

  // get direction and mesh size along polygon segments, from file
  sscanf(temp, " %lf %lf %lf %lf ", direction, direction + 1, direction + 2,
      &gridSize);

  //std::vector<double> xyz;
  while (!datafile.eof()) {
    datafile.getline(temp, 100);
    //int id = 0;
    double x, y, z;
    int nr = sscanf(temp, "%lf %lf %lf", &x, &y, &z);
    if (nr == 3) {
      points.push_back(x);
      points.push_back(y);
      points.push_back(z);
    }
  }
  return 0;
}
#endif // READPOLYLINE_H
