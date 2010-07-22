
#ifndef CamalPaveDriver_H
#define CamalPaveDriver_H

#include <vector>

#include "MBInterface.hpp"
#include <map>

//#include "iMesh.h"

class SmoothFaceEval;
class SmoothCurveEval;
class SmoothVertex;

class CamalPaveDriver
{
public:
  CamalPaveDriver( MBInterface * mb, MBEntityHandle set,
        MBInterface * out, double angle);
  
  bool remesh(double mesh_size, int mesh_intervals,
                     const bool force_intervals);
  
  virtual ~CamalPaveDriver();
private :

  bool initializeSmoothing();

  bool prepareCGMEvaluator();

  void mesh_vertices();

  MBErrorCode find_face_loops();

  void establish_mesh_curve_count();

  void mesh_curves();

  void mesh_surfaces();

  MBInterface * _mb;
  MBInterface * _mbo;
  MBEntityHandle _set;
  double _angle;
  double _mesh_size;
  int _mesh_intervals;
  bool _force_intervals;

  SmoothFaceEval ** _smthFace;
  SmoothCurveEval ** _smthCurve;
  SmoothVertex ** _smthVertex;
  
  std::map<MBEntityHandle, SmoothFaceEval*> _mapSurfaces;
  std::map<MBEntityHandle, SmoothCurveEval*> _mapCurves;
  std::map<MBEntityHandle, SmoothVertex*> _mapVertices;
};

#endif  // CamalPaveDriver_H

