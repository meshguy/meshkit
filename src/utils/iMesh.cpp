#include "meshkit/iMesh.hh"

// map from MB's entity type to TSTT's entity topology
const moab::EntityType iMesh::mb_topology_table[] =
{
    moab::MBVERTEX,
    moab::MBEDGE,
    moab::MBPOLYGON,
    moab::MBTRI,
    moab::MBQUAD,
    moab::MBPOLYHEDRON,
    moab::MBTET,
    moab::MBHEX,
    moab::MBPRISM,
    moab::MBPYRAMID,
    moab::MBMAXTYPE,
    moab::MBMAXTYPE,
};
