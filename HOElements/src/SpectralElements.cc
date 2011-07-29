#include "SpectralElements.h"


////////////////////////////////////////////////////////////////////////////////

void saveFace(int norder, vector<double> &xyz)
{
     static int findex = 0;

     ostringstream oss;
     oss << "/tmp/face";
     if (findex < 1000) oss << "0";
     if (findex < 100) oss << "0";
     oss << findex << ".vtk";
     string filename = oss.str();

     cout << " Filename " << filename << endl;

     ofstream ofile(filename.c_str(), ios::out);

     ofile << "# vtk DataFile Version 2.0 " << endl;
     ofile << "Cubit Face  " << endl;
     ofile << "ASCII  " << endl;
     ofile << "DATASET STRUCTURED_GRID  " << endl;
     ofile << "DIMENSIONS " << norder << " " << norder << " 1 " << endl;
     ofile << "POINTS " << norder * norder << " float " << endl;

     for (int i = 0; i < norder * norder; i++)
          ofile << xyz[3 * i + 0] << " " << xyz[3 * i + 1] << " " << xyz[3 * i + 2] << endl;

     findex++;
}

///////////////////////////////////////////////////////////////////////////////

SpectralElements::SpectralElements()
{
     subElements = 0;
     hasSubStructuredGrid = 0;
}

///////////////////////////////////////////////////////////////////////////////

int SpectralElements::canonical_quad_coords(const vector<Point2D> &uvCorners,
          const Point2D &qPoint, Point2D &xy)
{

     //****************************************************************************
     // Input :
     //    uvCorners: Four corner values on any Quad element in either clockwise
     //                or anti-clockwise order.
     //    qPoint   : A point within the Quad element, but no checking is done
     //                inside to ensure that the point is within the Quad Element.
     //                but may be done using CGAL.
     // Output:
     //     xy     : Canonical Coordinates of "qPoint".
     //
     // Return Value:
     //            0:  No Error:
     //            1:  Convergence not reached.
     //            2:  Error in x coordinate calculations.
     //            3:  Error in y coordinate calculations;
     //
     // Defination: A Canonical Quad Elements has corners values
     //             ( -1, -1) , (1, -1)  (1,1) and (-1,1)
     //
     // Programmed By:  Chaman Singh Verma
     //                 Argonne National Lab. Chicago.
     // Date : 29th July 2009
     //
     //
     // Suggestions: It could perform better with Newton Raphson method.
     //
     //****************************************************************************

     assert(uvCorners.size() == 4);

     double xi, eta, du, dv, dl, mindist, prev_dist, eps = 1.0E-15;
     vector<double> weights(4), u(4), v(4);

     for (int i = 0; i < 4; i++) {
          u[i] = uvCorners[i][0];
          v[i] = uvCorners[i][1];
     }
     double umin = *min_element(u.begin(), u.end());
     double umax = *max_element(u.begin(), u.end());
     double vmin = *min_element(v.begin(), v.end());
     double vmax = *max_element(v.begin(), v.end());

     double dx = 1.0E-06;
     if (fabs(umax - umin) > 1.0E-06) dx = 1.0 / fabs(umax - umin);

     double dy = 1.0E-06;
     if (fabs(vmax - vmin) > 1.0E-06) dy = 1.0 / fabs(vmax - vmin);

     mindist = numeric_limits<double>::max();
     prev_dist = mindist;
     xi = 0.0;
     eta = 0.0;

     int iter = 0, maxiter = 100;
     while (1) {
          bilinear_weights(xi, eta, weights);
          xy[0] = linear_interpolation(u, weights);
          xy[1] = linear_interpolation(v, weights);

          du = qPoint[0] - xy[0];
          dv = qPoint[1] - xy[1];
          dl = du * du + dv*dv;
          if (dl < eps * eps) break;

          mindist = min(mindist, dl);
          xi = xi + dx*du;
          eta = eta + dy*dv;

          if (iter++ == maxiter) {
               cout << "Warning: Searching not converged " << endl;
               return 1;
          }
     }

     if (xy[0] < -1.0 || xy[0] > 1.0) return 2;
     if (xy[1] < -1.0 || xy[1] > 1.0) return 3;

     return 0;

}
///////////////////////////////////////////////////////////////////////////////

void SpectralElements::init()
{
     int err;

     const char *tagname = "HO_POINTS";
     err = mesh->createTag(tagname, sizeof (struct HO_Points), iBase_BYTES, horder_tag);
     assert(!err);

     err = mesh->getTagHandle("GLOBAL_ID", idtag);
     assert(!err);

     meshRootSet = mesh->getRootSet();
}

///////////////////////////////////////////////////////////////////////////////

void SpectralElements::gauss_linear_nodes(int N, vector<double> &x)
{
     x.resize(N);
     double du = 2.0 / (double) (N - 1);

     for (int i = 0; i < N; i++)
          x[i] = -1.0 + 2.0 * du;

     // Just to ensure no rounding errors;
     x[0] = -1.0;
     x[N - 1] = +1.0;
}

///////////////////////////////////////////////////////////////////////////////

void SpectralElements::gauss_lobatto_nodes(int N, vector<double> &x)
{

     x.resize(N);

     for (int i = 0; i < N; i++)
          x[i] = cos(M_PI * i / (double) (N - 1));

     double P[N][N + 1];

     vector<double> xold(N);
     for (int i = 0; i < N; i++) xold[i] = 2.0;

     double eps = 1.0E-10;
     while (1) {
          double maxerror = 0.0;
          for (int i = 0; i < N; i++)
               maxerror = max(maxerror, fabs(x[i] - xold[i]));

          if (maxerror < eps) break;

          for (int i = 0; i < N; i++) xold[i] = x[i];

          for (int i = 0; i < N; i++) {
               P[i][1] = 1.0;
               P[i][2] = x[i];
          }

          for (int k = 2; k <= N - 1; k++) {
               for (int i = 0; i < N; i++) {
                    P[i][k + 1] = ((2 * k - 1) * x[i] * P[i][k] - (k - 1) * P[i][k - 1]) / k;
               }
          }

          for (int i = 0; i < N; i++) {
               x[i] = xold[i] - (x[i] * P[i][N] - P[i][N - 1]) / (N * P[i][N]);
          }
     }

     std::reverse(x.begin(), x.end());
}

///////////////////////////////////////////////////////////////////////////////

bool SpectralElements::hasStructuredMesh2D(iBase_EntityHandle gFace) const
{
     int err;
     iBase_EntitySetHandle meshSet;

     std::set<iBase_EntityHandle> boundNodes, domainNodes;

     vector<iBase_EntityHandle> gEdges;
     err = geom->getEntAdj(gFace, iBase_EDGE, gEdges);

     vector<iBase_EntityHandle> mEdges, edgenodes, facenodes, faceedges;

     iBase_TagHandle dim_tag;
     const char *tag2 = "GEOM_DIMENSION";
     err = mesh->getTagHandle(tag2, dim_tag);
     assert(!err);

     int geom_dim;
     for (size_t i = 0; i < gEdges.size(); i++) {
          err = relPair->getEntSetRelation(gEdges[i], 0, meshSet);
          assert(!err);

          err = mesh->getEntSetIntData(meshSet, dim_tag, geom_dim);
          assert(!err);

          if (geom_dim == 1) {
               mEdges.clear();

               err = mesh->getEntities(meshSet, iBase_EDGE, iMesh_ALL_TOPOLOGIES, mEdges);
               assert(!err);

               for (size_t j = 0; j < mEdges.size(); j++) {
                    edgenodes.clear();
                    err = mesh->getEntAdj(mEdges[j], iBase_VERTEX, edgenodes );
                    for (size_t k = 0; k < edgenodes.size(); k++)
                         boundNodes.insert(edgenodes[k]);
               }
          }
     }

     vector<iBase_EntityHandle> mFaces;
     err = relPair->getEntSetRelation(gFace, 0, meshSet);
     assert( !err );

     err = mesh->getEntities(meshSet, iBase_FACE, iMesh_ALL_TOPOLOGIES, mFaces);
     assert( !err );

     map<iBase_EntityHandle, set<iBase_EntityHandle> > relation00Map; // Vertex->Vertex Relations

     for (size_t i = 0; i < mFaces.size(); i++) {
          facenodes.clear();
          err = mesh->getEntAdj(mFaces[i], iBase_VERTEX, facenodes);
          if (facenodes.size() != 4) return 0;
          for (size_t j = 0; j < facenodes.size(); j++)
               domainNodes.insert(facenodes[j]);

          err = mesh->getEntAdj(mFaces[i], iBase_EDGE, faceedges);
          for (size_t j = 0; j < faceedges.size(); j++) {
               edgenodes.clear();
               mesh->getEntAdj(faceedges[j], iBase_VERTEX, edgenodes);
               relation00Map[edgenodes[0]].insert(edgenodes[1]);
               relation00Map[edgenodes[1]].insert(edgenodes[0]);
          }
     }

     std::vector<iBase_EntityHandle> internalNodes;
     std::set_difference(domainNodes.begin(), domainNodes.end(),
                         boundNodes.begin(), boundNodes.end(),
                         back_inserter(internalNodes));

     // All the internal nodes must have four neighboring faces.
     for (size_t i = 0; i < internalNodes.size(); i++)
          if (relation00Map[internalNodes[i]].size() != 4) return 0;

     return 1;
}

///////////////////////////////////////////////////////////////////////////////

bool SpectralElements::hasStructuredMesh3D(iBase_EntityHandle gRegion) const
{
     int err;
     iBase_EntitySetHandle meshSet;

     std::set<iBase_EntityHandle> boundNodes, domainNodes;

     vector<iBase_EntityHandle> gFaces;
     geom->getEntAdj(gRegion, iBase_FACE, gFaces);

     vector<iBase_EntityHandle> mFaces, edgenodes, facenodes, cellnodes, celledges;

     iBase_TagHandle dim_tag;
     const char *tag2 = "GEOM_DIMENSION";
     mesh->getTagHandle(tag2, dim_tag);
     assert(!err);

     int geom_dim;
     for (size_t i = 0; i < gFaces.size(); i++) {
          err = relPair->getEntSetRelation(gFaces[i], 0, meshSet);
          assert(!err);

          err = mesh->getEntSetIntData(meshSet, dim_tag, geom_dim);
          assert(!err);

          if (geom_dim == 2) {
               mFaces.clear();

               err = mesh->getEntities(meshSet, iBase_FACE, iMesh_ALL_TOPOLOGIES, mFaces);
               assert(!err);

               for (size_t j = 0; j < mFaces.size(); j++) {
                    facenodes.clear();
                    mesh->getEntAdj(mFaces[j], iBase_VERTEX, facenodes);
                    for (int k = 0; k < facenodes.size(); k++)
                         boundNodes.insert(facenodes[k]);
               }
          }
     }

     vector<iBase_EntityHandle> mCells;
     err = relPair->getEntSetRelation(gRegion, 0, meshSet);
     mesh->getEntities(meshSet, iBase_REGION, iMesh_ALL_TOPOLOGIES, mCells);

     map<iBase_EntityHandle, set<iBase_EntityHandle> > relation00Map; // Vertex->Vertex Relations

     for (int i = 0; i < mCells.size(); i++) {
          cellnodes.clear();
          mesh->getEntAdj(mCells[i], iBase_VERTEX, cellnodes);
          if (cellnodes.size() != 8) return 0;
          for (int j = 0; j < cellnodes.size(); j++)
               domainNodes.insert(cellnodes[j]);

          mesh->getEntAdj(mCells[i], iBase_EDGE, celledges);
          for (int j = 0; j < celledges.size(); j++) {
               edgenodes.clear();
               mesh->getEntAdj(celledges[j], iBase_VERTEX, edgenodes);
               relation00Map[edgenodes[0]].insert(edgenodes[1]);
               relation00Map[edgenodes[1]].insert(edgenodes[0]);
          }
     }

     std::vector<iBase_EntityHandle> internalNodes;
     std::set_difference(domainNodes.begin(), domainNodes.end(),
                         boundNodes.begin(), boundNodes.end(),
                         back_inserter(internalNodes));

     // All the internal nodes must have four neighboring faces.
     for (int i = 0; i < internalNodes.size(); i++)
          if (relation00Map[internalNodes[i]].size() != 6) return 0;

     for (int i = 0; i < gFaces.size(); i++)
          if( !hasStructuredMesh2D( gFaces[i] ) )   return 0;

     return 1;
}

///////////////////////////////////////////////////////////////////////////////

void SpectralElements::bilinear_weights(double xi, double eta, vector<double> &weight)
{
     weight.resize(4);

     assert(xi >= -1.0 && xi <= 1.0);
     assert(eta >= -1.0 && eta <= 1.0);

     double coeff = 1.0 / 4.0;

     weight[0] = coeff * (1.0 - xi)*(1.0 - eta);
     weight[1] = coeff * (1.0 + xi)*(1.0 - eta);
     weight[2] = coeff * (1.0 + xi)*(1.0 + eta);
     weight[3] = coeff * (1.0 - xi)*(1.0 + eta);
}

///////////////////////////////////////////////////////////////////////////////

void SpectralElements::trilinear_weights(double xi, double eta, double zeta, vector<double> &weight)
{
     weight.resize(8);

     assert(xi >= -1.0 && xi <= 1.0);
     assert(eta >= -1.0 && eta <= 1.0);
     assert(zeta >= -1.0 && zeta <= 1.0);

     double coeff = 1.0 / 8.0;

     weight[0] = coeff * (1.0 - xi)*(1.0 - eta)*(1.0 - zeta);
     weight[1] = coeff * (1.0 + xi)*(1.0 - eta)*(1.0 - zeta);
     weight[2] = coeff * (1.0 + xi)*(1.0 + eta)*(1.0 - zeta);
     weight[3] = coeff * (1.0 - xi)*(1.0 + eta)*(1.0 - zeta);

     weight[4] = coeff * (1.0 - xi)*(1.0 - eta)*(1.0 + zeta);
     weight[5] = coeff * (1.0 + xi)*(1.0 - eta)*(1.0 + zeta);
     weight[6] = coeff * (1.0 + xi)*(1.0 + eta)*(1.0 + zeta);
     weight[7] = coeff * (1.0 - xi)*(1.0 + eta)*(1.0 + zeta);
}
///////////////////////////////////////////////////////////////////////////////

double SpectralElements::linear_interpolation(const vector<double> &x,
          const vector<double> &w)
{
     int numNodes = x.size();
     assert(w.size() == numNodes);

     double sum = 0.0;
     for (int i = 0; i < numNodes; i++) sum += x[i] * w[i];

     return sum;
}

///////////////////////////////////////////////////////////////////////////////

void SpectralElements::linear_interpolation(iBase_EntityHandle mCellHandle,
          double xi, double eta, double zeta,
          double &x, double &y, double &z)
{
     int err;
     static vector<double> weight(8);

     trilinear_weights(xi, eta, zeta, weight);

     vector<iBase_EntityHandle> nodeHandles;
     mesh->getEntAdj(mCellHandle, iBase_VERTEX, nodeHandles);

     int numNodes = nodeHandles.size();

     double xsum = 0.0;
     double ysum = 0.0;
     double zsum = 0.0;
     for (int i = 0; i < numNodes; i++) {
          mesh->getVtxCoord(nodeHandles[i], x, y, z);
          xsum += weight[i] * x;
          ysum += weight[i] * y;
          zsum += weight[i] * z;
     }
     x = xsum;
     y = ysum;
     z = zsum;
}

///////////////////////////////////////////////////////////////////////////////

void SpectralElements::get_subedges(iBase_EntityHandle mEdgeHandle, vector<int> &connect)
{
     int err;
     vector<iBase_EntityHandle> edgenodes;
     mesh->getEntAdj(mEdgeHandle, iBase_VERTEX, edgenodes);

     connect.clear();

     int gid;
     mesh->getIntData(edgenodes[0], idtag, gid);
     connect.push_back(gid);

     vector<iBase_EntityHandle> nodesOnEdge = entityNodesMap[mEdgeHandle].nodes;

     for (int i = 0; i < nodesOnEdge.size(); i++) {
          mesh->getIntData(nodesOnEdge[i], idtag, gid);
          connect.push_back(gid);
     }

     mesh->getIntData(edgenodes[1], idtag, gid);
     connect.push_back(gid);
}

///////////////////////////////////////////////////////////////////////////////

void SpectralElements::linear_edge(iBase_EntityHandle mEdgeHandle, int numPoints)
{
     int err;
     iBase_EntityHandle newhandle;

     vector<iBase_EntityHandle> edgenodes;
     mesh->getEntAdj(mEdgeHandle, iBase_VERTEX, edgenodes);

     assert(edgenodes.size() == 2);

     double x0, x1, y0, y1, z0, z1, xm, ym, zm;
     mesh->getVtxCoord( edgenodes[0], x0, y0, z0 );
     mesh->getVtxCoord( edgenodes[1], x1, y1, z1 );

     entityNodesMap[mEdgeHandle].clear();

     vector<iBase_EntityHandle> arranged_nodes(numPoints);
     arranged_nodes[0] = edgenodes[0];

     double u;
     for (int i = 1; i < numPoints - 1; i++) {
          u = gllnodes[i];
          xm = 0.5 * (1.0 - u) * x0 + 0.5 * (1.0 + u) * x1;
          ym = 0.5 * (1.0 - u) * y0 + 0.5 * (1.0 + u) * y1;
          zm = 0.5 * (1.0 - u) * z0 + 0.5 * (1.0 + u) * z1;
          mesh->createVtx(xm, ym, zm, newhandle);
          arranged_nodes[i] = newhandle;
          entityNodesMap[mEdgeHandle].add(newhandle);
     }
     arranged_nodes[numPoints - 1] = edgenodes[1];

     int numHPoints = numPoints;

     HO_Points hopoints;
     hopoints.nodeHandles = new iBase_EntityHandle[numHPoints];
     hopoints.nx = numHPoints;

     for (int i = 0; i < numHPoints; i++)
          hopoints.nodeHandles[i] = arranged_nodes[i];

     err = mesh->setData(mEdgeHandle, horder_tag, (const char *) & hopoints);
     assert(!err);
}

///////////////////////////////////////////////////////////////////////////////

void SpectralElements::bilinear_face(iBase_EntityHandle faceHandle, int numPoints)
{
     int i, j, err;
     iBase_EntityHandle newhandle;

     vector<iBase_EntityHandle> facenodes;
     mesh->getEntAdj(faceHandle, iBase_VERTEX, facenodes);

     assert(facenodes.size() == 4);

     int nx = numPoints;
     int ny = numPoints;

     double u, v, xm, ym, zm;
     for (j = 1; j < ny - 1; j++) {
          v = gllnodes[j];
          for (i = 1; i < nx - 1; i++) {
               u = gllnodes[i];
               linear_interpolation(faceHandle, u, v, -1.0, xm, ym, zm);
               mesh->createVtx(xm, ym, zm, newhandle);
               entityNodesMap[faceHandle].add(newhandle);
          }
     }

     // Now arrange the nodes on the face...

     vector<iBase_EntityHandle> arranged_nodes, nodesOnEdge;

     vector<iBase_EntityHandle> faceedges, edgenodes;
     mesh->getEntAdj(faceHandle, iBase_EDGE, faceedges);
     int numEdges = faceedges.size();

     int offset, numHPoints = nx*ny;
     arranged_nodes.resize(numHPoints);

     for (i = 0; i < arranged_nodes.size(); i++) arranged_nodes[i] = 0;

     i = 0;
     j = 0;
     offset = j * nx + i;
     arranged_nodes[offset] = facenodes[0];

     i = nx - 1;
     j = 0;
     offset = j * nx + i;
     arranged_nodes[offset] = facenodes[1];

     i = nx - 1;
     j = ny - 1;
     offset = j * nx + i;
     arranged_nodes[offset] = facenodes[2];

     i = 0;
     j = ny - 1;
     offset = j * nx + i;
     arranged_nodes[offset] = facenodes[3];

     int side_no, sense;

     for (int iedge = 0; iedge < numEdges; iedge++) {
          mesh->getEntAdj(faceedges[iedge], iBase_VERTEX, edgenodes);

          GEdge::get_quad_edge_number(&facenodes[0], &edgenodes[0], side_no, sense);

          nodesOnEdge = entityNodesMap[ faceedges[iedge] ].nodes;

          switch (side_no) {
          case 0:
               // Edge 0-1
               j = 0;
               if (sense == +1) {
                    for (i = 0; i < nx - 2; i++) {
                         offset = j * nx + i + 1;
                         arranged_nodes[offset] = nodesOnEdge[i];
                    }
               } else {
                    for (i = 0; i < nx - 2; i++) {
                         offset = j * nx + i + 1;
                         arranged_nodes[offset] = nodesOnEdge[nx - i - 3];
                    }
               }
               break;
          case 1:
               // Edge 1-2
               i = nx - 1;
               if (sense == +1) {
                    for (j = 0; j < ny - 2; j++) {
                         offset = (j + 1) * nx + i;
                         arranged_nodes[offset] = nodesOnEdge[j];
                    }
               } else {
                    for (j = 0; j < ny - 2; j++) {
                         offset = (j + 1) * nx + i;
                         arranged_nodes[offset] = nodesOnEdge[ny - j - 3];
                    }
               }
               break;
          case 2:
               // Edge 2-3
               j = ny - 1;
               if (sense == +1) {
                    for (i = 0; i < nx - 2; i++) {
                         offset = j * nx + i + 1;
                         arranged_nodes[offset] = nodesOnEdge[i];
                    }
               } else {
                    for (i = 0; i < nx - 2; i++) {
                         offset = j * nx + i + 1;
                         arranged_nodes[offset] = nodesOnEdge[nx - i - 3];
                    }
               }
               break;
          case 3:
               // Edge 0-3
               i = 0;
               if (sense == +1) {
                    for (j = 0; j < ny - 2; j++) {
                         offset = (j + 1) * nx + i;
                         arranged_nodes[offset] = nodesOnEdge[j];
                    }
               } else {
                    for (j = 0; j < ny - 2; j++) {
                         offset = (j + 1) * nx + i;
                         arranged_nodes[offset] = nodesOnEdge[ny - j - 3];
                    }
               }
               break;
          }
     }

     vector<iBase_EntityHandle> nodesOnFace = entityNodesMap[faceHandle].nodes;

     int index = 0;
     for (int j = 1; j < ny - 1; j++) {
          for (int i = 1; i < nx - 1; i++) {
               offset = j * nx + i;
               arranged_nodes[offset] = nodesOnFace[index++];
          }
     }

     for (i = 0; i < arranged_nodes.size(); i++)
          assert(arranged_nodes[i]);

     HO_Points hopoints;
     hopoints.nodeHandles = new iBase_EntityHandle[numHPoints];
     hopoints.nx = nx;
     hopoints.ny = ny;

     for (i = 0; i < numHPoints; i++)
          hopoints.nodeHandles[i] = arranged_nodes[i];

     err = mesh->setData(faceHandle, horder_tag, (const char *) & hopoints );
     assert(!err);
}

///////////////////////////////////////////////////////////////////////////////

void SpectralElements::trilinear_cell(iBase_EntityHandle cellHandle, int numPoints)
{
     int err;
     iBase_EntityHandle newhandle;

     vector<iBase_EntityHandle> cellnodes;
     mesh->getEntAdj(cellHandle, iBase_VERTEX, cellnodes);

     assert(cellnodes.size() == 8);

     double u, v, w;

     int nx = numPoints;
     int ny = numPoints;
     int nz = numPoints;

     double xm, ym, zm;
     for (int k = 1; k < nz - 1; k++) {
          w = gllnodes[k];
          for (int j = 1; j < ny - 1; j++) {
               v = gllnodes[j];
               for (int i = 1; i < nx - 1; i++) {
                    u = gllnodes[i];
                    linear_interpolation(cellHandle, u, v, w, xm, ym, zm);
                    mesh->createVtx( xm, ym, zm, newhandle);
                    entityNodesMap[cellHandle].add(newhandle);
               }
          }
     }

     // Now Arranged the nodes in the cell

     int numHPoints = nx * ny*nz;

     vector<iBase_EntityHandle> arranged_nodes;

     arrange_subcell_nodes(cellHandle, arranged_nodes);

     assert(arranged_nodes.size() == numHPoints);

     HO_Points hopoints;
     hopoints.nx = nx;
     hopoints.ny = ny;
     hopoints.nz = nz;
     hopoints.nodeHandles = new iBase_EntityHandle[numHPoints];

     for (int i = 0; i < numHPoints; i++)
          hopoints.nodeHandles[i] = arranged_nodes[i];

     err = mesh->setData(cellHandle, horder_tag, (const char *) & hopoints);
     assert(!err);

}

///////////////////////////////////////////////////////////////////////////////

void SpectralElements::arrange_brick20_nodes(iBase_EntityHandle cellHandle, vector<iBase_EntityHandle> &nodes)
{

     int err, side_no, sense, offset;

     vector<iBase_EntityHandle> hexVertexHandles;
     mesh->getEntAdj(cellHandle, iBase_VERTEX, hexVertexHandles);
     assert(hexVertexHandles.size() == 8);

     nodes.resize(20);
     nodes[0] = hexVertexHandles[0];
     nodes[2] = hexVertexHandles[1];
     nodes[5] = hexVertexHandles[3];
     nodes[7] = hexVertexHandles[2];
     nodes[12] = hexVertexHandles[4];
     nodes[14] = hexVertexHandles[5];
     nodes[17] = hexVertexHandles[7];
     nodes[19] = hexVertexHandles[6];

     std::map<iBase_EntityHandle, int> handleMap;
     for (int i = 0; i < 8; i++)
          handleMap[hexVertexHandles[i]] = i;

     iBase_EntityHandle vhandle;

     //
     // Now search vertices on the edges;
     //
     vector<iBase_EntityHandle> edgeHandles, nodeHandles;
     mesh->getEntAdj(cellHandle, iBase_EDGE, edgeHandles);
     assert(edgeHandles.size() == 12);

     vector<int> eConnect(2), hexConnect(8);
     for (int i = 0; i < 8; i++) hexConnect[i] = i;

     for (int i = 0; i < edgeHandles.size(); i++) {
          nodeHandles.clear();
          mesh->getEntAdj(edgeHandles[i], iBase_VERTEX, nodeHandles);
          eConnect[0] = handleMap[nodeHandles[0]];
          eConnect[1] = handleMap[nodeHandles[1]];
          MBCN::SideNumber(MBHEX, &hexConnect[0], &eConnect[0], 2, 1, side_no, sense, offset);
//        mesh->getEHData(edgeHandles[i], h1order_tag, vhandle);
          switch (side_no) {
          case 0:
               nodes[1] = vhandle;
               break;
          case 1:
               nodes[4] = vhandle;
               break;
          case 2:
               nodes[6] = vhandle;
               break;
          case 3:
               nodes[3] = vhandle;
               break;
          case 4:
               nodes[8] = vhandle;
               break;
          case 5:
               nodes[9] = vhandle;
               break;
          case 6:
               nodes[11] = vhandle;
               break;
          case 7:
               nodes[10] = vhandle;
               break;
          case 8:
               nodes[13] = vhandle;
               break;
          case 9:
               nodes[16] = vhandle;
               break;
          case 10:
               nodes[18] = vhandle;
               break;
          case 11:
               nodes[15] = vhandle;
               break;
          default:
               cout << "Fatal Error: Invalid side number " << endl;
               exit(0);
          }
     }
}

///////////////////////////////////////////////////////////////////////////////

void SpectralElements::brick20()
{
     init();

     subCellDim[0] = 2;
     subCellDim[1] = 2;
     subCellDim[2] = 2;

     hasSubStructuredGrid = 0;

     int err, result;

     iBase_TagHandle h20_tag;

     const char *tagname = "HO_20NODES";
     mesh->createTag(tagname, 20, iBase_ENTITY_HANDLE, h20_tag);
     assert(!err);

     vector<iBase_EntityHandle> edgeHandles;
     mesh->getEntities(meshRootSet, iBase_EDGE, iMesh_ALL_TOPOLOGIES, edgeHandles);
     int numEdges = edgeHandles.size();
     for (int i = 0; i < numEdges; i++) linear_edge(edgeHandles[i], 3);

     vector<iBase_EntityHandle> cellHandles;
     mesh->getEntities(meshRootSet, iBase_REGION, iMesh_HEXAHEDRON, cellHandles);
     int numCells = cellHandles.size();

     /*
         vector<iBase_EntityHandle> nodes(20);
         for (int i = 0; i < numCells; i++)
         {
             arrange_brick20_nodes(cellHandles[i], nodes);
             err = mesh->setEHArrData(&cellHandles[i], 1, h20_tag, &nodes[0], 20);
         }
     */
}

////////////////////////////////////////////////////////////////////////////////

void SpectralElements::brick20(iMesh *m)
{
     mesh = m;
     hasGeometry = 0;

     brick20();
}

////////////////////////////////////////////////////////////////////////////////

void SpectralElements::brick20(iMesh *m, iGeom *g,
                               iRel *r, iRel::PairHandle *a)
{
     mesh = m;
     geom  = g;
     rel   = r;
     relPair = a;
     hasGeometry = 1;

     brick20();
}

////////////////////////////////////////////////////////////////////////////////

int SpectralElements::getIndex(int *dim, int i, int j, int k)
{
     return k * dim[1] * dim[2] + j * dim[0] + i;
}

////////////////////////////////////////////////////////////////////////////////

void SpectralElements::get_subcells(iBase_EntityHandle cellHandle, vector<int> &connect)
{
     int err;
     int nodedim[3];
     iBase_EntityHandle *nodeHandles;

     nodedim[0] = subCellDim[0] + 1;
     nodedim[1] = subCellDim[1] + 1;
     nodedim[2] = subCellDim[2] + 1;

     int numNodes = nodedim[0] * nodedim[1] * nodedim[2];
     int numCells = subCellDim[0] * subCellDim[1] * subCellDim[2];

     char *tag_val = NULL;
     err = mesh->getData(cellHandle, horder_tag, &tag_val);
     assert(!err);

     HO_Points *hopoints = (HO_Points *) tag_val;
     nodeHandles = hopoints->nodeHandles;

     connect.resize(8 * numCells);

     int lid, gid, index = 0;

     for (int k = 0; k < subCellDim[2]; k++)
          for (int j = 0; j < subCellDim[1]; j++)
               for (int i = 0; i < subCellDim[0]; i++) {

                    lid = getIndex(nodedim, i, j, k);
                    err = mesh->getIntData(nodeHandles[lid], idtag, gid);
                    assert(!err);
                    connect[index++] = gid;

                    lid = getIndex(nodedim, i + 1, j, k);
                    err = mesh->getIntData(nodeHandles[lid], idtag, gid);
                    assert(!err);
                    connect[index++] = gid;

                    lid = getIndex(nodedim, i, j + 1, k);
                    err = mesh->getIntData(nodeHandles[lid], idtag, gid);
                    assert(!err);
                    connect[index++] = gid;

                    lid = getIndex(nodedim, i + 1, j + 1, k);
                    err = mesh->getIntData(nodeHandles[lid], idtag, gid);
                    assert(!err);
                    connect[index++] = gid;

                    lid = getIndex(nodedim, i, j, k + 1);
                    err = mesh->getIntData(nodeHandles[lid], idtag, gid);
                    assert(!err);
                    connect[index++] = gid;

                    lid = getIndex(nodedim, i + 1, j, k + 1);
                    err = mesh->getIntData(nodeHandles[lid], idtag, gid);
                    assert(!err);
                    connect[index++] = gid;

                    lid = getIndex(nodedim, i, j + 1, k + 1);
                    err = mesh->getIntData(nodeHandles[lid], idtag, gid);
                    assert(!err);
                    connect[index++] = gid;

                    lid = getIndex(nodedim, i + 1, j + 1, k + 1);
                    err = mesh->getIntData(nodeHandles[lid], idtag, gid);
                    assert(!err);
                    connect[index++] = gid;

               }
}

////////////////////////////////////////////////////////////////////////////////

int SpectralElements::matching_node(const vector<Point3D> &points, const Point3D &qPoint)
{
     double mindist = std::numeric_limits<double>::max();
     int minindex = 0;

     for (int i = 0; i < points.size(); i++) {
          double dx = points[i][0] - qPoint[0];
          double dy = points[i][1] - qPoint[1];
          double dz = points[i][2] - qPoint[2];
          double dl = dx * dx + dy * dy + dz*dz;
          if (dl < mindist) {
               minindex = i;
               mindist = dl;
          }
     }
     mindist = sqrt(mindist);

     if (mindist > 1.0E-15)
          cout << "Warning: Maximum Matching error : " << mindist << endl;

     return minindex;
}

////////////////////////////////////////////////////////////////////////////////////////////////

void SpectralElements::arrange_subcell_nodes(iBase_EntityHandle cellHandle,
          vector<iBase_EntityHandle> &arranged_nodes)
{

     int err;

     vector<iBase_EntityHandle> cellnodes;
     mesh->getEntAdj(cellHandle, iBase_VERTEX, cellnodes);

     vector<iBase_EntityHandle> celledges;
     mesh->getEntAdj(cellHandle, iBase_EDGE, celledges);

     int i, j, k, nx, ny, nz;
     int side_no, sense, offset;

     vector<iBase_EntityHandle> edgenodes;
     for (size_t i = 0; i < celledges.size(); i++) {
          mesh->getEntAdj(celledges[i], iBase_VERTEX, edgenodes);

          GEdge::get_hex_edge_number(&cellnodes[0], &edgenodes[0], side_no, sense);

          switch (side_no) {
          case 0:
               nx = entityNodesMap[celledges[i]].nodes.size() + 2;
               break;
          case 3:
               ny = entityNodesMap[celledges[i]].nodes.size() + 2;
               break;
          case 4:
               nz = entityNodesMap[celledges[i]].nodes.size() + 2;
               break;
          }
     }

     arranged_nodes.resize(nx * ny * nz);
     for (size_t i = 0; i < arranged_nodes.size(); i++)
          arranged_nodes[i] = 0;

     i = 0;
     j = 0;
     k = 0;
     offset = k * nx * ny + j * nx + i;
     arranged_nodes[offset] = cellnodes[0];

     i = nx - 1;
     j = 0;
     k = 0;
     offset = k * nx * ny + j * nx + i;
     arranged_nodes[offset] = cellnodes[1];

     i = nx - 1;
     j = ny - 1;
     k = 0;
     offset = k * nx * ny + j * nx + i;
     arranged_nodes[offset] = cellnodes[2];

     i = 0;
     j = ny - 1;
     k = 0;
     offset = k * nx * ny + j * nx + i;
     arranged_nodes[offset] = cellnodes[3];

     i = 0;
     j = 0;
     k = nz - 1;
     offset = k * nx * ny + j * nx + i;
     arranged_nodes[offset] = cellnodes[4];

     i = nx - 1;
     j = 0;
     k = nz - 1;
     offset = k * nx * ny + j * nx + i;
     arranged_nodes[offset] = cellnodes[5];

     i = nx - 1;
     j = ny - 1;
     k = nz - 1;
     offset = k * nx * ny + j * nx + i;
     arranged_nodes[offset] = cellnodes[6];

     i = 0;
     j = ny - 1;
     k = nz - 1;
     offset = k * nx * ny + j * nx + i;
     arranged_nodes[offset] = cellnodes[7];

     for (int iedge = 0; iedge < celledges.size(); iedge++) {
          mesh->getEntAdj(celledges[iedge], iBase_VERTEX, edgenodes);
          GEdge::get_hex_edge_number(&cellnodes[0], &edgenodes[0], side_no, sense);

          vector<iBase_EntityHandle> nodesOnEdge = entityNodesMap[celledges[iedge]].nodes;
          switch (side_no) {
          case 0:
               // Edge 0-1
               j = 0;
               k = 0;
               if (sense == +1) {
                    for (i = 0; i < nx - 2; i++) {
                         offset = k * nx * ny + j * nx + i + 1;
                         arranged_nodes[offset] = nodesOnEdge[i];
                    }
               } else {
                    for (i = 0; i < nx - 2; i++) {
                         offset = k * nx * ny + j * nx + i + 1;
                         arranged_nodes[offset] = nodesOnEdge[nx - i - 3];
                    }
               }
               break;
          case 1:
               // Edge 1-2
               i = nx - 1;
               k = 0;
               if (sense == +1) {
                    for (j = 0; j < ny - 2; j++) {
                         offset = k * nx * ny + (j + 1) * nx + i;
                         arranged_nodes[offset] = nodesOnEdge[j];
                    }
               } else {
                    for (j = 0; j < ny - 2; j++) {
                         offset = k * nx * ny + (j + 1) * nx + i;
                         arranged_nodes[offset] = nodesOnEdge[ny - j - 3];
                    }
               }
               break;
          case 2:
               // Edge 2-3
               j = ny - 1;
               k = 0;
               if (sense == +1) {
                    for (i = 0; i < nx - 2; i++) {
                         offset = k * nx * ny + j * nx + i + 1;
                         arranged_nodes[offset] = nodesOnEdge[i];
                    }
               } else {
                    for (i = 0; i < nx - 2; i++) {
                         offset = k * nx * ny + j * nx + i + 1;
                         arranged_nodes[offset] = nodesOnEdge[nx - i - 3];
                    }
               }
               break;
          case 3:
               // Edge 0-3
               i = 0;
               k = 0;
               if (sense == +1) {
                    for (j = 0; j < ny - 2; j++) {
                         offset = k * nx * ny + (j + 1) * nx + i;
                         arranged_nodes[offset] = nodesOnEdge[j];
                    }
               } else {
                    for (j = 0; j < ny - 2; j++) {
                         offset = k * nx * ny + (j + 1) * nx + i;
                         arranged_nodes[offset] = nodesOnEdge[ny - j - 3];
                    }
               }
               break;
          case 4:
               // Edge 0-4
               i = 0;
               j = 0;
               if (sense == +1) {
                    for (int k = 0; k < nz - 2; k++) {
                         offset = (k + 1) * nx * ny + j * nx + i;
                         arranged_nodes[offset] = nodesOnEdge[k];
                    }
               } else {
                    for (k = 0; k < nz - 2; k++) {
                         offset = (k + 1) * nx * ny + j * nx + i;
                         arranged_nodes[offset] = nodesOnEdge[nz - k - 3];
                    }
               }
               break;
          case 5:
               // Edge 1-5
               i = nx - 1;
               j = 0;
               if (sense == +1) {
                    for (k = 0; k < nz - 2; k++) {
                         offset = (k + 1) * nx * ny + j * nx + i;
                         arranged_nodes[offset] = nodesOnEdge[k];
                    }
               } else {
                    for (k = 0; k < nz - 2; k++) {
                         offset = (k + 1) * nx * ny + j * nx + i;
                         arranged_nodes[offset] = nodesOnEdge[nz - k - 3];
                    }
               }
               break;
          case 6:
               // Edge 2-6
               i = nx - 1;
               j = ny - 1;
               if (sense == +1) {
                    for (k = 0; k < nz - 2; k++) {
                         offset = (k + 1) * nx * ny + j * nx + i;
                         arranged_nodes[offset] = nodesOnEdge[k];
                    }
               } else {
                    for (k = 0; k < nz - 2; k++) {
                         offset = (k + 1) * nx * ny + j * nx + i;
                         arranged_nodes[offset] = nodesOnEdge[nz - k - 3];
                    }
               }
               break;
          case 7:
               // Edge 3-7
               i = 0;
               j = ny - 1;
               if (sense == +1) {
                    for (k = 0; k < nz - 2; k++) {
                         offset = (k + 1) * nx * ny + j * nx + i;
                         arranged_nodes[offset] = nodesOnEdge[k];
                    }
               } else {
                    for (k = 0; k < nz - 2; k++) {
                         offset = (k + 1) * nx * ny + j * nx + i;
                         arranged_nodes[offset] = nodesOnEdge[nz - k - 3];
                    }
               }
               break;
          case 8:
               // Edge 4-5
               j = 0;
               k = nz - 1;
               if (sense == +1) {
                    for (i = 0; i < nx - 2; i++) {
                         offset = k * nx * ny + j * nx + i + 1;
                         arranged_nodes[offset] = nodesOnEdge[i];
                    }
               } else {
                    for (i = 0; i < nx - 2; i++) {
                         offset = k * nx * ny + j * nx + i + 1;
                         arranged_nodes[offset] = nodesOnEdge[nx - i - 3];
                    }
               }
               break;
          case 9:
               // Edge 5-6
               i = nx - 1;
               k = nz - 1;
               if (sense == +1) {
                    for (j = 0; j < ny - 2; j++) {
                         offset = k * nx * ny + (j + 1) * nx + i;
                         arranged_nodes[offset] = nodesOnEdge[j];
                    }
               } else {
                    for (j = 0; j < ny - 2; j++) {
                         offset = k * nx * ny + (j + 1) * nx + i;
                         arranged_nodes[offset] = nodesOnEdge[ny - j - 3];
                    }
               }
               break;
          case 10:
               // Edge 6-7
               j = ny - 1;
               k = nz - 1;
               if (sense == +1) {
                    for (i = 0; i < nx - 2; i++) {
                         offset = k * nx * ny + j * nx + i + 1;
                         arranged_nodes[offset] = nodesOnEdge[i];
                    }
               } else {
                    for (i = 0; i < nx - 2; i++) {
                         offset = k * nx * ny + j * nx + i + 1;
                         arranged_nodes[offset] = nodesOnEdge[nx - i - 3];
                    }
               }
               break;
          case 11:
               // Edge 4-7
               i = 0;
               k = nz - 1;
               if (sense == +1) {
                    for (j = 0; j < ny - 2; j++) {
                         offset = k * nx * ny + (j + 1) * nx + i;
                         arranged_nodes[offset] = nodesOnEdge[j];
                    }
               } else {
                    for (j = 0; j < ny - 2; j++) {
                         offset = k * nx * ny + (j + 1) * nx + i;
                         arranged_nodes[offset] = nodesOnEdge[ny - j - 3];
                    }
               }
               break;
          default:
               cout << "Fatal Error: Invalid Side number " << side_no << endl;
               exit(0);
          }
     }

     int pos;
     double x, y, z, xi, eta, zeta;

     vector<int> qnodes(4);
     vector<double> xcellCorner(8), ycellCorner(8), zcellCorner(8);
     vector<double> xfaceCorner(4), yfaceCorner(4), zfaceCorner(4);
     vector<double> weights(4);

     vector<iBase_EntityHandle> cellfaces, facenodes;
     mesh->getEntAdj(cellHandle, iBase_FACE, cellfaces);
     vector<iBase_EntityHandle> nodesOnFace;
     vector<Point3D> hopoints;
     Point3D qPoint;

     for (i = 0; i < 8; i++)
          mesh->getVtxCoord(cellnodes[i], xcellCorner[i], ycellCorner[i], zcellCorner[i]);

     for (int iface = 0; iface < cellfaces.size(); iface++) {
          mesh->getEntAdj(cellfaces[iface], iBase_VERTEX, facenodes);
          GFace::get_hex_face_number(cellnodes, facenodes, side_no, sense);

          nodesOnFace = entityNodesMap[cellfaces[iface]].nodes;
          int numHPoints = nodesOnFace.size();
          hopoints.resize(numHPoints);
          for (i = 0; i < numHPoints; i++)
               mesh->getVtxCoord(nodesOnFace[i], hopoints[i][0], hopoints[i][1], hopoints[i][2]);

          switch (side_no) {
          case 0:
               qnodes[0] = 0;
               qnodes[1] = 3;
               qnodes[2] = 7;
               qnodes[3] = 4;
               for (j = 0; j < 4; j++) {
                    xfaceCorner[j] = xcellCorner[ qnodes[j] ];
                    yfaceCorner[j] = ycellCorner[ qnodes[j] ];
                    zfaceCorner[j] = zcellCorner[ qnodes[j] ];
               }
               i = 0;
               for (k = 1; k < nz - 1; k++) {
                    zeta = gllnodes[k];
                    for (j = 1; j < ny - 1; j++) {
                         offset = k * nx * ny + j * nx + i;
                         eta = gllnodes[j];
                         bilinear_weights(eta, zeta, weights);
                         qPoint[0] = linear_interpolation(xfaceCorner, weights);
                         qPoint[1] = linear_interpolation(yfaceCorner, weights);
                         qPoint[2] = linear_interpolation(zfaceCorner, weights);
                         pos = matching_node(hopoints, qPoint);
                         if (arranged_nodes[offset]) cout << "Error: Wrong matching occured " << endl;
                         arranged_nodes[offset] = nodesOnFace[pos];

                    }
               }
               break;
          case 1:
               qnodes[0] = 1;
               qnodes[1] = 2;
               qnodes[2] = 6;
               qnodes[3] = 5;
               for (j = 0; j < 4; j++) {
                    xfaceCorner[j] = xcellCorner[ qnodes[j] ];
                    yfaceCorner[j] = ycellCorner[ qnodes[j] ];
                    zfaceCorner[j] = zcellCorner[ qnodes[j] ];
               }
               i = nx - 1;
               for (k = 1; k < nz - 1; k++) {
                    zeta = gllnodes[k];
                    for (j = 1; j < ny - 1; j++) {
                         eta = gllnodes[j];
                         offset = k * nx * ny + j * nx + i;
                         bilinear_weights(eta, zeta, weights);
                         qPoint[0] = linear_interpolation(xfaceCorner, weights);
                         qPoint[1] = linear_interpolation(yfaceCorner, weights);
                         qPoint[2] = linear_interpolation(zfaceCorner, weights);
                         pos = matching_node(hopoints, qPoint);
                         if (arranged_nodes[offset]) cout << "Error: Wrong matching occured " << endl;
                         arranged_nodes[offset] = nodesOnFace[pos];
                    }
               }
               break;
          case 2:
               qnodes[0] = 0;
               qnodes[1] = 1;
               qnodes[2] = 5;
               qnodes[3] = 4;
               for (j = 0; j < 4; j++) {
                    xfaceCorner[j] = xcellCorner[ qnodes[j] ];
                    yfaceCorner[j] = ycellCorner[ qnodes[j] ];
                    zfaceCorner[j] = zcellCorner[ qnodes[j] ];
               }
               j = 0;
               for (k = 1; k < nz - 1; k++) {
                    zeta = gllnodes[k];
                    for (i = 1; i < nx - 1; i++) {
                         xi = gllnodes[i];
                         offset = k * nx * ny + j * nx + i;
                         bilinear_weights(xi, zeta, weights);
                         qPoint[0] = linear_interpolation(xfaceCorner, weights);
                         qPoint[1] = linear_interpolation(yfaceCorner, weights);
                         qPoint[2] = linear_interpolation(zfaceCorner, weights);
                         pos = matching_node(hopoints, qPoint);
                         if (arranged_nodes[offset]) cout << "Error: Wrong matching occured " << endl;
                         arranged_nodes[offset] = nodesOnFace[pos];
                    }
               }
               break;
          case 3:
               qnodes[0] = 3;
               qnodes[1] = 2;
               qnodes[2] = 6;
               qnodes[3] = 7;
               for (j = 0; j < 4; j++) {
                    xfaceCorner[j] = xcellCorner[ qnodes[j] ];
                    yfaceCorner[j] = ycellCorner[ qnodes[j] ];
                    zfaceCorner[j] = zcellCorner[ qnodes[j] ];
               }
               j = ny - 1;
               for (k = 1; k < nz - 1; k++) {
                    zeta = gllnodes[k];
                    for (i = 1; i < nx - 1; i++) {
                         xi = gllnodes[i];
                         offset = k * nx * ny + j * nx + i;
                         bilinear_weights(xi, zeta, weights);
                         qPoint[0] = linear_interpolation(xfaceCorner, weights);
                         qPoint[1] = linear_interpolation(yfaceCorner, weights);
                         qPoint[2] = linear_interpolation(zfaceCorner, weights);
                         pos = matching_node(hopoints, qPoint);
                         if (arranged_nodes[offset]) cout << "Error: Wrong matching occured " << endl;
                         arranged_nodes[offset] = nodesOnFace[pos];
                    }
               }
               break;
          case 4:
               qnodes[0] = 0;
               qnodes[1] = 1;
               qnodes[2] = 2;
               qnodes[3] = 3;
               for (j = 0; j < 4; j++) {
                    xfaceCorner[j] = xcellCorner[ qnodes[j] ];
                    yfaceCorner[j] = ycellCorner[ qnodes[j] ];
                    zfaceCorner[j] = zcellCorner[ qnodes[j] ];
               }
               k = 0;
               for (j = 1; j < ny - 1; j++) {
                    eta = gllnodes[j];
                    for (i = 1; i < nx - 1; i++) {
                         xi = gllnodes[i];
                         offset = k * nx * ny + j * nx + i;
                         bilinear_weights(xi, eta, weights);
                         qPoint[0] = linear_interpolation(xfaceCorner, weights);
                         qPoint[1] = linear_interpolation(yfaceCorner, weights);
                         qPoint[2] = linear_interpolation(zfaceCorner, weights);
                         pos = matching_node(hopoints, qPoint);
                         if (arranged_nodes[offset]) cout << "Error: Wrong matching occured " << endl;
                         arranged_nodes[offset] = nodesOnFace[pos];
                    }
               }
               break;
          case 5:
               qnodes[0] = 4;
               qnodes[1] = 5;
               qnodes[2] = 6;
               qnodes[3] = 7;
               for (j = 0; j < 4; j++) {
                    xfaceCorner[j] = xcellCorner[ qnodes[j] ];
                    yfaceCorner[j] = ycellCorner[ qnodes[j] ];
                    zfaceCorner[j] = zcellCorner[ qnodes[j] ];
               }
               k = nz - 1;
               for (j = 1; j < ny - 1; j++) {
                    eta = gllnodes[j];
                    for (i = 1; i < nx - 1; i++) {
                         xi = gllnodes[i];
                         offset = k * nx * ny + j * nx + i;
                         bilinear_weights(xi, eta, weights);
                         qPoint[0] = linear_interpolation(xfaceCorner, weights);
                         qPoint[1] = linear_interpolation(yfaceCorner, weights);
                         qPoint[2] = linear_interpolation(zfaceCorner, weights);
                         pos = matching_node(hopoints, qPoint);
                         if (arranged_nodes[offset]) cout << "Error: Wrong matching occured " << endl;
                         arranged_nodes[offset] = nodesOnFace[pos];
                    }
               }
               break;
          }
     }

     vector<iBase_EntityHandle> nodesInCell;
     nodesInCell = entityNodesMap[cellHandle].nodes;
     int index = 0;
     for (int k = 1; k < nz - 1; k++) {
          for (int j = 1; j < ny - 1; j++) {
               for (int i = 1; i < nx - 1; i++) {
                    offset = k * nx * ny + j * nx + i;
                    arranged_nodes[offset] = nodesInCell[index++];
               }
          }
     }

     for (int i = 0; i < arranged_nodes.size(); i++)
          assert(arranged_nodes[i]);

}

////////////////////////////////////////////////////////////////////////////////

int SpectralElements::check_edge_discretization()
{

     int numnodes = gllnodes.size();

     vector<double> arclength_ratio(numnodes);

     double u;
     for (int i = 0; i < numnodes; i++) {
          u = 0.5 * (1.0 + gllnodes[i]);
          arclength_ratio[i] = u;
     }

     std::map<iBase_EntityHandle, EntityNodes> ::const_iterator miter, mstart, mend;
     mstart = entityNodesMap.begin();
     mend = entityNodesMap.end();

     vector<Point3D> vpoints;
     vector<iBase_EntityHandle> edgeNodes;

     int err;
     double x, y, z;

     iMesh::EntityType entity_type;

     for (miter = mstart; miter != mend; ++miter) {
          iBase_EntityHandle mhandle = miter->first;
          mesh->getEntType( mhandle, entity_type);
          if (entity_type == iBase_EDGE) {
               vector<iBase_EntityHandle> nodesOnEdge = entityNodesMap[mhandle].nodes;
               mesh->getEntAdj( mhandle, iBase_VERTEX, edgeNodes);
               numnodes = nodesOnEdge.size() + 2;

               vpoints.resize(numnodes);

               err = mesh->getVtxCoord(edgeNodes[0], x, y, z);
               assert(!err);
               vpoints[0][0] = x;
               vpoints[0][1] = y;
               vpoints[0][2] = z;

               err = mesh->getVtxCoord(edgeNodes[1], x, y, z);
               assert(!err);
               vpoints[numnodes - 1][0] = x;
               vpoints[numnodes - 1][1] = y;
               vpoints[numnodes - 1][2] = z;

               for (int i = 0; i < nodesOnEdge.size(); i++) {
                    err = mesh->getVtxCoord(nodesOnEdge[i], x, y, z);
                    assert(!err);
                    vpoints[i + 1][0] = x;
                    vpoints[i + 1][1] = y;
                    vpoints[i + 1][2] = z;
               }

               double len, lensum = 0.0;
               for (int i = 0; i < numnodes - 1; i++) {
                    double dx = vpoints[i + 1][0] - vpoints[i][0];
                    double dy = vpoints[i + 1][1] - vpoints[i][1];
                    double dz = vpoints[i + 1][2] - vpoints[i][2];
                    double dl = sqrt(dx * dx + dy * dy + dz * dz);
                    lensum += dl;
               }
               double arclength = lensum;

               lensum = 0.0;
               for (int i = 0; i < numnodes - 1; i++) {
                    double dx = vpoints[i + 1][0] - vpoints[i][0];
                    double dy = vpoints[i + 1][1] - vpoints[i][1];
                    double dz = vpoints[i + 1][2] - vpoints[i][2];
                    double dl = sqrt(dx * dx + dy * dy + dz * dz);
                    lensum += dl;
                    double ar = lensum / arclength;
                    if (fabs(ar - arclength_ratio[i + 1]) > 1.0E-06) {
                         cout << "Warning: Edge vertex is not properly distributed: " << i << endl;
                         cout << "Currnet Arc Position " << ar << " must be " << arclength_ratio[i + 1] << endl;
                    }
               }
          }
     }
     cout << " Testing done " << endl;
}
////////////////////////////////////////////////////////////////////////////////

void SpectralElements::generate(iMesh *m, int order)
{
     mesh = m;
     hasGeometry = 0;
     hasSubStructuredGrid = 1;

     subCellDim[0] = order - 1;
     subCellDim[1] = order - 1;
     subCellDim[2] = order - 1;

     gauss_lobatto_nodes(order, gllnodes);

     init();

     entityNodesMap.clear();

     vector<iBase_EntityHandle> edgeHandles;
     mesh->getEntities(meshRootSet, iBase_EDGE, iMesh_ALL_TOPOLOGIES, edgeHandles);
     int numEdges = edgeHandles.size();
     for (int i = 0; i < numEdges; i++) linear_edge(edgeHandles[i], order);

     vector<iBase_EntityHandle> faceHandles;
     mesh->getEntities(meshRootSet, iBase_FACE, iMesh_ALL_TOPOLOGIES, faceHandles);
     int numFaces = faceHandles.size();
     for (int i = 0; i < numFaces; i++) bilinear_face(faceHandles[i], order);

     vector<iBase_EntityHandle> cellHandles;
     mesh->getEntities(meshRootSet, iBase_REGION, iMesh_ALL_TOPOLOGIES, cellHandles);
     int numCells = cellHandles.size();
     for (int i = 0; i < numCells; i++) trilinear_cell(cellHandles[i], order);

}

////////////////////////////////////////////////////////////////////////////////

void SpectralElements::project_on_edges()
{
     int err, namelen;

     iBase_EntitySetHandle geomRootSet;
     geomRootSet = geom->getRootSet();

     vector<iBase_EntityHandle> gEdges, mEdges;
     geom->getEntities(geomRootSet, iBase_EDGE, gEdges);

     /*
         for (int i = 0; i < gEdges.size(); i++)
         {
             GEdge curredge( gEdges[i], mesh, geom, rel, assoc );
             curredge.projectHigherOrderNodes( gllnodes );
         }
     */
}

///////////////////////////////////////////////////////////////////////////////

void SpectralElements::project_on_quad_face(GFace &currface,
          iBase_EntityHandle &mFaceHandle)
{
     int err;

     int offset, nx, ny, numHPoints;

     vector<iBase_EntityHandle> mEdges, faceNodes;
     mesh->getEntAdj(mFaceHandle, iBase_EDGE, mEdges);
     mesh->getEntAdj(mFaceHandle, iBase_VERTEX, faceNodes );

     Point3D p3d, pCentroid;
     Point2D uv, uvCentroid;

     pCentroid[0] = 0.0;
     pCentroid[1] = 0.0;
     pCentroid[2] = 0.0;
     for (int i = 0; i < 4; i++) {
          mesh->getVtxCoord(faceNodes[i], p3d[0], p3d[1], p3d[2]);
          pCentroid[0] += p3d[0];
          pCentroid[1] += p3d[1];
          pCentroid[2] += p3d[2];
     }

     pCentroid[0] /= 4.0;
     pCentroid[1] /= 4.0;
     pCentroid[2] /= 4.0;
     uvCentroid = currface.getUVCoords(pCentroid);

     char *tag_val = NULL;
     err = mesh->getData(mFaceHandle, horder_tag, &tag_val);
     assert(!err);
     HO_Points *hopoints = (HO_Points *) tag_val;
     nx = hopoints->nx;
     ny = hopoints->ny;
     numHPoints = nx*ny;

     iBase_EntityHandle *nodeHandles = hopoints->nodeHandles;

     vector<double> u(numHPoints);
     vector<double> v(numHPoints);

     iBase_EntityHandle currvertex;

     offset = 0;
     currvertex = nodeHandles[offset];
     mesh->getVtxCoord(currvertex, p3d[0], p3d[1], p3d[2]);
     uv = currface.getUVCoords(p3d, uvCentroid);
     u[offset] = uv[0];
     v[offset] = uv[1];

     offset = nx - 1;
     currvertex = nodeHandles[offset];
     mesh->getVtxCoord(currvertex, p3d[0], p3d[1], p3d[2]);
     uv = currface.getUVCoords(p3d, uvCentroid);
     u[offset] = uv[0];
     v[offset] = uv[1];

     offset = (ny - 1) * nx;
     currvertex = nodeHandles[offset];
     mesh->getVtxCoord(currvertex, p3d[0], p3d[1], p3d[2]);
     uv = currface.getUVCoords(p3d, uvCentroid);
     u[offset] = uv[0];
     v[offset] = uv[1];

     offset = nx * ny - 1;
     currvertex = nodeHandles[offset];
     mesh->getVtxCoord(currvertex, p3d[0], p3d[1], p3d[2]);
     uv = currface.getUVCoords(p3d, uvCentroid);
     u[offset] = uv[0];
     v[offset] = uv[1];

     for (int i = 1; i < nx - 1; i++) {
          offset = i;
          currvertex = nodeHandles[offset];
          mesh->getVtxCoord(currvertex, p3d[0], p3d[1], p3d[2]);
          uv = currface.getUVCoords(p3d, uvCentroid);
          u[offset] = uv[0];
          v[offset] = uv[1];

          offset = i + (ny - 1) * nx;
          currvertex = nodeHandles[offset];
          mesh->getVtxCoord(currvertex, p3d[0], p3d[1], p3d[2] );
          uv = currface.getUVCoords(p3d, uvCentroid);
          u[offset] = uv[0];
          v[offset] = uv[1];
     }

     for (int j = 1; j < ny - 1; j++) {
          offset = j*nx;
          currvertex = nodeHandles[offset];
          mesh->getVtxCoord(currvertex, p3d[0], p3d[1], p3d[2] );
          uv = currface.getUVCoords(p3d, uvCentroid);
          u[offset] = uv[0];
          v[offset] = uv[1];

          offset = j * nx + (nx - 1);
          currvertex = nodeHandles[offset];
          mesh->getVtxCoord(currvertex, p3d[0], p3d[1], p3d[2]);
          uv = currface.getUVCoords(p3d, uvCentroid);
          u[offset] = uv[0];
          v[offset] = uv[1];
     }

     TFIMap::blend_from_edges(u, gllnodes, gllnodes);
     TFIMap::blend_from_edges(v, gllnodes, gllnodes);

     for (int j = 1; j < ny - 1; j++) {
          for (int i = 1; i < nx - 1; i++) {
               offset = j * nx + i;
               uv[0] = u[offset];
               uv[1] = v[offset];
               Point3D pon = currface.getXYZCoords(uv);
               mesh->setVtxCoord(nodeHandles[offset], pon[0], pon[1], pon[2]);
          }
     }
}

////////////////////////////////////////////////////////////////////////////////

void SpectralElements::project_on_faces()
{
     int err, namelen;
     int geom_dim, index = 0;

     iBase_EntitySetHandle geomRootSet, meshSet;
     geomRootSet = geom->getRootSet();

     vector<iBase_EntityHandle> gFaces, mFaces, mEdges;
     geom->getEntities(geomRootSet, iBase_FACE, gFaces);

     iBase_TagHandle id_tag, dim_tag;
     const char *tag1 = "GLOBAL_ID";
     mesh->getTagHandle(tag1, id_tag);

     const char *tag2 = "GEOM_DIMENSION";
     mesh->getTagHandle(tag2, dim_tag);
     assert(!err);

     vector<int> gFaceProcessed(gFaces.size());

     TFIMap tfimap(mesh, geom, rel, relPair);
//  for (int i = 0; i < gFaces.size(); i++)
     for (int i = 3; i < gFaces.size(); i++) {
          gFaceProcessed[i] = 0;
          if (tfimap.getTFI2D(gFaces[i], gllnodes) == 0) {
               gFaceProcessed[i] = 1;
               ostringstream oss;
               oss << "tfi" << i << ".dat";
               tfimap.saveAs( oss.str() );
               cout << " CSV Return " << endl;
               return;
          }
     }

     // To Visualize TFI Mapping, store the mesh
     iBase_EntitySetHandle rootSet;
     rootSet = mesh->getRootSet();

     const char *outfile = "meshcubit1.vtk";
     namelen = strlen(outfile);
//    err = mesh->save(rootSet, outfile, NULL, &err, namelen, 0);

     vector<iBase_EntityHandle> edgeNodes;
     for (int i = 0; i < gFaces.size(); i++) {
          if( !gFaceProcessed[i] ) {
               GFace currface(gFaces[i], geom, mesh, rel, relPair);
               currface.projectEdgeHigherOrderNodes(gllnodes);
          }
     }

     // Only at this stage all edges have been projected and discretized. After this
     // stage, no modification on the edges will be made. Check it now.

     for (int i = 0; i < gFaces.size(); i++) {
          if( !gFaceProcessed[i] ) {
               GFace currface(gFaces[i], geom, mesh, rel, relPair);
               currface.projectFaceHigherOrderNodes(gllnodes);
          }
     }
}

////////////////////////////////////////////////////////////////////////////////

void SpectralElements::generate(iMesh *m, int order,
                                iGeom *g, iRel *r, iRel::PairHandle *p)
{
     generate(m, order);

     hasGeometry = 1;
     geom    = g;
     rel     = r;
     relPair = p;

     int err, namelen, geom_dim;

     iBase_EntitySetHandle geomRootSet;
     geomRootSet = geom->getRootSet();

     vector<iBase_EntityHandle> gFaces, mFaces, gRegions;
     geom->getEntities(geomRootSet, iBase_FACE, gFaces);

     // Check wthere we have structured mesh on geometric faces and cells. If we have
     // structured mesh on them, we can apply TFI on them.

     for( size_t i = 0; i < gFaces.size(); i++) {
          if( hasStructuredMesh2D( gFaces[i] )  )
               structured_mesh[gFaces[i]] = 1;
          else
               structured_mesh[gFaces[i]] = 0;
     }

     geom->getEntities(geomRootSet, iBase_REGION, gRegions);
     for( size_t i = 0; i < gRegions.size(); i++) {
          if( hasStructuredMesh3D( gRegions[i] )  )
               structured_mesh[gRegions[i]] = 1;
          else
               structured_mesh[gRegions[i]] = 0;
     }

     project_on_edges();
     project_on_faces();

     // So after projection, even linear faces, adjoining the geometry got
     // changed. It is necessary to remesh both the surface and volume mesh;

     iBase_TagHandle dim_tag;
     const char *tag = "GEOM_DIMENSION";
     err = mesh->getTagHandle(tag, dim_tag);
     assert(!err);

     vector<iBase_EntityHandle> facecells;
     std::set<iBase_EntityHandle> boundFaces, boundCells;

     /*
         iBase_EntitySetHandle meshSet;
         for (int i = 0; i < gFaces.size(); i++)
         {
     //      iRel_getEntSetAssociation(assoc, rel, gFaces[i], 0, &meshSet, &err);
             mesh->getEntSetIntData(meshSet, dim_tag, geom_dim);

             if (geom_dim == 2)
             {
                 mFaces.clear();
                 mesh->getEntities(meshSet, iBase_FACE, iMesh_ALL_TOPOLOGIES, mFaces);

                 for (int j = 0; j < mFaces.size(); j++)
                 {
                     facecells.clear();
                     mesh->getEntAdj(mFaces[j], iBase_REGION, facecells);
                     assert(facecells.size() == 1);
                     boundFaces.insert(mFaces[j]);
                     boundCells.insert(facecells[0]);
                 }
             }
         }

         iBase_EntityHandle mFaceHandle, mCellHandle;

         vector<iBase_EntityHandle> cellfaces;
         std::set<iBase_EntityHandle> planarFaces;

         BOOST_FOREACH( mCellHandle, boundCells)
         {
             mesh->getEntAdj(mCellHandle, iBase_FACE, cellfaces);
             for (int i = 0; i < cellfaces.size(); i++)
                 planarFaces.insert(cellfaces[i]);
         }

         BOOST_FOREACH(mFaceHandle, boundFaces)
               planarFaces.erase(mFaceHandle);

         // Do remapping for both bounding faces and cells.
         TFIMap  tfimap(mesh);
         tfimap.setType( TFIMap::PHYSICAL_TFI );
         BOOST_FOREACH(mFaceHandle, planarFaces)
              tfimap.projectHigherOrderNodes2D( mFaceHandle, gllnodes );

         BOOST_FOREACH( mCellHandle, boundCells)
              tfimap.projectHigherOrderNodes3D( mCellHandle, gllnodes);
     */
}

////////////////////////////////////////////////////////////////////////////////

void SpectralElements::saveAs(const string &filename)
{
     ofstream ofile(filename.c_str(), ios::out);

     if (ofile.fail()) {
          cout << "Warning: Cann't open file " << filename << endl;
          return;
     }

     int err;
     vector<iBase_EntityHandle> nodeHandles;
     mesh->getEntities(meshRootSet, iBase_VERTEX, iMesh_ALL_TOPOLOGIES, nodeHandles);
     int numNodes = nodeHandles.size();

     ofile << "#Nodes  " << numNodes << endl;

     double x, y, z;

     int index = 0;
     for (int i = 0; i < numNodes; i++) {
          err = mesh->getVtxCoord(nodeHandles[i], x, y, z);
          err = mesh->setIntData(nodeHandles[i], idtag, index);
          ofile << index << "  " << fixed << x << " " << y << " " << z << endl;
          index++;
     }

     if (!hasSubStructuredGrid) return;

     vector<int> connect;

     vector<iBase_EntityHandle> cellHandles;
     err = mesh->getEntities(meshRootSet, iBase_REGION, iMesh_ALL_TOPOLOGIES, cellHandles);
     int numCells = cellHandles.size();

     int numSubCells = subCellDim[0] * subCellDim[1] * subCellDim[2];

     ofile << "#Hex  " << numSubCells * numCells << endl;

     for (int k = 0; k < numCells; k++) {
          get_subcells(cellHandles[k], connect);
          index = 0;
          for (int j = 0; j < numSubCells; j++) {
               for (int i = 0; i < 8; i++) ofile << connect[index++] << " ";
               ofile << endl;
          }
     }
     cout << " Data Saved " << endl;
}

////////////////////////////////////////////////////////////////////////////////
