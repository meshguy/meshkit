/*
 * MBVolOp.cpp
 *
 *  Created on: Jan 13, 2012
 */

#include "meshkit/MBVolOp.hpp"
#include "meshkit/MBSplitOp.hpp"
#include "meshkit/ModelEnt.hpp"
#include "moab/CartVect.hpp"

//#include "moab/Core.hpp"
//#include "meshkit/FBiGeom.hpp"


namespace MeshKit {


//Entity Type initialization for splitting; no mesh output
moab::EntityType MBVolOp_tps[] = { moab::MBMAXTYPE }; // no mesh, really
const moab::EntityType* MBVolOp::output_types()
{
  return MBVolOp_tps;
}

MBVolOp::MBVolOp(MKCore *mk_core, const MEntVector &me_vec) :
  MeshScheme(mk_core, me_vec)
{
  _direction[0]=_direction[1]=0.;
  _direction[2]=1.0;
  _pGTT = NULL;
  _rootSet = 0; // it means no set yet, although 0 mean root set in moab
  _fbe = NULL;
}

MBVolOp::~MBVolOp() {
  // TODO Auto-generated destructor stub
}
void MBVolOp::setup_this()
{
  // construct an FBEngine object, used more like a container, and to trigger the weaving
  // it is not really required; this object is based on gtt object build from top and bottom faces
  // ( which are extracted from model entities of dimension 2)
  // so, involve FBEngine just because we have something we need there
  // collect the top and bottom faces from ment vector
  int nTotSurf = this->mentSelection.size();
  std::cout << " total number of faces:" << nTotSurf << "\n";

  // grab all surfaces, and duplicate model using gtt, to get new gsets, that will be
  // continued with volume sets in the end
  // establish the loops from faces
  MEntSelection::iterator mit;

  moab::Interface * MBI = mk_core()->moab_instance();
  moab::GeomTopoTool gtt(MBI, true);// to retrieve the gsets
  moab::EntityHandle mset;
  // these should all be of dimension 2, faces
  std::vector<moab::EntityHandle> vSurfaces;

  moab::ErrorCode rval;

  for (mit = mentSelection.begin(); mit != mentSelection.end(); mit++) {
    mset = (mit->first)->mesh_handle();
    // get the globalId tag
    vSurfaces.push_back(mset);
  }

  rval = gtt.duplicate_model(_pGTT, &vSurfaces);
  MBERRCHK(rval, MBI);
  // now, _pGTT will have the new sets, which will form the basis for the volume creation
  // new gsets will be added to this _pGTT, and this will (unfortunately) create an OBB tree too,
  // although we do not really need it (actually, 2 obb trees that are not needed!) bummer!
  // the result will be the model root set of the new GTT, stored in result range
  // that corresponds to the first model ent
  ModelEnt * firstMe = (*(mentSelection.begin()) ).first;
  moab::Range & resultRange = mentSelection[firstMe];
  _rootSet = _pGTT->get_root_model_set();
  resultRange.insert(_rootSet);
  _fbe = new moab::FBEngine(MBI, _pGTT, false);// smooth=false, not needed here (are you sure?)
  //we will use FBEngine for querying; although ModelEnt would have helped
  // still, the model ent works with geometry adjacency; in this case, there is no geometry
  // we need to work directly with MOAB;
  // not pretty :(
  // we build all this infrastructure and I am not using it


  establish_mapping();
  return;
}
void MBVolOp::establish_mapping()
{
  // the first half of the surfaces are bottom, next are top
  // establish the connection between top and bottom entities
  // use the direction for vertices, and tangent direction for edges

  // these are all vertices from current gtt
  moab::Interface * MBI = mk_core()->moab_instance();
  moab::Range verticeSets= _pGTT->geoRanges()[0];
  int num_vertices = verticeSets.size();
  std::vector<moab::CartVect> coordVert;
  coordVert.resize(num_vertices);
  std::vector<int> corrV;
  corrV.resize(num_vertices);
  moab::ErrorCode rval = moab::MB_SUCCESS;
  std::map<moab::EntityHandle, int> indexVertex;
  for (int i=0; i<num_vertices; i++)
  {
    rval = _fbe->getVtxCoord(verticeSets[i], &(coordVert[i][0]),
        &(coordVert[i][1]), &(coordVert[i][2]));
    MBERRCHK(rval, MBI);
    corrV[i]=-1;// no correspondence yet
    indexVertex[verticeSets[i]] = i;
  }

  // decide with an expensive linear search, what vertices correspond to what
  moab::CartVect dirr(_direction);
  dirr.normalize();
  int j=-1;
  for ( j=0; j<num_vertices; j++)
  {
    if (corrV[j]!=-1)
      continue; // we already know about this one
    moab::CartVect & v1 = coordVert[j];
    int minIndex = -1;
    double minVal = 1.e38; // HUGE
    for (int k=0; k<num_vertices; k++)
    {
      if (j==k || corrV[k]>-1)
        continue;
      moab::CartVect product = (v1-coordVert[k])*dirr;
      double valProd = product.length_squared();
      if (valProd<minVal)
      {
        minVal = valProd;
        minIndex=k;
      }
    }
    std::cout<<"j: "<< j <<" val min:" << minVal << " min index:" << minIndex << "\n";
    if (minIndex==-1)
    {
      std::cout<< "Error Stop it. j: " << j << "\n";
      continue;
    }
    // make sure that the lower index is on bottom
    double dotProd = (coordVert[minIndex]-coordVert[j])%dirr;
    if (dotProd < 0)
    {
      std::cout<<"wrong orientation, bottom vertices should be of lower index\n";
      MBERRCHK(moab::MB_FAILURE, MBI);// bail out from here, stop it!
    }

    corrV[j] = minIndex;
    corrV[minIndex] = j;
  }
  // check that we do not have any left vertices
  for (j=0; j<num_vertices; j++)
  {
    if (corrV[j]==-1)
    {
      std::cout<<"Error in finding a corresponding vertex to j:"<<j <<"\n";
      MBERRCHK(moab::MB_FAILURE, MBI);// bail out from here, stop it!
    }
    vertexMap[verticeSets[j]]= verticeSets[corrV[j]];
  }

  // so now we have mappings between vertices
  // we need to build mappings between edges and faces.
  // we usually start with bottom and continue to the top
  moab::Range edgeSets= _pGTT->geoRanges()[1];
  int numEdges = edgeSets.size();
  std::vector<moab::CartVect> edgeTangents;
  edgeTangents.resize(numEdges*2);// we need 2 tangents per edge
  // for each vertex, compute the tangents at ends, projected on a plan perpendicular to the direction
  // (usually, xy plan...)
  // tangents at both ends!// we know that the u range is from 0 to 1
  std::vector<int> corrE;
  corrE.resize(numEdges);
  std::map<moab::EntityHandle, int> indexEdge;
  for (j=0; j<numEdges; j++)
  {
    corrE[j] = -1;// no correspondence yet
    indexEdge[edgeSets[j]] = j;
    // careful about the FBEngine, it is not smooth
    std::vector<moab::EntityHandle> meshEdges;
    rval = MBI->get_entities_by_type(edgeSets[j], moab::MBEDGE,  meshEdges);
    MBERRCHK(rval, MBI);
    // get first and last edges, connectivity and compute tangent
    // it will be used to map edges
    // first edge, last edge
    moab::EntityHandle firstMeshEdge = meshEdges[0];
    moab::EntityHandle lastMeshEdge = meshEdges[meshEdges.size()-1];
    const moab::EntityHandle  * conn=NULL;
    int nnodes;
    rval = MBI->get_connectivity(firstMeshEdge, conn, nnodes);
    MBERRCHK(rval, MBI);
    moab::CartVect posi[2];
    rval = MBI->get_coords(conn, 2, &(posi[0][0]));
    MBERRCHK(rval, MBI);
    posi[0] = posi[1]-posi[0]; //
    posi[0] = posi[0]*dirr;
    edgeTangents[2*j] = posi[0]*dirr;// this is in the plane
    edgeTangents[2*j].normalize();// it should be non null, but accidents happen :)
    rval = MBI->get_connectivity(lastMeshEdge, conn, nnodes);
    MBERRCHK(rval, MBI);
    rval = MBI->get_coords(conn, 2, &(posi[0][0]));
    MBERRCHK(rval, MBI);
    posi[0] = posi[1]-posi[0]; //
    posi[0] = posi[0]*dirr;
    edgeTangents[2*j+1] = posi[0]*dirr;// this is in the plane
    edgeTangents[2*j+1].normalize();
  }
  // now try to match edges based on their vertices matching, and start and end tangents matching
  for (j = 0; j<numEdges; j++)
  {
    if (corrE[j]>=0)
      continue; // we already have a correspondent edge for it
    moab::Range adjVertices;
    rval = _fbe-> getEntAdj(edgeSets[j], /*vertex type*/0, adjVertices);
    MBERRCHK(rval, MBI);
    if (adjVertices.size()==2)
    {
      int sense=0;
      rval = _fbe->getEgVtxSense(edgeSets[j], adjVertices[0], adjVertices[1], sense);
      MBERRCHK(rval, MBI);
      // find the edges that are adjacent to map(v0) and map(v1)
      moab::EntityHandle mapV0 = vertexMap [adjVertices[0]];
      moab::EntityHandle mapV1 = vertexMap [adjVertices[1]];
      // get all edges adjacent to both vertices
      moab::Range adjEdges0, adjEdges1;
      rval = _fbe-> getEntAdj(mapV0, /*edge type*/1, adjEdges0);
      MBERRCHK(rval, MBI);
      rval = _fbe-> getEntAdj(mapV1, /*edge type*/1, adjEdges1);
      adjEdges0 = intersect(adjEdges0, adjEdges1);
      // the mapped edge should be in this range
      moab::CartVect & t0=edgeTangents[2*j];
      moab::CartVect & t1=edgeTangents[2*j+1];
      for (moab::Range::iterator rit= adjEdges0.begin(); rit!=adjEdges0.end(); rit++)
      {
        moab::EntityHandle candidateMapEdge = *rit;
        int sense1=0;
        int indexCandEdge = indexEdge[candidateMapEdge];
        if (indexCandEdge == j)
          // error
          MBERRCHK(moab::MB_FAILURE, MBI);
        moab::CartVect & tm0=edgeTangents[2*indexCandEdge];
        moab::CartVect & tm1=edgeTangents[2*indexCandEdge+1];
        rval = _fbe->getEgVtxSense(candidateMapEdge, mapV0, mapV1, sense1);
        MBERRCHK(rval, MBI);

        if (sense*sense1>0)// same sense
        {
          if ( (t0%tm0 >= 0.99 ) && (t1%tm1>=0.99))// same sense, almost..
          {
            corrE[j]=indexCandEdge;
            corrE[indexCandEdge] = j;
          }
        }
        else
        {
          // different senses, check the opposite tangents
          if ( (t0%tm0 <=-0.99) && (t1%tm1<=-0.99) )
          {
            // how do we store the opposite senses?
            corrE[j]=indexCandEdge;
            corrE[indexCandEdge] = j;
          }
        }
      }
      if (corrE[j]==-1)
        MBERRCHK(moab::MB_FAILURE, MBI);// can't find a corresponding edge

    }
  }
  // so we have edge matching for every edge
  for (j=0; j<numEdges; j++)
  {
    if (corrE[j]==-1)
    {
      std::cout<<"Error in finding a corresponding edge: "<<j <<"\n";
      MBERRCHK(moab::MB_FAILURE, MBI);// bail out from here, stop it!
    }
    std::cout<< "edge j=" << j << "  mapped to edge " << corrE[j] << "\n";
    edgeMap[edgeSets[j]]= edgeSets[corrE[j]];
  }
  // now, we still need to separate top and bottom faces, somehow
  // we assume faces were stenciled correctly with the same polylines
  // start from bottom vertices (first half of the vertices is on bottom)
  moab::Range bottomFaces;
  for (j=0; j<num_vertices/2; j++)// we know that half are on bottom, we verified that already
  {
    moab::Range faces;
    rval = _fbe-> getEntAdj(verticeSets[j], /*face type*/2, faces);
    MBERRCHK(rval, MBI);
    bottomFaces.merge(faces);
  }
  // now find the mapped face for each of these, based on mapped edges
  for (j=0; j<(int)bottomFaces.size(); j++)
  {
    moab::Range edges;
    rval = _fbe-> getEntAdj(bottomFaces[j], /*edge type*/1, edges);
    MBERRCHK(rval, MBI);
    moab::Range mappedEdges;
    // get now the face connected to the mapped edges
    int i=0;
    // find the face with all those edges among the other faces

    moab::Range mapFaces;
    rval = _fbe-> getEntAdj(edgeMap[edges[0]], /*edge type*/2, mapFaces);
    MBERRCHK(rval, MBI);
    for (i=1; i<(int)edges.size(); i++)
    {
      moab::Range faces;
      rval = _fbe-> getEntAdj(edgeMap[edges[i]], /*edge type*/2, faces);
      mapFaces = intersect(mapFaces, faces);
    }
    if (mapFaces.size()!=1)
    {
      std::cout<<"Can't find unique mapped face to face index " << j << "\n";
      MBERRCHK(moab::MB_FAILURE, MBI);// bail out from here, stop it!
    }
    faceMap[bottomFaces[j]] = mapFaces[0];
    faceMap[mapFaces[0]]=bottomFaces[j];
    std::cout << "face j:" << j << " set:" << MBI->id_from_handle(bottomFaces[j]) << " mapped to "
        << MBI->id_from_handle(mapFaces[0]) << "\n";
  }

}

void MBVolOp::execute_this()
{
  // now, we have the maps between top and bottom faces established, also for edges and vertices
  // first build edges between corresponding vertices, then faces between edges, then
  // finally, volumes
  // we start from bottom towards top
  moab::ErrorCode rval = moab::MB_SUCCESS;
  moab::Interface * MBI = mk_core()->moab_instance();

  moab::Range verticeSets=_pGTT->geoRanges()[0];
  int num_vertices = verticeSets.size();
  moab::Range bottomEdges;
  int j=0;
  for (j=0; j<num_vertices/2; j++)// we know that half are on bottom, we verified that already
  {
    moab::Range edges;
    rval = _fbe-> getEntAdj(verticeSets[j], /*edge type*/1, edges);
    MBERRCHK(rval, MBI);
    bottomEdges.merge(edges);
  }
  // we need to decide bottom faces before we create new faces, edges, etc
  moab::Range bottomFaces;
  for (j=0; j<num_vertices/2; j++)// we know that half are on bottom, we verified that already
  {
    moab::Range faces;
    rval = _fbe-> getEntAdj(verticeSets[j], /*face type*/2, faces);
    MBERRCHK(rval, MBI);
    bottomFaces.merge(faces);
  }

  // we have verified that the first half are from bottom surfaces
  // for each bottom face we will create a volume;
  // start with the edges, create weaving faces between edges (no new vertices!)
  // now, for each bottom edge, create a weaving face
  std::vector<moab::EntityHandle> newFaces;
  newFaces.resize(bottomEdges.size());
  for (j=0; j<(int)bottomEdges.size(); j++)
  {
    rval = _fbe->weave_lateral_face_from_edges(bottomEdges[j], edgeMap[bottomEdges[j]],  _direction, newFaces[j]);
    MBERRCHK(rval, MBI);
  }
  // now, to create volumes, we need to loop over bottom faces and add volumes one by one, and the
  // orientation of faces in them

  int volumeMatId = 0;// just give a mat id, for debugging, mostly
  moab::Tag matTag;
  rval = MBI->tag_get_handle("MATERIAL_SET", 1, moab::MB_TYPE_INTEGER, matTag);
  MBERRCHK(rval, MBI);

  for (j = 0; j<(int)bottomFaces.size(); j++)
  {
    // create volume here , finally!!!
    moab::EntityHandle volume;
    rval = MBI->create_meshset(moab::MESHSET_SET, volume);
    MBERRCHK(rval, MBI);
    volumeMatId++;
    rval= MBI->tag_set_data(matTag, &volume, 1, &volumeMatId);
    MBERRCHK(rval, MBI);
    moab::EntityHandle botFace = bottomFaces[j];
    moab::EntityHandle topFace = faceMap[botFace];
    // start copy
    // get the edges of bot face
    rval= MBI->add_parent_child(volume, botFace);
    MBERRCHK(rval, MBI);

    rval= MBI->add_parent_child(volume, topFace);
    MBERRCHK(rval, MBI);

    rval = _pGTT->add_geo_set(volume, 3);
    MBERRCHK(rval, MBI);

    // set senses
    // bottom face is negatively oriented, its normal is toward interior of the volume
    rval = _pGTT->set_sense(botFace, volume, -1);
    MBERRCHK(rval, MBI);

    // the top face is positively oriented
    rval = _pGTT->set_sense(topFace, volume, 1);
    MBERRCHK(rval, MBI);

    // the children should be in the same direction
    //   get the side edges of each face, and form lateral faces, along direction
    std::vector<moab::EntityHandle> edges1;

    rval = MBI->get_child_meshsets(botFace, edges1); // no hops
    MBERRCHK(rval, MBI);

    for (unsigned int i = 0; i < edges1.size(); ++i)
    {
      // the orientation of edge in face  will give the sense of face in volume
      int indexB = bottomEdges.index(edges1[i]);
      if (indexB<=-1)
        MBERRCHK(moab::MB_FAILURE, MBI);
      int sense = 1;
      rval = _pGTT->get_sense(edges1[i], botFace, sense);
      MBERRCHK(rval, MBI);
      rval=MBI->add_parent_child(volume, newFaces[indexB]);
      MBERRCHK(rval, MBI);

      // set sense as the one decided from the bottom edge sense within the bottom face
      rval = _pGTT->set_sense(newFaces[indexB], volume, sense);
      MBERRCHK(rval, MBI);
    }
    // end copy
  }


  delete _fbe;
  delete _pGTT; // when we are done, remove the _pGTT;
  // at this point, the result will be the model root set of the first ment of the model
  // ( (*(mentSelection.begin()) ).second )
}

}
