#include <meshkit/DijkstraShortestPath.h>

///////////////////////////////////////////////////////////////////////////////

void DijkstraShortestPath::initialize()
{
    if (!mesh->isPruned()) mesh->prune();

    size_t numNodes = mesh->getSize(0);
    vnodes.resize(numNodes);

    for (size_t i = 0; i < numNodes; i++)
    {
        Vertex *vertex = mesh->getNodeAt(i);
        vertex->setVisitMark(0);
    }

    vsrc->setVisitMark(1);

    for (size_t i = 0; i < numNodes; i++)
    {
        Vertex *vertex = mesh->getNodeAt(i);
        vnodes[i].vertex = vertex;
        vnodes[i].previous = NULL;
        if (vertex->isVisited())
            vnodes[i].distance = 0.0;
        else
            vnodes[i].distance = MAXDOUBLE;
    }

}

///////////////////////////////////////////////////////////////////////////////

int DijkstraShortestPath::atomicOp(LVertex &currnode)
{
    Vertex *vi = currnode.vertex;
    if (vi == vdst) return 0;

    if( filter )  {
        if( !vi->isVisited() )  {           // Avoids the sources ...
             if( !filter->pass( vi ) ) {
                  vdst = vi;
                  return 0;
             }
        }
    }

    NodeSequence vneighs = vi->getRelations0();
    std::random_shuffle( vneighs.begin(), vneighs.end() ); // For randomization ...

    int nid;
    double vcost;
    for (size_t i = 0; i < vneighs.size(); i++)
    {
        Vertex *vj = vneighs[i];
        nid = vj->getID();
        vcost = getCost(currnode, vnodes[nid]);
        if (vcost < vnodes[nid].distance)
        {
            vnodes[nid].previous = vi;
            vnodes[nid].distance = vcost;
            vertexQ.push(vnodes[nid]);
        }
    }

    return 1;
}

///////////////////////////////////////////////////////////////////////////////

void DijkstraShortestPath::traceback()
{
    nodepath.clear();

#ifdef SEQUENCE_IS_VECTOR
    nodepath.reserve(vnodes.size());
#endif
    //
    // If the destination is set, then return the sequence from source to the
    // destination, otherwise return all the nodes with the increasing distance
    // from the source ...
    //
    if (vdst)
    {
        for (size_t i = 0; i < vnodes.size(); i++)
        {
            if (vnodes[i].vertex == vdst)
            {
                LVertex currnode = vnodes[i];
                while (1)
                {
                    nodepath.push_back(currnode.vertex);
                    Vertex *v = currnode.previous;
                    if (v == NULL) break;
                    currnode = vnodes[v->getID()];
                }
                break;
            }
        }
    }
    else
    {
        sort(vnodes.begin(), vnodes.end());
        for (size_t i = 0; i < vnodes.size() - 1; i++)
            assert(vnodes[i].distance <= vnodes[i + 1].distance);

        nodepath.resize(vnodes.size());
        for (size_t i = 0; i < vnodes.size(); i++)
            nodepath[i] = vnodes[i].vertex;
    }

}
///////////////////////////////////////////////////////////////////////////////

void DijkstraShortestPath::fastmarching()
{
    int relexist0 = mesh->build_relations(0, 0);

    initialize();

    while (!vertexQ.empty()) vertexQ.pop();

    vertexQ.push(vnodes[vsrc->getID() ]);

    int progress;
    while (!vertexQ.empty())
    {
        LVertex currVertex = vertexQ.top();
        vertexQ.pop();
        progress = atomicOp(currVertex);
        if (!progress) break;
    }

    if( !relexist0) mesh->clear_relations(0, 0);

}
///////////////////////////////////////////////////////////////////////////////

int Jaal::dijkstra_shortest_path_test()
{
  double origin[]   = { 0.0, 0.0, 0.0};
  double length[]   = { 1.0, 1.0, 1.0};
  int    gridim[]   = { 5, 5, 5};

  Jaal::Mesh *mesh = Jaal::create_structured_mesh(origin, length, gridim, 2 );

  struct MyFilter: public MeshFilter 
  {
      MyFilter( int i ) { id = i; }
      int id;
      bool pass( const Vertex *v) { 
         return v->getID() != id; 
      }
  };

  MeshFilter *filter = new MyFilter(18);

  DijkstraShortestPath djk(mesh);

  NodeSequence sq = djk.getPath(mesh->getNodeAt(0), filter);

  cout << " Final Path " << endl;
  for( int i = 0; i < sq.size(); i++)  
       cout << sq[i]->getID() << endl;
  
  delete filter;
}

///////////////////////////////////////////////////////////////////////////////


#ifdef JAAL_EXAMPLE
void Example::shortest_path(const string &filename)
{
    OFFMeshImporter mimp;
    SharedFaceMesh facemesh = mimp.read_faces(filename);

    string attribname = "Geodesic";
    DijkstraShortestPath2D djk(facemesh);
    djk.setAttributeName(attribname);
    djk.execute();

    VTKMeshExporter vtk;
    vtk.addAttribute(attribname, 0);
    vtk.saveAs(facemesh, "geodesic.vtk");
}

///////////////////////////////////////////////////////////////////////////////
#endif
