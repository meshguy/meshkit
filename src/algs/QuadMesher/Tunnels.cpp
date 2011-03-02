#include "QuadCleanUp.hpp"

#ifdef NOT_IMPLEMENTED

vector<Edge>
QuadCleanUp::search_tunnels()
{
    size_t numnodes = mesh->getSize(0);

    int relexist = mesh->build_relations(0, 0);

    mesh->search_boundary();

    vector<PVertex> vneighs;
    vTunnels.clear();

    int ncount;
    for (size_t i = 0; i < numnodes; i++)
    {
        Vertex *vertex = mesh->getNodeAt(i);
        ncount = 0;
        if (!vertex->isBoundary())
        {
            vneighs = vertex->getRelations0();
            for (size_t j = 0; j < vneighs.size(); j++)
                if (vneighs[j]->isBoundary()) {
                    Edge edge(vertex, vneighs[j] );
                    if( vertex < vneighs[j] ) {
                        if( isTunnel( &edge) ) 
                            vTunnels.push_back( Edge( vertex, vneighs[j] ) );
                   }
                }
        }
    }

    if (!relexist)
        mesh->clear_relations(0, 0);

    return vTunnels;
}


int
QuadCleanUp::remove_tunnels_once()
{
    /*
      vector<Vertex*> vrestrict = search_restricted_nodes ();

      if (vrestrict.size () == 0) return 0;

      int relexist0 = mesh->build_relations (0, 0);
      int relexist2 = mesh->build_relations (0, 2);

      vector<Vertex*> vneighs;
      vector<RestrictedEdge> redges;

      for (size_t i = 0; i < vrestrict.size (); i++)
        {
          Vertex *vertex = vrestrict[i];
          vneighs = vertex->getRelations0 ();
          for (size_t j = 0; j < vneighs.size (); j++)
            {
              if (vneighs[j]->isBoundary ())
                {
                  RestrictedEdge edge (mesh, vertex, vneighs[j]);
                  int err = edge.build ();
                  if (!err)
                    {
                      redges.push_back (edge);
                    }
                }
            }
        }

      if (redges.size ())
        cout << "Info:  # of valid restricted edges : " << redges.size () << endl;

      int ncount = 0;
      for (size_t i = 0; i < redges.size (); i++)
        {
          int err = redges[i].commit ();
          if (!err) ncount++;
        }

      if (ncount)
        {
          mesh->prune ();
          mesh->enumerate (0);
          mesh->enumerate (2);
          cout << "Info: # of restricted edges committed : " << ncount << endl;
        }

      if (!relexist0)
        mesh->clear_relations (0, 0);

      if (!relexist2)
        mesh->clear_relations (0, 2);

      lapsmooth->execute ();

      return ncount;
     */
    return 1;
}

////////////////////////////////////////////////////////////////////////////

void
QuadCleanUp::remove_tunnels()
{
    int ncount;
    while (1)
    {
        ncount = free_restricted_nodes_once();
        if (ncount == 0) break;
    }

    if (!mesh->isConsistentlyOriented())
        mesh->makeConsistentlyOriented();
}

////////////////////////////////////////////////////////////////////////////

#endif


