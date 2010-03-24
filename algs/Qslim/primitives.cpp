/*
 * primitives.cpp
 *
 *  Created on: Mar 22, 2010
 *      Author: iulian
 */

#include "primitives.h"
#include "MBRange.hpp"

MBErrorCode contractionRegion(MBInterface * mb, MBEntityHandle v1, MBEntityHandle v2,
		std::vector<MBEntityHandle>& changed)
{
	MBEntityHandle vlist[2];
	vlist[0] = v1;
	vlist[1] = v2;
        // it makes more sense to use vector, than range
        //  the entities could be very disjoint at some point
	// MBRange adj_ents;
	MBErrorCode rval = mb->get_adjacencies(vlist, 2, 2, false, changed, MBInterface::UNION);
	return rval;
}

//
int classifyVertex(MBInterface * mb, MBEntityHandle v1)
{
   // return 0, 1, or 2 if the vertex is interior, on the border, or
   // border only (corner)
   // to do
   std::vector<MBEntityHandle> adjEdges;
   MBErrorCode rval = mb->get_adjacencies(&v1, 1, 1, false,  adjEdges,
                                      MBInterface::UNION);
   if (MB_SUCCESS!=rval)
       return 0; // interior??
   int nBorder = 0;
   for (int i=0; i<adjEdges.size(); i++)
   {
      MBEntityHandle edg=adjEdges[i];
      std::vector<MBEntityHandle> adjFaces;
      rval = mb->get_adjacencies(&edg, 1, 2, false,  adjFaces,
                                      MBInterface::UNION);
      if (adjFaces.size()==1)
         nBorder++;
   }
   if (nBorder==0)
      return 0;// everything interior
   else
      if (nBorder == adjEdges.size() )
	  return 2;
      else
          return 1; // some edges are interior
}

Vec3 convertMBVertexToVec3(MBInterface * mbi, MBEntityHandle v)
{
	double c[3];
	mbi->get_coords(&v, 1, c);
	return Vec3(c[0], c[1], c[2]);
}
// every time we are getting the normal, we compute a new plane
// maybe we should store it??
// this is debatable...
Plane trianglePlane(MBInterface * mb, MBEntityHandle tri)
{
	// get connectivity of triangle
	const MBEntityHandle * conn;
	int num_nodes;
	MBErrorCode rval = mb->get_connectivity(tri, conn, num_nodes);
	assert(3==num_nodes && rval == MB_SUCCESS);
	Vec3 ve1=convertMBVertexToVec3( mb, conn[0]);
	Vec3 ve2=convertMBVertexToVec3( mb, conn[1]);
	Vec3 ve3=convertMBVertexToVec3( mb, conn[2]);
	return Plane(ve1, ve2, ve3);

}

MBErrorCode contract (MBInterface * mb, MBEntityHandle v0, MBEntityHandle v1, Vec3 & vnew, std::vector<MBEntityHandle> & changed )
{

//
//// Collect all the faces that are going to be changed
////
     std::vector<MBEntityHandle> adj_entities;
     contractionRegion(mb, v0, v1, adj_entities);	
     // those are all triangles that are affected 
     // find also all edges that are affect
     MBEntityHandle vlist[2];
     vlist[0] = v0;
     vlist[1] = v1;
     // it makes more sense to use vector, than range
     //  the entities could be very disjoint at some point
     // MBRange adj_ents;
     std::vector<MBEntityHandle> edges;
     MBErrorCode rval = mb->get_adjacencies(vlist, 2, 1, false, edges, MBInterface::UNION);     
    // we have the edges and the triangles that are affected
    // 2 situations
    //   1) edge v0 v1 is existing 
    //      we will delete edge (v0, v1), and triangles formed 
    //       with edge (v0, v1)
    //   2) edge v0 v1 is not existing, but due to proximity
    //      only edges v2 v1 , v2, v0 will need to merge
    // more important is case 1)

    // first, find edge v0, v1
    MBEntityHandle ev0v1;
    int foundEdge =0;
    for (int i=0; i<edges.size(); i++)
    {
       MBEntityHandle e = edges[i];
       int nnodes;
       const MBEntityHandle * conn2;
       mb->get_connectivity(e, conn2, nnodes);
       if ( (conn2[0] == v0 && conn2[1]== v1) || 
            (conn2[0] == v1 && conn2[1]== v0) )
       {
           foundEdge = 1;
           ev0v1 = e; // could be ev1v0, but who cares?
           break;
       } 
           
    }
    // set the position of new vertices in vnew
    double newCoords[3];
    newCoords[0] = vnew[0];
    newCoords[1] = vnew[1];
    newCoords[2] = vnew[2];
    mb->set_coords(&v0, 1, newCoords);
    mb->set_coords(&v1, 1, newCoords);
// although, vertex v1 will be deleted!
//
   
    MBRange entitiesToDelete;
    // entitiesToDelete.insert(v1);// the vertex for sure
    // although vertex v1 will be merged!!
    std::vector<MBEntityHandle> edgePairsToMerge; // the one that has v0 will 
    // be kept
    if (foundEdge)
    {
       // this is case 1, the most complicated
       // get triangles connected to edge ev0v1
       std::vector<MBEntityHandle> tris;
       rval = mb->get_adjacencies(&ev0v1, 1, 2, false, tris, MBInterface::UNION);
       // find all edges that will be merged ( xv0, xv1, etc)
       for (int i=0; i<tris.size(); i++)
       {
	  MBEntityHandle triangleThatCollapses= tris[i];
          std::vector<MBEntityHandle> localEdges;
          rval = mb->get_adjacencies(&triangleThatCollapses, 1, 1, false, 
             localEdges, MBInterface::UNION);
          // find the edges that contains v0 
          MBEntityHandle e[2];// the 2 edges e0, e1, that are not ev0v1;
          if (localEdges.size()!=3)
              return MB_FAILURE; // failure
          int index=0;
          for (int k=0; k<3; k++)
             if (localEdges[k]!=ev0v1)
                 e[index++]=localEdges[k]; 
          // among those 2 edges, find out which one has v0, and which one v1
          if (index!=2)
              return MB_FAILURE; // failure
          for (int j=0; j<2; j++)
          {
              int nn;
              const MBEntityHandle * conn2;
	      mb->get_connectivity(e[j], conn2, nn);
              if (conn2[0] == v0 || conn2[1] == v0)
              {
		 // this is the edge that will be kept, the other one collapsed
		 edgePairsToMerge.push_back(e[j]);
                 j=(j+1)%2;// the other one 
                 edgePairsToMerge.push_back(e[j]);
                 break; // no need to check the other one. it 
                        // will contain v1  
              }
          }
       }
       // first merge vertices v0 and v1 : will also NOT delete v1 (yet)
       // the tag on v1 will be deleted too, and we do not want that, at least until
       // after the merging of edges, and deleting the pair
       rval = mb->merge_entities(v0, v1, false, false);
       // merge edgePairsToMerge // now, v0 and v1 should be collapsed!
       for (int j=0; j<edgePairsToMerge.size(); j+=2)
       {
           // will also delete edges that contained v1 before
	       mb->merge_entities(edgePairsToMerge[j], edgePairsToMerge[j+1],
              false, true);
       }
       // the only things that need deleted are triangles
       rval = mb->delete_entities(&(tris[0]), tris.size() );
       validFaceCount -= tris.size();
       // hopefully, all adjacencies are preserved
       // delete now the edge ev0v1
       rval=mb->delete_entities(&ev0v1 , 1);
       // among adj_entities, remove the tris triangles, and return them
       for (int j=0; j<adj_entities.size(); j++)
       {
    	   MBEntityHandle F = adj_entities[j];
    	   int inTris=0;
    	   for (int k=0; k<tris.size(); k++)
    		   if (F==tris[k])
    		   {
    			   inTris=1;break;
    		   }
    	   if (!inTris)
    		   changed.push_back(F);
       }
      

    }
    /*
 *  merge_entities(EntityHandle entity_to_keep, 
 *                          EntityHandle entity_to_remove,
 *                       bool auto_merge,
 *                   bool delete_removed_entity)
 */
    
     return MB_SUCCESS;
  
}

