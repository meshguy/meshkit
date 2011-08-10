#include "QuadCleanUp.hpp"

using namespace Jaal;

size_t OneDefectPatch :: MAX_FACES_ALLOWED = 500;

size_t OneDefectPatch :: num_boundaries  = 0;
size_t OneDefectPatch :: num_3_patches   = 0;
size_t OneDefectPatch :: num_4_patches   = 0;
size_t OneDefectPatch :: num_5_patches   = 0;
double OneDefectPatch :: exec_time       = 0.0;

///////////////////////////////////////////////////////////////////////////////
int TriRemeshTemplate::remesh(Mesh *msh,
                              NodeSequence &anodes,
                              NodeSequence &bnodes,
                              NodeSequence &cnodes, int *partition)
{
     mesh = msh;

     // First thing to do is to clear the existing record.
     newnodes.clear();
     newfaces.clear();

     // We need atleast three nodes on each side ...
     if (anodes.size() < 3) return 1;
     if (bnodes.size() < 3) return 1;
     if (cnodes.size() < 3) return 1;

     if (partition == NULL) {
          segments[0] = anodes.size() - 1;
          segments[1] = bnodes.size() - 1;
          segments[2] = cnodes.size() - 1;

          if (!Face::is_3_sided_convex_loop_quad_meshable(segments, partSegments))
               return 1;
     } else {
          for (int i = 0; i < 6; i++)
               partSegments[i] = partition[i];
     }


     int err;
     if (anodes.back() == bnodes.back()) {
          reverse(bnodes.begin(), bnodes.end());
          swap(partSegments[2], partSegments[3]);
     }

     if (anodes.front() == cnodes.front()) {
          reverse(cnodes.begin(), cnodes.end());
          swap(partSegments[4], partSegments[5]);
     }

     // Ensure that input is closed loop of three sides ...
     assert(anodes.front() == cnodes.back());
     assert(bnodes.front() == anodes.back());
     assert(cnodes.front() == bnodes.back());

     // Ensure that segments are valid ...
     assert((int)anodes.size() == partSegments[0] + partSegments[1] + 1);
     assert((int)bnodes.size() == partSegments[2] + partSegments[3] + 1);
     assert((int)cnodes.size() == partSegments[4] + partSegments[5] + 1);

     // Split Side "A" nodes into a1nodes and b1nodes.
     int a1 = partSegments[0];
     err = split_stl_vector(anodes, a1 + 1, a1nodes, b1nodes);
     if (err) return 1;

     // Split Side "B" nodes into a2nodes and b2nodes.
     int a2 = partSegments[2];
     err = split_stl_vector(bnodes, a2 + 1, a2nodes, b2nodes);
     if (err) return 1;

     // Split Side "C" nodes into a3nodes and b3nodes.
     int a3 = partSegments[4];
     err = split_stl_vector(cnodes, a3 + 1, a3nodes, b3nodes);
     if (err) return 1;

     // Splitting nodes on each side
     Vertex *ca = a1nodes.back();
     Vertex *cb = a2nodes.back();
     Vertex *cc = a3nodes.back();

     // We are not if the centroid will be in the patch, take the risk.
     Vertex *co = Face::centroid(ca, cb, cc);

     mesh->addNode(co);
     newnodes.push_back(co);

     linear_interpolation(mesh, co, ca, a2 + 1, oa_nodes);
     size_t nSize = oa_nodes.size();
     for (size_t i = 1; i < nSize - 1; i++)
          newnodes.push_back(oa_nodes[i]);

     linear_interpolation(mesh, co, cb, a3 + 1, ob_nodes);
     nSize = ob_nodes.size();
     for (size_t i = 1; i < nSize - 1; i++)
          newnodes.push_back(ob_nodes[i]);

     linear_interpolation(mesh, co, cc, a1 + 1, oc_nodes);
     nSize = oc_nodes.size();
     for (size_t i = 1; i < nSize - 1; i++)
          newnodes.push_back(oc_nodes[i]);

     err = quadtemplate.remesh(mesh, b2nodes, a3nodes, oc_nodes, ob_nodes);
     if (!err) addNewElements( quadtemplate.newnodes, quadtemplate.newfaces);

     if (!err) {
          err = quadtemplate.remesh(mesh, a1nodes, oa_nodes, oc_nodes, b3nodes);
          if (!err) addNewElements( quadtemplate.newnodes, quadtemplate.newfaces);
     }

     if (!err) {
          err = quadtemplate.remesh(mesh, a2nodes, ob_nodes, oa_nodes, b1nodes);
          if (!err) addNewElements( quadtemplate.newnodes, quadtemplate.newfaces);
     }

     if (err) {
          discard();
          return 1;
     }

     // Check for any inversion of the element, if there is inversion,
     // undo everthing (i.e. remove new nodes and faces).
     //
     nSize = newfaces.size();
     for (size_t i = 0; i < nSize; i++) {
          if (newfaces[i]->concaveAt() >= 0) {
               discard();
               return 2;
          }
     }

     return 0;
}

///////////////////////////////////////////////////////////////////////////////

int
QuadRemeshTemplate :: remesh()
{
     size_t nSize = boundnodes.size();
     assert( nSize == size_t(2 * nx + 2 * ny - 4));

     qnodes.resize(nx * ny);

     //
     // Put the boundary nodes on the structured mesh: The orientation is
     // Counter clockwise ( south->east->north->west );
     //
     int offset, index = 0;

     // South Side ...
     index = 0;
     for (int i = 0; i < nx; i++) {
          offset = i;
          qnodes[i] = boundnodes[index++];
     }

     // East Side
     for (int j = 1; j < ny; j++) {
          offset = j * nx + (nx - 1);
          qnodes[offset] = boundnodes[index++];
     }

     // North Side
     for (int i = nx - 2; i >= 0; i--) {
          offset = (ny - 1) * nx + i;
          qnodes[offset] = boundnodes[index++];
     }

     // West Side
     for (int j = ny - 2; j >= 1; j--) {
          offset = j*nx;
          qnodes[j * nx] = boundnodes[index++];
     }

     // New internal nodes must be generated.(either new or from  the trash bin)
     mesh->objects_from_pool( (nx-2)*(ny-2), newnodes );

     // Assign the coordinates using Transfinite-interpolation method. For 3D
     // Coordinates must be projected on the surface.
     //
     index = 0;
     for (int j = 1; j < ny - 1; j++) {
          for (int i = 1; i < nx - 1; i++) {
               offset = j * nx + i;
               qnodes[offset] = newnodes[index++];
               set_tfi_coords(i, j, nx, ny, qnodes); // Coordinates values
          }
     }

     // All new faces must be generated. ( Objects are either new or from the trash bin)
     mesh->objects_from_pool( (nx-1)*(ny-1), newfaces );

     NodeSequence qc(4);

     // Create new faces ...
     index = 0;
     for (int j = 0; j < ny - 1; j++) {
          for (int i = 0; i < nx - 1; i++) {
               int offset = j * nx + i;
               qc[0] = qnodes[offset];
               qc[1] = qnodes[offset + 1];
               qc[2] = qnodes[offset + 1 + nx];
               qc[3] = qnodes[offset + nx];
               Face *face = newfaces[index++];
               face->setNodes(qc);
          }
     }

     // Update the mesh ...
     nSize = newnodes.size();
     for (size_t i = 0; i < nSize; i++)
          mesh->reactivate( newnodes[i] );

     nSize = newfaces.size();
     for (size_t i = 0; i < nSize; i++)
          mesh->reactivate( newfaces[i] );

     // Check for any inversion of the element, if there is inversion,
     // undo everthing (i.e. remove new nodes and faces).
     //
     nSize = newfaces.size();
     for (size_t i = 0; i < nSize; i++) {
          if (newfaces[i]->concaveAt() >= 0) {
               discard();
               return 1;
          }
     }

     return 0;
}

////////////////////////////////////////////////////////////////////////////////

int
QuadRemeshTemplate ::remesh(Mesh *msh,
                            NodeSequence &anodes,
                            NodeSequence &bnodes,
                            NodeSequence &cnodes,
                            NodeSequence &dnodes )
{
     mesh = msh;

     newnodes.clear();
     newfaces.clear();

     assert(anodes.size() == cnodes.size());
     assert(bnodes.size() == dnodes.size());

     nx = anodes.size();
     ny = bnodes.size();

     boundnodes.resize(2 * nx + 2 * ny - 4);

     int index = 0;

     // First edge is always the reference edge which decides the orientation.
     // Append anodes..
     size_t nSize = anodes.size();
     for (size_t i = 0; i < nSize; i++)
          boundnodes[index++] = anodes[i];

     // Append bnodes ...
     if (bnodes.front() != anodes.back())
          reverse(bnodes.begin(), bnodes.end());

     assert(anodes.back() == bnodes.front());

     nSize = bnodes.size();
     for (size_t i = 1; i < nSize; i++)
          boundnodes[index++] = bnodes[i];

     // Append cnodes ...
     if (cnodes.front() != bnodes.back())
          reverse(cnodes.begin(), cnodes.end());

     assert(bnodes.back() == cnodes.front());

     nSize = cnodes.size();
     for (size_t i = 1; i < nSize; i++)
          boundnodes[index++] = cnodes[i];

     // Append dnodes ...
     if (dnodes.front() != cnodes.back())
          reverse(dnodes.begin(), dnodes.end());

     assert(cnodes.back() == dnodes.front());
     nSize = dnodes.size();
     for (size_t i = 1; i < nSize; i++)
          boundnodes[index++] = dnodes[i];

     // Ensure that loop is closed ...
     assert(anodes.front() == dnodes.back());

     return remesh();
}

////////////////////////////////////////////////////////////////////////////////

int
PentaRemeshTemplate :: remesh(Mesh *msh,
                              NodeSequence &anodes,
                              NodeSequence &bnodes,
                              NodeSequence &cnodes,
                              NodeSequence &dnodes,
                              NodeSequence &enodes,
                              int *partition)
{
     //
     // Implicit in this implementation is the orientation of the newly formed
     // faces which is decided by the orientation of the boundary. User need
     // not check for the consistency of the mesh. It will be, if the initial
     // mesh is consistent.
     //
     mesh = msh;

     newnodes.clear();
     newfaces.clear();

     if (partition == NULL) {
          segments[0] = anodes.size() - 1;
          segments[1] = bnodes.size() - 1;
          segments[2] = cnodes.size() - 1;
          segments[3] = dnodes.size() - 1;
          segments[4] = enodes.size() - 1;
          if (!Face::is_5_sided_convex_loop_quad_meshable(segments, partSegments)) {
               cout << "Warning: Penta patch not quad meshable " << endl;
               return 1;
          }
     } else {
          for (int i = 0; i < 10; i++)
               partSegments[i] = partition[i];
     }

     int err;
     // Ensure that input is closed loop of three sides ...
     assert(anodes.front() == enodes.back());
     assert(bnodes.front() == anodes.back());
     assert(cnodes.front() == bnodes.back());
     assert(dnodes.front() == cnodes.back());
     assert(enodes.front() == dnodes.back());

     // Ensure that segments are valid ...
     assert((int)anodes.size() == partSegments[0] + partSegments[1] + 1);
     assert((int)bnodes.size() == partSegments[2] + partSegments[3] + 1);
     assert((int)cnodes.size() == partSegments[4] + partSegments[5] + 1);
     assert((int)dnodes.size() == partSegments[6] + partSegments[7] + 1);
     assert((int)enodes.size() == partSegments[8] + partSegments[9] + 1);

     // Split Side "A" nodes into a1nodes and b1nodes.
     int a1 = partSegments[0];
     err = split_stl_vector(anodes, a1 + 1, a1nodes, b1nodes);
     if (err) return 1;

     // Split Side "B" nodes into a2nodes and b2nodes.
     int a2 = partSegments[2];
     err = split_stl_vector(bnodes, a2 + 1, a2nodes, b2nodes);
     if (err) return 1;

     // Split Side "C" nodes into a3nodes and b3nodes.
     int a3 = partSegments[4];
     err = split_stl_vector(cnodes, a3 + 1, a3nodes, b3nodes);
     if (err) return 1;

     // Split Side "C" nodes into a3nodes and b3nodes.
     int a4 = partSegments[6];
     split_stl_vector(dnodes, a4 + 1, a4nodes, b4nodes);
     if (err) return 1;

     // Split Side "C" nodes into a3nodes and b3nodes.
     int a5 = partSegments[8];
     err = split_stl_vector(enodes, a5 + 1, a5nodes, b5nodes);
     if (err) return 1;

     // Splitting nodes on each side
     Vertex *ca = a1nodes.back();
     Vertex *cb = a2nodes.back();
     Vertex *cc = a3nodes.back();
     Vertex *cd = a4nodes.back();
     Vertex *ce = a5nodes.back();

     // Centoid may not be within the patch. But take the risk.
     Vertex *co = Face::centroid(ca, cb, cc, cd, ce);

     mesh->addNode(co);
     newnodes.push_back(co);

     size_t nSize;
     linear_interpolation(mesh, co, ca, a2 + 1, oa_nodes);
     nSize = oa_nodes.size();
     for (size_t i = 1; i < nSize - 1; i++) {
          newnodes.push_back(oa_nodes[i]);
     }

     linear_interpolation(mesh, co, cb, a3 + 1, ob_nodes);
     nSize = ob_nodes.size();
     for (size_t i = 1; i < nSize - 1; i++) {
          newnodes.push_back(ob_nodes[i]);
     }

     linear_interpolation(mesh, co, cc, a4 + 1, oc_nodes);
     nSize = oc_nodes.size();
     for (size_t i = 1; i < nSize - 1; i++) {
          newnodes.push_back(oc_nodes[i]);
     }

     linear_interpolation(mesh, co, cd, a5 + 1, od_nodes);
     nSize = od_nodes.size();
     for (size_t i = 1; i < nSize - 1; i++) {
          newnodes.push_back(od_nodes[i]);
     }

     linear_interpolation(mesh, co, ce, a1 + 1, oe_nodes);
     nSize = oe_nodes.size();
     for (size_t i = 1; i < nSize - 1; i++) {
          newnodes.push_back(oe_nodes[i]);
     }

     err = quadtemplate.remesh(mesh, a1nodes, oa_nodes, oe_nodes, b5nodes);
     if (!err) addNewElements( quadtemplate.newnodes, quadtemplate.newfaces);

     if (!err) {
          err = quadtemplate.remesh(mesh, a2nodes, ob_nodes, oa_nodes, b1nodes);
          if (!err) addNewElements( quadtemplate.newnodes, quadtemplate.newfaces);
     }

     if (!err) {
          err = quadtemplate.remesh(mesh, a3nodes, oc_nodes, ob_nodes, b2nodes);
          if (!err) addNewElements( quadtemplate.newnodes, quadtemplate.newfaces);
     }

     if (!err) {
          err = quadtemplate.remesh(mesh, a4nodes, od_nodes, oc_nodes, b3nodes );
          if (!err) addNewElements( quadtemplate.newnodes, quadtemplate.newfaces);
     }

     if (!err) {
          err = quadtemplate.remesh(mesh, a5nodes, oe_nodes, od_nodes, b4nodes);
          if (!err) addNewElements( quadtemplate.newnodes, quadtemplate.newfaces);
     }

     if (err) {
          discard();
          return 1;
     }

     // Check for any inversion of the element, if there is inversion,
     // undo everthing (i.e. remove new nodes and faces).
     //
     nSize = newfaces.size();
     for (size_t i = 0; i < nSize; i++) {
          if (newfaces[i]->concaveAt() >= 0) {
               discard();
               return 2;
          }
     }

     return 0;
}

//////////////////////////////////////////////////////////////////////////////


void OneDefectPatch::setTags()
{
     size_t nSize;

     nSize = mesh->getSize(0);
     for (size_t i = 0; i < nSize; i++) {
          Vertex *v = mesh->getNodeAt(i);
          v->setTag(0);
     }

     nSize = mesh->getSize(2);
     for (size_t i = 0; i < nSize; i++) {
          Face *f = mesh->getFaceAt(i);
          f->setTag(0);
     }

     nSize = nodepath.size();
     for (size_t i = 0; i < nSize; i++)
          nodepath[i]->setTag(1);

     nSize = bound_nodes.size();
     for (size_t i = 0; i <  nSize; i++)
          bound_nodes[i]->setTag(3);

     NodeSet::const_iterator vit;
     for (vit = corners.begin(); vit != corners.end(); ++vit)
          (*vit)->setTag(2);

     FaceSet::const_iterator it;
     for (it = faces.begin(); it != faces.end(); ++it)
          (*it)->setTag(1);

}

///////////////////////////////////////////////////////////////////////////////

bool OneDefectPatch::isSafe()
{
     FaceSet::const_iterator it;
     for (it = faces.begin(); it != faces.end(); ++it)
          if ((*it)->isRemoved()) return 0;

     // Check if the previous updates modified the splitting node degree.
     if (quad_splitting_node) {
          int nSize = quad_splitting_node->getNumRelations(2);
          if ( nSize != quad_splitting_node_degree) return 0;
     }

     return 1;
}
///////////////////////////////////////////////////////////////////////////////

void OneDefectPatch::start_boundary_loop_from(Vertex *vmid)
{
     assert(corners.find(vmid) != corners.end());

     NodeSequence::iterator middle;
     middle = find(bound_nodes.begin(), bound_nodes.end(), vmid);
     assert(middle != bound_nodes.end());

     std::rotate(bound_nodes.begin(), middle, bound_nodes.end());
     assert(bound_nodes[0] == vmid);

     set_boundary_segments();
}
///////////////////////////////////////////////////////////////////////////////

void OneDefectPatch::get_bound_nodes(const Vertex *src, const Vertex *dst, NodeSequence &seq)
{
     int start_pos = getPosOf(src);
     int end_pos = getPosOf(dst);
     int nsize = bound_nodes.size();

     assert(nsize > 1);

     if (end_pos == 0) end_pos = nsize;
     assert(end_pos > start_pos);

     seq.resize(end_pos - start_pos + 1);
     int index = 0;
     for (int i = start_pos; i <= end_pos; i++)
          seq[index++] = bound_nodes[i % nsize];
}

///////////////////////////////////////////////////////////////////////////////

bool OneDefectPatch::has_irregular_node_on_first_segment() const
{
     int start_pos = cornerPos[0];
     int end_pos = start_pos + segSize[0] - 1;

     for (int i = 1; i < end_pos - 1; i++)
          if( bound_nodes[i]->getNumRelations(2) < 4 ) return 1;

     return 0;
}

///////////////////////////////////////////////////////////////////////////////

bool OneDefectPatch::is_simply_connected()
{
     size_t V = inner_nodes.size() + bound_nodes.size();
     size_t F = faces.size();

     vector<Edge> edges;
     edges.reserve(4 * F);

     FaceSet::const_iterator it;
     for (it = faces.begin(); it != faces.end(); ++it) {
          Face *face = *it;
          for (int i = 0; i < 4; i++) {
               Vertex *v0 = face->getNodeAt(i);
               Vertex *v1 = face->getNodeAt(i + 1);
               Edge edge(v0, v1);
               int found = 0;
               size_t nSize =  edges.size();
               for (size_t j = 0; j < nSize; j++) {
                    if (edges[j].isSameAs(edge)) {
                         found = 1;
                         break;
                    }
               }
               if (!found) edges.push_back(edge);
          }
     }

     size_t E = edges.size();

     if (V - E + F == 1) return 1;
     return 0;
}
///////////////////////////////////////////////////////////////////////////////

int OneDefectPatch::reorient_4_sided_loop()
{
     //Always remeshable, nothing has to be done...
     if ((segSize[0] == segSize[2]) && (segSize[1] == segSize[3])) return 0;

     //////////////////////////////////////////////////////////////////////////
     // Defination:  A four sided convex loop has four segments.
     // Objectives:  A four sided convex loop must be orietned such that
     //   1.  First segment must be smaller than the third one, because
     //       we need to create a triangle patch based at segment#3.
     //
     //   2.  If there are two choices, then the side having irregular
     //       node must be given higher priority. ( Here irregulaty means
     //       that vertex valency < 4 ).
     //
     //   3.  A side having less number of nodes on the first segment than
     //       the third is given preference.
     //
     // Pre-Conditions  :  A loop must be oriented ( CW or CCW ).
     //
     // Date: 17th Nov. 2010.
     //////////////////////////////////////////////////////////////////////////

     Vertex *start_corner = NULL;

     if (segSize[0] == segSize[2]) {
          if (min(segSize[1], segSize[3]) == 2) return 1;
          //  Either Segment 2 or 3 must be starting node
          if (segSize[1] < segSize[3])
               start_corner = bound_nodes[ cornerPos[1] ];
          else
               start_corner = bound_nodes[ cornerPos[3] ];
          start_boundary_loop_from(start_corner);
     }

     if (min(segSize[0], segSize[2]) == 2) return 1;

     if (segSize[2] < segSize[0]) {
          start_corner = bound_nodes[ cornerPos[2] ];
          start_boundary_loop_from(start_corner);
     }

     // By this stage, the loop must be reoriented correctly.
     assert(segSize[0] < segSize[2]);

     // Great, we found one irregular node on the first boundary...
     if (has_irregular_node_on_first_segment()) return 1;

     // If the segment 2 and 3 have same size, Alas, nothing can be done.
     if (segSize[1] == segSize[3]) return 1;

     if (min(segSize[1], segSize[3]) == 2) return 1;

     if (segSize[3] < segSize[1]) {
          start_corner = bound_nodes[ cornerPos[3] ];
          start_boundary_loop_from(start_corner);
     } else {
          start_corner = bound_nodes[ cornerPos[1] ];
          start_boundary_loop_from(start_corner);
     }

     //
     // Note that we didn't check for irregular node here. So if this segment
     // has at least one irregular node, then we are lucky. Otherwise decision
     // to remesh it done based wthere remeshing will result in the reduction
     // of irregular nodes in patch.
     //
     assert(segSize[0] < segSize[2]);
     return 0;
}

//////////////////////////////////////////////////////////////////////////////

int OneDefectPatch::remesh_3_sided_patch()
{
     Vertex *c0, *c1, *c2;

     c0 = bound_nodes[ cornerPos[0] ];
     c1 = bound_nodes[ cornerPos[1] ];
     c2 = bound_nodes[ cornerPos[2] ];

     get_bound_nodes(c0, c1, anodes);
     get_bound_nodes(c1, c2, bnodes);
     get_bound_nodes(c2, c0, cnodes);

     int err = template3.remesh(mesh, anodes, bnodes, cnodes, partSegments);
     if( !err ) num_3_patches++;
     return err;
}
//////////////////////////////////////////////////////////////////////////////

int OneDefectPatch::remesh_4_sided_patch()
{
     int err;
     size_t nSize;
     Vertex *c0, *c1, *c2, *c3, *ta, *tb, *tc;

     // Collect landmark nodes on the patch..
     c0 = bound_nodes[ cornerPos[0] ];
     c1 = bound_nodes[ cornerPos[1] ];
     c2 = bound_nodes[ cornerPos[2] ];
     c3 = bound_nodes[ cornerPos[3] ];
     ta = quad_splitting_node;

     if (ta == NULL) {
          if ((segSize[0] == segSize[2]) && (segSize[1] == segSize[3])) {
               get_bound_nodes(c0, c1, anodes);
               get_bound_nodes(c1, c2, bnodes);
               get_bound_nodes(c2, c3, cnodes);
               get_bound_nodes(c3, c0, dnodes);
               err = template4.remesh(mesh, anodes, bnodes, cnodes, dnodes);
               return err;
          }
     }

     assert(ta != NULL);

     get_bound_nodes(c0, ta, a1nodes);
     get_bound_nodes(ta, c1, a2nodes);

     tb = bound_nodes[cornerPos[2] + a2nodes.size() - 1 ];
     tc = bound_nodes[cornerPos[3] - a1nodes.size() + 1];
     get_bound_nodes(tb, tc, c0nodes);

     linear_interpolation(mesh, ta, tb, segSize[1], abnodes);

     nSize = abnodes.size();
     for (size_t  i = 1; i < nSize - 1; i++)
          newnodes.push_back(abnodes[i]);

     linear_interpolation(mesh, tc, ta, segSize[3], canodes);

     nSize = canodes.size();
     for (size_t i = 1; i < nSize - 1; i++)
          newnodes.push_back(canodes[i]);

     get_bound_nodes(tb, tc, bcnodes);
     get_bound_nodes(c1, c2, b1nodes);
     get_bound_nodes(c2, tb, c2nodes);
     get_bound_nodes(tc, c3, c1nodes);
     get_bound_nodes(c3, c0, d1nodes);

     err = template3.remesh(mesh, abnodes, bcnodes, canodes, partSegments);
     if (!err) {
          for( size_t i = 0; i < template3.newnodes.size(); i++)
               newnodes.push_back( template3.newnodes[i] );
          for( size_t i = 0; i < template3.newfaces.size(); i++)
               newfaces.push_back( template3.newfaces[i] );
     }

     if (!err) {
          err = template4.remesh(mesh, a1nodes, canodes, c1nodes, d1nodes);
          if (!err) {
               for( size_t i = 0; i < template4.newnodes.size(); i++)
                    newnodes.push_back( template4.newnodes[i] );
               for( size_t i = 0; i < template4.newfaces.size(); i++)
                    newfaces.push_back( template4.newfaces[i] );
          }
     }

     if (!err) {
          err = template4.remesh(mesh, a2nodes, b1nodes, c2nodes, abnodes);
          if (!err) {
               for( size_t i = 0; i < template4.newnodes.size(); i++)
                    newnodes.push_back( template4.newnodes[i] );
               for( size_t i = 0; i < template4.newfaces.size(); i++)
                    newfaces.push_back( template4.newfaces[i] );
          }
     }

     if (err) {
          nSize = newfaces.size();
          for (size_t i = 0; i < nSize; i++)
               mesh->remove(newfaces[i]);

          nSize = newnodes.size();
          for (size_t i = 0; i < nSize; i++)
               mesh->remove(newnodes[i]);

          newnodes.clear();
          newfaces.clear();
          return 2;
     }

     if( !err ) num_4_patches++;

     return err;
}

//////////////////////////////////////////////////////////////////////////////

int OneDefectPatch::remesh_5_sided_patch()
{
     Vertex *c0, *c1, *c2, *c3, *c4;
     c0 = bound_nodes[ cornerPos[0] ];
     c1 = bound_nodes[ cornerPos[1] ];
     c2 = bound_nodes[ cornerPos[2] ];
     c3 = bound_nodes[ cornerPos[3] ];
     c4 = bound_nodes[ cornerPos[4] ];

     get_bound_nodes(c0, c1, anodes);
     get_bound_nodes(c1, c2, bnodes);
     get_bound_nodes(c2, c3, cnodes);
     get_bound_nodes(c3, c4, dnodes);
     get_bound_nodes(c4, c0, enodes);

     int err = template5.remesh(mesh, anodes, bnodes, cnodes, dnodes, enodes, partSegments);

     if( !err) num_5_patches++;

     return err;
}

//////////////////////////////////////////////////////////////////////////////

bool OneDefectPatch::is_quad_breakable_at(const Vertex *ta)
{
     // Landmark nodes on the boundary. There are seven such nodes.
     Vertex *c0, *c1, *c2, *c3, *tb, *tc;
     c0 = bound_nodes[ cornerPos[0] ];
     c1 = bound_nodes[ cornerPos[1] ];
     c2 = bound_nodes[ cornerPos[2] ];
     c3 = bound_nodes[ cornerPos[3] ];

     get_bound_nodes(c0, ta, a1nodes);
     get_bound_nodes(ta, c1, a2nodes);

     tb = bound_nodes[cornerPos[2] + a2nodes.size() - 1 ];
     tc = bound_nodes[cornerPos[3] - a1nodes.size() + 1];
     get_bound_nodes(tb, tc, c0nodes);

     int segments[3];
     segments[0] = c0nodes.size() - 1;
     segments[1] = segSize[3] - 1;
     segments[2] = segSize[2] - 1;

     if (Face::is_3_sided_convex_loop_quad_meshable(segments, partSegments))
          return 0;

     return 1;
}
//////////////////////////////////////////////////////////////////////////////

bool OneDefectPatch::is_4_sided_convex_loop_quad_meshable()
{
     assert(corners.size() == 4);

     // It is always quadrangulable.
     if ((segSize[0] == segSize[2]) && (segSize[1] == segSize[3])) return 1;

     int nr0 = count_irregular_nodes(0);
     int nr1 = count_irregular_nodes(1);

     //
     // If there are only two irregular nodes in the patch ( inner + boundary)
     // then it may be quadrangable, but it will not reduce number of irregular
     // nodes ( one irregular node and one irregular node on the boundary
     // will be created. Since there is no reduction in #of irregular nodes,
     // it is not profitable, and we just return it.
     //
     // If there are more than two irrgular nodes, then we can ensure that
     // that there will at most two irregular nodes in the patch ( one inside
     // and probably one on the boundary.
     //

     if (nr0 + nr1 == 2) return 0;

     // The loop is reoriented so that the first edge is smaller than the
     // opposite edge ....
     if (reorient_4_sided_loop()) return 0;

     // We need at least 3 points on the first seqment ...
     if (segSize[0] < 3) return 0;

     // Landmark nodes on the boundary. There are seven such nodes.
     Vertex *c0, *c1, *c2, *c3, *ta, *tb, *tc;
     c0 = bound_nodes[ cornerPos[0] ];
     c1 = bound_nodes[ cornerPos[1] ];
     c2 = bound_nodes[ cornerPos[2] ];
     c3 = bound_nodes[ cornerPos[3] ];

     // First let us prioritieze vertices on the first segment that will
     // be tried for splitting the region. First of all, if #of irregular
     // nodes inside the region is 2, then we must have one irregular
     // node on the boundary. But even if the #irregular nodes inside
     // the region > 2, we still prefer irregular node in the boundary
     // as it will reduce irregularity by one.

     // Both end nodes on the segment are not included in splitting.

     // Method: Start node(1,...n-1) and check it could lead to acceptable
     // splitting. A proper splitting means that the triangle sandwitched
     // between two quads regions is "quad-meshable" with one defect.
     //
     // Future Work: It is possible that the triangle region is not "one
     // defect" meshable, but with more than one defect. We will extend
     // this case later.
     //

     NodeSequence trynodes;
     for (int i = 1; i < segSize[0] - 1; i++) {
          if (bound_nodes[i]->getNumRelations(2) == 3)
               trynodes.push_back(bound_nodes[i]);
     }

     if (trynodes.empty() && nr0 == 2) {
          // Patch not remeshable because there are no irregular node on the boundary
          return 0;
     }

     if (nr0 > 2) {
          for (int i = 1; i < segSize[0] - 1; i++) {
               if (bound_nodes[i]->getNumRelations(2) == 4)
                    trynodes.push_back(bound_nodes[i]);
          }
     }

     quad_splitting_node = NULL;

     int segments[3];
     int nSize = trynodes.size();
     for (int i = 0; i < nSize; i++) {
          ta = trynodes[i];
          get_bound_nodes(c0, ta, a1nodes);
          get_bound_nodes(ta, c1, a2nodes);

          tb = bound_nodes[cornerPos[2] + a2nodes.size() - 1 ];
          tc = bound_nodes[cornerPos[3] - a1nodes.size() + 1];

          get_bound_nodes(tb, tc, bcnodes);
          segments[0] = segSize[1] - 1;
          segments[1] = bcnodes.size() - 1;
          segments[2] = segSize[3] - 1;

          if (Face::is_3_sided_convex_loop_quad_meshable(segments, partSegments)) {
               quad_splitting_node = ta;
               quad_splitting_node_degree = ta->getNumRelations(2);
               return 1;
          }
     }

     return 0;
}

///////////////////////////////////////////////////////////////////////////////

void OneDefectPatch::init_blob()
{
     nodes.clear();
     inner_nodes.clear();
     bound_nodes.clear();

     faces.clear();
     inner_faces.clear();

     relations02.clear();

     size_t nSize;
     if (apex) {
          expand_blob( apex );
     } else {
          nSize = nodepath.size();
          for (size_t i = 0; i < nSize; i++) {
               Vertex *v = nodepath[i];
               expand_blob(v );
          }
     }
     create_boundary();
}

////////////////////////////////////////////////////////////////////

size_t OneDefectPatch::count_irregular_nodes(int where)
{
     NodeSequence::const_iterator it;
     assert(mesh->getAdjTable(0, 2));

     size_t ncount = 0;

     if (where == 0) {
          for (it = inner_nodes.begin(); it != inner_nodes.end(); ++it) {
               Vertex *v = *it;
               if( !v->isRemoved() && (v->getNumRelations(2) != 4)) ncount++;
          }
     } else {
          for (size_t i = 0; i < bound_nodes.size(); i++) {
               Vertex *v = bound_nodes[i];
               if( !v->isRemoved() )
                    if (!v->isBoundary() && v->getNumRelations(2) == 3) ncount++;
          }
     }

     return ncount;
}
////////////////////////////////////////////////////////////////////

int
OneDefectPatch::create_boundary()
{
     num_boundaries++;
     StopWatch wc;

     corners.clear();
     boundary.clear();

     wc.start();

     // We need to rebuild relations locally to identfy corners and boundary.
     FaceSet::const_iterator fiter;

     // A boundary edge must have exactly one face neighbor...
     FaceSequence faceneighs;
     for (fiter = faces.begin(); fiter != faces.end(); ++fiter) {
          Face *face = *fiter;
          if( inner_faces.find(face) == inner_faces.end() ) {
               int ncount = 0;
               for (int j = 0; j < 4; j++) {
                    Vertex *v0 = face->getNodeAt(j + 0);
                    Vertex *v1 = face->getNodeAt(j + 1);
                    faceneighs.clear();
                    assert(relations02[v0].size() > 0);
                    assert(relations02[v1].size() > 0);
                    set_intersection(relations02[v0].begin(), relations02[v0].end(),
                                     relations02[v1].begin(), relations02[v1].end(),
                                     back_inserter(faceneighs));
                    if (faceneighs.size() == 1) {
                         ncount = 1;
                         Edge newedge(v0, v1);
                         boundary.push_back(newedge);
                    }
               }
               if( ncount == 0) inner_faces.insert(face);
          }
     }

     wc.stop();
     exec_time += wc.getSeconds();

     // Sequence the chain and start from one of the corner...
     int err = Mesh::make_chain(boundary);
     if (err) return 2;

     size_t nSize = boundary.size();

     Vertex *vertex;
     for (size_t k = 0; k < nSize; k++) {
          vertex = boundary[k].getNodeAt(0);
          if (relations02[vertex].size() == 1) corners.insert(vertex);

          vertex = boundary[k].getNodeAt(1);
          if (relations02[vertex].size() == 1) corners.insert(vertex);
     }

     if( corners.size() < 3) return 3;

     // Start the chain from one of the corners.
     err = Mesh::rotate_chain(boundary, *corners.begin());
     if (err) return 4;

     bound_nodes.resize( nSize );
     for (size_t k = 0; k < nSize; k++)
          bound_nodes[k] = boundary[k].getNodeAt(0);
     sort( bound_nodes.begin(), bound_nodes.end() );
     inner_nodes.clear();
     set_difference(nodes.begin(), nodes.end(),
                    bound_nodes.begin(), bound_nodes.end(), back_inserter(inner_nodes));

     // We need nodes in cycle, so collect the sequence of boundary nodes...,
     bound_nodes.resize( nSize );
     for (size_t k = 0; k < nSize; k++)
          bound_nodes[k] = boundary[k].getNodeAt(0); // Only the first node.

     // Split the boundary loop into segments.
     // (i.e. End of the segments are the corners identified earlier )
     set_boundary_segments();

     return 0;
}

////////////////////////////////////////////////////////////////////

void OneDefectPatch::set_boundary_segments()
{
     size_t nSize = corners.size();
     // Although this stage will not come in this algorithm...
     if ( nSize == 0) return;

     cornerPos.resize( nSize + 1);

     NodeSet::const_iterator it;
     int index = 0;
     for (it = corners.begin(); it != corners.end(); ++it) {
          cornerPos[index++] = getPosOf(*it);
     }

     cornerPos[corners.size()] = bound_nodes.size();

     sort(cornerPos.begin(), cornerPos.end());

     segSize.resize( nSize );

     for (size_t i = 0; i < nSize; i++)
          segSize[i] = cornerPos[(i + 1)] - cornerPos[i] + 1;
}

////////////////////////////////////////////////////////////////////

int OneDefectPatch::get_topological_outer_angle(Vertex *vertex)
{
     const FaceSequence &vfaces = vertex->getRelations2();

     // How many faces are outside the regions.
     int ncount = 0;
     int nSize  = vfaces.size();
     for (int i = 0; i < nSize; i++)
          if (faces.find(vfaces[i]) == faces.end()) ncount++;

     return ncount - 2;
}
////////////////////////////////////////////////////////////////////////////////

void OneDefectPatch::expand_blob(Vertex *vertex)
{
     const FaceSequence &vfaces = vertex->getRelations2();

     int nSize = vfaces.size();
     for (int i = 0; i < nSize; i++) {
          Face *face = vfaces[i];
          assert(!face->isRemoved() );
          faces.insert(face);
          for (int j = 0; j < 4; j++) {
               Vertex *vf = face->getNodeAt(j);
               relations02[vf].insert(face);
               nodes.insert( vf );
          }
     }
}
////////////////////////////////////////////////////////////////////////////////

int
OneDefectPatch::build_remeshable_boundary()
{
     assert(mesh->getAdjTable(0, 2));

     /////////////////////////////////////////////////////////////////////////////
     // A valid boundary patch has the following propeties:
     //
     // 1) The loop is closed.
     // 2) The region is simply connected. The Euler characteristic must be
     //    one, which mean there are no holes in the region.
     // 3) The loop has 3, 4 or 5 corner nodes.
     // 4) The loop form almost convex boundary.
     //
     // If these conditions are met, we may try to remesh the region with
     // quadrilateral elements. Sometimes it may not be possible to remesh
     // a reason, or the resulting elements have unacceptable quality that
     // we would like to avoid.
     //
     // A resulting mesh MUST decrease the irregular nodes count.
     //////////////////////////////////////////////////////////////////////////
     init_blob();

     // There are less than two irregular nodes in the patch, so return..
     if (count_irregular_nodes(0) < 2) return 1;

     int err;
     size_t  nSize;
     while (1) {
          err = create_boundary();
          if (err) return 2;

          if (faces.size() > MAX_FACES_ALLOWED ) return 3;

          // Expand the concave blob till convex patch is formed...
          int topo_convex_region = 1;
          nSize = bound_nodes.size();
          for (size_t i = 0; i < nSize; i++) {
               if (!bound_nodes[i]->isBoundary()) {
                    int topo_angle = get_topological_outer_angle(bound_nodes[i]);
                    if (topo_angle < 0) {
                         expand_blob(bound_nodes[i]);
                         topo_convex_region = 0;
                    }
               }
          }

          // If convex patch is formed, check it is remeshable ...
          if (topo_convex_region) {
               int nsides = corners.size();

               if (nsides > 5) return 6;

               int segments[5];

               int meshable = 0;
               switch (nsides) {
               case 3:
                    segments[0] = segSize[0] - 1;
                    segments[1] = segSize[1] - 1;
                    segments[2] = segSize[2] - 1;
                    meshable = Face::is_3_sided_convex_loop_quad_meshable(segments, partSegments);
                    break;

               case 4:
                    meshable = is_4_sided_convex_loop_quad_meshable();
                    break;

               case 5:
                    segments[0] = segSize[0] - 1;
                    segments[1] = segSize[1] - 1;
                    segments[2] = segSize[2] - 1;
                    segments[3] = segSize[3] - 1;
                    segments[4] = segSize[4] - 1;
                    meshable = Face::is_5_sided_convex_loop_quad_meshable(segments, partSegments);
                    break;
               }

               // If the region is meshable, it is good, it will be remeshed in the next step...
               if (meshable) {
                    if (!is_simply_connected()) return 4;
                    return 0;
               }

               // Otherwise, expand the convex region further from the boundary...
               size_t before_expansion = faces.size();

               nSize = bound_nodes.size();
               for (size_t i = 0; i < nSize; i++)
                    expand_blob(bound_nodes[i]);

               size_t after_expansion = faces.size();

               // bad luck; no expansion took place....
               if( after_expansion == before_expansion ) return 1;
          }
     }

     // You should never come at this stage...
     return 1;
}

///////////////////////////////////////////////////////////////////////////////

void OneDefectPatch::rollback()
{
     NodeSequence::const_iterator it;
     for (it = inner_nodes.begin(); it != inner_nodes.end(); ++it) 
          mesh->reactivate(*it);

     FaceSet::const_iterator fiter;
     for (fiter = faces.begin(); fiter != faces.end(); ++fiter)
          mesh->reactivate(*fiter);

     size_t nSize = newfaces.size();
     for (size_t i = 0; i < nSize; i++)
          mesh->remove(newfaces[i]);
     newfaces.clear();

     nSize = newnodes.size();
     for (size_t i = 0; i < nSize; i++)
          mesh->remove(newnodes[i]);
     newnodes.clear();

}

////////////////////////////////////////////////////////////////////////////////

void OneDefectPatch::pre_remesh()
{
     //
     // Before remeshing deactivate all the internal faces and nodes so that all
     // the relations are cleared. They are not removed because if something goes
     // wrong, the data structures are restored using rollback operations.
     //
     FaceSet::const_iterator fiter;
     for (fiter = faces.begin(); fiter != faces.end(); ++fiter)
          mesh->deactivate(*fiter);

     irregular_nodes_removed.clear();
     NodeSequence::const_iterator niter;
     for (niter = inner_nodes.begin(); niter != inner_nodes.end(); ++niter) {
          if (!QuadCleanUp::isRegular(*niter)) irregular_nodes_removed.push_back(*niter);
          mesh->deactivate(*niter);
     }

     // At this stage, we will have an empty patch ...
}

////////////////////////////////////////////////////////////////////////////////

void OneDefectPatch::post_remesh()
{
     size_t nSize;

     // Remove the old faces ( will go to trash bin );
     FaceSet::const_iterator fiter;
     for (fiter = faces.begin(); fiter != faces.end(); ++fiter)
          mesh->remove(*fiter);

     // Remove the old nodes.. (will go to trash bin );
     NodeSequence::const_iterator niter;
     for (niter = inner_nodes.begin(); niter != inner_nodes.end(); ++niter)
          mesh->remove(*niter);

     // Patch contains new nodes new faces now ...

     faces.clear();
     nSize = newfaces.size();
     for (size_t i = 0; i < nSize; i++) {
          faces.insert(newfaces[i]);
     }
     newfaces.clear();

     nSize = newnodes.size();
     inner_nodes.clear();
     for (size_t i = 0; i < nSize; i++)
          inner_nodes.push_back(newnodes[i]);
     newnodes.clear();

     int ncount = 0;
     new_defective_node = NULL;
     for (niter = inner_nodes.begin(); niter != inner_nodes.end(); ++niter) {
          Vertex *v = *niter;
          if (!QuadCleanUp::isRegular(v)) {
               ncount++;
               new_defective_node = v;
          }
     }

     if( ncount > 2) {
          cout << "Fatal Error:There must be at the most one defective node in the patch" << endl;
          exit(0);
     }
}

////////////////////////////////////////////////////////////////////////////////

int OneDefectPatch::remesh()
{
     // Check if any node or element is modifed by previous operations.
     if (!isSafe()) return 1;

     newnodes.clear();
     newfaces.clear();

     if (corners.size() > 5) cout << "Warning: Corners more than 5 " << endl;

#ifdef DEBUG
     int nirregular0 = count_irregular_nodes(0) + count_irregular_nodes(1);
#endif

     // Do some pre-processing before remeshing the patch
     pre_remesh();

     int err = 1;
     switch (corners.size()) {
     case 3:
          err = remesh_3_sided_patch();
          break;
     case 4:
          err = remesh_4_sided_patch();
          break;
     case 5:
          err = remesh_5_sided_patch();
          break;
     }

     if (err) {
          rollback();
          return 1;
     }

     // Do some post-processing of the patch
     post_remesh();

#ifdef DEBUG
     int nirregular1 = count_irregular_nodes(0) + count_irregular_nodes(1);

     // Our method must stricly decrease the number of irregular nodes,
     // otherwise it will be go into infinite loop.

     assert(nirregular1 < nirregular0);
#endif

     return err;
}

////////////////////////////////////////////////////////////////////////////////

OneDefectPatch *QuadCleanUp::build_one_defect_patch(Vertex *vertex)
{
     if( defective_patch == NULL )
          defective_patch = new OneDefectPatch(mesh);

     defective_patch->clear();

     if( vertex->isRemoved() || isRegular(vertex) ) return NULL;

     if( djkpath == NULL ) {
          djkpath = new DijkstraShortestPath(mesh, 1);
     }

     MeshFilter *filter = new FirstIrregularNode();
     const NodeSequence &nodepath = djkpath->getPath(vertex, filter);
     delete filter;

     if (isRegular(nodepath.front())) return NULL;
     if (isRegular(nodepath.back()))  return NULL;

     defective_patch->set_initial_path( nodepath );

     int err = defective_patch->build_remeshable_boundary();

     if (err) return NULL;

     return defective_patch;
}

////////////////////////////////////////////////////////////////////////////////

int QuadCleanUp::remesh_defective_patches()
{
     int relexist0 = mesh->build_relations(0, 0);
     int relexist2 = mesh->build_relations(0, 2);

     mesh->search_boundary();

     if( !mesh->is_consistently_oriented() ) 
          mesh->make_consistently_oriented();

     djkpath = new DijkstraShortestPath(mesh, 1);

     assert( djkpath);

     assert(mesh->getAdjTable(0, 2));
     OneDefectPatch *patch;

     Jaal::MeshOptimization mopt;

     NodeSequence nextSeq;
     irregular_nodes = mesh->get_irregular_nodes(4);

     // Record execution time ...
     StopWatch wc[4];

     int ncount1 = 0;
     for( int i = 4; i < 14; i++) {
          OneDefectPatch :: MAX_FACES_ALLOWED = ( int)pow(2,i);
          cout << "Remeshing with maximum patch size :  " << (int)pow(2,i) << endl;
          cout << "#of irregular nodes in the mesh " << irregular_nodes.size() << endl;

          wc[3].start();

          while (1) {
               nextSeq.clear();
               size_t ncount2 = 0;
               size_t nSize = irregular_nodes.size();
               for( size_t i = 0; i < nSize; i++) {
                    Vertex *vertex = irregular_nodes[i];
                    if( !vertex->isRemoved() ) {
                         wc[0].start();
                         patch = build_one_defect_patch(vertex);
                         wc[0].stop();
                         if (patch) {
                              wc[1].start();
                              int err = patch->remesh();
                              wc[1].stop();
                              if (!err) {
                                   if( ncount1% 100 == 0) cout << "patch count : " << ncount1 << endl;
                                   Vertex *vtx = patch->get_new_defective_node();
                                   if (vtx) nextSeq.push_back(vtx);
                                   ncount2++;
                                   ncount1++;
                              } else
                                   nextSeq.push_back( vertex );  // Because patch was not remeshed
                         } else
                              nextSeq.push_back( vertex ); // Because invalid patch occured.
                    }
               }
               irregular_nodes.swap(nextSeq);

               if (ncount2) {
                    wc[2].start();
//                mopt.shape_optimize(mesh);
                    wc[2].stop();
               }
               if( ncount2 == 0) break;
          }
          wc[3].stop();
          cout << "Execution Summary: " << endl;
          cout << "    Time for searching patches  : " << wc[0].getSeconds()  << endl;
          cout << "    Time for remeshing patches  : " << wc[1].getSeconds()  << endl;
          cout << "    Time for shape optimization : " << wc[2].getSeconds()  << endl;
          cout << "    Total Execution time        : " << wc[3].getSeconds() << endl;
     }

     cout << "#of times boundaries created : " << OneDefectPatch::num_boundaries << endl;
     cout << "Exec time boundaries created : " << OneDefectPatch::exec_time      << endl;
     cout << "3 Sided Patches : " << OneDefectPatch::num_3_patches << endl;
     cout << "4 Sided Patches : " << OneDefectPatch::num_4_patches << endl;
     cout << "5 Sided Patches : " << OneDefectPatch::num_5_patches << endl;

     cout << "Execution Summary: " << endl;
     cout << "    Time for searching patches  : " << wc[0].getSeconds()  << endl;
     cout << "    Time for remeshing patches  : " << wc[1].getSeconds()  << endl;
     cout << "    Time for shape optimization : " << wc[2].getSeconds()  << endl;

     if (!relexist0)
          mesh->clear_relations(0, 0);

     if (!relexist2)
          mesh->clear_relations(0, 2);

     return 0;
}

/////////////////////////////////////////////////////////////////////////
