/*! 
    \page thegraph A Graph-Based Model for Mesh Generation
 
%MeshKit is designed around a graph-based view of the meshing process.  What is a graph-based view?
If you're used to a meshing process built around BREP-based geometry (e.g. like the one used
in CUBIT and other solid modeler-based meshing tools), you're probably used to the idea of meshing BREP
entities in increasing dimension, i.e. the 
geometric vertices, then edges, then faces, then regions.  In this case, the graph is the topology
graph, with graph nodes corresponding to geometric entities, and graph edges the adjacency relationships
between entities.  When the BREP is the graph, each meshing operation is synonymous with the geometric 
tity it operates on.  However, this is only a special case of a graph-based process; more generally, each 
meshing operation can be a node in the graph, with 
graph edges representing dependencies between operations.  For example, in the case of meshing
then refinement, the graph would have two nodes, one for generating the initial
mesh, and the other for refining that mesh.  The (single) edge linking
the two nodes represents the dependency between these operations.  In other words, you can't perform
the refinement operation until you have something to refine.  A graphical representation of this graph is:
 
 \image html meshnrefine.png "The graph for a 2-step mesh-then-refine operation, along with the root and leaf nodes of the graph."

%MeshKit represents the meshing process as a single directed, acyclic graph (DAG).  The graph
always has a single root and leaf node; this is for convenience, so that traversals 
always begin at one node and end at another.  When a meshing operation is instantiated in %MeshKit, it
is automatically made a child of the root node and a parent of the leaf node.  The
node can then be inserted elsewhere in the graph, making it a child/parent of other nodes (often disconnecting
the node from the root and/or leaf nodes).  Any meshing operation (implemented in a class derived from MeshOp) 
can be a node in the graph.

Using the graph, the meshing process is implemented as two phases: a "setup" phase, 
when the graph is traversed in the reverse direction
from the leaf to the root node; and an "execute" phase, when the graph is traversed in the forward direction,
from the root to the leaf node.  During the setup phase, a meshing operation examines the model entity(ies) it
is assigned to, possibly creating other meshing operations on which it depends.  These new meshing operations
are inserted in the graph, upstream of the meshing operation / graph node being examined.  Since the new nodes
are inserted upstream, that means the current node depends on those operations.  The Setup phase traverses
up the graph in a "topological sort" order, in the process visiting any new graph nodes (and any new nodes those 
nodes produce), until the root node is visited.  After the Setup phase is completed, the graph is traversed 
in the forward direction, again in a topological sort order, starting at the root node and ending at the leaf
node.  This phase is called the Execute phase, where mesh is generated (or operated on, depending on the
type of meshing operation in each graph node).  Traversing in a topological sorted order ensures that all
dependencies for a given node have been executed before that node is visited.

Using a graph-based approach for the meshing process has several advantages over the more traditional 
BREP-based method: 
- It is easy to make local changes to the meshing process and to determine exactly which other steps of the
meshing procedure need to be re-executed.  Namely, if one graph node is changed (usually because input affecting
the meshing operation is changed), that node and all descendents in the graph must be re-executed.
- Parallelism is represented explicitly in the graph, and can be exploited by executing independent nodes
concurrently.
- Conflicting input in the meshing specification appears explicitly as cycles in the (directed) graph, making
them easy to detect and avoid.
- The resulting graphs are much less complex than those representing the BREP, while still capturing the
relevant dependencies between meshing operations.

A more detailed example of graph traversal is given \subpage detailedtraversal "here".

\page detailedtraversal A more detailed example of graph traversal

In this example, an example of meshing graph traversal is given where graph nodes create other graph nodes
dynamically during the setup phase.  Consider the generation of a mesh for a surface bounded by two geometric 
edges and vertices:

 \image html detailedmesh.png "Geometry and mesh for a more detailed example; a surface bounded by 2 geometric vertices and 2 edges."

For the starting mesh graph, the user specifies a "trimesh" meshing operation, assigned to the surface:

\image html detailedgraph1.png "User-specified meshing graph, with a single \"trimesh\" node."

Starting the setup traversal, the trimesh operation determines that it needs a mesh for the bounding edges, and
creates a graph node representing that operation (an EqualEdgeMesh meshing operation applied to edges e1 and e2).  
Continuing the traversal, the setup operation on the EqualEdgeMesh node determines that the two bounding vertices
must be meshed; this is accomplished by creating a VertexMesher graph node, applied to vertices v1 and v2.
The setup operation on the VertexMesher node creates no new graph nodes, so the setup traversal terminates at
the root node.  The graph resulting from the setup traversal appears as:

\image html detailedgraph2.png "The mesh graph after setup traversal is completed."

with user-specified nodes colored cyan, and automatically-created nodes colored magenta.  Note that a given 
node can represent a meshing operation applied to multiple entities in the BREP model.  In practice, this greatly
simplifies the mesh graph, with no loss of detail in the overall meshing process.

The execute phase traverses the graph in the forward direction, generating mesh for the geometric vertices 
in the first (non-root) node, for geometric edges in the second node, and finally the surface mesh in the
third node:

\image html detailedgraph3.png "The mesh graph and resulting mesh during execute traversal."

Top: \ref index Prev: \ref datamodel Next: \ref interfaces

*/
