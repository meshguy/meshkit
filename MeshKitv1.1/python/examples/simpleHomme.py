from MeshKit import *
from itaps import iBase
import math

# this is the length of the side of the enclosed cube
length = 10.
# each edge of the cube will be divided using this meshcount
meshcount = 9

# circumscribed sphere radius
radius = length * math.sqrt(3) /2

mk = MKCore()
geom = mk.igeom_instance()

geom.createBrick([length, length, length])
mk.populate_model_ents(0,0,0)

vols = MEntVector()
mk.get_entities_by_dimension(3, vols)

surfs = MEntVector()
mk.get_entities_by_dimension(2, surfs)

edges= MEntVector()
mk.get_entities_by_dimension(1, edges)
esize = SizingFunction(mk, meshcount, -1)
edm = mk.construct_meshop("EdgeMesher", edges)
# scheme number 5 is equignomonic
edm.set_edge_scheme(5)
tfi=mk.construct_meshop("TFIMapping", surfs)
tfi.skip_improve()

vols[0].sizing_function_index(esize.core_index())

mk.setup_and_execute()

mesh_file="../data/cubeSurface.h5m"
mk.save_mesh(mesh_file)

mesh = mk.imesh_instance()

nodes = mesh.getEntities(iBase.Type.vertex)

# project each node in the mesh on the sphere
for node in nodes:
  x, y, z = mesh.getVtxCoords(node)
  dist1 = math.sqrt( x*x + y*y + z*z )
  ratio = radius/dist1
  x1 = x*ratio
  y1 = y*ratio
  z1 = z*ratio
  mesh.setVtxCoords(node, [x1, y1, z1])

mesh_file2="../data/Homme1.h5m"
mk.save_mesh(mesh_file2)
  

