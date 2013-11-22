from MeshKit import *
from itaps import iBase
import math
import random

def cart2sph(x,y,z):
    XsqPlusYsq = x**2 + y**2
    r = math.sqrt(XsqPlusYsq + z**2)               # r
    elev = math.atan2(z,math.sqrt(XsqPlusYsq))     # theta
    az = math.atan2(y,x)                           # phi
    return r, elev, az

# this is the length of the side of the enclosed cube
length = 6.
# each edge of the cube will be divided using this meshcount
meshcount = 11

# circumscribed sphere radius
radius = length * math.sqrt(3) /2

mk = MKCore()
geom = mk.igeom_instance()

brick = geom.createBrick([length, length, length])
# geom.rotateEnt(brick, 10, [1., 1., 0])

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
#compute also the spherical coordinates, and some u, v velocities
utag=mesh.createTag('Uvel', 1, 'd')
vtag=mesh.createTag('Vvel', 1, 'd')
vectag=mesh.createTag('Velo', 3, 'd')
T=5
t=0.1

Cos=math.cos
Sin=math.sin
for node in nodes:
  x, y, z = mesh.getVtxCoords(node)
  dist1 = math.sqrt( x*x + y*y + z*z )
  ratio = radius/dist1
  x1 = x*ratio
  y1 = y*ratio
  z1 = z*ratio
  mesh.setVtxCoords(node, [x1, y1, z1])
  [r, elev, az] = cart2sph(x,y,z)
  elev1=elev-2*math.pi*t/T
  uu=3*radius/T*math.pow(Sin(elev1),2)*Sin(2*az)*Cos(math.pi*t/T)+2*math.pi*radius*Cos(az)/T
  vv=3*radius/T*(Sin(2*elev1))*Cos(az)*Cos(math.pi*t/T)
  
  utag[node]=uu
  vtag[node]=vv
  vx=-uu*Sin(az)-vv*Sin(elev)*Cos(az)
  vy= uu*Cos(az)-vv*Sin(elev)*Sin(az)
  vz= vv*Cos(elev)
  vectag[node] = [vx, vy, vz]

mesh_file2="../data/Homme3.h5m"
mk.save_mesh(mesh_file2)
mesh_file3="../data/eulerHomme.vtk"
mk.save_mesh(mesh_file3)

time=0.05
# rotate an extra 360/(4*meshcount) around z, so the deformation will be bigger 
rA = math.pi/(4*meshcount)*2.5
# something random
# mesh size is about 2*radius*PI/(4*meshcount)
meshSize = 2*radius*math.pi/(4*meshcount)
for node in nodes:
  x, y, z = mesh.getVtxCoords(node)
  [vx, vy, vz] = vectag[node]
  x = x + vx*time
  y = y + vy*time
  z = z + vz*time
  x2 = x * Cos(rA) - y*Sin(rA)
  y2 = x * Sin(rA) + y*Cos(rA) 
  x = x2
  y = y2
  x=x+0.2*random.random()
  y=y+0.2*random.random()
  z=z+0.2*random.random()
  dist1 = math.sqrt( x*x + y*y + z*z )
  ratio = radius/dist1
  x1 = x*ratio
  y1 = y*ratio
  z1 = z*ratio
  mesh.setVtxCoords(node, [x1, y1, z1])


mesh_file4="../data/lagrangeHomme.vtk"
mk.save_mesh(mesh_file4)

  

