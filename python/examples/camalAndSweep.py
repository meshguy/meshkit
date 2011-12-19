from MeshKit import *

gfile_name = "../data/brick.stp"

mk = MKCore()
# this will create also relations
mk.load_geometry(gfile_name)

vols = MEntVector()
mk.get_entities_by_dimension(3, vols)

surfs = MEntVector()
vols[0].get_adjacencies(2, surfs)

curves = MEntVector()
vols[0].get_adjacencies(1, curves)

vertices = MEntVector()
vols[0].get_adjacencies(0, vertices)

s1= MEntVector()
s1.push_back(surfs[0])

cm = mk.construct_meshop("CAMALPaver", s1)

esize = SizingFunction(mk, -1, 0.1)

# we know only one surface
s1[0].sizing_function_index(esize.core_index())

first_vol = MEntVector()
first_vol.push_back(vols.front())
sw = mk.construct_meshop("OneToOneSwept", vols)

sw.SetSourceSurface(0)
sw.SetTargetSurface(1)

swSize = SizingFunction(mk, 6, -1)
vols[0].sizing_function_index(swSize.core_index())

mk.add_arc(cm, sw)

mk.setup_and_execute()

mesh_file = "meshSwept.h5m"
mk.save_mesh(mesh_file)
