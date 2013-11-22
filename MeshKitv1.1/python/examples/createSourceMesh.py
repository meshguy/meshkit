from MeshKit import *

file_name = "../data/brick.stp"


mk = MKCore()
mk.load_geometry(file_name)

vols = MEntVector()
mk.get_entities_by_dimension(3, vols)

surfs = MEntVector()
vols[0].get_adjacencies(2, surfs)

firstFace = MEntVector()
firstFace.push_back(surfs[0])
esize = SizingFunction(mk, -1, 0.1)

curves = MEntVector()
vols[0].get_adjacencies(1, curves)

mk.construct_meshop("CAMALPaver", firstFace)
firstFace[0].sizing_function_index(esize.core_index())

mk.setup_and_execute()

mesh_file="../data/sf.h5m"
mk.save_mesh(mesh_file)
