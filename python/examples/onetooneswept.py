from MeshKit import *

file_name = "../data/BrickWithSrcMeshed.cub"
mk = MKCore()
mk.load_geometry_mesh(file_name, file_name)

vols = MEntVector()
mk.get_entities_by_dimension(3, vols)

this_vol = vols[0]

surfs = MEntVector()
this_vol.get_adjacencies(2, surfs)

curves = MEntVector()
this_vol.get_adjacencies(1, curves)

vertices = MEntVector()
this_vol.get_adjacencies(0, vertices)

sw = mk.construct_meshop("OneToOneSwept", vols)

sw.SetSourceSurface(1)
sw.SetTargetSurface(0)

swSize = SizingFunction(mk, 6, -1)
this_vol.sizing_function_index(swSize.core_index())

mk.setup_and_execute()

mk.save_mesh("OneToOneSwept.vtk")
