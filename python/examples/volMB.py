from MeshKit import *

mk = MKCore()

filename="../data/volIce.h5m"
mk.load_mesh(filename, "", 0, 0, 0, False, False);

mk.initialize_mesh_based_geometry()

vols = MEntVector()
mk.get_entities_by_dimension(3, vols)

mk.construct_meshop("CAMALTetMesher", vols)

esize = SizingFunction(mk, -1, 300)

vols[0].sizing_function_index(esize.core_index())


mk.setup_and_execute()

mk.save_mesh_from_model_ents("out.h5m", vols)


