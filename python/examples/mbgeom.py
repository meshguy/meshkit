from MeshKit import *

mk = MKCore()

filename="../data/shell.h5m"
mk.load_mesh(filename, "", 0, 0, 0, False, False);

indx = mk.initialize_mesh_based_geometry()

surfs = MEntVector()
mk.get_entities_by_dimension(2, surfs)

mk.construct_meshop("CAMALPaver", surfs)

esize = SizingFunction(mk, -1, 0.25)
surfs[0].sizing_function_index(esize.core_index())
surfs[1].sizing_function_index(esize.core_index())

mk.setup_and_execute()

mk.remove_mesh_based_geometry(indx)

mk.save_mesh_from_model_ents("out.h5m", surfs)


