from MeshKit import *

file_name = "../data/simpletet.sat"
mk = MKCore()
mk.load_geometry(file_name)

vols = MEntVector()
mk.get_entities_by_dimension(3, vols)

first_vol = MEntVector()
first_vol.push_back(vols.front())

esize = SizingFunction(mk, -1, 0.25)
first_vol[0].sizing_function_index(esize.core_index())
mk.construct_meshop("CAMALTetMesher", first_vol)
