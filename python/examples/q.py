from MeshKit import *

qopts = QslimOptions()

file_name = "../data/partBed.smf"
outfile = "out.h5m"
qopts.will_constrain_boundaries = 1
qopts.face_target = 4500
qopts.boundary_constraint_weight = 100
qopts.height_fields = 1
print file_name 

mk = MKCore()

mk.load_mesh(file_name)

surfs = MEntVector()
mk.get_entities_by_dimension(2, surfs)

qm = mk.construct_meshop("QslimMesher", surfs)
qm.set_options(qopts)

mk.setup_and_execute()

mk.save_mesh(outfile)




