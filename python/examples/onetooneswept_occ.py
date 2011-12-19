from MeshKit import *
from optparse import OptionParser

parser = OptionParser(usage="Usage: %prog [options] [FILE_GEOM] [FILE_MESH]")
parser.add_option("-s", "--save", dest="save",
                  help="(optional) filename to save")

options, args = parser.parse_args()

if len(args) == 0:
    file_name = "../data/brick.stp"
    file_name_mesh = "../data/sf.h5m"
elif len(args) == 2:
    file_name = args[0]
    file_name_mesh = args[1]
else:
    print parser.get_usage()


mk = MKCore()
mk.load_geometry_mesh(file_name, file_name_mesh)

vols = MEntVector()
mk.get_entities_by_dimension(3, vols)

surfs = MEntVector()
vols[0].get_adjacencies(2, surfs)

curves = MEntVector()
vols[0].get_adjacencies(1, curves)

vertices = MEntVector()
vols[0].get_adjacencies(0, vertices)

first_vol = MEntVector()
first_vol.push_back(vols.front())
sw = mk.construct_meshop("OneToOneSwept", vols)

sw.SetSourceSurface(0)
sw.SetTargetSurface(1)

swSize = SizingFunction(mk, 6, -1)
vols[0].sizing_function_index(swSize.core_index())

mk.setup_and_execute()

if options.save:
    mk.save_mesh(options.save)
