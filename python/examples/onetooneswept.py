from MeshKit import *
from optparse import OptionParser

parser = OptionParser(usage="Usage: %prog [options] [FILE]")
parser.add_option("-s", "--save", dest="save",
                  help="(optional) filename to save")

options, args = parser.parse_args()

if len(args) == 0:
    file_name = "../data/BrickWithSrcMeshed.cub"
elif len(args) == 1:
    file_name = args[0]
else:
    print parser.get_usage()


mk = MKCore()
mk.load_geometry_mesh(file_name, file_name)

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

sw.SetSourceSurface(1)
sw.SetTargetSurface(0)

swSize = SizingFunction(mk, 6, -1)
vols[0].sizing_function_index(swSize.core_index())

mk.setup_and_execute()

if options.save:
    mk.save_mesh(options.save)
