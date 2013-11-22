from MeshKit import *
from optparse import OptionParser

parser = OptionParser(usage="Usage: %prog [options] [FILE]")
parser.add_option("-s", "--save", dest="save",
                  help="(optional) filename to save")

options, args = parser.parse_args()

if len(args) == 0:
    file_name = "../data/simpletet.sat"
elif len(args) == 1:
    file_name = args[0]
else:
    print parser.get_usage()

mk = MKCore()
mk.load_geometry(file_name)

vols = MEntVector()
mk.get_entities_by_dimension(3, vols)

first_vol = MEntVector()
first_vol.push_back(vols.front())
mk.construct_meshop("CAMALTetMesher", vols)

esize = SizingFunction(mk, -1, 0.25)
first_vol[0].sizing_function_index(esize.core_index())

mk.setup_and_execute()

if options.save:
    mk.save_mesh(options.save)
