import unittest
from MeshKit import *
from itaps import iBase

file_name = "/home/jvporter/dev/meshkit/data/BrickWithSrcMeshed.cub"

class TestPytaps(unittest.TestCase):
    def setUp(self):
        self.mk = MKCore()

    def testBrick(self):
        self.mk.load_geometry_mesh(file_name, file_name)

        vols = MEntVector()
        self.mk.get_entities_by_dimension(3, vols)

        surfs = MEntVector()
        vols[0].get_adjacencies(2, surfs)
        self.assertEqual(surfs.size(), 6)

        curves = MEntVector()
        vols[0].get_adjacencies(1, curves)
        self.assertEqual(curves.size(), 12)

        vertices = MEntVector()
        vols[0].get_adjacencies(0, vertices)
        self.assertEqual(vertices.size(), 8)

        sw = self.mk.construct_meshop("OneToOneSwept", vols)

        sw.SetSourceSurface(1)
        sw.SetTargetSurface(0)

        swSize = SizingFunction(self.mk, 6, -1)
        vols[0].sizing_function_index(swSize.core_index())

        self.mk.setup_and_execute()
        self.assertEqual(self.mk.imesh_instance().getNumOfType(
                iBase.Type.region), 600)

if __name__ == '__main__':
    unittest.main()
