import unittest
from MeshKit import *
from itaps import iGeom, iMesh, iRel

class TestPytaps(unittest.TestCase):
    def setUp(self):
        self.mk = MKCore()

    def testInstance(self):
        self.assertTrue(isinstance(self.mk.igeom_instance(), iGeom.Geom))
        self.assertTrue(isinstance(self.mk.imesh_instance(), iMesh.Mesh))
        self.assertTrue(isinstance(self.mk.irel_instance(),  iRel.Rel))

if __name__ == '__main__':
    unittest.main()
