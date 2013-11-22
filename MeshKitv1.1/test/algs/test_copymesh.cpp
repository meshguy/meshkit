/** \file test_copymesh.cpp \test
 *
 * Test CopyMesh & MergeMesh
 *
 */

#include "meshkit/MKCore.hpp"
#include "meshkit/MeshOp.hpp"
#include "meshkit/CopyMesh.hpp"
#include "meshkit/ModelEnt.hpp"
#include "meshkit/MergeMesh.hpp"

using namespace MeshKit;

#include "TestUtil.hpp"
#define DEFAULT_TEST_FILE1 "c.exo"
#define DEFAULT_TEST_FILE2 "d.exo"

MKCore *mk;

void test_copy_quad();
void test_load_and_copymove();
void test_set_imesh();

int main(int argc, char **argv)
{
    mk = new MKCore();
    int num_fail = 0;
    num_fail += RUN_TEST(test_set_imesh);
    num_fail += RUN_TEST(test_load_and_copymove);
    num_fail += RUN_TEST(test_copy_quad);

    delete mk;
    return num_fail;
}

void test_set_imesh()
{
    std::vector <CopyMesh*> cm;
    iMesh *imesh = mk->imesh_instance();
    int err =0;
    std::vector<std::string> files;
    std::string name1 = TestDir + "/" + DEFAULT_TEST_FILE1;
    std::string name2 = TestDir + "/" + DEFAULT_TEST_FILE2;

    files.push_back(name1);
    files.push_back(name2);
    cm.resize(files.size());
    iBase_EntitySetHandle orig_set;

    for (unsigned int i = 0; i < files.size(); i++) {
        iMesh_createEntSet(imesh->instance(), 0, &orig_set, &err);
        //ERRORR("Couldn't create file set.", err);
        std::cout << "Loading File: " << files[i].c_str() << std::endl;
        iMesh_load(imesh->instance(), orig_set, files[i].c_str(), NULL, &err, strlen(files[i].c_str()), 0);

        ModelEnt *me;
        me = NULL;
        me = new ModelEnt(mk, iBase_EntitySetHandle(0), /*igeom instance*/0, (moab::EntityHandle)orig_set, 0);
        MEntVector assm_set;
        //assm_set.clear();
        assm_set.push_back(me);
        cm[i] = (CopyMesh*) mk->construct_meshop("CopyMesh", assm_set);
        cm[i]->set_name("copy_move_mesh");
        cm[i]->copy_sets().add_set(orig_set);
        // some entity tag types are always copy or expand
        cm[i]->expand_sets().add_tag("MATERIAL_SET");
        cm[i]->expand_sets().add_tag("DIRICHLET_SET");
        cm[i]->expand_sets().add_tag("NEUMANN_SET");
        cm[i]->copy_sets().add_tag("GEOM_DIMENSION");

        Vector<3> dx; dx[0] = 10; dx[1] = 0; dx[2] = 0;
        cm[i]->set_transform(Copy::Translate(dx));
    }
    // copy and merge
    mk->setup_and_execute();
     for (unsigned int i = 0; i < files.size(); i++)
         delete cm[i];
}


void test_load_and_copymove()
{
    //clear_mesh(mk);
    mk->delete_all();
    std::string filename = TestDir + "/" + DEFAULT_TEST_FILE1;
    mk->load_mesh(filename.c_str());

    // get the hexes
    MEntVector vols, v1;
    //mk->get_entities_by_dimension(3, vols);
    moab::EntityHandle root_set = mk->moab_instance()->get_root_set();
    ModelEnt me(mk, iBase_EntitySetHandle(0), /*igeom instance*/0, root_set);
    vols.push_back(&me);

    CopyMesh *cm = (CopyMesh*) mk->construct_meshop("CopyMesh", vols);
    cm->set_name("copy_move_mesh");
    // some entity tag types are always copy or expand
    cm->expand_sets().add_tag("MATERIAL_SET");
    cm->expand_sets().add_tag("DIRICHLET_SET");
    cm->expand_sets().add_tag("NEUMANN_SET");

    Vector<3> dx; dx[0] = 10; dx[1] = 0; dx[2] = 0;
    cm->set_transform(Copy::Translate(dx));

    //accesses entities for merging directly from moab instance, vols are used for finding the dimension
    MergeMesh *mm = (MergeMesh*) mk->construct_meshop("MergeMesh", v1);
    mm->set_name("merge_mesh");

    double merge_tol = 1e-3; int updatesets = 0, domerge = 1; iBase_TagHandle merge_tag = NULL;
    mm->set_merge_params(merge_tol, updatesets, domerge, merge_tag);

    // put them in the graph
    mk->get_graph().addArc(mk->root_node()->get_node(), cm->get_node());
    mk->get_graph().addArc(cm->get_node(), mm->get_node());
    mk->get_graph().addArc(mm->get_node(), mk->leaf_node()->get_node());

    // copy and merge
    mk->setup_and_execute();
    delete cm;
}


void test_copy_quad()
{
    mk->delete_all();
    iMesh *mesh = mk->imesh_instance();

    double coords[] = {
        0, 0, 0,
        1, 0, 0,
        1, 1, 0,
        0, 1, 0,
    };

    iMesh::EntityHandle verts[4];
    iMesh::EntityHandle quad;
    iMesh::EntitySetHandle set;

    mesh->createVtxArr(4, iBase_INTERLEAVED, coords, verts);
    mesh->createEnt(iMesh_QUADRILATERAL, verts, 4, quad);
    mesh->createEntSet(true, set);
    mesh->addEntToSet(quad, set);

    ModelEnt me(mk, iBase_EntitySetHandle(0), /*igeom instance*/0, (moab::EntityHandle)set);
    MEntVector selection;
    selection.push_back(&me);

    CopyMesh *cm = (CopyMesh*) mk->construct_meshop("CopyMesh", selection);
    cm->set_name("copy_move_mesh");

    Vector<3> dx; dx[0] = 0; dx[1] = 0; dx[2] = 10;
    cm->set_transform(Copy::Translate(dx));

    // put them in the graph
    mk->get_graph().addArc(mk->root_node()->get_node(), cm->get_node());
    mk->get_graph().addArc(cm->get_node(), mk->leaf_node()->get_node());

    cm->copy_sets().add_set(set);

    // mesh embedded boundary mesh, by calling execute
    mk->setup_and_execute();

    std::vector<iMesh::EntityHandle> new_verts;
    double new_coords[4*3];
    std::vector<iMesh::EntityHandle> new_quad;
    iMesh::EntitySetHandle new_set;

    mesh->getEntSetEHData(set, cm->copy_tag(), (iMesh::EntityHandle&)new_set);

    mesh->getEntities(new_set, iBase_ALL_TYPES, iMesh_ALL_TOPOLOGIES, new_quad);
    CHECK_EQUAL(new_quad.size(), size_t(1));

    mesh->getEntAdj(new_quad[0], iBase_ALL_TYPES, new_verts);
    CHECK_EQUAL(new_verts.size(), size_t(4));

    mesh->getVtxArrCoords(&new_verts[0], new_verts.size(), iBase_INTERLEAVED,
            new_coords);

    for(int i=0; i<4; i++) {
        CHECK_REAL_EQUAL(coords[i*3+0],    new_coords[i*3+0], 0.00001);
        CHECK_REAL_EQUAL(coords[i*3+1],    new_coords[i*3+1], 0.00001);
        CHECK_REAL_EQUAL(coords[i*3+2]+10, new_coords[i*3+2], 0.00001);
    }
    delete cm;
}
