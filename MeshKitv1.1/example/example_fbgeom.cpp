/*!
\example example_fbgeom.cpp

\section example_fbgeom_cpp_title FBiGeom Example Using RayTracing

\subsection example_fbgeom_cpp_in Input
It reads top (cropSE.h5m) and bottom (quadsOnBed.vtk) mesh file from data directory.

\subsection example_fbgeom_cpp_out Output
\image html fbgeom_out.jpg "Resulting hex mesh from input top and bottom surfaces."

\subsection example_fbgeom_cpp_src Source Code
*/
#include "meshkit/MKCore.hpp"
#include "meshkit/FBiGeom.hpp"
#include "meshkit/SizingFunction.hpp"
#include "meshkit/ModelEnt.hpp"

using namespace MeshKit;

#include "stdlib.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <time.h>
#include "moab/CartVect.hpp"
#include "moab/Interface.hpp"
#include "moab/Skinner.hpp"
//#include "MBiMesh.hpp"
#include <vector>

MKCore * mk;
bool save_mesh = true;
//std::string file_name;
std::string bottom_mesh; //="quadsOnBed.vtk";
std::string top_filename; // "cropSE.h5m", to be smaller
std::string out_mesh; // deleted by default
int nbLayers;
std::vector<double> grades;

bool create_a_set_with_tag_value(iMesh_Instance meshIface,
                                 iBase_EntityHandle * ents, int num_ents, char * tag_name,
                                 int size_name_tag, int value, iBase_EntitySetHandle & mesh_set);

int main(int argc, char *argv[])
{
    // Check command line arg

    nbLayers = 5; // a good number of layers
    if (argc < 4) {
        std::cout << "Usage: " << argv[0]
                  << " <bottom_mesh> <top_filename> <out_mesh> [-n layers]  <ratios> "

                  << std::endl;
        bottom_mesh = (std::string) MESH_DIR + "/" + "quadsOnBed.vtk";
        top_filename = std::string (MESH_DIR) + "/" + "cropSE.h5m"; //, to be smaller
        out_mesh = "hexas.h5m"; // deleted by default
        std::cout << " Using : " << bottom_mesh << " " << top_filename << " "
                  << out_mesh << " -n 5" << "\n";

        std::cout << " the default output file is deleted \n";
    } else {
        bottom_mesh = argv[1];
        top_filename = argv[2];
        out_mesh = argv[3];
        save_mesh = true;
        int argno = 4;
        while (argno < argc) {
            if (!strcmp(argv[argno], "-n")) {
                argno++;
                nbLayers = atoi(argv[argno]);
                argno++;
            } else {
                //everything that is not interpreted must be a ratio
                grades.push_back(atof(argv[argno]));
                argno++;
            }
        }
    }
    int nbGrades = (int) grades.size();
    if (nbGrades > 0) {
        // add all grades, then divide by total to get the new accumulated ratios
        double total = 0.;
        int i;
        for (i = 0; i < nbGrades; i++) {
            total += grades[i];
            //grades[i] = total; // accumulated so far
        }
        for (i = 0; i < nbGrades; i++) {
            grades[i] /= total;
        }
        nbLayers = nbGrades;
    } else {
        nbGrades = nbLayers;
        for (int i = 0; i < nbLayers; i++)
            grades.push_back(1. / nbLayers);
    }

    clock_t start_time = clock();
    mk = new MKCore;
    iBase_ErrorType err;
    FBiGeom * fbiGeom = new FBiGeom(); // true for smooth, false for linear
    mk->add_igeom_instance(fbiGeom);
    // read initial mesh (quad mesh surface bottom)
    char opts[]="SMOOTH;";
    err = fbiGeom->load(top_filename.c_str(), opts);// loading the top surface , as geometry
    IBERRCHK(err, "Failure to load smooth geometry.");
    clock_t load_time1 = clock();// load fb file (expensive),
    // load here the quads from bottom
    iBase_EntitySetHandle setForQuads;

    err = mk->imesh_instance()->createEntSet(false, setForQuads);
    IBERRCHK(err, "Failure to create set for loading the quads");
    // load the quads in this set
    err = mk->imesh_instance()->load(setForQuads, bottom_mesh.c_str());
    IBERRCHK(err, "Failure to load quad file.");

    clock_t load_time2 = clock();
    // and quads file (cheap)
    //
    double direction[3] = { 0., 0., 1. }; // normalized

    // can we get to iMesh from iMesh.hpp?
    iMesh_Instance mesh1 = mk->imesh_instance()->instance();

    // first get all nodes from bottom mesh
    // get the quads and the vertices in one shot
    iBase_EntityHandle *quads = NULL;
    int quads_alloc = 0;
    iBase_EntityHandle *vert_adj = NULL;
    int vert_adj_alloc = 0, vert_adj_size, numQuads;
    int * offsets = NULL, offsets_alloc = 0, indices_size;
    int * indices = NULL, indices_alloc = 0, offsets_size;
    int ierr=0;
    iMesh_getAdjEntIndices(mesh1, setForQuads, iBase_FACE, iMesh_QUADRILATERAL,
                           iBase_VERTEX, &quads, &quads_alloc, &numQuads, &vert_adj,
                           &vert_adj_alloc, &vert_adj_size, &indices, &indices_alloc, &indices_size,
                           &offsets, &offsets_alloc, &offsets_size, &ierr);
    // create one set with desired tag (1)
    iBase_EntitySetHandle mesh_set1;
    char nsname[]= "NEUMANN_SET";
    create_a_set_with_tag_value(mesh1, quads, numQuads, nsname , 11, 1, mesh_set1);
    /* get the coordinates in one array */

    // delete now the original set?
    int vert_coords_alloc = 0;
    double * xyz = 0; // not allocated
    int vertex_coord_size = 0;

    iMesh_getVtxArrCoords(mesh1, vert_adj, vert_adj_size, iBase_INTERLEAVED,
                          &xyz, &vert_coords_alloc, &vertex_coord_size, &ierr);
    //IBERRCHK(err, "failed to get vertex coordinates of entities.");

    // then, go to shoot rays to decide the sweeping thickness
    int & numNodes = vert_adj_size; // to not change too much
    double * dArr = new double[numNodes];
    // then, go to shoot rays
    int j = 0;
    int numRaysIntersected = 0;
    // double factorFloating = (1.-937./1026.);
    for (j = 0; j < numNodes; j++) {
        //dArr[j] = xyz[j * 3 + 2];
        dArr[j] = 0; // no intersection situation, marked by 0
        // for a point, see if it is inside the polygon, with winding number

        iBase_StorageOrder order=iBase_INTERLEAVED;
        std::vector<iGeom::EntityHandle> entities_out;
        std::vector<double> points_out;
        std::vector<double> params_out;

        double pos[3] = { xyz[3 * j], xyz[3 * j + 1], xyz[3 * j + 2] };
        err= fbiGeom->getPntRayIntsct(pos[0], pos[1], pos[2], direction[0],
                direction[1], direction[2], order, entities_out, points_out,
                params_out);
        // get the first coordinate
        if (err != 0 || entities_out.size() < 1)
            continue;// we should bail out
        //double zTop = points_out[2]; // the z of the top, intersection computation

        // consider only the first intersection point
        dArr[j] = params_out[0]; // the first intersection only
        numRaysIntersected++;
    }
    // to xyz coordinates, add some thickness, and (thick/layers), to create
    // a next layer, and the hexas
    iBase_EntityHandle * layer1 = vert_adj;
    iBase_EntitySetHandle setForHexas;

    err = mk->imesh_instance()->createEntSet(false, setForHexas);
    IBERRCHK(err, "Failure to create set for loading the hexas");
    // create in a loop layer 2 , and vertices in layer 2, and hexas
    // between layer 1 and layer2
    iBase_EntitySetHandle mesh_set2;
    for (int i = 0; i < nbLayers; i++) {
        for (int j = 0; j < numNodes; j++) {
            for (int k = 0; k < 3; k++)
                xyz[3 * j + k] += direction[k] * dArr[j] * grades[i];
        }
        // now create new vertices at this position, layer 2
        iBase_EntityHandle * newVerts = NULL; // no vertices yet
        int size1, size2;
        iMesh_createVtxArr(mesh1,
                           /*in*/numNodes,
                           /*in*/iBase_INTERLEAVED,
                           /*in*/xyz,
                           /*in*/numNodes * 3,
                           /*inout*/&newVerts,
                           /*inout*/&size1,
                           /*inout*/&size2,
                           /*out*/&ierr);

        int numHexa = numQuads;
        long int * adjacency = new long int[8 * numHexa];
        iBase_EntityHandle * conn = (iBase_EntityHandle *) adjacency;
        for (int L = 0; L < numHexa; L++) {

            for (int k = 0; k < 4; k++) {
                conn[8 * L + k] = layer1[indices[offsets[L] + k]];
                conn[8 * L + k + 4] = newVerts[indices[offsets[L] + k]];
            }
        }

        int n = numHexa;
        int junk1 = n, junk2 = n, junk3 = n, junk4 = n;
        int * stat = new int[numHexa];
        int* ptr2 = stat;
        //int ierr;
        iBase_EntityHandle * new_entity_handles = NULL;
        iMesh_createEntArr(mesh1, iMesh_HEXAHEDRON, conn, 8 * n,
                           &new_entity_handles, &junk1, &junk2, &ptr2, &junk3, &junk4, &ierr);
        // add hexas to the setForHexas set
        iMesh_addEntArrToSet(mesh1,
                             new_entity_handles, numHexa,
                             setForHexas,
                             &ierr );
        // end copy
        // at the end, layer 1 becomes layer 2
        layer1 = newVerts;
        // if last layer, create top quads too
        if (i == nbLayers - 1) {
            // last layer, create top quads too
            for (int L = 0; L < numQuads; L++) {

                for (int k = 0; k < 4; k++) {
                    //conn[8*L+k]   = layer1[ indices[offsets[L]+k ] ];
                    conn[4 * L + k] = newVerts[indices[offsets[L] + k]];
                }
            }

            iMesh_createEntArr(mesh1, iMesh_QUADRILATERAL, conn, 4 * n,
                               &new_entity_handles, &junk1, &junk2, &ptr2, &junk3, &junk4, &ierr);
            // create another set
            create_a_set_with_tag_value(mesh1, new_entity_handles, numQuads,
                                        nsname, 11, 2, mesh_set2); // top
        }
    }

    clock_t compute_time = clock();

    // now get MOAB , to use skin on the hexas, to get the rest of quads
    moab::Interface * mb = mk->moab_instance();
    // get all hexas, and get the skin
    moab::Range hexas, iniQuads;
    mb->get_entities_by_type(0, moab::MBHEX, hexas);
    mb->get_entities_by_type(0, moab::MBQUAD, iniQuads);
    moab::Skinner skinner(mb);
    moab::Range allQuads;
    skinner.find_skin(0, hexas, 2, allQuads);

    moab::Range latQuads = subtract(allQuads, iniQuads);

    // add this one too to the hexa set
    moab::EntityHandle latSet;
    mb->create_meshset(moab::MESHSET_SET, latSet);
    mb->add_entities(latSet, latQuads);

    // add range and tag
    moab::Tag ntag;
    mb->tag_get_handle("NEUMANN_SET", 1,
                       moab::MB_TYPE_INTEGER,  ntag);
    int val = 3;
    mb->tag_set_data(ntag, &latSet, 1, (void*) &val);

    // to the hexa set, add chlidren sets, mesh_set1, etc
    mb->add_parent_child( (moab::EntityHandle)setForHexas,
                          (moab::EntityHandle) mesh_set1);
    mb->add_parent_child( (moab::EntityHandle)setForHexas,
                          (moab::EntityHandle) mesh_set2);
    mb->add_parent_child( (moab::EntityHandle)setForHexas,
                          latSet);

    delete fbiGeom;
    assert(iBase_SUCCESS == err);
    iMesh_save(mesh1, setForHexas, out_mesh.c_str(), 0, &ierr, out_mesh.length(), 0);
    if (iBase_SUCCESS != err) {
        std::cerr << "ERROR saving mesh to " << out_mesh << std::endl;
        return 0;
    }
    clock_t out_time = clock();
    std::cout << "Total time is " << (double) (out_time - start_time)
                 / CLOCKS_PER_SEC << " s\n  load top : " << (double) (load_time1
                                                                      - start_time) / CLOCKS_PER_SEC << " s\n  load bottom : "
              << (double) (load_time2 - load_time1) / CLOCKS_PER_SEC
              << " s\n  compute time : " << (double) (compute_time - load_time2)
                 / CLOCKS_PER_SEC << " s\n  write time : " << (double) (out_time
                                                                        - compute_time) / CLOCKS_PER_SEC << std::endl;
    std::cout << "num rays intersected : " << numRaysIntersected << "\n";

    free(xyz);
    delete[] dArr;

    if (!save_mesh)
    {
        remove(out_mesh.c_str());
    }
    return 0;
}

// Helper function: usually a NEUMANN set
bool create_a_set_with_tag_value(iMesh_Instance meshIface,
                                 iBase_EntityHandle * ents, int num_ents, char * tag_name,
                                 int size_name_tag, int value, iBase_EntitySetHandle & mesh_set)
{
    // will first create a new entity set in the mesh; will add the
    // entities;
    // will associate the tag with the given value (create the tag if not existent yet)
    // NEUMANN_SET
    int result = iBase_SUCCESS;

    bool isList = false;
    //iBase_EntitySetHandle mesh_set;

    iMesh_createEntSet(meshIface, isList, &mesh_set, &result);
    if (iBase_SUCCESS != result)
        return false;

    iBase_TagHandle tag_handle;
    iMesh_getTagHandle(meshIface, tag_name, &tag_handle, &result, size_name_tag);
    assert (0 != tag_handle);
    iMesh_setEntSetIntData(meshIface, mesh_set, tag_handle, value, &result);
    if (iBase_SUCCESS != result)
        return false;
    // add the entities to the set
    //
    iMesh_addEntArrToSet(meshIface, ents, num_ents, mesh_set, &result);
    if (iBase_SUCCESS != result)
        return false;

    return true;
}
