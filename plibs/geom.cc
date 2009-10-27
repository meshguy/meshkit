#include "itaps.h"

///////////////////////////////////////////////////////////////////////////////

double get_geom_edge_length(iGeom_Instance geom, iBase_EntityHandle gEdge)
{
    int err;
    SimpleArray<double> edgeLength;
    iGeom_measure( geom, &gEdge, 1, ARRAY_INOUT(edgeLength), &err);
    return edgeLength[0];
}
  
///////////////////////////////////////////////////////////////////////////////

int readGeometry(string &filename, iGeom_Instance &geom)
{
    string engine_opt = ";engine=OCC";

    int err;
    iGeom_newGeom(engine_opt.c_str(), &geom, &err, engine_opt.length());
    iGeom_load(geom, &filename[0], 0, &err, filename.length(), 0);

    iBase_EntitySetHandle rootSet;
    iGeom_getRootSet(geom, &rootSet, &err);

    cout << "Model Contents " << endl;
    const char *gtype[] = {"vertices: ", "edges: ", "faces: ", "regions: "};

    int count;
    for (int i = 0; i <= 3; ++i)
    {
        iGeom_getNumOfType(geom, rootSet, i, &count, &err);
        std::cout << gtype[i] << count << std::endl;
    }
}

///////////////////////////////////////////////////////////////////////////////

void Discretize_Geometric_Model(GeomMesh &geomesh)
{
    char *options = NULL;
    int optlen = 0;
    int err;

    iMesh_Instance mesh = geomesh.mesh;
    iGeom_Instance geom = geomesh.geom;
    iRel_Instance assoc = geomesh.assoc;
    iRel_RelationHandle relation00 = geomesh.relation00;

    iBase_EntitySetHandle geomRootSet = geomesh.geomRootSet;

    ///////////////////////////////////////////////////////////////////////////////
    // Step I: Collect all the feature nodes in Geometry:
    ///////////////////////////////////////////////////////////////////////////////
    SimpleArray<iBase_EntityHandle> geoNodes;
    iGeom_getEntities(geom, geomRootSet, iBase_VERTEX, ARRAY_INOUT(geoNodes), &err);

    double x, y, z;
    iBase_EntityHandle newHandle;
    for (int i = 0; i < geoNodes.size(); i++)
    {
        iGeom_getVtxCoord(geom, geoNodes[i], &x, &y, &z, &err);
        iMesh_createVtx(mesh, x, y, z, &newHandle, &err);
        iRel_setEntEntAssociation(assoc, relation00, geoNodes[i], newHandle, &err);
    }

    ///////////////////////////////////////////////////////////////////////////////
    // Step II: Collect all the geometric edges and discretize them.
    ///////////////////////////////////////////////////////////////////////////////
    SimpleArray<iBase_EntityHandle> geoEdges;
    iGeom_getEntities(geom, geomRootSet, iBase_EDGE, ARRAY_INOUT(geoEdges), &err);

    int numGeoEdges = geoEdges.size();

    vector<double> elen(numGeoEdges);
    for (int i = 0; i < numGeoEdges; i++) {
        elen[i] = get_geom_edge_length(geom, geoEdges[i]);
    }
    vector<double> tmpelen;
    tmpelen = elen;
    sort( tmpelen.begin(), tmpelen.end() );
    double mean_edge_length = tmpelen[numGeoEdges/2];

    cout << " Minimum Edge Length : " << *min_element( elen.begin(), elen.end() ) << endl;
    cout << " Maximum Edge Length : " << *max_element( elen.begin(), elen.end() ) << endl;
    cout << " Mean Length         : " <<  mean_edge_length << endl;

    double spacing = mean_edge_length /(double) 10.0;

    int numEdges;
    EdgeMesher *edgeMesher = geomesh.edgeMesher;
    for (int i = 0; i < numGeoEdges; i++) {
        if( elen[i] < 0.5*mean_edge_length) 
            numEdges = 1;
        else
            numEdges = elen[i]/ spacing;

        edgeMesher->discretize(geomesh, geoEdges[i], numEdges);
    }
    return;
   
    ///////////////////////////////////////////////////////////////////////////////
    //Step III: Collect all the geometric faces and discretize them.
    ///////////////////////////////////////////////////////////////////////////////
    SimpleArray<iBase_EntityHandle> geoFaces;
    iGeom_getEntities(geom, geomRootSet, iBase_FACE, ARRAY_INOUT(geoFaces), &err);

    int numGeomFaces = geoFaces.size();

    SurfaceMesher *surfMesher = geomesh.surfMesher;
    for (int i = 2; i < numGeomFaces; i++) {
        surfMesher->discretize(geomesh, geoFaces[i]);
        exit(0);
    }

    return ;

    ///////////////////////////////////////////////////////////////////////////////
    //Step IV: Collect all the geometric regions and discretize them.
    ///////////////////////////////////////////////////////////////////////////////
    SimpleArray<iBase_EntityHandle> geoCells;
    iGeom_getEntities(geom, geomRootSet, iBase_REGION, ARRAY_INOUT(geoCells), &err);

    VolumeMesher *volMesher = geomesh.volMesher;
    for (int i = 0; i < geoCells.size(); i++)
        volMesher->discretize(geomesh, geoCells[i]);

    ///////////////////////////////////////////////////////////////////////////////
    //Step V: Generate Spectral (Higher order) elements:
    ///////////////////////////////////////////////////////////////////////////////
    generate_spectral_elements(geomesh, 10); // 10 noded tetrahedra.
}

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
    example_volume_mesher( "cube.off");
    exit(0);
    int err, result;

    if (argc != 2)
    {
        cout << "Usage: executable geomfile (*.brep) " << endl;
        return 1;
    }

    ////////////////////////////////////////////////////////////////////////////
    //Read the BRep geometric file
    ////////////////////////////////////////////////////////////////////////////

    iGeom_Instance geom;
    string geomfile = argv[1];
    readGeometry(geomfile, geom);

    ////////////////////////////////////////////////////////////////////////////
    // Generate the mesh.
    ////////////////////////////////////////////////////////////////////////////
    EdgeMesher    *edgeMesher = new EdgeMesher;
    SurfaceMesher *surfMesher = SurfaceMesherFactory::getProduct("NetGen");
    VolumeMesher   *volMesher = VolumeMesherFactory::getProduct("NetGen");

    GeomMesh geomesh;
    geomesh.setGeometry(geom);

    geomesh.edgeMesher = edgeMesher;
    geomesh.surfMesher = surfMesher;
    geomesh.volMesher  = volMesher;

    Discretize_Geometric_Model(geomesh);

    ////////////////////////////////////////////////////////////////////////////
    //Save the mesh
    ////////////////////////////////////////////////////////////////////////////

    const char *outfile = "model.vtk";
    int namelen = strlen(outfile);
    const char *options = NULL;
    int optlen = 0;

    iMesh_Instance mesh = geomesh.mesh;
    iBase_EntitySetHandle meshRootSet = geomesh.meshRootSet;

    iMesh_save(mesh, meshRootSet, outfile, options, &err, namelen, optlen);

    return 0;
}

////////////////////////////////////////////////////////////////////////////////
