int VTKMeshExporter::saveAs(const SharedFaceMesh &fmesh, const string &fname)
{
    if (fmesh == NULL) return FAIL;

    ofstream ofile(fname.c_str(), ios::out);

    nodeset = fmesh->getNodes(1);
    nodeset->getRenumbered();

    int numNodes = nodeset->getSize();

    ofile << "# vtk DataFile Version 2.0" << endl;
    ofile << " ANL BioMesh Generator " << endl;

    switch (fileFormat) {
    case ASCII_FILE:
        ofile << "ASCII " << endl;
        break;
    case BINARY_FILE:
        ofile << "BINARY " << endl;
        break;
    default:
        cout << " Error: No file format selected " << endl;
        break;
    }
    ofile << endl;

    ofile << "DATASET POLYDATA " << endl;
    ofile << "POINTS " << numNodes << " float " << endl;

    ofile << setiosflags(ios::fixed);

    size_t index1 = 0;

    SharedVertex vertex;
    if (fileFormat == ASCII_FILE) {

        FOR_EACH_VERTEX( vertex, nodeset)
        {
            assert(vertex->getLocalID() == index1);
            index1++;
            Point3D &xyz = vertex->getCoords();
            ofile << xyz[0] << " " << xyz[1] << " " << xyz[2] << endl;
        }
    }

    if (fileFormat == BINARY_FILE) {
        float fbuf[3];

        FOR_EACH_VERTEX( vertex, nodeset)
        {
            assert(vertex->getLocalID() == index1);
            index1++;
            Point3D &xyz = vertex->getCoords();
            fbuf[0] = xyz[0];
            fbuf[1] = xyz[1];
            fbuf[2] = xyz[2];
            if (enableByteSwap) {
                swap_bytes(&fbuf[0]);
                swap_bytes(&fbuf[1]);
                swap_bytes(&fbuf[2]);
            }
            ofile.write(reinterpret_cast<char*> (fbuf), 3 * sizeof(float));
        }
    }

    ofile << endl;

    int numCells = fmesh->getSize();

    int numConn = 0;

    SharedFace face;
    FOR_EACH_FACE(face, fmesh)
    {
        numConn += face->getNumNodes();
    }

    ofile << "POLYGONS " << numCells << "  " << numCells + numConn << endl;

    vector<int> celltype;
    celltype.resize(fmesh->getSize());

    if (fileFormat == ASCII_FILE) {
        index1 = 0;

        FOR_EACH_FACE( face, fmesh)
        {
            int nsize = face->getNumNodes();
            ofile << nsize << " ";
            if (face->getShapeID() == 23) {
                ofile << face->getVertex(0)->getLocalID() << " "
                        << face->getVertex(1)->getLocalID() << " "
                        << face->getVertex(2)->getLocalID() << endl;
                celltype[index1++] = 5;
            } else if (face->getShapeID() == 24) {
                ofile << face->getVertex(0)->getLocalID() << " "
                        << face->getVertex(1)->getLocalID() << " "
                        << face->getVertex(3)->getLocalID() << " "
                        << face->getVertex(2)->getLocalID() << endl;
                celltype[index1++] = 9;
            } else {
                for (int i = 0; i < nsize; i++)
                    ofile << face->getVertex(i)->getLocalID() << " ";
                celltype[index1++] = 7;
                ofile << endl;
            }
        }
    }

    if (fileFormat == BINARY_FILE) {
        vector<int> ibuf;
        index1 = 0;

        FOR_EACH_FACE( face, fmesh)
        {
            int nsize = face->getNumNodes();
            ibuf.resize(nsize + 1);
            ibuf[0] = nsize;
            if (enableByteSwap) swap_bytes(&ibuf[0]);
            if (face->getShapeID() == 23) {
                celltype[index1++] = 5;
                ibuf[1] = face->getVertex(0)->getLocalID();
                ibuf[2] = face->getVertex(1)->getLocalID();
                ibuf[3] = face->getVertex(2)->getLocalID();
            } else if (face->getShapeID() == 24) {
                celltype[index1++] = 9;
                ibuf[1] = face->getVertex(0)->getLocalID();
                ibuf[2] = face->getVertex(1)->getLocalID();
                ibuf[3] = face->getVertex(3)->getLocalID();
                ibuf[4] = face->getVertex(2)->getLocalID();
            } else {
                celltype[index1++] = 7;
                for (int i = 0; i < nsize; i++)
                    ibuf[i + 1] = face->getVertex(i)->getLocalID();
            }
            if (enableByteSwap)
                for (int i = 1; i < nsize + 1; i++) swap_bytes(&ibuf[i]);
            ofile.write(reinterpret_cast<char*> (&ibuf[0]), (nsize + 1) * sizeof(int));
        }
    }

    ofile << endl;

    /*
    ofile << "CELL_TYPES " << numCells << endl;
    if( fileFormat == ASCII_FILE ) {
      for( int i = 0; i < celltype.size(); i++) 
        ofile << celltype[i] << endl;
    }

    if( fileFormat == BINARY_FILE ) {
      int nsize = celltype.size();
      if( enableByteSwap ) 
        for( int i = 0; i < nsize; i++) 
          swap_bytes( &celltype[i] );
      ofile.write( reinterpret_cast<char*>( &celltype[0] ), nsize*sizeof(int));
    }
    ofile << endl;
     */


    if (attribs[0].size()) {
        ofile << "CELL_DATA " << fmesh->getSize() << endl;
        ofile << "POINT_DATA " << nodeset->getSize() << endl;
        for (int i = 0; i < attribs[0].size(); i++) {
            string aname = attribs[0][i];
            if (!aname.empty()) {
                ofile << "SCALARS " << aname << " float " << endl;
                ofile << "LOOKUP_TABLE  default " << endl;
                write_node_attrib(nodeset, attribs[0][i], ofile);
            }
        }
    }

    for (int i = 0; i < attribs[2].size(); i++)
        write_face_attrib(fmesh, attribs[2][i], ofile);

    return PASS;
}

