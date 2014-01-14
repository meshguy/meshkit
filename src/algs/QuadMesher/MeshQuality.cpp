#include "meshkit/Mesh.hpp"

using namespace Jaal;

void quality_values( const vector<double> &quality)
{
   cout << "   MinValue     : " << *min_element( quality.begin(), quality.end() ) << endl;
   cout << "   MaxValue     : " << *max_element( quality.begin(), quality.end() ) << endl;
   cout << "   Mean Value   : " << quality[quality.size()/2] << endl;
   size_t nSize = quality.size();
   double sum = 0.0;
   for( unsigned int i = 0; i < nSize; i++)
        sum += quality[i];
   cout << "   Avg Value    : " << sum/(double)quality.size() << endl;
}

void print_histogram( const vector<double> &quality, const string &header, int n )
{
    string filename;
    switch( n ) {
    case 3:
        filename ="./QData/T" +  header + ".gnu";
        break;
    case 4:
        filename = "./QData/Q" +  header + ".gnu";
        break;
    }

    ofstream ofile( filename.c_str(), ios::out);
    if( ofile.fail() ) {
        cout << "Warning: Can not open file " << filename << endl;
        return;
    }

    ofile << "set title \" QuadMesh Quality:  " <<  header << "\"" << "  font \"Helvetica,20\"" << endl;
    ofile << "set xlabel \"Quad Elements\"  font \"Helvetica,15 \"" << endl;
    ofile << "set ylabel \"" << header  << "\" font \"Helvetica,15\""  << endl;
    ofile << "set format y \"%11.3E\"" << endl;
    ofile << "set grid" << endl;
    ofile << "set terminal png" << endl;

    switch( n ) {
    case 3:
        filename = "./QData/T" +  header + ".data";
        ofile << "set output \"" << "T" << header << ".png\"" <<  endl;
        ofile << "plot \"" << "T" << header << ".data\"" <<  "  notitle with impulses" << endl;
        break;
    case 4:
        filename = "./QData/Q" +  header + ".data";
        ofile << "set output \"" << "Q" << header << ".png\"" <<  endl;
        ofile << "plot \"" << "Q" << header << ".data\"" <<  "  notitle with impulses" << endl;
        break;
    }
    ofile.close();

    ofile.open( filename.c_str(), ios::out);
    if( ofile.fail() ) return;

    for( unsigned int i = 0; i < quality.size(); i++)
        ofile << i << " " << quality[i] << endl;

}
///////////////////////////////////////////////////////////////////////////////

void plot_all_quad_quality_measures( Jaal::Mesh *mesh )
{
#ifdef HAVE_VERDICT
    VERDICT_REAL  coords[4][3];
    Point3D xyz;
    int num_nodes = 4;

    size_t numfaces = mesh->getSize(2);

    vector<double> quality( numfaces );

    ///////////////////////////////////////////////////////////////////////////

    for( size_t i = 0; i < numfaces; i++) {
        Face *face = mesh->getFaceAt(i);
        for( int j = 0; j < 4; j++) {
            Vertex *vertex = face->getNodeAt(j);
            xyz  = vertex->getXYZCoords();
            coords[j][0] = xyz[0];
            coords[j][1] = xyz[1];
            coords[j][2] = xyz[2];
        }
        quality[i] =  v_quad_aspect( num_nodes, coords);
    }
    sort( quality.begin(), quality.end() );
    print_histogram( quality, "AspectRatio", 4  );
    cout << "Quality: AspectRatio " << endl;
    quality_values( quality );

    ///////////////////////////////////////////////////////////////////////////

    for( size_t i = 0; i < numfaces; i++) {
        Face *face = mesh->getFaceAt(i);
        for( int j = 0; j < 4; j++) {
            Vertex *vertex = face->getNodeAt(j);
            xyz  = vertex->getXYZCoords();
            coords[j][0] = xyz[0];
            coords[j][1] = xyz[1];
            coords[j][2] = xyz[2];
        }
        quality[i] = v_quad_skew( num_nodes, coords );
    }
    sort( quality.begin(), quality.end() );
    print_histogram( quality, "Skewness", 4 );

    ///////////////////////////////////////////////////////////////////////////

    for( size_t i = 0; i < numfaces; i++) {
        Face *face = mesh->getFaceAt(i);
        for( int j = 0; j < 4; j++) {
            Vertex *vertex = face->getNodeAt(j);
            xyz  = vertex->getXYZCoords();
            coords[j][0] = xyz[0];
            coords[j][1] = xyz[1];
            coords[j][2] = xyz[2];
        }
        quality[i] =  v_quad_taper( num_nodes, coords );
    }
    sort( quality.begin(), quality.end() );
    print_histogram( quality, "Taper", 4 );

    ///////////////////////////////////////////////////////////////////////////

    for( size_t i = 0; i < numfaces; i++) {
        Face *face = mesh->getFaceAt(i);
        for( int j = 0; j < 4; j++) {
            Vertex *vertex = face->getNodeAt(j);
            xyz  = vertex->getXYZCoords();
            coords[j][0] = xyz[0];
            coords[j][1] = xyz[1];
            coords[j][2] = xyz[2];
        }
        quality[i]  = v_quad_warpage( num_nodes, coords );
    }
    sort( quality.begin(), quality.end() );
    print_histogram( quality, "Warpage", 4 );


    ///////////////////////////////////////////////////////////////////////////

    for( size_t i = 0; i < numfaces; i++) {
        Face *face = mesh->getFaceAt(i);
        for( int j = 0; j < 4; j++) {
            Vertex *vertex = face->getNodeAt(j);
            xyz  = vertex->getXYZCoords();
            coords[j][0] = xyz[0];
            coords[j][1] = xyz[1];
            coords[j][2] = xyz[2];
        }
        quality[i] = v_quad_area( num_nodes, coords );
    }
    sort( quality.begin(), quality.end() );
    print_histogram( quality, "Area", 4 );
    cout << "Quality: Area " << endl;
    quality_values( quality );


    ///////////////////////////////////////////////////////////////////////////

    for( size_t i = 0; i < numfaces; i++) {
        Face *face = mesh->getFaceAt(i);
        for( int j = 0; j < 4; j++) {
            Vertex *vertex = face->getNodeAt(j);
            xyz  = vertex->getXYZCoords();
            coords[j][0] = xyz[0];
            coords[j][1] = xyz[1];
            coords[j][2] = xyz[2];
        }
        quality[i] = v_quad_stretch( num_nodes, coords );
    }
    sort( quality.begin(), quality.end() );
    print_histogram( quality, "Stretch", 4 );

    ///////////////////////////////////////////////////////////////////////////

    for( size_t i = 0; i < numfaces; i++) {
        Face *face = mesh->getFaceAt(i);
        for( int j = 0; j < 4; j++) {
            Vertex *vertex = face->getNodeAt(j);
            xyz  = vertex->getXYZCoords();
            coords[j][0] = xyz[0];
            coords[j][1] = xyz[1];
            coords[j][2] = xyz[2];
        }
        quality[i] = v_quad_minimum_angle( num_nodes, coords );
    }
    sort( quality.begin(), quality.end() );
    print_histogram( quality, "MinAngle", 4 );
    cout << "Quality: MinAngle (in Degrees) " << endl;
    quality_values( quality );

    ///////////////////////////////////////////////////////////////////////////

    for( size_t i = 0; i < numfaces; i++) {
        Face *face = mesh->getFaceAt(i);
        for( int j = 0; j < 4; j++) {
            Vertex *vertex = face->getNodeAt(j);
            xyz  = vertex->getXYZCoords();
            coords[j][0] = xyz[0];
            coords[j][1] = xyz[1];
            coords[j][2] = xyz[2];
        }
        quality[i] = v_quad_maximum_angle( num_nodes, coords );
    }
    sort( quality.begin(), quality.end() );
    print_histogram( quality, "MaxAngle", 4 );
    cout << "Quality: MaxAngle (in Degrees) " << endl;
    quality_values( quality );

    ///////////////////////////////////////////////////////////////////////////

    for( size_t i = 0; i < numfaces; i++) {
        Face *face = mesh->getFaceAt(i);
        for( int j = 0; j < 4; j++) {
            Vertex *vertex = face->getNodeAt(j);
            xyz  = vertex->getXYZCoords();
            coords[j][0] = xyz[0];
            coords[j][1] = xyz[1];
            coords[j][2] = xyz[2];
        }
        quality[i] = v_quad_oddy( num_nodes, coords );
    }
    sort( quality.begin(), quality.end() );
    print_histogram( quality, "Oddy", 4 );

    ///////////////////////////////////////////////////////////////////////////

    for( size_t i = 0; i < numfaces; i++) {
        Face *face = mesh->getFaceAt(i);
        for( int j = 0; j < 4; j++) {
            Vertex *vertex = face->getNodeAt(j);
            xyz  = vertex->getXYZCoords();
            coords[j][0] = xyz[0];
            coords[j][1] = xyz[1];
            coords[j][2] = xyz[2];
        }
        quality[i] = v_quad_condition( num_nodes, coords );
    }
    sort( quality.begin(), quality.end() );
    print_histogram( quality, "ConditionNumber", 4 );

    ///////////////////////////////////////////////////////////////////////////

    for( size_t i = 0; i < numfaces; i++) {
        Face *face = mesh->getFaceAt(i);
        for( int j = 0; j < 4; j++) {
            Vertex *vertex = face->getNodeAt(j);
            xyz  = vertex->getXYZCoords();
            coords[j][0] = xyz[0];
            coords[j][1] = xyz[1];
            coords[j][2] = xyz[2];
        }
        quality[i] = v_quad_jacobian( num_nodes, coords );
    }
    sort( quality.begin(), quality.end() );
    print_histogram( quality, "Jacobian", 4 );

    ///////////////////////////////////////////////////////////////////////////

    for( size_t i = 0; i < numfaces; i++) {
        Face *face = mesh->getFaceAt(i);
        for( int j = 0; j < 4; j++) {
            Vertex *vertex = face->getNodeAt(j);
            xyz  = vertex->getXYZCoords();
            coords[j][0] = xyz[0];
            coords[j][1] = xyz[1];
            coords[j][2] = xyz[2];
        }
        quality[i] = v_quad_scaled_jacobian( num_nodes, coords );
    }
    sort( quality.begin(), quality.end() );
    print_histogram( quality, "ScaledJacobian", 4 );

    ///////////////////////////////////////////////////////////////////////////

    for( size_t i = 0; i < numfaces; i++) {
        Face *face = mesh->getFaceAt(i);
        for( int j = 0; j < 4; j++) {
            Vertex *vertex = face->getNodeAt(j);
            xyz  = vertex->getXYZCoords();
            coords[j][0] = xyz[0];
            coords[j][1] = xyz[1];
            coords[j][2] = xyz[2];
        }
        quality[i] = v_quad_shear( num_nodes, coords );
    }
    sort( quality.begin(), quality.end() );
    print_histogram( quality, "Shear", 4 );

    ///////////////////////////////////////////////////////////////////////////

    for( size_t i = 0; i < numfaces; i++) {
        Face *face = mesh->getFaceAt(i);
        for( int j = 0; j < 4; j++) {
            Vertex *vertex = face->getNodeAt(j);
            xyz  = vertex->getXYZCoords();
            coords[j][0] = xyz[0];
            coords[j][1] = xyz[1];
            coords[j][2] = xyz[2];
        }
        quality[i] = v_quad_shape( num_nodes, coords );
    }
    sort( quality.begin(), quality.end() );
    print_histogram( quality, "Shape", 4 );

    ///////////////////////////////////////////////////////////////////////////

    for( size_t i = 0; i < numfaces; i++) {
        Face *face = mesh->getFaceAt(i);
        for( int j = 0; j < 4; j++) {
            Vertex *vertex = face->getNodeAt(j);
            xyz  = vertex->getXYZCoords();
            coords[j][0] = xyz[0];
            coords[j][1] = xyz[1];
            coords[j][2] = xyz[2];
        }
        quality[i] = v_quad_relative_size_squared( num_nodes, coords );
    }
    sort( quality.begin(), quality.end() );
    print_histogram( quality, "RelativeSizeSquared", 4 );

    ///////////////////////////////////////////////////////////////////////////

    for( size_t i = 0; i < numfaces; i++) {
        Face *face = mesh->getFaceAt(i);
        for( int j = 0; j < 4; j++) {
            Vertex *vertex = face->getNodeAt(j);
            xyz  = vertex->getXYZCoords();
            coords[j][0] = xyz[0];
            coords[j][1] = xyz[1];
            coords[j][2] = xyz[2];
        }
        quality[i] = v_quad_shape_and_size( num_nodes, coords );
    }
    sort( quality.begin(), quality.end() );
    print_histogram( quality, "ShapeSize", 4 );

    ///////////////////////////////////////////////////////////////////////////

    for( size_t i = 0; i < numfaces; i++) {
        Face *face = mesh->getFaceAt(i);
        for( int j = 0; j < 4; j++) {
            Vertex *vertex = face->getNodeAt(j);
            xyz  = vertex->getXYZCoords();
            coords[j][0] = xyz[0];
            coords[j][1] = xyz[1];
            coords[j][2] = xyz[2];
        }
        quality[i] = v_quad_shear_and_size( num_nodes, coords );
    }
    sort( quality.begin(), quality.end() );
    print_histogram( quality, "ShearSize", 4 );

    ///////////////////////////////////////////////////////////////////////////

    for( size_t i = 0; i < numfaces; i++) {
        Face *face = mesh->getFaceAt(i);
        for( int j = 0; j < 4; j++) {
            Vertex *vertex = face->getNodeAt(j);
            xyz  = vertex->getXYZCoords();
            coords[j][0] = xyz[0];
            coords[j][1] = xyz[1];
            coords[j][2] = xyz[2];
        }
        quality[i] = v_quad_distortion( num_nodes, coords );
    }
    sort( quality.begin(), quality.end() );
    print_histogram( quality, "Distortion", 4 );
    cout << "Quality: Distortion " << endl;
    quality_values( quality );


    ///////////////////////////////////////////////////////////////////////////
#endif
}

void plot_all_tri_quality_measures( Jaal::Mesh *mesh )
{
#ifdef HAVE_VERDICT
    VERDICT_REAL  coords[3][3];
    Point3D xyz;
    int num_nodes = 3;

    size_t numfaces = mesh->getSize(2);

    vector<double> quality( numfaces );

    ///////////////////////////////////////////////////////////////////////////

    for( size_t i = 0; i < numfaces; i++) {
        Face *face = mesh->getFaceAt(i);
        for( int j = 0; j < num_nodes; j++) {
            Vertex *vertex = face->getNodeAt(j);
            xyz  = vertex->getXYZCoords();
            coords[j][0] = xyz[0];
            coords[j][1] = xyz[1];
            coords[j][2] = xyz[2];
        }
        quality[i] =  v_tri_aspect( num_nodes, coords);
    }
    sort( quality.begin(), quality.end() );
    print_histogram( quality, "AspectRatio",3  );

    ///////////////////////////////////////////////////////////////////////////

    for( size_t i = 0; i < numfaces; i++) {
        Face *face = mesh->getFaceAt(i);
        for( int j = 0; j < 3; j++) {
            Vertex *vertex = face->getNodeAt(j);
            xyz  = vertex->getXYZCoords();
            coords[j][0] = xyz[0];
            coords[j][1] = xyz[1];
            coords[j][2] = xyz[2];
        }
        quality[i] = v_tri_area( num_nodes, coords );
    }
    sort( quality.begin(), quality.end() );
    print_histogram( quality, "Area" ,3);

    ///////////////////////////////////////////////////////////////////////////

    for( size_t i = 0; i < numfaces; i++) {
        Face *face = mesh->getFaceAt(i);
        for( int j = 0; j < 3; j++) {
            Vertex *vertex = face->getNodeAt(j);
            xyz  = vertex->getXYZCoords();
            coords[j][0] = xyz[0];
            coords[j][1] = xyz[1];
            coords[j][2] = xyz[2];
        }
        quality[i] = v_tri_minimum_angle( num_nodes, coords );
    }
    sort( quality.begin(), quality.end() );
    print_histogram( quality, "MinAngle" ,3);


    ///////////////////////////////////////////////////////////////////////////

    for( size_t i = 0; i < numfaces; i++) {
        Face *face = mesh->getFaceAt(i);
        for( int j = 0; j < 3; j++) {
            Vertex *vertex = face->getNodeAt(j);
            xyz  = vertex->getXYZCoords();
            coords[j][0] = xyz[0];
            coords[j][1] = xyz[1];
            coords[j][2] = xyz[2];
        }
        quality[i] = v_tri_maximum_angle( num_nodes, coords );
    }
    sort( quality.begin(), quality.end() );
    print_histogram( quality, "MaxAngle" ,3);

    ///////////////////////////////////////////////////////////////////////////

    for( size_t i = 0; i < numfaces; i++) {
        Face *face = mesh->getFaceAt(i);
        for( int j = 0; j < 3; j++) {
            Vertex *vertex = face->getNodeAt(j);
            xyz  = vertex->getXYZCoords();
            coords[j][0] = xyz[0];
            coords[j][1] = xyz[1];
            coords[j][2] = xyz[2];
        }
        quality[i] = v_tri_condition( num_nodes, coords );
    }
    sort( quality.begin(), quality.end() );
    print_histogram( quality, "ConditionNumber" ,3);

    ///////////////////////////////////////////////////////////////////////////

    for( size_t i = 0; i < numfaces; i++) {
        Face *face = mesh->getFaceAt(i);
        for( int j = 0; j < 3; j++) {
            Vertex *vertex = face->getNodeAt(j);
            xyz  = vertex->getXYZCoords();
            coords[j][0] = xyz[0];
            coords[j][1] = xyz[1];
            coords[j][2] = xyz[2];
        }
        quality[i] = v_tri_scaled_jacobian( num_nodes, coords );
    }
    sort( quality.begin(), quality.end() );
    print_histogram( quality, "ScaledJacobian" ,3);

    ///////////////////////////////////////////////////////////////////////////

    for( size_t i = 0; i < numfaces; i++) {
        Face *face = mesh->getFaceAt(i);
        for( int j = 0; j < 3; j++) {
            Vertex *vertex = face->getNodeAt(j);
            xyz  = vertex->getXYZCoords();
            coords[j][0] = xyz[0];
            coords[j][1] = xyz[1];
            coords[j][2] = xyz[2];
        }
        quality[i] = v_tri_shape( num_nodes, coords );
    }
    sort( quality.begin(), quality.end() );
    print_histogram( quality, "Shape" ,3);

    ///////////////////////////////////////////////////////////////////////////

    for( size_t i = 0; i < numfaces; i++) {
        Face *face = mesh->getFaceAt(i);
        for( int j = 0; j < 3; j++) {
            Vertex *vertex = face->getNodeAt(j);
            xyz  = vertex->getXYZCoords();
            coords[j][0] = xyz[0];
            coords[j][1] = xyz[1];
            coords[j][2] = xyz[2];
        }
        quality[i] = v_tri_shape_and_size( num_nodes, coords );
    }
    sort( quality.begin(), quality.end() );
    print_histogram( quality, "ShapeSize" ,3);

    ///////////////////////////////////////////////////////////////////////////

    for( size_t i = 0; i < numfaces; i++) {
        Face *face = mesh->getFaceAt(i);
        for( int j = 0; j < 3; j++) {
            Vertex *vertex = face->getNodeAt(j);
            xyz  = vertex->getXYZCoords();
            coords[j][0] = xyz[0];
            coords[j][1] = xyz[1];
            coords[j][2] = xyz[2];
        }
        quality[i] = v_tri_distortion( num_nodes, coords );
    }
    sort( quality.begin(), quality.end() );
    print_histogram( quality, "Distortion" ,3);

    ///////////////////////////////////////////////////////////////////////////
#endif
}
