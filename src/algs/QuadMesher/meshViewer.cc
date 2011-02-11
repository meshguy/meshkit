#include "meshViewer.h"
#ifdef CSV
#include <assert.h>
#include <fstream>

using namespace std;

void Viewer::normalize()
{
  size_t numndes = nodes.size();

  float xmin, xmax, ymin, ymax, zmin, zmax;
  xmin = nodes[0].xyz[0];
  xmax = nodes[0].xyz[0];
  ymin = nodes[0].xyz[1];
  ymax = nodes[0].xyz[1];
  zmin = nodes[0].xyz[1];
  zmax = nodes[0].xyz[1];

  for( int i = 0; i < numnodes; i++)  {
       xmin = min( xmin, nodes[i].xyz[0] );
       xmax = max( xmax, nodes[i].xyz[0] );
       ymin = min( ymin, nodes[i].xyz[1] );
       ymax = max( ymax, nodes[i].xyz[1] );
       zmin = min( zmin, nodes[i].xyz[2] );
       zmax = max( zmax, nodes[i].xyz[2] );
  }

  double xlen  = fabs(xmax-xmin);
  double ylen  = fabs(ymax-ymin);
  double zlen  = fabs(zmax-zmin);
  double scale = max( max( xlen, ylen ), zlen);

  for( int i = 0; i < numnodes; i++)  
  {
       nodes[i].xyz[0] /= scale;
       nodes[i].xyz[1] /= scale;
       nodes[i].xyz[2] /= scale;
  }
}

///////////////////////////////////////////////////////////////////////////////

void Viewer::readData( const string &fname)
{
  ifstream infile( fname.c_str(), ios::in);
  if( infile.fail() )  {
      cout << "Warning: cann't open node file " << fname << endl;
      return;
  }

  int numnodes, numfaces;
  infile >> numnodes >> numfaces;

  nodes.resize( numnodes );

  double x, y, z;
  for( int i = 0; i < numnodes; i++)  {
       infile >> x >> y >> z ;
       nodes[i].xyz[0]  = x;
       nodes[i].xyz[1]  = y;
       nodes[i].xyz[2]  = z;
  }
 
  normalize();

  faces.resize( numfaces);

  int nnodes;
  for( int i = 0; i < numfaces; i++) {
       infile >> nnodes;
       faces[i].connect.resize( nnodes);
       for( int j = 0; j < nnodes; j++) 
            infile >> faces[i].connect[j];
  }

  int itag;
  float r, g, b;
  Color clr;

  clr.rgb[0] = 1.0;
  clr.rgb[1] = 1.0;
  clr.rgb[2] = 1.0;
  colormap[0] = clr;

  clr.rgb[0] = 1.0;
  clr.rgb[1] = 0.0;
  clr.rgb[2] = 0.0;
  colormap[1] = clr;

  clr.rgb[0] = 0.0;
  clr.rgb[1] = 1.0;
  clr.rgb[2] = 0.0;
  colormap[2] = clr;

  clr.rgb[0] = 0.0;
  clr.rgb[1] = 0.0;
  clr.rgb[2] = 1.0;
  colormap[3] = clr;

  for( int i = 0; i < numnodes; i++) {
       infile >> itag;
       nodes[i].iTag = itag;
  }

  for( int i = 0; i < numfaces; i++) {
       infile >> itag;

       if( colormap.find(itag) == colormap.end() ) {
	   r  = drand48() ;  if( r < 0.2) r = 0.2;
	   g  = drand48() ;  if( g < 0.2) g = 0.2;
	   b  = drand48() ;  if( b < 0.2) b = 0.2;
	   clr.rgb[0] = r;
	   clr.rgb[1] = g;
	   clr.rgb[2] = b;
	   colormap[itag] = clr;
       }
       faces[i].iTag = itag;
  }

}
///////////////////////////////////////////////////////////////////////////////

void Viewer:: draw_mesh()
{
 size_t numfaces = faces.size();

 glBegin(GL_TRIANGLES);
  for (size_t i=0; i < numfaces; i++) 
  {
    int nnodes = faces[i].connect.size();
    if( nnodes == 3 ) {
        for( int j = 0; j < 3; j++) 
        {
             int n0 = faces[i].connect[j];
             glVertex3fv( nodes[n0].xyz );
        }
     }
  }
  glEnd();

 glBegin(GL_QUADS);
  for (size_t i = 0; i < numfaces; i++) 
  {
    int nnodes = faces[i].connect.size();
    if( nnodes == 4 ) {
        for( int j = 0; j < 4; j++) 
        {
             int n0 = faces[i].connect[j];
             glVertex3fv( nodes[n0].xyz );
        }
    }
  }
 glEnd();

 glBegin(GL_POLYGON);
  for (int i=0; i < numfaces; i++) 
  {
    int nnodes = faces[i].connect.size();
    if( nnodes > 4 ) {
    for( int j = 0; j < nnodes; j++) 
      {
         int n0 = faces[i].connect[j];
         glVertex3fv( nodes[n0].xyz );
      }
    }
  }
  glEnd();
}

///////////////////////////////////////////////////////////////////////////////

void Viewer::draw()
{
  qglviewer::Vec wc, sc;
  glClearColor( 1.0, 1.0, 1.0, 0.0);

  int numfaces = faces.size();

  glDisable(GL_LIGHTING);

  glPointSize(2.0);
  glColor3f( 0.0, 0.0, 0.0);
  for( int j = 0; j < nodes.size(); j++) {
       glBegin(GL_POINTS);
       if( nodes[j].iTag == 0) 
          glColor3f( 0.0, 1.0, 0.0);
       else
          glColor3f( 1.0, 0.0, 0.0);
       glVertex3fv( nodes[j].xyz );
       glEnd();
  }

  glColor3f( 0.0, 0.0, 0.0);
  glLineWidth(1.0);
  glPolygonMode( GL_FRONT_AND_BACK, GL_LINE);

  for (int i=0; i < numfaces; i++) 
  {
      glBegin(GL_POLYGON);
         for( int j = 0; j < faces[i].connect.size(); j++) 
         {
              int n0 = faces[i].connect[j];
              glVertex3fv( nodes[n0].xyz );
         }
       glEnd();
  }

  glPolygonMode( GL_FRONT_AND_BACK, GL_FILL);

/*

  for (int i=0; i < numfaces; i++) 
  {
      int fnodes = faces[i].connect.size();
      
      if( fnodes == 4)
          glColor3f( 0.0, 1.0, 0.0);
      else
          glColor3f( 1.0, 0.0, 0.0);
      glBegin(GL_POLYGON);
         for( int j = 0; j < fnodes; j++) 
         {
              int n0 = faces[i].connect[j];
              glVertex3fv( nodes[n0].xyz );
         }
       glEnd();
  }
  */

  for (int i=0; i < numfaces; i++) 
  {
      int fnodes = faces[i].connect.size();
      int itag   = faces[i].iTag;
      Color clr  = colormap[itag];
      
      glColor3f( clr.rgb[0], 
                 clr.rgb[1], 
		 clr.rgb[2] );

      glBegin(GL_POLYGON);
         for( int j = 0; j < fnodes; j++) 
         {
              int n0 = faces[i].connect[j];
              glVertex3fv( nodes[n0].xyz );
         }
       glEnd();
  }

  glColor3f( 1.0, 0.0, 0.0);
  /*
  for (int i=0; i < numfaces; i++) 
  {
         wc.x = 0.0;
         wc.y = 0.0;
         wc.z = 0.0;
	 int nsize = faces[i].connect.size();
         for( int j = 0; j < nsize; j++) 
         {
              int n0 = faces[i].connect[j];
              wc.x += nodes[n0].xyz[0];
              wc.y += nodes[n0].xyz[1];
              wc.z += nodes[n0].xyz[2];
         }
         wc.x /= (double)nsize;
         wc.y /= (double)nsize;
         wc.z /= (double)nsize;
         sc   = camera()->projectedCoordinatesOf(wc);
         drawText(int(sc.x), int(sc.y), QString::number(i));
  }


  glColor3f( 0.0, 0.0, 0.0);
  for (int i=0; i < nodes.size(); i++) 
  {
       wc.x = nodes[i].xyz[0];
       wc.y = nodes[i].xyz[1];
       wc.z = nodes[i].xyz[2];
       sc   = camera()->projectedCoordinatesOf(wc);
       drawText(int(sc.x), int(sc.y), QString::number(i));
  }
  */
}

void Viewer::init()
{
    int type = 2;
    if (type < 3)
    {   
      // Move camera according to viewer type (on X, Y or Z axis)
      camera()->setPosition(Vec((type==0)? 1.0 : 0.0, (type==1)? 1.0 : 0.0, (type==2)? 1.0 : 0.0));
      camera()->lookAt(sceneCenter());

      camera()->setType(Camera::ORTHOGRAPHIC);
      camera()->showEntireScene();

      // Forbid rotation
      WorldConstraint* constraint = new WorldConstraint();
      constraint->setRotationConstraintType(AxisPlaneConstraint::FORBIDDEN);
      camera()->frame()->setConstraint(constraint);
    }   


  readData( filename );
  // Restore previous viewer state.
//  restoreStateFromFile();
  
  // Opens help window
  help();
}

QString Viewer::helpString() const
{
  QString text("<h2>S i m p l e V i e w e r</h2>");
  text += "Use the mouse to move the camera around the object. ";
  text += "You can respectively revolve around, zoom and translate with the three mouse buttons. ";
  text += "Left and middle buttons pressed together rotate around the camera view direction axis<br><br>";
  text += "Pressing <b>Alt</b> and one of the function keys (<b>F1</b>..<b>F12</b>) defines a camera keyFrame. ";
  text += "Simply press the function key again to restore it. Several keyFrames define a ";
  text += "camera path. Paths are saved when you quit the application and restored at next start.<br><br>";
  text += "Press <b>F</b> to display the frame rate, <b>A</b> for the world axis, ";
  text += "<b>Alt+Return</b> for full screen mode and <b>Control+S</b> to save a snapshot. ";
  text += "See the <b>Keyboard</b> tab in this window for a complete shortcut list.<br><br>";
  text += "Double clicks automates single click actions: A left button double click aligns the closer axis with the camera (if close enough). ";
  text += "A middle button double click fits the zoom of the camera and the right button re-centers the scene.<br><br>";
  text += "A left button double click while holding right button pressed defines the camera <i>Revolve Around Point</i>. ";
  text += "See the <b>Mouse</b> tab and the documentation web pages for details.<br><br>";
  text += "Press <b>Escape</b> to exit the viewer.";
  return text;
}
#endif

int main(int argc, char** argv)
{
#ifdef CSV
   assert( argc == 2 );
  // Read command lines arguments.
  QApplication application(argc,argv);

/*
  // Instantiate the viewer.
  Viewer viewer;
  viewer.setDataFile( argv[1] );

#if QT_VERSION < 0x040000
  // Set the viewer as the application main widget.
  application.setMainWidget(&viewer);
#else
  viewer.setWindowTitle("simpleViewer");
#endif

  // Make the viewer window visible on screen.
  viewer.show();

  // Run main loop.
  return application.exec();
*/
#endif
}

