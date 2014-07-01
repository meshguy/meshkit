#ifndef JAAL_MESH
#define JAAL_MESH


#include <meshkit/qglviewer.h>
#include <meshkit/qapplication.h>

#ifdef CSV
#include <string>
#include <vector>
#include <map>

using namespace std;
using namespace qglviewer;

class Viewer : public QGLViewer
{
public:
     void setDataFile( const string &f)
     { 
       filename = f; 
     }
protected :
  struct Color
  {
     float rgb[3];
  };
  struct Vertex
  {
    int  iTag;
     float xyz[3];
  };

  struct Edge
  {
     int boundid;
     int connect[2];
  };

  struct Face
  {
    int  iTag;
    vector<int>  connect;
  };

  string filename;
  void readData( const string &f);
  int facetype;
  vector<Vertex> nodes;
  vector<Edge>   edges;
  vector<Face>   faces;
  std::map<int,int>   global2local;
  std::map<int,Color> colormap;

  virtual void draw();
  virtual void init();
  virtual QString helpString() const;
private:
   void readNodes( const string &f);
   void readEdges( const string &f);
   void readFaces( const string &f);
   void normalize();
   void draw_mesh();
};

#endif

#endif
