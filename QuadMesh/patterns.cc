#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>

#include <vector>

using namespace std;

struct Vertex
{
  int  id;
  double coords[2];
};

struct Edge
{
   Vertex *connect[2];
};

struct Shape2D
{
  double center[2];
  vector<Vertex*> nodes;
  vector<Edge> lineSegments;
};

struct Circle : public Shape2D
{
  double radius;
  int    numSegments;

  void discretize()
  {
     nodes.resize(numSegments);

     double dtheta = 2*M_PI/(double)numSegments;
     double angle  = 0.0;
     for( int i = 0; i < numSegments; i++) {
          double x = center[0] + radius*cos( angle );
          double y = center[1] + radius*sin( angle );
	  angle   += dtheta;
	  Vertex *v= new Vertex;
	  v->coords[0] = x;
	  v->coords[1] = y;
	  nodes[i]  = v;
     }

     lineSegments.resize(numSegments);
     for( int i = 0; i < numSegments; i++) {
          lineSegments[i].connect[0] =  nodes[(i+0)%numSegments];
          lineSegments[i].connect[1] =  nodes[(i+1)%numSegments];
     }
  }

};

struct Square : public Shape2D
{
  double length;
  double center[2];
  int    numSegments;
};

void write_pslg( const string &fname, vector<Shape2D*> &shapes)
{
    string filename = fname  + ".poly";
    ofstream ofile( filename.c_str(), ios::out);

    size_t index = 0;
    for( int i = 0; i < shapes.size(); i++) {
         int nsize = shapes[i]->nodes.size();
	 for( int j = 0; j < nsize; j++) { 
	      shapes[i]->nodes[j]->id = index+1;
	      index++;
         }
    }

    size_t numNodes = index;

    ofile << numNodes << " 2 0 0 " << endl;

    index = 1;
    for( int i = 0; i < shapes.size(); i++) {
         int nsize = shapes[i]->nodes.size();
	 for( int j = 0; j < nsize; j++) 
	      ofile << index++ << " " << shapes[i]->nodes[j]->coords[0] << " " 
	                              << shapes[i]->nodes[j]->coords[1] << endl;
   }

   int numSegments = 0;
   for( int i = 0; i < shapes.size(); i++) 
        numSegments += shapes[i]->lineSegments.size();

   ofile << numSegments << " 0 " << endl;

   index = 1;
   for( int i = 0; i < shapes.size(); i++) {
        int nsize = shapes[i]->lineSegments.size();
        for( int j = 0; j < nsize; j++) 
	      ofile << index++ << " " 
	            << shapes[i]->lineSegments[j].connect[0]->id << " " 
	            << shapes[i]->lineSegments[j].connect[1]->id << endl;
   }

   int numHoles = shapes.size()-1;

   ofile << numHoles << endl;

   index = 1;
   for( int i = 0; i < shapes.size()-1; i++) 
        ofile << index++ << " " << shapes[i]->center[0] << " " << shapes[i]->center[1] << endl; 

}

void gen_circle_pattern( )
{
   Circle *circle;
   vector<Shape2D*> shapes;

/*
   circle = new Circle;
   circle->radius    = 1.2;
   circle->center[0] = 0.0;
   circle->center[1] = 0.0;
   circle->numSegments = 20;
   circle->discretize();
   shapes.push_back(circle);
*/

   int nCount = 8;
   double radius1 = 10.00, radius2  = 50.0, holeRadius = 2.0;
   int nRings = 5;

   double dr  = (radius2-radius1)/(double)nRings;

   for( int i = 0; i <= nRings; i++) {
      double angle = 0;
      double radius = radius1 + i*dr;
      cout << radius << endl;
      for( int j = 0; j < nCount; j++) {
           circle = new Circle;
           circle->radius    = holeRadius;
           circle->center[0] = radius*cos(angle);
           circle->center[1] = radius*sin(angle);
	   angle += 2.0*M_PI/(double)nCount;
           circle->numSegments = 32;
           circle->discretize();
           shapes.push_back(circle);
      }
      nCount += 6;
   }

   circle = new Circle;
   circle->radius    = radius2 + dr;
   circle->center[0] = 0.0;
   circle->center[1] = 0.0;
   circle->numSegments = 200;
   circle->discretize();
   shapes.push_back(circle);

   write_pslg( "circles" , shapes );
}


int main()
{
  gen_circle_pattern();
}

