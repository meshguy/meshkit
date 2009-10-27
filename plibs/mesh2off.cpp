#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#include <fstream>
#include <vector>
#include <set>

using namespace std;


int main( int argc, char **argv)
{
   if( argc < 3) {
       cout << "Usage: executable infile outfile " << endl;
       return 1;
   }

   ifstream ifile( argv[1], ios::in);
   if( ifile.fail() ) {
       cout << "Cann't open file " << argv[1] << endl;
       return 1;
   }

   ofstream ofile( argv[2], ios::out);
   if( ofile.fail() ) {
       cout << "Cann't open file " << argv[2] << endl;
       return 1;
   }

   ofstream partfile( "model.part", ios::out);
   if( partfile.fail() ) {
       cout << "Cann't open file " << endl;
       return 1;
   }

   int numNodes;
   vector<double> xCoord, yCoord, zCoord;

   ifile >> numNodes;
   xCoord.resize(numNodes);
   yCoord.resize(numNodes);
   zCoord.resize(numNodes);

   for( int i = 0; i < numNodes; i++)
        ifile >> xCoord[i] >> yCoord[i] >> zCoord[i];

   int numCells;
   int n0, n1, n2, n3;
   int partid;

   ifile >> numCells; 
   for( int i = 0; i < numCells; i++)
        ifile >> partid >> n0 >> n1 >> n2 >> n3;

   int numFaces;
   ifile >> numFaces; 
   cout << numFaces << endl;

   ofile << "OFF" << endl;
   ofile << numNodes << "  " << numFaces <<  " 0 " << endl;
   for( int i = 0; i < numNodes; i++)
        ofile << xCoord[i] << " " << yCoord[i] << "  " << zCoord[i] << endl;

   vector<int> part;
   set<int>    partSet;
   part.resize( numFaces );
   for( int i = 0; i < numFaces; i++) {
        ifile >> partid >>  n0 >> n1 >> n2;
        ofile << "3 " << n0-1 << "  " << n1-1 << "  " << n2-1 << endl;
        part[i] = partid;
        partSet.insert(partid);
   }

   partfile << numFaces <<  "  " << partSet.size() << endl;
   for( int i = 0; i < numFaces; i++) 
        partfile << part[i] << endl;

}

   


