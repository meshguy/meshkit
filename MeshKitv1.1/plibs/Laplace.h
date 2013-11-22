#ifndef LAPLACE_H
#define LAPLACE_H

#include <vector>
#include <math.h>

using namespace std;

class Laplacian
{
public:
    int solve2d(int nx, int ny, vector<double> &phi);

    int solve2d(int nx, int ny, const vector<double> &x, 
              const vector<double> &y, vector<double> &phi);

    int solve3d(int nx, int ny, int nz, vector<double> &phi);
    int solve3d(int nx, int ny, int nz, const vector<double> &x, 
              const vector<double> &y, const vector<double> &z, 
              vector<double> &phi);


private:

};


#endif



 
