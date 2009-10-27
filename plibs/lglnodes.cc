#include <stdio.h>
#include <stdlib.h>

#include <math.h>

#include <iostream>

using namespace std;


int main(int argc, char **argv)
{
   int N = atoi( argv[1] );

   int N1 = N + 1;

   double *x = new double[N1];

   for(int i = 0; i < N1; i++) 
       x[i] = cos( M_PI*i /(double)N );

   double P[N1][N1+1];

   double *xold = new double[N1];
   for(int i = 0; i < N1; i++) 
       xold[i] = 2.0;

   double eps = 1.0E-10;
   while(1 ) 
   {
       double maxerror = 0.0;
       for(int i = 0; i < N1; i++) 
           maxerror = max( maxerror, fabs( x[i] - xold[i] ) );

       if( maxerror < eps ) break;

       for(int i = 0; i < N1; i++) xold[i] = x[i];

       for(int i = 0; i < N1; i++) {
           P[i][1] = 1.0;
           P[i][2] = x[i];
       }

       for(int k = 2; k <= N; k++)  {
          for(int i = 0; i < N1; i++) {
              P[i][k+1] = ( ( 2*k-1)*x[i]*P[i][k] - (k-1)*P[i][k-1])/k;
          }
       }

      for(int i = 0; i < N1; i++)  {
           x[i] = xold[i] - (x[i]*P[i][N1]-P[i][N])/(N1*P[i][N1]);
      }
   }

   for(int i = 0; i < N1; i++)  
         cout << fixed << x[i] << endl;

   

}

   



  


