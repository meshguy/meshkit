#include <stdio.h>
#include <iostream>
#include <algorithm>
#include <vector>
#include <list>

#include <vector>

using namespace std;

int main()

{
  vector<int>  connect(4);

  connect[0] = 1;
  connect[1] = 2;
  connect[2] = 3;
  connect[3] = 4;

  std::reverse( connect.begin(), connect.end() );

  for( int i = 0; i < 4; i++)
       cout << connect[i] << endl;



}

  
