#include "EquipotentialSmooth.hpp"
#include <algorithm>


//---------------------------------------------------------------------------//
// construction function for EquipotentialSmooth class
EquipotentialSmooth::EquipotentialSmooth()
{
	

}



//---------------------------------------------------------------------------//
// deconstruction function for EquipotentialSmooth class
EquipotentialSmooth::~EquipotentialSmooth()
{

}

//setup the data for EquipotentialSmooth
void EquipotentialSmooth::SetupData(std::vector<std::set<int> > AdjElements, std::vector<std::vector<int> > Quads, std::vector<std::vector<double> > coords, std::vector<bool> isBnd, std::vector<std::vector<int> > conn)
{
    //pass the node coordinates to the local variables
    coordinates.resize(coords.size());
    for (unsigned int i = 0; i < coords.size(); i++){
        coordinates[i].resize(coords[i].size());
        for (unsigned int j = 0; j < coords[i].size(); j++)
	    coordinates[i][j] = coords[i][j];
    }

    //pass the connectivity information to the local variables
    adjElements.resize(AdjElements.size());
    for (unsigned int i = 0; i < AdjElements.size(); i++){
	for (std::set<int>::iterator it = AdjElements[i].begin(); it != AdjElements[i].end(); it++)
	    adjElements[i].insert(*it);
    }

    quads.resize(Quads.size());
    for (unsigned int i = 0; i < Quads.size(); i++){
	quads[i].resize(Quads[i].size());
	for (unsigned int j = 0; j < Quads[i].size(); j++)
	    quads[i][j] = Quads[i][j];
    }

    //pass the boundary info to the local variables
    isBoundary.resize(isBnd.size());
    for (unsigned int i = 0; i < isBnd.size(); i++)
	isBoundary[i] = isBnd[i];

	//pass the connectivity to the local variables
	connect.resize(conn.size());
	for (unsigned int i = 0; i < conn.size(); i++){
		connect[i].resize(conn[i].size());
		for (unsigned int j = 0; j < conn[i].size(); j++)
			connect[i][j] = conn[i][j];
	}
    //done
}

//return the results
void EquipotentialSmooth::GetCoords(std::vector<std::vector<double> > &coords)
{
    //return the results to the user
    coords.resize(coordinates.size());
    for (unsigned int i = 0; i < coords.size(); i++){
	if (!isBoundary[i]){
           coords[i].resize(coordinates[i].size());
	   for (unsigned int j = 0; j < coordinates[i].size(); j++)
	       coords[i][j] = coordinates[i][j];
	}
    }
}

//Execute function
void EquipotentialSmooth::Execute()
{
    /*   Algorithm --- Iterative smoothing   */
    /*   p'(i) = p(i) + sum{w(n)*(p(n)-p(i)), n = 1...N}     elements 1...N are the neighbors of node i  */
    /*   w(1) = -belta/2		w(2) = alpha		w(3) = belta/2		w(4) = gamma     				 */
    /*   w(5) = -belta/2		w(6) = alpha		w(7) = belta/2		w(8) = gamma                     */
    /*   where alpha = xp^2 + yp^2 + zp^2 							         							 */
    /*	       belta = xp*xq + yp*yq + zp*zq							         						 */
    /*         gamma = xq^2 + yq^2 + zq^2								 								 */
    /*         xp = (x2-x6)/2	yp = (y2-y6)/2	 zp = (z2-z6)/2						  					 */
    /*         xq = (x8-x4)/2	yq = (y8-y4)/2	 zq = (z8-z4)/2						 				     */

    //		3------------2------------1
    //      |            |            |
    //		|            |            |
    //		|            |            |
    //		4------------i------------8
    //		|            |            |	
    //		|            |            |
    //		|            |            |
    //      5------------6------------7

    double epilson = 0.01;
    double error = 0.0;
    int step = 0;
    while(true){
		error = 0.0;
		//loop over all the nodes
		for (unsigned int i = 0; i < coordinates.size(); i++){
			if (!(isBoundary[i])){
				double weight[9];
				double alpha = 0.0, belta = 0.0, gamma = 0.0;				
				double p[3] = {0.0, 0.0, 0.0};
				double q[3] = {0.0, 0.0, 0.0};
				for (int j = 0; j < 3; j++){
					p[j] = (coordinates[connect[i][2]][j] - coordinates[connect[i][6]][j])/2.0;
					q[j] = (coordinates[connect[i][4]][j] - coordinates[connect[i][8]][j])/2.0;
				}
				alpha = pow(p[0],2) + pow(p[1],2) + pow(p[2],2);
				gamma = pow(q[0],2) + pow(q[1],2) + pow(q[2],2);
				belta = p[0]*q[0] + p[1]*q[1] + p[2]*q[2];

				weight[1] = belta/2.0;
				weight[2] = gamma;
				weight[3] = -belta/2.0;
				weight[4] = alpha;
				weight[5] = belta/2.0;
				weight[6] = gamma;
				weight[7] = -belta/2.0;
				weight[8] = alpha;

				double sum_dis[3] = {0.0, 0.0, 0.0};
				for (int j = 0; j < 3; j++){
					for (int k = 1; k < 9; k++){
						sum_dis[j]+=weight[k]*coordinates[connect[i][k]][j]/(2.0*(alpha+gamma));
					}
				}
				double e = fabs(sum_dis[0] - coordinates[i][0]) + fabs(sum_dis[1] - coordinates[i][1]) + fabs(sum_dis[2]- coordinates[i][2]);
				for (int j = 0; j < 3; j++)
					coordinates[i][j] = sum_dis[j];
				if(e > error) error = e;	
			}
		}
		step++;

		std::cout << "EquipotentialSmooth smoothing step = " << step << "\tError = " << error << std::endl;
		if (error  < epilson) break;
    }	
}

