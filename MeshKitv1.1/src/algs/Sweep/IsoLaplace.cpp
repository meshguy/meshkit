#include "IsoLaplace.hpp"
#include <algorithm>


//---------------------------------------------------------------------------//
// construction function for IsoLaplace class
IsoLaplace::IsoLaplace()
{
	

}



//---------------------------------------------------------------------------//
// deconstruction function for IsoLaplace class
IsoLaplace::~IsoLaplace()
{

}

//setup the data for IsoLaplace
void IsoLaplace::SetupData(std::vector<std::set<int> > AdjElements, std::vector<std::vector<int> > Quads, std::vector<std::vector<double> > coords, std::vector<bool> isBnd, std::vector<double> w)
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

    //pass the element's weight to the local variables
    weight.resize(w.size());
    for (unsigned int i = 0; i < w.size(); i++)
	weight[i] = w[i];
    //done
}

//return the results
void IsoLaplace::GetCoords(std::vector<std::vector<double> > &coords)
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
void IsoLaplace::Execute()
{

    //test the input
    for (unsigned int i = 0; i < coordinates.size(); i++){
		if (!isBoundary[i]){
			std::cout << "IsoLaplace index = " << i << "\t[" << coordinates[i][0] << ", " << coordinates[i][1] << ", " << coordinates[i][2] << "]\n";
			std::cout << "\n\n";
		}
    }
    std::cout << "test the adjacent element info\n";
    for (unsigned int i = 0; i < adjElements.size(); i++){
		std::set<int>::iterator it = adjElements[i].begin();
		std::cout << "node index = " << i << "\t";	
		for (; it != adjElements[i].end(); it++)
	    	std::cout << *it << "\t";
		std::cout << "\n";
    }
    
    std::cout << "test the quad info\n";
    for (unsigned int i = 0; i < quads.size(); i++){
		std::cout << "quad index = " << i << "\nthe adjacent nodes are ";
		for (unsigned int j = 0; j < quads[i].size(); j++){
	    	std::cout << quads[i][j] << "\t";	
		}
		std::cout << std::endl;
    }
    std::cout << "test the weights\n";
    for (unsigned int i = 0; i < weight.size(); i++){
        std::cout << "index = " << i << "\tweight is " << weight[i] << std::endl;
    }
	



    double epilson = 0.01;
    double error = 0.0;
    int step = 0;
    double r = 0.1;
    while(true){
		error = 0.0;
		//loop over all the nodes
		for (unsigned int i = 0; i < coordinates.size(); i++){
	    	if (!(isBoundary[i])){
				double sum_neighbors[3] = {0.0, 0.0, 0.0};
				double sum_opposite[3] = {0.0, 0.0, 0.0};
				double sum_weight = 0.0;
				//loop over the adjacent quads
				for (std::set<int>::iterator it = adjElements[i].begin(); it != adjElements[i].end(); it++){
		    		int stNodeIndex = -1;
		    		if (int(i) == quads[*it][0])		
		    			stNodeIndex = 0;
		    		else if (int(i) == quads[*it][1])
		    			stNodeIndex = 1;
		    		else if (int(i) == quads[*it][2])
		    			stNodeIndex = 2;
		    		else
		    			stNodeIndex = 3;
		    		//sum up the neighbor nodes and opposite nodes
					for (int k = 0; k < 3; k++){
						sum_neighbors[k] += weight[*it]*coordinates[quads[*it][(stNodeIndex+1)%4]][k];
						sum_neighbors[k] += weight[*it]*coordinates[quads[*it][(stNodeIndex+3)%4]][k];
						sum_opposite[k] += weight[*it]*r*coordinates[quads[*it][(stNodeIndex+2)%4]][k];
					}
		    		sum_weight += weight[*it];
				}
				double tmp[3];
				double e = 0.0;
				//average the sum_neighbors and sum_opposite
				for (int k = 0; k < 3; k++){
					sum_neighbors[k] = sum_neighbors[k]/sum_weight;
					sum_opposite[k] = sum_opposite[k]/sum_weight;
				}
				//update the node position and calculate the error
				for (int j = 0; j < 3; j++){
					tmp[j] = (sum_neighbors[j] - sum_opposite[j])/(2-r);
					e += fabs(tmp[j] - coordinates[i][j]);
					coordinates[i][j] = tmp[j];
				}
				if(e > error) 
					error = e;
	    	}
		}
		step++;
		//UpdateWeight();

		std::cout << "IsoLaplace smoothing step = " << step << "\tError = " << error << std::endl;
		if (error  < epilson) break;
    }


	for (unsigned int i = 0; i < coordinates.size(); i++){
		if (!isBoundary[i]){
			std::cout << "IsoLaplace index = " << i << "\t[" << coordinates[i][0] << ", " << coordinates[i][1] << ", " << coordinates[i][2] << "]\n";
			std::cout << "\n\n";
		}
    }

	std::cout << std::endl;
	
}

void IsoLaplace::UpdateWeight(){
    for (unsigned int i = 0; i < quads.size(); i++){
	//calculate the area
	double a, b, c, s;
	weight[i] = 0.0;
	//calculate one half triangle
	a = sqrt(pow(coordinates[quads[i][0]][0] - coordinates[quads[i][1]][0], 2) + pow(coordinates[quads[i][0]][1] - coordinates[quads[i][1]][1],2) + pow(coordinates[quads[i][0]][2] - coordinates[quads[i][1]][2], 2));
	b = sqrt(pow(coordinates[quads[i][0]][0] - coordinates[quads[i][2]][0], 2) + pow(coordinates[quads[i][0]][1] - coordinates[quads[i][2]][1],2) + pow(coordinates[quads[i][0]][2] - coordinates[quads[i][2]][2], 2));
	c = sqrt(pow(coordinates[quads[i][1]][0] - coordinates[quads[i][2]][0], 2) + pow(coordinates[quads[i][1]][1] - coordinates[quads[i][2]][1],2) + pow(coordinates[quads[i][1]][2] - coordinates[quads[i][2]][2], 2));
	s = 0.5*(a + b + c);
	weight[i] += sqrt(fabs(s*(s-a)*(s-b)*(s-c)));

	//calculate the other half triangle
	a = sqrt(pow(coordinates[quads[i][0]][0] - coordinates[quads[i][3]][0], 2) + pow(coordinates[quads[i][0]][1] - coordinates[quads[i][3]][1],2) + pow(coordinates[quads[i][0]][2] - coordinates[quads[i][3]][2], 2));
	b = sqrt(pow(coordinates[quads[i][3]][0] - coordinates[quads[i][2]][0], 2) + pow(coordinates[quads[i][3]][1] - coordinates[quads[i][2]][1],2) + pow(coordinates[quads[i][3]][2] - coordinates[quads[i][2]][2], 2));
	c = sqrt(pow(coordinates[quads[i][1]][0] - coordinates[quads[i][2]][0], 2) + pow(coordinates[quads[i][1]][1] - coordinates[quads[i][2]][1],2) + pow(coordinates[quads[i][1]][2] - coordinates[quads[i][2]][2], 2));
	s = 0.5*(a + b + c);
	weight[i] += sqrt(fabs(s*(s-a)*(s-b)*(s-c)));
     }
}

