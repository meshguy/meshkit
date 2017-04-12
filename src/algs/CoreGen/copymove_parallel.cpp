/*********************************************
 June,10
 Reactor Assembly Mesh Assembler
 Argonne National Laboratory

 CCrgen class functions for cop/moving meshes
 based on symmetry and geometry type
*********************************************/

#include "meshkit/CoreGen.hpp"
using namespace MeshKit;
int CoreGen::copymove(const int nrank, const int numprocs)
// ---------------------------------------------------------------------------
// Function: copy/move the assemblies based on the geometrytype and symmetry - Assume 1 meshfile in each instance
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
    if (nrank < ((int) core_alias.size() + ((int) files.size() - nassys))){
        err = set_copymove_coords();
        // ERRORR("Failed to set cm coords.", err);

        // now copy/move
        err = copymove_all(nrank, numprocs);
        // ERRORR("Failed to cm hexflat.", err);
    }
    return iBase_SUCCESS;
}

int CoreGen::set_copymove_coords()
// ---------------------------------------------------------------------------
// Function:
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
    x_coord.resize(tot_assys);
    y_coord.resize(tot_assys);
    int i = 0;
    int assm_index;
    double dx[3] = {0.0, 0.0, 0.0};
    // move the assys based on the geometry type
    if (!strcmp(geom_type.c_str(), "hexflat") && symm == 6) {
        for (int n1 = 0; n1 < nrings; n1++) {
            for (int n2 = 0; n2 < n1 + 1; n2++) {
                err = find_assm(i, assm_index);
                y_coord[i] = n1 * pitch * sin(PII / 3.0) - n2 * pitch * sin(PII / 3.0);
                x_coord[i] = n1 * pitch * cos(PII / 3.0) + n2 * pitch * cos(PII / 3.0);
                i++;
            }
        }
    }
    if (!strcmp(geom_type.c_str(), "rectangular") && symm == 1) {
        double dx = 0.0, dy = 0.0;
        int i = 0;
        for (int m = 1; m <= nringsy; m++) {
            for (int n = 1; n <= nringsx; n++) {
                if (n==1){
                    dx = 0;
                    if(m==1)
                      dy = 0;
                  }
                else{
                    dx+= pitchx;
                  }
                if (m > 1 && n==1){
                    dy+= pitchy;
                  }
                y_coord[i] = -dy;
                x_coord[i] = dx;
                i++;
              }
          }
      }
    if (!strcmp(geom_type.c_str(), "hexflat") && symm == 1) {
        int t, width = 2 * nrings - 1;
        for (int n1 = 1; n1 <= width; n1++) {
            if(n1 > nrings)
                t = 2 * nrings - n1;
            else
                t = n1;

            for (int n2 = 1; n2 <= (nrings + t - 1); n2++) {
                err = find_assm(i,assm_index);

                if (n1 < nrings){
                    x_coord[i] = (nrings - n2 + 1) * pitch / 2.0 + n2 * pitch / 2.0 +
                            (n2 - 1) * pitch - (n1 - 1) * pitch / 2.0;
                    y_coord[i] = -((n1 - 1) * (0.5 * pitch / sin(PII/3.0) + 0.5 * pitch * sin(PII/6.0) / sin(PII/3.0)));
                }
                else{
                    x_coord[i] = (nrings - n2 + 1) * pitch / 2.0 + n2 * pitch / 2.0 + (n2 - 1) * pitch -
                            (2 * nrings - n1 -1) * pitch / 2.0;
                    y_coord[i] = -((n1 -1) * (0.5 * pitch / sin(PII/3.0) + 0.5 * pitch * sin(PII/6.0) / sin(PII/3.0)));
                }
                i++;
            }
        }
    }

    if (!strcmp(geom_type.c_str(), "hexflat") && symm == 12) {
        int flag = 0;
        for (int n1 = 0; n1 < nrings; n1++) {
            int loc = (n1 + 2)/2;

            if( flag == 0 ){
                dx[0] = (n1 + loc - 1) * pitch / 2.0;
                dx[1] = (n1 - loc + 1) * pitch * sin(PII/3.0);
                flag = 1;
            }
            else{
                dx[0] = (n1 + loc) * pitch / 2.0;
                dx[1] = (n1 - loc) * pitch * sin(PII/3.0);
                flag = 0;
            }

            for (int n2 = 0; n2 < loc; n2++) {
                err = find_assm(i,assm_index);
                y_coord[i] = dx[1] - n2 * pitch * sin(PII/3.0);
                x_coord[i] = dx[0] + n2 * pitch * cos(PII/3.0);
                i++;
            }
        }
    }

    if (!strcmp(geom_type.c_str(), "hexvertex") && symm == 6) {
        int bd = 0;
        for (int n1 = 0; n1 < nrings; n1++) {
            if(n1%2==0){//check if n1 is even
                for (int n2 = 0; n2 < n1+1; n2++) {

                    err = find_assm(i, assm_index);
                    if (-1 == assm_index){
                        i++;
                        if(n2 > (n1+1)/2)
                            ++bd; // index for assemblies below diagonal needs updatation
                        continue;
                    }
                    if(n2 <= n1/2){// before or equal to diagonal
                        dx[0] = n2 * pitch;
                        dx[1] = n1 * pitch * sin(PII/3.0);
                    }
                    else{//below the diagonal
                        dx[0] = (n1 + 1 + bd) * pitch / 2.0;
                        dx[1] = (n1 - 1 - bd) * pitch * sin(PII/3.0);
                        ++bd;
                    }
                    x_coord[i] = (dx[0] * cos(PII/6.0) + dx[1] * sin(PII/6.0));
                    y_coord[i] = (dx[1] * cos(PII/6.0) - dx[0] * sin(PII/6.0));
                    i++;
                }
            }
            else{//n1 is odd
                for (int n2 = 0; n2 < n1; n2++) {
                    err = find_assm(i, assm_index);
                    if (-1 == assm_index){
                        i++;
                        if(n2 > (n1+1)/2)
                            ++bd; // index for assemblies below diagonal needs updatation
                        continue;

                    }
                    if(n2 <= (n1-1)/2){// before or equal to diagonal
                        dx[0] = (2 * n2 + 1) * pitch / 2.0;
                        dx[1] = n1 * pitch * sin(PII/3.0);

                    }
                    else{//below the diagonal
                        dx[0] = (n1 + 1 + bd) * pitch / 2.0;
                        if (bd == 0) // first n2 = 1 assembly
                            dx[1] = pitch * sin(PII/3.0);
                        dx[1] = (n1 - 1 - bd) * pitch * sin(PII/3.0);
                        ++bd;
                    }


                    // starting from x-axis
                    x_coord[i] = (dx[0] * cos(PII/6.0) + dx[1] * sin(PII/6.0));
                    y_coord[i] = (dx[1] * cos(PII/6.0) - dx[0] * sin(PII/6.0));
                    i++;
                }
            }
            bd = 0;
        }
    }
    return iBase_SUCCESS;
}

int CoreGen::copymove_all(const int nrank, const int numprocs)
// ---------------------------------------------------------------------------
// Function:
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
    if(prob_type == "mesh"){
        // get the copy/expand sets
        const char *ctag_names[] = { "GEOM_DIMENSION" };
        const char *ctag_vals[] = { (const char*) &set_DIM };

        int flag = 1;
        int assm_index = -1;
        double dx_orig[3], dx[3];

        if(numprocs <= (int) files.size()){
            // no distribution of task for copy/move; each file loaded only once
            int *flags = new int[assys.size()];
            int move_index = -1;
            int *run_count = new int[assys.size()];
            double dx[3] = { 0.0, 0.0, 0.0 };
            double dx_move[3] = { 0.0, 0.0, 0.0 };
            CMatrix<double> dx_orig(assys.size(), 3);
            dx_orig.Set(0.0);

            for (unsigned int i = 0; i < assys.size(); i++) {
                flags[i]=0;
                run_count[i]=0;
                cm[i]->expand_sets().add_tag("MATERIAL_SET");
                cm[i]->expand_sets().add_tag("DIRICHLET_SET");
                cm[i]->expand_sets().add_tag("NEUMANN_SET");

                cm[i]->copy_sets().add_tag("GEOM_DIMENSION");

                cm[i]->update_sets();
              }
            for (int k=0; k< tot_assys; k++){
                err = find_assm(k, assm_index);
                if(assm_index >= 0){
                    // check if this file is with the proc
                    int out = 0;
                    for(int c=0; c < (int) assys.size(); c++){
                        if(assys_index[c] == assm_index){
                            //found assembly
                            move_index=c;
                            out = 0;
                            break;
                        }
                        else{
                            out = 1;
                        }
                    }
                    if(out == 1){
                        continue;
                    }
                    if (-1 == assm_index) {
                        continue;
                    }

                    if (flags[move_index] == 0) {
                        dx_move[0] = x_coord[k];
                        dx_move[1] = y_coord[k];
                        dx_move[2] = 0;
                        dx_orig(move_index + 1, 1) = dx_move[0];
                        dx_orig(move_index + 1, 2) = dx_move[1];
                        dx_orig(move_index + 1, 3) = dx_move[2];

                        move_verts(assys[move_index], dx_move);
                        logfile << "Moved Assembly: " << assm_index  << " dX = " << dx_move[0] << " dY = "
                                  << dx_move[1] << " in proc with rank " << nrank << std::endl;
                        if(strcmp(info.c_str(),"on") == 0)
                            info_file << assm_index << " \t" << k  << " \t" << dx_move[0] << " \t" << dx_move[1]  << " \t" << dx_move[2]  << " \t" << nrank << std::endl;

                        flags[move_index]=1;
                    }
                    else{

                        dx[0] =  x_coord[k] - dx_orig(move_index+1, 1);
                        dx[1] =  y_coord[k] - dx_orig(move_index+1, 2);
                        dx[2] =  0.0;
                        cm[move_index]->set_transform(Copy::Translate(dx));
                        cm[move_index]->execute_this();
                        ++run_count[move_index];

                        logfile << "Copy/moved A: " << assm_index
                                  <<" dX = " <<dx[0]<< " dY = " << dx[1] << " rank " << nrank << std::endl;
                        if(strcmp(info.c_str(),"on") == 0)
                            info_file << assm_index << " " << k  << "  " << dx[0] << " \t" << dx[1]  << " \t" << dx[2]  << " \t" << nrank << std::endl;
                        cm[move_index]->tag_copied_sets(ctag_names, ctag_vals, 1);
                    }

                }
            }
            delete[] flags;
            delete[] run_count;
        }
        else{
            for(int i =0; i < (int) position_core[nrank].size(); i++){
                assm_index = position_core[nrank][i];
                if(assm_index >= 0){
                    // some files are loaded in two or more processors, copy/move task distribution takes place
                    int j = 0;
                    cm[j]->expand_sets().add_tag("MATERIAL_SET");
                    cm[j]->expand_sets().add_tag("DIRICHLET_SET");
                    cm[j]->expand_sets().add_tag("NEUMANN_SET");

                    cm[j]->copy_sets().add_tag("GEOM_DIMENSION");

                    cm[j]->update_sets();

                    if(flag == 0){
                        dx[0] = x_coord[assm_index] - dx_orig[0];
                        dx[1] = y_coord[assm_index] - dx_orig[1];
                        dx[2] = 0.0;

                        cm[0]->set_transform(Copy::Translate(dx));
                        cm[0]->execute_this();
                        logfile << "Copy/moved Assm: " << assm_index << " dX = " << dx[0] << " dY = "
                                  << dx[1]  << " rank " << nrank << std::endl;
                        if(strcmp(info.c_str(),"on") == 0)
                            info_file << assm_index << " \t" << i  << " \t" << dx[0] << " \t" << dx[1]  << " \t" << dx[2]  << " \t" << nrank << std::endl;
                      } else {
                        flag = 0;
                        dx_orig[0] = x_coord[assm_index];
                        dx_orig[1] = y_coord[assm_index];
                        dx_orig[2] = 0;
                        move_verts(assys[0], dx_orig);
                        logfile << "Moved Assm: " << assm_index << " dX = " << dx_orig[0] << " dY = "
                                  << dx_orig[1] << " rank " << nrank << std::endl;
                        if(strcmp(info.c_str(),"on") == 0)
                            info_file << assm_index << " \t" << i  << " \t" << dx_orig[0] << " \t" << dx_orig[1]  << " \t" << dx_orig[2]  << " \t" << nrank << std::endl;

                    }
                }
            }
        }
    }
    else{ // prob type is geometry
        int assm_index = -1;
        // no distribution of task for copy/move; each file loaded only once
        int *flags = new int[assys.size()];
        int move_index = -1;
        double dx[3] = { 0.0, 0.0, 0.0 };
        double dx_move[3] = { 0.0, 0.0, 0.0 };
        CMatrix<double> dx_orig(assys.size(), 3);
        dx_orig.Set(0.0);

        for (unsigned int i = 0; i < assys.size(); i++) {
            flags[i]=0;
        }

        for (int k=0; k< tot_assys; k++){
            err = find_assm(k, assm_index);
            if(assm_index >= 0){
                // check if this file is with the proc
                int out = 0;
                for(int c=0; c < (int) assys.size(); c++){
                    if(assys_index[c] == assm_index){
                        //found assembly
                        move_index=c;
                        out = 0;
                        break;
                    }
                    else{
                        out = 1;
                    }
                }
                if(out == 1){
                    continue;
                }
                if (-1 == assm_index) {
                    continue;
                }

                if (flags[move_index] == 0) {
                    dx_move[0] = x_coord[k];
                    dx_move[1] = y_coord[k];
                    dx_move[2] = 0;
                    dx_orig(move_index + 1, 1) = dx_move[0];
                    dx_orig(move_index + 1, 2) = dx_move[1];
                    dx_orig(move_index + 1, 3) = dx_move[2];

                    move_geoms(assys[move_index], dx_move);
                    logfile << "Moved Assembly: " << assm_index  << " dX = " << dx_move[0] << " dY = "
                              << dx_move[1] << " in proc with rank " << nrank << std::endl;
                    flags[move_index]=1;
                }
                else{

                    dx[0] =  x_coord[k] - dx_orig(move_index+1, 1);
                    dx[1] =  y_coord[k] - dx_orig(move_index+1, 2);
                    dx[2] =  0.0;

                    cg[move_index]->set_location(dx);
                    cg[move_index]->execute_this();
                    logfile << "Copy/moved A: " << assm_index
                              <<" dX = " <<dx[0]<< " dY = " << dx[1] << " rank " << nrank << std::endl;
                }

            }
        }
        delete[] flags;
    }
    return iBase_SUCCESS;
}
