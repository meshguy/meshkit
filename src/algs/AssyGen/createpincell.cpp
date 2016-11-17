
#include "meshkit/AssyGen.hpp"

namespace MeshKit
{

  void AssyGen::CreatePinCell(int i, double dX, double dY, double dZ)
  //---------------------------------------------------------------------------
  //Function: Create pincell i in location dX dY and dZ
  //Input:    none
  //Output:   none
  //---------------------------------------------------------------------------
  {
    int nRadii=0, nCyl=0, nCells = 0;
    double dCylMoveX = 0.0, dCylMoveY = 0.0, dHeightTotal = 0.0;
    double dHeight =0.0,dZMove = 0.0, PX = 0.0,PY = 0.0,PZ = 0.0, dP=0.0;
    CVector<double> dVCylZPos(2), dVCylXYPos(2), dVStartZ, dVEndZ;;
    CVector<std::string> szVMatName, szVMatAlias, szVCellMat;
    iBase_EntityHandle cell = NULL, cyl= NULL, tmp_vol= NULL,tmp_vol1= NULL, tmp_new= NULL;
    std::vector<iBase_EntityHandle> cp_in;
    // name tag handle
    iBase_TagHandle this_tag= NULL;
    char* tag_name = (char*)"NAME";

    std::string sMatName = "";
    std::string sMatName1 = "";
    int nDuctIndex = -1;

    if(strcmp(m_szInfo.c_str(),"on") == 0){
        std::ostringstream os;
        pin_name = "_xp";
        os << (m_nTotalPincells + m_nStartpinid - 1);
        os << "_";
        std::string pid = os.str(); //retrieve as a string
        pin_name+=pid;
      }

    // get tag handle for 'NAME' tag, already created as iGeom instance is created
#if defined (HAVE_ACIS) || defined (HAVE_OCC)
    iGeom_getTagHandle(igeomImpl->instance(), tag_name, &this_tag, &err, 4);
#endif

    // get cell material
    m_Pincell(i).GetCellMatSize(nCells);
    SimpleArray<iBase_EntityHandle> cells(nCells);

    // branch when cells are present
    if(nCells > 0){
        dVStartZ.SetSize(nCells);
        dVEndZ.SetSize(nCells);
        szVCellMat.SetSize(nCells);
        m_Pincell(i).GetCellMat(dVStartZ, dVEndZ, szVCellMat);

        // get cylinder data
        m_Pincell(i).GetNumCyl(nCyl);

        for(int n=1;n<=nCells; n++){
            // get cylinder locations
            if(nCyl > 0){
                m_Pincell(i).GetCylZPos(n, dVCylZPos);
                nDuctIndex = -1;

                // get the index for cp_inpins based on Z-heights
                for (int dd = 1; dd <= m_nDuct; dd++){
                    if((m_dMZAssm(dd, 2)) >= (dVCylZPos(2)) && (m_dMZAssm(dd, 1)) >= (dVCylZPos(1)))
                      nDuctIndex = dd;
                    if (nDuctIndex != -1)
                      break;
                  }
              }
            dHeight = fabs(dVEndZ(n) - dVStartZ(n));
            if(m_szGeomType =="hexagonal"){

                m_Pincell(i).GetPitch(dP, dHeightTotal); // this dHeight is not used in creation
                double dSide = dP/(sqrt(3));

                if(nCells >0){
                    // create prism
#if defined (HAVE_ACIS) || defined (HAVE_OCC)
                    iGeom_createPrism(igeomImpl->instance(), dHeight, 6,
                                      dSide, dSide,
                                      &cell, &err);
#endif
                    m_PyCubGeomFile << "cell = cubit.prism( " << dHeight << ", 6, " << dSide << ", " << dSide << ")" << std::endl;
                  }
              }
            // if rectangular geometry
            if(m_szGeomType =="rectangular"){

                m_Pincell(i).GetPitch(PX, PY, PZ);
                m_PyCubGeomFile << "cells = []" << std::endl;
                m_PyCubGeomFile << "sub1 = [] \nsub2 = []" << std::endl;

                if(nCells >0){
#if defined (HAVE_ACIS) || defined (HAVE_OCC)
                    // create brick
                    iGeom_createBrick( igeomImpl->instance(),PX,PY,dHeight,&cell,&err );
#endif
                    m_PyCubGeomFile << "cell = cubit.brick( " << PX << ", " << PY << ", " << dHeight << ")" << std::endl;
                  }
              }

            dZMove = (dVStartZ(n)+dVEndZ(n))/2.0;
            if(nCells > 0){
                // position the brick in assembly
#if defined (HAVE_ACIS) || defined (HAVE_OCC)
                iGeom_moveEnt(igeomImpl->instance(), cell, dX, dY, dZMove, &err);
#endif
                m_PyCubGeomFile << "vector = [" << dX << ", " << dY << ", " << dZMove << "]" << std::endl;
                m_PyCubGeomFile << "cubit.move( cell, vector)" << std::endl;

                m_PyCubGeomFile << "cells.append(cell)" << std::endl;
                cells[n-1]=cell;

                //search for the full name of the abbreviated Cell Mat and set name
                for(int p=1;p<= m_szAssmMatAlias.GetSize();p++){
                    if(strcmp (szVCellMat(n).c_str(), m_szAssmMatAlias(p).c_str()) == 0){
                        sMatName = m_szAssmMat(p);
                      }
                  }
                std::cout << "created: " << sMatName << std::endl;
#if defined (HAVE_ACIS) || defined (HAVE_OCC)
                iGeom_setData(igeomImpl->instance(), cell, this_tag,
                              sMatName.c_str(), sMatName.size(), &err);
#endif
                m_PyCubGeomFile  << "lid = cells[0].id()" << std::endl;
                m_PyCubGeomFile << "cubit.set_entity_name(\"body\", lid, \""  << sMatName <<  "\" )" << std::endl;

                                   if(strcmp(m_szInfo.c_str(),"on") == 0){
#if defined (HAVE_ACIS) || defined (HAVE_OCC)
                    iGeom_setData(igeomImpl->instance(), cell, this_tag,
                                  pin_name.c_str(), pin_name.size(), &err);
#endif
                    m_PyCubGeomFile << "cubit.set_entity_name(\"body\", lid, \""  << pin_name <<  "\" )" << std::endl;
                    std::cout << "Naming pin body :" <<  pin_name << std::endl;
                  }


                Name_Faces(sMatName, cell, this_tag);
                m_PyCubGeomFile << "name_faces(\"" << sMatName << "\", cell) " << std::endl;
              }
            // loop and create cylinders
            if(nCyl > 0){
                m_Pincell(i).GetCylSizes(n, nRadii);
                SimpleArray<iBase_EntityHandle> cyls(nRadii);

                //declare variables
                CVector<double> dVCylRadii(2*nRadii);
                CVector<std::string> szVMat(nRadii);
                CVector<std::string> szVCylMat(nRadii);
                int nType = 0;
                //get values
                m_Pincell(i).GetCylRadii(n, dVCylRadii);
                m_Pincell(i).GetCylPos(n, dVCylXYPos);
                m_Pincell(i).GetCylMat(n, szVCylMat);
                m_Pincell(i).GetCylZPos(n, dVCylZPos);
                m_Pincell(i).GetCellType(n, nType);

                dHeight = dVCylZPos(2)-dVCylZPos(1);

                for (int m=1; m<=nRadii; m++){

                    if (nType == 0){
#if defined (HAVE_ACIS) || defined (HAVE_OCC)
                        iGeom_createCylinder(igeomImpl->instance(), dHeight, dVCylRadii(m), dVCylRadii(m),
                                             &cyl, &err);
#endif
                        m_PyCubGeomFile << "#\n#\ncyl = cubit.cylinder(" << dHeight << ", " << dVCylRadii(m) << ", " << dVCylRadii(m) << ", " << ", " << dVCylRadii(m) << ")" << std::endl;
                        std::cout << m << ": Creating cylinder with radii " << dVCylRadii(m) << std::endl;
                      }
                    else{
#if defined (HAVE_ACIS) || defined (HAVE_OCC)
                        iGeom_createCone(igeomImpl->instance(), dHeight, dVCylRadii(2*m-1), dVCylRadii(2*m-1), dVCylRadii(2*m),
                                         &cyl, &err);
#endif
                        m_PyCubGeomFile << "#\n#\ncyl = cubit.cylinder(" << dHeight << ", " << dVCylRadii(2*m-1) << ", " << dVCylRadii(2*m-1) << ", " << dVCylRadii(2*m) << ")" << std::endl;
                      }
                    // move their centers and also move to the assembly location  ! Modify if cyl is outside brick
                    dCylMoveX = dVCylXYPos(1)+dX;
                    dCylMoveY = dVCylXYPos(2)+dY;
                    dZMove = (dVCylZPos(1)+dVCylZPos(2))/2.0;
#if defined (HAVE_ACIS) || defined (HAVE_OCC)
                    iGeom_moveEnt(igeomImpl->instance(), cyl, dCylMoveX,dCylMoveY,dZMove, &err);
#endif
                    m_PyCubGeomFile << "cyls.append(cyl)" << std::endl;

                    m_PyCubGeomFile << "vector = [" << dCylMoveX << ", " << dCylMoveY << ", " << dZMove << "]" << std::endl;
                    m_PyCubGeomFile << "cubit.move(cyl, vector)" << std::endl;
                    cyls[m-1] = cyl;
                  }

                m_PyCubGeomFile << "sub1[:] = []\nsub2[:] = []" << std::endl;

                if(nCells > 0){
                    // copy cyl before subtract
#if defined (HAVE_ACIS) || defined (HAVE_OCC)
                    iGeom_copyEnt(igeomImpl->instance(), cyls[nRadii-1], &tmp_vol, &err);
#endif
                    m_PyCubGeomFile << "tmp_vol = cubit.copy_body(cyls[" << nRadii-1 << "])" << std::endl;

                    // subtract outer most cyl from brick
#if defined (HAVE_ACIS) || defined (HAVE_OCC)
                    iGeom_subtractEnts(igeomImpl->instance(), cells[n-1], tmp_vol, &tmp_new, &err);
#endif
                    m_PyCubGeomFile << "sub1.append(cells[" << n-1 << "])\nsub2.append(tmp_vol)" << std::endl;
                    m_PyCubGeomFile << "sub1[:] = []\nsub2[:] = []" << std::endl;

                    m_PyCubGeomFile << "cells[n-1] = tmp_new \ncell = tmp_new" << std::endl;

                    // copy the new into the cyl array
                    cells[n-1] = tmp_new; cell = tmp_new;

                  }
                cp_in.push_back(tmp_new);
                m_PyCubGeomFile << "cp_in.append(tmp_new[0])" << std::endl;

                //set tag on inner most cylinder, search for the full name of the abbreviated Cell Mat
                for(int p=1;p<=m_szAssmMatAlias.GetSize();p++){
                    if(strcmp (szVCylMat(1).c_str(), m_szAssmMatAlias(p).c_str()) == 0){
                        sMatName = m_szAssmMat(p);
                      }
                  }
                tmp_vol1=cyls[0]; //inner most cyl
                cp_in.push_back(tmp_vol1);
                m_PyCubGeomFile << "tmpvol1 = cyls[0]\ncp_in.append(tmp_vol1)" << std::endl;
#if defined (HAVE_ACIS) || defined (HAVE_OCC)
                iGeom_setData(igeomImpl->instance(), tmp_vol1, this_tag,
                              sMatName.c_str(), 10, &err);
#endif
                m_PyCubGeomFile  << "lid = tmp_vol1.id()" << std::endl;
                m_PyCubGeomFile << "cubit.set_entity_name(\"body\", lid, \""  << sMatName <<  "\" )" << std::endl;

                 if(strcmp(m_szInfo.c_str(),"on") == 0){
#if defined (HAVE_ACIS) || defined (HAVE_OCC)
                    iGeom_setData(igeomImpl->instance(), tmp_vol1, this_tag,
                                  pin_name.c_str(), pin_name.size(), &err);
#endif
                    std::cout << "Naming pin body :" <<  pin_name << std::endl;
                  }

                Name_Faces(sMatName, tmp_vol1, this_tag);
                m_PyCubGeomFile << "name_faces(\"" << sMatName << "\", tmp_vol1) " << std::endl;

                // other cyl annulus after substraction
                for (int b=nRadii; b>1; b--){
#if defined (HAVE_ACIS) || defined (HAVE_OCC)
                    iGeom_copyEnt(igeomImpl->instance(), cyls[b-2], &tmp_vol, &err);
#endif
                    m_PyCubGeomFile << "tmp_vol = cubit.copy_body(cyls[" << b-2 << "])" << std::endl;

                    //subtract tmp vol from the outer most
#if defined (HAVE_ACIS) || defined (HAVE_OCC)
                    iGeom_subtractEnts(igeomImpl->instance(), cyls[b-1], tmp_vol, &tmp_new, &err);
#endif
                    m_PyCubGeomFile << "sub1.append(cyls[" << b-1 << "])\nsub2.append(tmp_vol)" << std::endl;
                    m_PyCubGeomFile << "tmp_new = cubit.subtract(sub2, sub1)" << std::endl;
                    m_PyCubGeomFile << "sub1[:] = []\nsub2[:] = []" << std::endl;

                    // now search for the full name of the abbreviated Cell Mat
                    //    int tag_no;
                    for(int p=1;p<=m_szAssmMatAlias.GetSize();p++){
                        if(strcmp (szVCylMat(b).c_str(), m_szAssmMatAlias(p).c_str()) == 0){
                            //        tag_no = p;
                            sMatName =  m_szAssmMat(p);
                          }
                      }
                    std::cout << "created: " << sMatName << std::endl;
                    cp_in.push_back(tmp_new);
                    m_PyCubGeomFile << "cp_in.append(tmp_new[0])" << std::endl;

                    // set the name of the annulus
#if defined (HAVE_ACIS) || defined (HAVE_OCC)
                    iGeom_setData(igeomImpl->instance(), tmp_new, this_tag,
                                  sMatName.c_str(),sMatName.size(), &err);
#endif
                    m_PyCubGeomFile  << "lid = tmp_new[0].id()" << std::endl;
                    m_PyCubGeomFile << "cubit.set_entity_name(\"body\", lid, \""  << sMatName <<  "\" )" << std::endl;

                    if(strcmp(m_szInfo.c_str(),"on") == 0){
#if defined (HAVE_ACIS) || defined (HAVE_OCC)
                        iGeom_setData(igeomImpl->instance(), tmp_new, this_tag,
                                      pin_name.c_str(), pin_name.size(), &err);
#endif
                        std::cout << "Naming pin body :" <<  pin_name<< std::endl;

                        m_PyCubGeomFile << "cubit.set_entity_name(\"body\", lid, \""  << pin_name <<  "\" )" << std::endl;
                      }
                    Name_Faces(sMatName, tmp_new, this_tag);
                    m_PyCubGeomFile << "name_faces(\"" << sMatName << "\", tmp_new1) " << std::endl;
                    m_PyCubGeomFile << "cyls[" << b-1 << "] = tmp_new" << std::endl;

                    // copy the new into the cyl array
                    cyls[b-1] = tmp_new;
                    tmp_vol=NULL;
                  }
              }
            if(nDuctIndex > 0){
                m_PyCubGeomFile << "cp_inpins.append([])" << std::endl;
                for (int count = 0; count < (int) cp_in.size(); count++){
                    cp_inpins[nDuctIndex-1].push_back(cp_in[count]);
                    m_PyCubGeomFile << "cp_inpins[" << nDuctIndex-1 << "].append(cp_in[" << count << "])" << std::endl;
                  }
              }
            cp_in.clear();
            m_PyCubGeomFile << "cp_in[:] =[]" << std::endl;
          }
      }
    // this branch of the routine is responsible for creating cylinders with '0' cells
    if(nCells == 0){

        // get cylinder data
        m_Pincell(i).GetNumCyl(nCyl);
        nCells = nCyl;

        for(int n=1;n<=nCells; n++){
            nDuctIndex = -1;
            if(m_szGeomType =="hexagonal"){

                m_Pincell(i).GetPitch(dP, dHeightTotal); // this dHeight is not used in creation
              }
            // if rectangular geometry
            if(m_szGeomType =="rectangular"){

                m_Pincell(i).GetPitch(PX, PY, PZ);
              }

            // loop and create cylinders
            if(nCyl > 0){
                m_Pincell(i).GetCylSizes(n, nRadii);
                SimpleArray<iBase_EntityHandle> cyls(nRadii);

                //declare variables
                CVector<double> dVCylRadii(2*nRadii);
                CVector<std::string> szVMat(nRadii);
                CVector<std::string> szVCylMat(nRadii);
                int nType = 0;
                //get values
                m_Pincell(i).GetCylRadii(n, dVCylRadii);
                m_Pincell(i).GetCylPos(n, dVCylXYPos);
                m_Pincell(i).GetCylMat(n, szVCylMat);
                m_Pincell(i).GetCylZPos(n, dVCylZPos);
                m_Pincell(i).GetCellType(n, nType);

                dHeight = dVCylZPos(2)-dVCylZPos(1);

                // get the index for cp_inpins based on Z-heights
                for (int dd = 1; dd <= m_nDuct; dd++){
                    if((m_dMZAssm(dd, 2)) >= (dVCylZPos(2)) && (m_dMZAssm(dd, 1)) >= (dVCylZPos(1)))
                      nDuctIndex = dd;
                    if (nDuctIndex != -1)
                      break;
                  }

                m_PyCubGeomFile << "cyls = [] \ncp_in = []" << std::endl;
                m_PyCubGeomFile << "sub1 = [] \nsub2 = []" << std::endl;

                for (int m=1; m<=nRadii; m++){
                    if (nType == 0){
#if defined (HAVE_ACIS) || defined (HAVE_OCC)
                        iGeom_createCylinder(igeomImpl->instance(), dHeight, dVCylRadii(m), dVCylRadii(m),
                                             &cyl, &err);
#endif
                        m_PyCubGeomFile << "#\n#\ncyl = cubit.cylinder(" << dHeight << ", " << dVCylRadii(m) << ", " << dVCylRadii(m) << ", " << dVCylRadii(m) << ")" << std::endl;
                      }
                    else{
#if defined (HAVE_ACIS) || defined (HAVE_OCC)
                        iGeom_createCone(igeomImpl->instance(), dHeight, dVCylRadii(2*m - 1), dVCylRadii(2*m - 1), dVCylRadii(2*m),
                                         &cyl, &err);
#endif
                        m_PyCubGeomFile << "#\n#\ncyl = cubit.cylinder(" << dHeight << ", " << dVCylRadii(2*m-1) << ", " << dVCylRadii(2*m-1) << ", " << dVCylRadii(2*m) << ")" << std::endl;
                      }

                    // move their centers and also move to the assembly location  ! Modify if cyl is outside brick
                    dCylMoveX = dVCylXYPos(1)+dX;
                    dCylMoveY = dVCylXYPos(2)+dY;
                    dZMove = (dVCylZPos(1)+dVCylZPos(2))/2.0;
#if defined (HAVE_ACIS) || defined (HAVE_OCC)
                    iGeom_moveEnt(igeomImpl->instance(), cyl, dCylMoveX, dCylMoveY, dZMove, &err);
#endif
                    m_PyCubGeomFile << "vector = [" << dCylMoveX << ", " << dCylMoveY << ", " << dZMove << "]" << std::endl;
                    m_PyCubGeomFile << "cubit.move( cyl, vector)" << std::endl;
                    m_PyCubGeomFile << "cyls.append(cyl)" << std::endl;
                    cyls[m-1] = cyl;
                  }

                //set tag on inner most cylinder, search for the full name of the abbreviated Cell Mat
                for(int p=1;p<=m_szAssmMatAlias.GetSize();p++){
                    if(strcmp (szVCylMat(1).c_str(), m_szAssmMatAlias(p).c_str()) == 0){
                        sMatName = m_szAssmMat(p);
                      }
                  }
                std::cout << "created: " << sMatName << std::endl;
                tmp_vol1=cyls[0]; //inner most cyl

                cp_in.push_back(tmp_vol1);
                m_PyCubGeomFile << "tmp_vol1 = cyls[0] \ncp_in.append(tmp_vol1)" << std::endl;
#if defined (HAVE_ACIS) || defined (HAVE_OCC)
                iGeom_setData(igeomImpl->instance(), tmp_vol1, this_tag,
                              sMatName.c_str(), 10, &err);
#endif
                m_PyCubGeomFile  << "lid =cyls[0].id()" << std::endl;
                m_PyCubGeomFile << "cubit.set_entity_name(\"body\", lid, \""  << sMatName <<  "\" )" << std::endl;


                if(strcmp(m_szInfo.c_str(),"on") == 0){
#if defined (HAVE_ACIS) || defined (HAVE_OCC)
                    iGeom_setData(igeomImpl->instance(), tmp_vol1, this_tag,
                                  pin_name.c_str(), pin_name.size(), &err);
#endif
                    m_PyCubGeomFile << "cubit.set_entity_name(\"body\", lid, \""  << pin_name <<  "\" )" << std::endl;

                    std::cout << "Naming pin body :" <<  pin_name<< std::endl;
                  }

                Name_Faces(sMatName, tmp_vol1, this_tag);
                m_PyCubGeomFile << "name_faces(\"" << sMatName << "\", tmp_vol1) " << std::endl;

                // other cyl annulus after substraction
                for (int b=nRadii; b>1; b--){
#if defined (HAVE_ACIS) || defined (HAVE_OCC)
                    iGeom_copyEnt(igeomImpl->instance(), cyls[b-2], &tmp_vol, &err);
#endif
                    m_PyCubGeomFile << "# SUBTRACTING ANNULUS ##\ntmp_vol = cubit.copy_body(cyls[" << b-2 << "])" << std::endl;

                    //subtract tmp vol from the outer most
#if defined (HAVE_ACIS) || defined (HAVE_OCC)
                    iGeom_subtractEnts(igeomImpl->instance(), cyls[b-1], tmp_vol, &tmp_new, &err);
#endif
                    m_PyCubGeomFile << "sub1.append(cyls[" << b-1 << "])\nsub2.append(tmp_vol)" << std::endl;
                    m_PyCubGeomFile << "tmp_new = cubit.subtract(sub2, sub1)" << std::endl;
                    m_PyCubGeomFile << "sub1[:] = []\nsub2[:] = []" << std::endl;

                    // now search for the full name of the abbreviated Cell Mat
                    for(int p=1;p<=m_szAssmMatAlias.GetSize();p++){
                        if(strcmp (szVCylMat(b).c_str(), m_szAssmMatAlias(p).c_str()) == 0){
                            sMatName =  m_szAssmMat(p);
                          }
                      }
                    std::cout <<"created: " << sMatName << std::endl;

                    cp_in.push_back(tmp_new);
                    m_PyCubGeomFile << "cp_in.append(tmp_new[0])" << std::endl;

                    // set the name of the annulus
#if defined (HAVE_ACIS) || defined (HAVE_OCC)
                    iGeom_setData(igeomImpl->instance(), tmp_new, this_tag,
                                  sMatName.c_str(),sMatName.size(), &err);
#endif
                    m_PyCubGeomFile  << "lid = tmp_new[0].id()" << std::endl;
                    m_PyCubGeomFile << "cubit.set_entity_name(\"body\", lid, \""  << sMatName <<  "\" )" << std::endl;
                    if(strcmp(m_szInfo.c_str(),"on") == 0){
#if defined (HAVE_ACIS) || defined (HAVE_OCC)
                        iGeom_setData(igeomImpl->instance(), tmp_new, this_tag,
                                      pin_name.c_str(), pin_name.size(), &err);
#endif
                        m_PyCubGeomFile << "cubit.set_entity_name(\"body\", lid, \""  << pin_name <<  "\" )" << std::endl;

                        std::cout << "Naming pin body :" <<  pin_name<< std::endl;
                      }

                    Name_Faces(sMatName, tmp_new, this_tag);
                    m_PyCubGeomFile << "name_faces(\"" << sMatName << "\", tmp_new[0]) " << std::endl;
                    m_PyCubGeomFile << "cyls[" << b-1 << "] = tmp_new" << std::endl;

                    // copy the new into the cyl array
                    cyls[b-1] = tmp_new;
                  }
              }
            if(nDuctIndex > 0){
                m_PyCubGeomFile << "cp_inpins.append([])" << std::endl;
                for (int count = 0; count < (int) cp_in.size(); count++){
                    cp_inpins[nDuctIndex-1].push_back(cp_in[count]);
                    m_PyCubGeomFile << "cp_inpins["<< nDuctIndex -1 << "].append(cp_in[" << count << "])" << std::endl;

                  }
              }
            cp_in.clear();
            m_PyCubGeomFile << "cp_in[:] =[]" << std::endl;

          }
      }

  }


  void AssyGen::CreatePinCell_Intersect(int i, double dX, double dY, double dZ)
  //---------------------------------------------------------------------------
  //Function: Create pincell i in location dX dY and dZ
  //Input:    none
  //Output:   none
  //---------------------------------------------------------------------------
  {
    int nRadii=0, nCyl=0, nCells = 0;
    double dCylMoveX = 0.0, dCylMoveY = 0.0, dHeightTotal = 0.0;
    double dHeight =0.0,dZMove = 0.0, PX = 0.0,PY = 0.0,PZ = 0.0, dP=0.0;
    CVector<double> dVCylZPos(2), dVCylXYPos(2), dVEndZ, dVStartZ;
    CVector<std::string> szVMatName, szVMatAlias, szVCellMat;
    iBase_EntityHandle cell = NULL, cyl= NULL, tmp_vol1= NULL, tmp_new= NULL, cell_copy = NULL, intersec = NULL;
    std::vector<iBase_EntityHandle> cp_in;

    // name tag handle
    iBase_TagHandle this_tag= NULL;
    char* tag_name = (char*)"NAME";

    std::string sMatName = "";
    std::string sMatName0 = "";
    std::string sMatName1 = "";
    int nDuctIndex = -1;

    if(strcmp(m_szInfo.c_str(),"on") == 0){
        std::ostringstream os;
        pin_name = "_xp";
        os << (m_nTotalPincells + m_nStartpinid - 1);
        os << "_";
        std::string pid = os.str(); //retrieve as a string
        pin_name+=pid;
      }

    // get tag handle for 'NAME' tag, already created as iGeom instance is created
#if defined (HAVE_ACIS) || defined (HAVE_OCC)
    iGeom_getTagHandle(igeomImpl->instance(), tag_name, &this_tag, &err, 4);
#endif

    // get cell material
    m_Pincell(i).GetCellMatSize(nCells);
    m_Pincell(i).GetNumCyl(nCyl);
    SimpleArray<iBase_EntityHandle> cells(nCells);

    // branch when cells are present
    if(nCells > 0){
        dVStartZ.SetSize(nCells);
        dVEndZ.SetSize(nCells);
        szVCellMat.SetSize(nCells);
        m_Pincell(i).GetCellMat(dVStartZ, dVEndZ, szVCellMat);

        // get cylinder data
        m_Pincell(i).GetNumCyl(nCyl);
        m_PyCubGeomFile << "cells = []" << std::endl;

        for(int n=1;n<=nCells; n++){

            dHeight = fabs(dVEndZ(n) - dVStartZ(n));
            // get the index for cp_inpins based on Z-heights
            for (int dd = 1; dd <= m_nDuct; dd++){
                if((m_dMZAssm(dd, 2)) >= (dVCylZPos(2)) && (m_dMZAssm(dd, 1)) >= (dVCylZPos(1)))
                  nDuctIndex = dd;
                if (nDuctIndex != -1)
                  break;
              }
            if(m_szGeomType =="hexagonal"){

                m_Pincell(i).GetPitch(dP, dHeightTotal); // this dHeight is not used in creation

                double dSide = dP/(sqrt(3));

                if(nCells >0){
                    // create prism
#if defined (HAVE_ACIS) || defined (HAVE_OCC)
                    iGeom_createPrism(igeomImpl->instance(), dHeight, 6,
                                      dSide, dSide,
                                      &cell, &err);
#endif
                    m_PyCubGeomFile << "cell = cubit.prism(' " << dHeight << ", 6, " << dSide << ", " << dSide << ")" << std::endl;                  }
              }
            // if rectangular geometry
            if(m_szGeomType =="rectangular"){

                m_Pincell(i).GetPitch(PX, PY, PZ);

                if(nCells >0){
                    // create brick
#if defined (HAVE_ACIS) || defined (HAVE_OCC)
                    iGeom_createBrick( igeomImpl->instance(),PX,PY,dHeight,&cell,&err );
#endif
                    m_PyCubGeomFile << "cell = cubit.brick(' " << PX << ", " << PY << ", " << dHeight << ")" << std::endl;                  }
              }

            dZMove = (dVStartZ(n)+dVEndZ(n))/2.0;

            if(nCells > 0){
                // position the brick in assembly
#if defined (HAVE_ACIS) || defined (HAVE_OCC)
                iGeom_moveEnt(igeomImpl->instance(), cell, dX, dY, dZMove, &err);
#endif
                m_PyCubGeomFile << "vector = [" << dX << ", " << dY << ", " << dZMove << "]" << std::endl;
                m_PyCubGeomFile << "cubit.move( cell, vector)" << std::endl;
                m_PyCubGeomFile << "cells.append(cell)" << std::endl;

                cells[n-1]=cell;
              }
            // loop and create cylinders
            if(nCyl > 0){
                m_Pincell(i).GetCylSizes(n, nRadii);
                SimpleArray<iBase_EntityHandle> cyls(nRadii);
                SimpleArray<iBase_EntityHandle> cell_copys(nRadii);
                SimpleArray<iBase_EntityHandle> intersec_main(nRadii);
                iBase_EntityHandle  tmp_intersec;
                //declare variables
                CVector<double> dVCylRadii(2*nRadii);
                CVector<std::string> szVMat(nRadii);
                CVector<std::string> szVCylMat(nRadii);
                int nType = 0;
                //get values
                m_Pincell(i).GetCylRadii(n, dVCylRadii);
                m_Pincell(i).GetCylPos(n, dVCylXYPos);
                m_Pincell(i).GetCylMat(n, szVCylMat);
                m_Pincell(i).GetCylZPos(n, dVCylZPos);
                m_Pincell(i).GetCellType(n, nType);

                dHeight = dVCylZPos(2)-dVCylZPos(1);

                m_PyCubGeomFile << "cells_copy=[]\nintersec_main=[]" << std::endl;
                m_PyCubGeomFile << "cyls = [] \ncp_in = []" << std::endl;
                m_PyCubGeomFile << "sub1 = [] \nsub2 = []" << std::endl;
                for (int m=1; m<=nRadii; m++){
                    if (nType == 0){
#if defined (HAVE_ACIS) || defined (HAVE_OCC)
                        iGeom_createCylinder(igeomImpl->instance(), dHeight, dVCylRadii(m), dVCylRadii(m),
                                             &cyl, &err);
#endif
                        m_PyCubGeomFile << "#\n#\ncyl = cubit.cylinder(" << dHeight << ", " << dVCylRadii(m) << ", " << dVCylRadii(m) << ", " << ", " << dVCylRadii(m) << ")" << std::endl;
                      }
                    else{
#if defined (HAVE_ACIS) || defined (HAVE_OCC)
                        iGeom_createCone(igeomImpl->instance(), dHeight, dVCylRadii(2*m-1), dVCylRadii(2*m-1), dVCylRadii(2*m),
                                         &cyl, &err);
#endif
                        m_PyCubGeomFile << "#\n#\ncyl = cubit.cylinder(" << dHeight << ", " << dVCylRadii(2*m-1) << ", " << dVCylRadii(2*m-1) << ", " << dVCylRadii(2*m) << ")" << std::endl;
                      }

                    // move their centers and also move to the assembly location  ! Modify if cyl is outside brick
                    dCylMoveX = dVCylXYPos(1)+dX;
                    dCylMoveY = dVCylXYPos(2)+dY;
                    dZMove = (dVCylZPos(1)+dVCylZPos(2))/2.0;
#if defined (HAVE_ACIS) || defined (HAVE_OCC)
                    iGeom_moveEnt(igeomImpl->instance(), cyl, dCylMoveX,dCylMoveY,dZMove, &err);
#endif
                    m_PyCubGeomFile << "vector = [" << dCylMoveX << ", " << dCylMoveY << ", " << dZMove << "]" << std::endl;
                    m_PyCubGeomFile << "cubit.move( cell, vector)" << std::endl;
                    m_PyCubGeomFile << "cyls.append(cyl)" << std::endl;
                    cyls[m-1] = cyl;


                    //copy cell nRadii  times for intersection with cylinders
#if defined (HAVE_ACIS) || defined (HAVE_OCC)
                    iGeom_copyEnt(igeomImpl->instance(), cells[n-1], &cell_copy, &err);
#endif
                    m_PyCubGeomFile << "cell_copy = cubit.copy_body(cells[" << n-1 << "])" << std::endl;
                    m_PyCubGeomFile << "cells_copy.append(cell_copy)" << std::endl;
                    cell_copys[m-1] = cell_copy;
#if defined (HAVE_ACIS) || defined (HAVE_OCC)
                    iGeom_intersectEnts(igeomImpl->instance(), cell_copys[m-1], cyls[m-1],&intersec,&err);
#endif
                    m_PyCubGeomFile << "tmpunite = cubit.unite(cells_copy[" << m-1 << "]\ncyls [" << m-1 << "])" << std::endl;
                    m_PyCubGeomFile << "intersec = cubit.subtract(tmpunite, cells_copy[" << m-1 << "])" << std::endl;
                    m_PyCubGeomFile << "intersec_main.append(intersec)" << std::endl;
                    intersec_main[m-1] = intersec;
                    intersec = NULL;
                  }

                //set tag on inner most cylinder, search for the full name of the abbreviated Cell Mat
                tmp_vol1=intersec_main[0];
                m_PyCubGeomFile << "tmp_vol1 = intersec_main[0]" << std::endl;
                for(int p=1;p<=m_szAssmMatAlias.GetSize();p++){
                    if(strcmp (szVCylMat(1).c_str(), m_szAssmMatAlias(p).c_str()) == 0){
                        sMatName = m_szAssmMat(p);
                      }
                  }
                m_PyCubGeomFile << "tmp_vol1 = cyls[0] \ncp_in.append(tmp_vol1)" << std::endl;
                cp_in.push_back(tmp_vol1);

#if defined (HAVE_ACIS) || defined (HAVE_OCC)
                iGeom_setData(igeomImpl->instance(), tmp_vol1, this_tag,
                              sMatName.c_str(), 10, &err);
#endif
                m_PyCubGeomFile  << "lid =tmp_vol1.id()" << std::endl;
                m_PyCubGeomFile << "cubit.set_entity_name(\"body\", lid, \""  << sMatName <<  "\" )" << std::endl;

                if(strcmp(m_szInfo.c_str(),"on") == 0){
#if defined (HAVE_ACIS) || defined (HAVE_OCC)
                    iGeom_setData(igeomImpl->instance(), tmp_vol1, this_tag,
                                  pin_name.c_str(), pin_name.size(), &err);
#endif
                    std::cout << "Naming pin body :" <<  pin_name<< std::endl;
                  }
                Name_Faces(sMatName, tmp_vol1, this_tag);
                m_PyCubGeomFile << "name_faces(\"" << sMatName << "\", tmp_vol1) " << std::endl;

                // copy the outermost cyl
#if defined (HAVE_ACIS) || defined (HAVE_OCC)
                iGeom_copyEnt(igeomImpl->instance(), intersec_main[nRadii-1], &tmp_intersec, &err);
#endif
                m_PyCubGeomFile << "tmp_intersec = cubit.copy_body(intersec_main[" << nRadii-1 << "])" << std::endl;

                // subtract the outermost cyl from the cell
#if defined (HAVE_ACIS) || defined (HAVE_OCC)
                iGeom_subtractEnts(igeomImpl->instance(), cells[n-1], tmp_intersec, &tmp_new, &err);
#endif
                m_PyCubGeomFile << "sub1.append(tmp_intersec\nsub2.append(cells["<< n-1 << "])" << std::endl;
                m_PyCubGeomFile << "tmp_new = cubit.subtract(sub2, sub1)" << std::endl;
                m_PyCubGeomFile << "sub1[:] = []\nsub2[:] = []" << std::endl;                
                // now search for the full name of the abbreviated Cell Mat
                //  int tag_no;
                for(int p=1;p<=m_szAssmMatAlias.GetSize();p++){
                    if(strcmp (szVCellMat(n).c_str(), m_szAssmMatAlias(p).c_str()) == 0){
                        //      tag_no = p;
                        sMatName =  m_szAssmMat(p);
                      }
                  }
                std::cout << "created: " << sMatName << std::endl;

                cp_in.push_back(tmp_new);
                m_PyCubGeomFile << "cp_in.append(tmp_new)" << std::endl;

                // set the name of the annulus
#if defined (HAVE_ACIS) || defined (HAVE_OCC)
                iGeom_setData(igeomImpl->instance(), tmp_new, this_tag,
                              sMatName.c_str(),sMatName.size(), &err);
#endif
                m_PyCubGeomFile  << "lid =tmp_new.id()" << std::endl;
                m_PyCubGeomFile << "cubit.set_entity_name(\"body\", lid, \""  << sMatName <<  "\" )" << std::endl;
#if defined (HAVE_ACIS) || defined (HAVE_OCC)
                iGeom_setData(igeomImpl->instance(), tmp_new, this_tag,
                              pin_name.c_str(), pin_name.size(), &err);
#endif
                std::cout << "Naming pin body :" <<  pin_name<< std::endl;
                m_PyCubGeomFile << "cubit.set_entity_name(\"body\", lid, \""  << pin_name <<  "\" )" << std::endl;

                Name_Faces(sMatName, tmp_new, this_tag);
                m_PyCubGeomFile << "name_faces(\"" << sMatName << "\", tmp_new) " << std::endl;

                for (int b=nRadii; b>1; b--){
#if defined (HAVE_ACIS) || defined (HAVE_OCC)
                    iGeom_copyEnt(igeomImpl->instance(), intersec_main[b-2], &tmp_intersec, &err);
#endif
                    m_PyCubGeomFile << "tmp_intersec = cubit.copy_body(intersec_main[" << b-1 << "])" << std::endl;

                    //subtract tmp vol from the outer most
#if defined (HAVE_ACIS) || defined (HAVE_OCC)
                    iGeom_subtractEnts(igeomImpl->instance(), intersec_main[b-1], tmp_intersec, &tmp_new, &err);
#endif
                    m_PyCubGeomFile << "sub1.append(tmp_intersec\nsub2.append(intersec_main["<< b-1 << "])" << std::endl;
                    m_PyCubGeomFile << "tmp_new = cubit.subtract(sub2, sub1)" << std::endl;
                    m_PyCubGeomFile << "sub1[:] = []\nsub2[:] = []" << std::endl; 
                    // now search for the full name of the abbreviated Cell Mat
                    //    int tag_no;
                    for(int p=1;p<=m_szAssmMatAlias.GetSize();p++){
                        if(strcmp (szVCylMat(b).c_str(), m_szAssmMatAlias(p).c_str()) == 0){
                            //        tag_no = p;
                            sMatName =  m_szAssmMat(p);
                          }
                      }
                    std::cout << "created: " << sMatName << std::endl;

                    cp_in.push_back(tmp_new);
                    m_PyCubGeomFile << "cp_in.append(tmp_new)" << std::endl;

                    // set the name of the annulus
#if defined (HAVE_ACIS) || defined (HAVE_OCC)
                    iGeom_setData(igeomImpl->instance(), tmp_new, this_tag,
                                  sMatName.c_str(),sMatName.size(), &err);
#endif
                    m_PyCubGeomFile  << "lid =tmp_new.id()" << std::endl;
                    m_PyCubGeomFile << "cubit.set_entity_name(\"body\", lid, \""  << sMatName <<  "\" )" << std::endl;

                    if(strcmp(m_szInfo.c_str(),"on") == 0){
#if defined (HAVE_ACIS) || defined (HAVE_OCC)
                        iGeom_setData(igeomImpl->instance(), tmp_new, this_tag,
                                      pin_name.c_str(), pin_name.size(), &err);
#endif
                        std::cout << "Naming pin body :" <<  pin_name<< std::endl;
                      }
                    Name_Faces(sMatName, tmp_new, this_tag);
                    m_PyCubGeomFile << "name_faces(\"" << sMatName << "\", tmp_new) " << std::endl;
                    m_PyCubGeomFile << "cyls[" << b-1 << "] = tmp_new" << std::endl;

                    // copy the new into the cyl array
                    cyls[b-1] = tmp_new;

                  }
              }
            if(nDuctIndex > 0){
                for (int count = 0; count < (int) cp_in.size(); count++){
                  cp_inpins[nDuctIndex-1].push_back(cp_in[count]);
                  m_PyCubGeomFile << "cp_inpins["<< nDuctIndex -1 << "].append(cp_in[" << count << "])" << std::endl;
                }

              }
            cp_in.clear();
            m_PyCubGeomFile << "cp_in[:] =[]" << std::endl;

          }
      }
    // this branch of the routine is responsible for creating cylinders with '0' cells
    if(nCells == 0){

        // get cylinder data
        m_Pincell(i).GetNumCyl(nCyl);
        nCells = nCyl;
        cells.resize(nCells);

        for(int n=1;n<=nCells; n++){

            // get some cylinder parameters to create the cell material for intersection
            m_Pincell(i).GetCylZPos(n, dVCylZPos);
            dHeight = dVCylZPos(2)-dVCylZPos(1);
            dZMove = (dVCylZPos(1)+dVCylZPos(2))/2.0;

            if(m_szGeomType =="hexagonal"){

                m_Pincell(i).GetPitch(dP, dHeightTotal); // this dHeight is not used in creation
                double dSide = dP/(sqrt(3));
#if defined (HAVE_ACIS) || defined (HAVE_OCC)
                iGeom_createPrism(igeomImpl->instance(), dHeight, 6,
                                  dSide, dSide,
                                  &cell, &err);
#endif
              }
            m_PyCubGeomFile << "cyls = [] \ncp_in = []\ncells = []" << std::endl;
            m_PyCubGeomFile << "sub1 = [] \nsub2 = []" << std::endl;
            // if rectangular geometry
            if(m_szGeomType =="rectangular"){

                m_Pincell(i).GetPitch(PX, PY, PZ);
                // create brick
#if defined (HAVE_ACIS) || defined (HAVE_OCC)
                iGeom_createBrick( igeomImpl->instance(),PX,PY,PZ, &cell,&err );
#endif
                m_PyCubGeomFile << "cell = cubit.brick(' " << PX << ", " << PY << ", " << PZ << ")" << std::endl;                 
              }
#if defined (HAVE_ACIS) || defined (HAVE_OCC)
            iGeom_moveEnt(igeomImpl->instance(), cell, dX, dY, dZMove, &err);
#endif
            m_PyCubGeomFile << "vector = [" << dX << ", " << dY << ", " << dZMove << "]" << std::endl;
            m_PyCubGeomFile << "cubit.move( cell, vector)" << std::endl;
            m_PyCubGeomFile << "cells.append(cell)" << std::endl;
            cells[n-1]=cell;
            // loop and create cylinders
            if(nCyl > 0){
                m_Pincell(i).GetCylSizes(n, nRadii);

                //declare variables
                SimpleArray<iBase_EntityHandle> cyls(nRadii), cell_copys(nRadii), intersec_main(nRadii), intersec_copy(nRadii);
                iBase_EntityHandle  tmp_intersec;
                CVector<double> dVCylRadii(2*nRadii);
                CVector<std::string> szVMat(nRadii), szVCylMat(nRadii);
                int nType = 0;
                //get values
                m_Pincell(i).GetCylRadii(n, dVCylRadii);
                m_Pincell(i).GetCylPos(n, dVCylXYPos);
                m_Pincell(i).GetCylMat(n, szVCylMat);
                m_Pincell(i).GetCellType(n, nType);


                // get the index for cp_inpins based on Z-heights
                for (int dd = 1; dd <= m_nDuct; dd++){
                    if((m_dMZAssm(dd, 2)) >= (dVCylZPos(2)) && (m_dMZAssm(dd, 1)) >= (dVCylZPos(1)))
                      nDuctIndex = dd;
                    if (nDuctIndex != -1)
                      break;
                  }

                for (int m=1; m<=nRadii; m++){
                    if (nType == 0){
#if defined (HAVE_ACIS) || defined (HAVE_OCC)
                        iGeom_createCylinder(igeomImpl->instance(), dHeight, dVCylRadii(m), dVCylRadii(m),
                                             &cyl, &err);
#endif
                        m_PyCubGeomFile << "#\n#\ncyl = cubit.cylinder(" << dHeight << ", " << dVCylRadii(m) << ", " << dVCylRadii(m) << ", " << ", " << dVCylRadii(m) << ")" << std::endl;
                      }
                    else{
#if defined (HAVE_ACIS) || defined (HAVE_OCC)
                        iGeom_createCone(igeomImpl->instance(), dHeight, dVCylRadii(2*m-1), dVCylRadii(2*m-1), dVCylRadii(2*m),
                                         &cyl, &err);
#endif
                        m_PyCubGeomFile << "#\n#\ncyl = cubit.cylinder(" << dHeight << ", " << dVCylRadii(2*m-1) << ", " << dVCylRadii(2*m-1) << ", " << dVCylRadii(2*m) << ")" << std::endl;
                      }

                    // move their centers and also move to the assembly location  ! Modify if cyl is outside brick
                    dCylMoveX = dVCylXYPos(1)+dX;
                    dCylMoveY = dVCylXYPos(2)+dY;
#if defined (HAVE_ACIS) || defined (HAVE_OCC)
                    iGeom_moveEnt(igeomImpl->instance(), cyl, dCylMoveX,dCylMoveY,dZMove, &err);
#endif
                    m_PyCubGeomFile << "vector = [" << dCylMoveX << ", " << dCylMoveY << ", " << dZMove << "]" << std::endl;
                    m_PyCubGeomFile << "cubit.move( cell, vector)" << std::endl;
                    m_PyCubGeomFile << "cyls.append(cyl)" << std::endl;

                    cyls[m-1] = cyl;

                    //copy cell nRadii  times for intersection with cylinders
#if defined (HAVE_ACIS) || defined (HAVE_OCC)
                    iGeom_copyEnt(igeomImpl->instance(), cells[n-1], &cell_copy, &err);
#endif
                    m_PyCubGeomFile << "cell_copy = cubit.copy_body(cells[" << n-1 << "])" << std::endl;
                    m_PyCubGeomFile << "cells_copy.append(cell_copy)" << std::endl;

                    cell_copys[m-1] = cell_copy;
#if defined (HAVE_ACIS) || defined (HAVE_OCC)
                    iGeom_intersectEnts(igeomImpl->instance(), cell_copy, cyls[m-1],&intersec,&err);
#endif
                    m_PyCubGeomFile << "tmpunite = cubit.unite(cells_copy[" << m-1 << "], cyls [" << m-1 << "])" << std::endl;
                    m_PyCubGeomFile << "intersec = cubit.subtract(tmpunite, cells_copy[" << m-1 << "])" << std::endl;
                    m_PyCubGeomFile << "intersec_main.append(intersec)" << std::endl;

                    intersec_main[m-1] = intersec;
                    intersec = NULL;
                  }

                //set tag on inner most cylinder, search for the full name of the abbreviated Cell Mat
                tmp_vol1=intersec_main[0];
                m_PyCubGeomFile << "tmp_vol1 = intersec_main[0]" << std::endl;

                for(int p=1;p<=m_szAssmMatAlias.GetSize();p++){
                    if(strcmp (szVCylMat(1).c_str(), m_szAssmMatAlias(p).c_str()) == 0){
                        sMatName = m_szAssmMat(p);
                      }
                  }

                m_PyCubGeomFile << "tmp_vol1 = cyls[0] \ncp_in.append(tmp_vol1)" << std::endl;
                cp_in.push_back(tmp_vol1);
#if defined (HAVE_ACIS) || defined (HAVE_OCC)
                iGeom_setData(igeomImpl->instance(), tmp_vol1, this_tag,
                              sMatName.c_str(), 10, &err);
#endif
                m_PyCubGeomFile  << "lid =tmp_vol1.id()" << std::endl;
                m_PyCubGeomFile << "cubit.set_entity_name(\"body\", lid, \""  << sMatName <<  "\" )" << std::endl;

                if(strcmp(m_szInfo.c_str(),"on") == 0){
#if defined (HAVE_ACIS) || defined (HAVE_OCC)
                    iGeom_setData(igeomImpl->instance(), tmp_vol1, this_tag,
                                  pin_name.c_str(), pin_name.size(), &err);
#endif
                    std::cout << "Naming pin body :" <<  pin_name<< std::endl;
                  }

                Name_Faces(sMatName, tmp_vol1, this_tag);
                m_PyCubGeomFile << "name_faces(\"" << sMatName << "\", tmp_vol1) " << std::endl;

                // delete the cell as this is the case when no. cell material is specified
#if defined (HAVE_ACIS) || defined (HAVE_OCC)
                iGeom_deleteEnt(igeomImpl->instance(), cells[n-1], &err);
#endif
                m_PyCubGeomFile << "cubit.cmd('delete vol << cells[" << n-1 << "].id()')" << std::endl;


                // other cyl annulus after substraction
                for (int b=nRadii; b>1; b--){
#if defined (HAVE_ACIS) || defined (HAVE_OCC)
                    iGeom_copyEnt(igeomImpl->instance(), intersec_main[b-2], &tmp_intersec, &err);
#endif
                    m_PyCubGeomFile << "tmp_intersec = cubit.copy_body(intersec_main[" << b-2 << "])" << std::endl;

                    //subtract tmp vol from the outer most
#if defined (HAVE_ACIS) || defined (HAVE_OCC)
                    iGeom_subtractEnts(igeomImpl->instance(), intersec_main[b-1], tmp_intersec, &tmp_new, &err);
#endif
                    m_PyCubGeomFile << "sub1.append(tmp_intersec\nsub2.append(intersec_main["<< b-1 << "])" << std::endl;
                    m_PyCubGeomFile << "tmp_new = cubit.subtract(sub2, sub1)" << std::endl;
                    m_PyCubGeomFile << "sub1[:] = []\nsub2[:] = []" << std::endl; 

                    // now search for the full name of the abbreviated Cell Mat
                    //    int tag_no;
                    for(int p=1;p<=m_szAssmMatAlias.GetSize();p++){
                        if(strcmp (szVCylMat(b).c_str(), m_szAssmMatAlias(p).c_str()) == 0){
                            //        tag_no = p;
                            sMatName =  m_szAssmMat(p);
                          }
                      }
                    std::cout << "created: " << sMatName << std::endl;

                    m_PyCubGeomFile << "cp_in.append(tmp_new)" << std::endl;
                    cp_in.push_back(tmp_new);

                    // set the name of the annulus
#if defined (HAVE_ACIS) || defined (HAVE_OCC)
                    iGeom_setData(igeomImpl->instance(), tmp_new, this_tag,
                                  sMatName.c_str(),sMatName.size(), &err);
#endif

                    m_PyCubGeomFile  << "lid =tmp_new.id()" << std::endl;
                    m_PyCubGeomFile << "cubit.set_entity_name(\"body\", lid, \""  << sMatName <<  "\" )" << std::endl;
                    
                    if(strcmp(m_szInfo.c_str(),"on") == 0){
#if defined (HAVE_ACIS) || defined (HAVE_OCC)
                        iGeom_setData(igeomImpl->instance(), tmp_new, this_tag,
                                      pin_name.c_str(), pin_name.size(), &err);
#endif
                        std::cout << "Naming pin body :" <<  pin_name<< std::endl;
                        m_PyCubGeomFile << "cubit.set_entity_name(\"body\", lid, \""  << pin_name <<  "\" )" << std::endl;
                      }

                    Name_Faces(sMatName, tmp_new, this_tag);
                    m_PyCubGeomFile << "name_faces(\"" << sMatName << "\", tmp_new) " << std::endl;
                    m_PyCubGeomFile << "cyls[" << b-1 << "] = tmp_new" << std::endl;

                    // copy the new into the cyl array
                    cyls[b-1] = tmp_new;

                  }
              }
            if(nDuctIndex > 0){
                for (int count = 0; count < (int) cp_in.size(); count++){
                  cp_inpins[nDuctIndex-1].push_back(cp_in[count]);
                  m_PyCubGeomFile << "cp_inpins["<< nDuctIndex -1 << "].append(cp_in[" << count << "])" << std::endl;
                }

              }
            cp_in.clear();
            m_PyCubGeomFile << "cp_in[:] =[]" << std::endl;

          }
      }

  }

}
