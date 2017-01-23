
#include "meshkit/AssyGen.hpp"

namespace MeshKit
{
void AssyGen::CreateCubitJournal()
//---------------------------------------------------------------------------
//Function: Create Cubit Journal File for generating mesh
//Input:    none
//Output:   none
//---------------------------------------------------------------------------
{
  if(m_szMeshScheme == "hole")
    m_FileOutput << "surf in group hole_surfaces scheme hole" << std::endl;

  if (m_nBLAssemblyMat !=0){
      // Also look for material name in BL material list
      for (int ll=1; ll<= m_nBLAssemblyMat; ll++){
          //if(szVCylMat(m) == m_szBLAssmMat(ll)) {
#ifdef HAVE_RGG16
            m_FileOutput << "group 'tmpgrp' equals surf with name '" <<  m_szBLAssmMat(ll)  << "*_top*'" << std::endl;
#else
            m_FileOutput << "group 'tmpgrp' equals surf with name '" <<  m_szBLAssmMat(ll)  << "_top'" << std::endl;
#endif

          m_FileOutput << "surf in tmpgrp size {RADIAL_MESH_SIZE}" << std::endl;
          m_FileOutput << "group '" << m_szBLAssmMat(ll) << "_hole_surfaces' equals surf in tmpgrp"<< std::endl;
          m_FileOutput << "surface in group " << m_szBLAssmMat(ll) << "_hole_surfaces scheme hole rad_interval " << m_nBLMatIntervals(ll) << " bias " << m_dBLMatBias(ll) << std::endl;
          if(strcmp(m_szSmooth.c_str(),"on") == 0)
            m_FileOutput << "surf in group " << m_szBLAssmMat(ll) << "_hole_surfaces" << " smooth scheme condition number beta 2.0 cpu 10" << std::endl;
          //         m_FileOutput << "mesh surf in group " << m_szBLAssmMat(ll) << "_hole_surfaces" << std::endl;
          // }
          m_FileOutput << "group 'bl_surfaces' add surf in tmpgrp" << std::endl;
        }
    }
  // variables
  int nColor;
  std::string color[21] = {" ", "thistle", "grey", "deepskyblue", "red", "purple",  "green",
                           "yellow", "royalblue", "magenta", "cyan", "lightsalmon", "springgreen",
                           "gold", "orange", "brown", "pink", "khaki", "black", "aquamarine", "mediumslateblue"};

  // if creating only journal file load the geometry file to compute bounding box for automatic size specification
  if(m_nJouFlag == 1){
#if defined (HAVE_ACIS) || defined (HAVE_OCC)
      iGeom_load(igeomImpl->instance(), m_szGeomFile.c_str(), NULL, &err, m_szGeomFile.length() , 0);
#endif
    }

  // get the max and min coordinates of the geometry
#if defined (HAVE_ACIS) || defined (HAVE_OCC)
  double x1, y1, z1, x2, y2, z2;
  iGeom_getBoundBox( igeomImpl->instance(), &x1, &y1, &z1, &x2, &y2, &z2, &err );
#endif

  int nSideset=m_nNeumannSetId;

  m_SchemesFile << "group \"gall\" add vol all\n#{gtempid = Id(\"group\")}\n" << std::endl;
  m_SchemesFile << "#{Zmax = BBox_ZMax(\"group\", gtempid)}" << std::endl;
  m_SchemesFile << "#{Zmin = BBox_ZMin(\"group\", gtempid)}" << std::endl;

  std::string szGrp, szBlock, szSurfTop, szSurfBot, szSize, szSurfSide;
  double dHeight = 0.0, dMid = 0.0;
  int nTemp = 1;
#if defined (HAVE_ACIS) || defined (HAVE_OCC)
  if(m_nDimensions > 0){
      dHeight= fabs(z2 - z1);
      dMid = z2 - dHeight/2.0;
    }
#endif

  // writing to template.jou
  m_SchemesFile << "## This file is created by rgg program in MeshKit ##\n";
  m_SchemesFile << "##Schemes " << std::endl  ;
  m_SchemesFile << "#{CIRCLE =\"circle interval 1 fraction 0.8\"}" << std::endl;
  m_SchemesFile << "#{HOLE = \"hole rad_interval 2 bias 0.0\"}" << std::endl;
  m_SchemesFile << "#{PAVE = \"pave\"}" << std::endl;
  m_SchemesFile << "#{MAP = \"map\"}" << std::endl;
  m_SchemesFile << "#{SWEEP = \"sweep\"}" << std::endl;
  m_SchemesFile << "#{TET = \"tetmesh\"}" << std::endl;
  m_SchemesFile << "#{TOP_EDGE_INTERVAL = " << m_edgeInterval << " }" << std::endl;
  m_SchemesFile << "## Dimensions" << std::endl;
  if(m_szGeomType == "hexagonal"){
      if(m_nDimensions > 0){
          m_SchemesFile << "#{PITCH =" << m_dMAssmPitch(nTemp, m_nDimensions) << "}" << std::endl;
        }
    }
  else if(m_szGeomType == "rectangular"){
      if(m_nDimensions > 0){
          m_SchemesFile << "#{PITCHX =" << m_dMAssmPitchX(nTemp, m_nDimensions)<< "}" << std::endl;
          m_SchemesFile << "#{PITCHY =" << m_dMAssmPitchY(nTemp, m_nDimensions) << "}" << std::endl;
        }
    }
  if( m_nPlanar ==0){
      m_SchemesFile << "#{Z_HEIGHT = Zmax - Zmin}" << std::endl;
      m_SchemesFile << "#{Z_MID = (Zmax + Zmin)/2.0}" << std::endl;
    }
  m_SchemesFile << "##Set Mesh Sizes" << std::endl;

  if (m_szMeshType == "hex"){
      // volume only
      if(m_nPlanar == 0 ){
          if (m_dAxialSize.GetSize() == 0){
              m_SchemesFile << "#{AXIAL_MESH_SIZE = 0.1*Z_HEIGHT}" << std::endl;
            }
          else {
              m_SchemesFile << "#{AXIAL_MESH_SIZE = " << m_dAxialSize(1) << "}" << std::endl;
            }

          // create templates for specifying block z intervals
          if (m_nDuct > 1){
              m_SchemesFile << "## Set interval along Z direction ## " << std::endl;

              for( int p=1; p<= m_nDuct; p++){
                  if (m_dAxialSize.GetSize() != 0)
                    m_SchemesFile << "#{AXIAL_MESH_SIZE" << p << "=" << m_dAxialSize(p) << "}" << std::endl;
                  else
                    m_SchemesFile << "#{AXIAL_MESH_SIZE" << p << "= 0.1*Z_HEIGHT}" << std::endl;
                  m_SchemesFile << "#{BLOCK" << p << "_Z_INTERVAL = AXIAL_MESH_SIZE" << p << "}" << std::endl;
                  m_SchemesFile << "#{BLOCK" << p << "_ZBOT = " << m_dMZAssm(p, 1) << "}" << std::endl;
                  m_SchemesFile << "#{BLOCK" << p << "_ZTOP = " << m_dMZAssm(p, 2) << "}" << std::endl;
                }
              m_SchemesFile << "##" << std::endl;
            }
        }
      if (-1.0 == m_dRadialSize) {
          if (m_szGeomType == "hexagonal")
            m_SchemesFile << "#{RADIAL_MESH_SIZE = 0.1*PITCH}" << std::endl;
          else
            m_SchemesFile << "#{RADIAL_MESH_SIZE = 0.02*0.5*(PITCHX+PITCHY)}" << std::endl;
        }
      else
        m_SchemesFile << "#{RADIAL_MESH_SIZE = " << m_dRadialSize << "}" << std::endl;
    }
  else if (m_szMeshType == "tet"){
      if (-1.0 == m_dTetMeshSize) {
          if (m_szGeomType == "hexagonal")
            m_SchemesFile << "#{TET_MESH_SIZE = 0.1*PITCH}" << std::endl;
          else
            m_SchemesFile << "#{TET_MESH_SIZE = 0.02*0.5*(PITCHX+PITCHY)}" << std::endl;
        }
      else {
          m_SchemesFile << "#{TET_MESH_SIZE = " << m_dTetMeshSize  << "}" << std::endl;
        }
    }

  if(m_nHblock == -1){ // if more blocks are needed axially, create'em using hexes and the end
      // block creation dumps
      m_FileOutput << "#Creating blocks, Note: you might need to combine some blocks" << std::endl;
      // group creation dumps. each material has a group
      m_FileOutput << "#Creating groups" << std::endl;
      for(int p=1;p<=(m_szAssmMatAlias.GetSize() - m_nBLAssemblyMat);p++){
          szGrp = "g_"+ m_szAssmMat(p);
          m_szAssmMat(p);
          if(m_nPlanar ==1){
              m_FileOutput << "group \"" << szGrp << "\" add surface name \"" << m_szAssmMat(p) <<"\"" << std::endl;
            }
          else{
              m_FileOutput << "group \"" << szGrp << "\" add body name \"" << m_szAssmMat(p) <<"\"" << std::endl;
            }
        }
      for(int p = 1; p <=  (m_szAssmMatAlias.GetSize() - m_nBLAssemblyMat);p++){
          szBlock = "b_"+ m_szAssmMat(p);
          szGrp = "g_"+ m_szAssmMat(p);
          m_FileOutput << "#{nb" << p << " =NumInGrp('" << szGrp << "')}" << std::endl;
          m_FileOutput << "#{Ifndef(nb" << p << ")}" << "\n" << "#{else}" << std::endl;
          if(m_nPlanar ==1){
              m_FileOutput << "block " << m_nMaterialSetId + p -1 << " surface in " << szGrp  << std::endl;
              m_FileOutput << "block " << m_nMaterialSetId + p -1 << " name \"" << szBlock <<"\""<< std::endl;
            }
          else{
              m_FileOutput << "block " << m_nMaterialSetId + p -1 << " body in " << szGrp  << std::endl;
              m_FileOutput << "block " << m_nMaterialSetId + p -1 << " name \"" << szBlock <<"\""<< std::endl;
            }
          m_FileOutput << "#{endif}" << std::endl;
        }
      m_FileOutput << "#" << std::endl;
    }
  if(m_szMeshType == "hex"){
      // imprint
      m_FileOutput << "#Imprint geometry" << std::endl;
      m_FileOutput << "imprint all" << std::endl;
      m_FileOutput << "#" << std::endl;
      // merge

      m_FileOutput << "Merge Tolerance " << m_dMergeTol << std::endl;
      m_FileOutput << "#" << std::endl;

      m_FileOutput << "#Merge geometry" << std::endl;
      m_FileOutput << "merge all" << std::endl;
      m_FileOutput << "#" << std::endl;
    }

  // for info keyword
  if(strcmp(m_szInfo.c_str(),"on") == 0){
      int temp = 9700;
      m_FileOutput << "# stuff for info keyword, remove if not desired " << std::endl;
      m_FileOutput << "# putting pins in seperate blocks " << std::endl;
      m_FileOutput << "#" << std::endl;
      for (int i=0; i<m_nTotalPincells; i++){
          m_FileOutput << "group 'g"<< i+m_nStartpinid << "' add body with name '_xp" << i+m_nStartpinid << "_'" << std::endl;

          m_FileOutput << "#{nbody" << i+1 << " =NumInGrp('g" <<i+m_nStartpinid << "')}" << std::endl;
          m_FileOutput << "#{Ifndef(nbody" << i+1 << ")}" << "\n" << "#{else}" << std::endl;
          m_FileOutput << "block " << temp+i << " body in group g" << i+m_nStartpinid << std::endl;
          m_FileOutput << "block " << temp+i << " name '_xp" << i+m_nStartpinid << "'" << std::endl;
          m_FileOutput << "#{endif}" << std::endl;
        }
    }

  //surface only
  if(m_nPlanar ==1){
      m_FileOutput << "# Pointing surface normals to 0.0, 0.0, -1.0 or -ve Z or correct STARCCM+ cell-face orientation" << std::endl;
      m_FileOutput << "surface all normal opposite" << std::endl;
      m_FileOutput << "#" << std::endl;
    }
  // volume only
  else{
      if(m_szSideset == "yes"){

          // rename the skin surfaces, so that they don't appear as sidesets
          for (int p=1; p<=m_nDuct; p++){
              for(int q=1;q<=m_nSides; q++){
                  m_FileOutput << "group 'edge" << (m_nSides*(p-1) + q ) <<"' equals curve with name 'side_edge"
                               << (m_nSides*(p-1) + q ) << "@'" << std::endl;

                  m_FileOutput << "group 'vt" <<  (m_nSides*(p-1) + q )  <<"' equals vertex with z_max == z_min in curve in edge"
                               <<  (m_nSides*(p-1) + q ) << std::endl;

                }
            }

          // creating groups for vertices on the top surface of the duct
          for (int p=1; p<=m_nDuct; p++){
              for(int q=1;q<=m_nSides; q++){

                  if(q != m_nSides){
                      m_FileOutput << "group 'v" << (m_nSides*(p-1) + q ) <<"' intersect group vt" << (m_nSides*(p-1) + q )
                                   << " with group vt" <<  (m_nSides*(p-1) + q + 1 )  << std::endl;
                    }
                  else {
                      m_FileOutput << "group 'v" << (m_nSides*(p-1) + q ) <<"' intersect group vt" << (m_nSides*(p-1) + q )
                                   << " with group vt" <<  (m_nSides*(p-1) + 1 )  << std::endl;
                    }
                }
            }
          // creating temp surfaces groups
          for (int p=1; p<=m_nDuct; p++){
              for(int q=1;q<=m_nSides; q++){
                  m_FileOutput << "group 'st" << (m_nSides*(p-1) + q ) <<"' equals surface with z_max <> z_min in vert in v"
                               << (m_nSides*(p-1) + q ) << "'" << std::endl;
                }
            }

          // creating surface groups for obtaining surfaces
          for (int p=1; p<=m_nDuct; p++){
              for(int q=1;q<=m_nSides; q++){
                  if(q != 1){
                      m_FileOutput << "group 's" << (m_nSides*(p-1) + q ) <<"' intersect group st"  << (m_nSides*(p-1) + q )
                                   << " with group st" <<  (m_nSides*(p-1) + q - 1 )  << std::endl;
                    }
                  else {
                      m_FileOutput << "group 's" << (m_nSides*(p-1) + q ) <<"' intersect group st" << (m_nSides*(p-1) + q )
                                   << " with group st" <<  (m_nSides*(p-1) + m_nSides )  << std::endl;
                    }
                }
            }

          // renaming the skin side surfaces
          for (int p=1; p<=m_nDuct; p++){
              for(int q=1;q<=m_nSides; q++){
                  m_FileOutput << "surface in group s" <<  (m_nSides*(p-1) + q ) << " rename 'side_surface"
                               <<  (m_nSides*(p-1) + q ) << "'" << std::endl;

                }
            }
        }

      if(m_szMeshType == "hex"){

          //now set the sizes
          m_FileOutput << "#Set Meshing Scheme and Sizes, use template.jou to specify sizes" << std::endl;

          for(int p=1;p<=(m_szAssmMatAlias.GetSize() - m_nBLAssemblyMat);p++){
              szGrp = "g_"+ m_szAssmMat(p);
              szSize =  m_szAssmMat(p) + "_size";
              szSurfBot = m_szAssmMat(p) + "_bot";
              szSize =  m_szAssmMat(p) + "_surf_size";
              m_FileOutput << "group 'tmpgrp' equals surface name \""  << szSurfBot  << "\"" << std::endl;
              m_FileOutput << "surface in tmpgrp  size {"  << szSize <<"}" << std::endl;
            }
          m_FileOutput << "#" << std::endl;
        }
    }
  if(m_szMeshType == "hex"){
      // some more common stuff meshing top surfaces set the sizes and mesh
      m_FileOutput << "#Surfaces mesh, use template.jou to specify sizes" << std::endl;
      for(int p=1;p<=(m_szAssmMatAlias.GetSize() - m_nBLAssemblyMat);p++){
          szSurfTop = m_szAssmMat(p) + "_top";
          szGrp = "g_"+ m_szAssmMat(p);
          szSize =  m_szAssmMat(p) + "_surf_size";
          if(m_szMeshScheme == "hole" && m_nBLAssemblyMat == 0){
              m_FileOutput << "group 'tmpgrp' equals surface name \""  << szSurfTop  << "\"" << std::endl;
              m_FileOutput << "group 'remove_hole' intersect group tmpgrp with group hole_surfaces" << std::endl;
              m_FileOutput << "#{nIntersect=NumInGrp('remove_hole')}" << std::endl;
              m_FileOutput << "#{If(nIntersect==0)}" << std::endl;
              m_FileOutput << "surface in tmpgrp  size {"  << szSize <<"}" << std::endl;
              m_FileOutput << "surface in tmpgrp scheme {" << "PAVE" << "}"  << std::endl;
              m_FileOutput << "#{endif}"  << std::endl;
            }
          else{
              m_FileOutput << "group 'tmpgrp' equals surface name \""  << szSurfTop  << "\"" << std::endl;
              m_FileOutput << "surface in tmpgrp  size {"  << szSize <<"}" << std::endl;
              m_FileOutput << "surface in tmpgrp scheme {" << "PAVE" << "}"  << std::endl;
            }

          if (p==1 && m_edgeInterval != 99){
              m_FileOutput << "group 'edges" <<"' equals curve with name 'side_edge'"<< std::endl;
              m_FileOutput << "curve in edges interval {TOP_EDGE_INTERVAL}" << std::endl;
            }

          //    m_FileOutput << "mesh surface in " << szGrp << "\n#" << std::endl;

          // dumping these sizes schemes.jou also
          m_SchemesFile << "#{"  << szSize <<" = RADIAL_MESH_SIZE}" << std::endl;
        }
      m_FileOutput << "#" << std::endl;

      // mesh all command after meshing surface
      if (m_nDuct <= 1 ){
          m_FileOutput << "group 'tmpgrp' add surface name '_top'" << std::endl;
          if (m_nBLAssemblyMat !=0){ // only if boundary layers are specified
              m_FileOutput << "group 'tmpgrp1' subtract innerduct from tmpgrp" << std::endl;
              m_FileOutput << "group 'tmpgrp2' subtract bl_surfaces from tmpgrp1" << std::endl;
              m_FileOutput << "mesh tmpgrp2" << std::endl;
            }
          else
            {
              m_FileOutput << "mesh tmpgrp" << std::endl;
            }
        }
      else {
          m_FileOutput << "#Meshing top surface" << std::endl;
          //m_FileOutput << "mesh surface with z_coord = " << z2 << std::endl;
          if(m_szMeshScheme == "hole" && m_nBLAssemblyMat == 0){
              m_FileOutput << "group 'tmpgrp' equals surface name \""  << szSurfTop  << "\"" << std::endl;
              m_FileOutput << "group 'remove_hole' intersect group tmpgrp with group hole_surfaces" << std::endl;
              m_FileOutput << "#{nIntersect=NumInGrp('remove_hole')}" << std::endl;
              m_FileOutput << "#{If(nIntersect==0)}" << std::endl;
              m_FileOutput << "surface in tmpgrp  size {"  << szSize <<"}" << std::endl;
              m_FileOutput << "surface in tmpgrp scheme {" << "PAVE" << "}"  << std::endl;
              m_FileOutput << "#{endif}"  << std::endl;
              m_FileOutput << "mesh surface with z_coord = {Zmax}" << std::endl;
            }
          else if (m_nBLAssemblyMat != 0){ // mesh by spefifying boundary layers or mesh partially
              m_FileOutput << "group 'tmpgrp' equals surface name '_top'" << std::endl;
              m_FileOutput << "group 'tmpgrp1' subtract innerduct from tmpgrp" << std::endl;
              m_FileOutput << "group 'tmpgrp2' subtract bl_surfaces from tmpgrp1" << std::endl;
              m_FileOutput << "group 'tmpgrp3' equals surface in tmpgrp2 with z_coord = {Zmax}" << std::endl;

              m_FileOutput << "surface in tmpgrp3  size {"  << szSize <<"}" << std::endl;
              m_FileOutput << "surface in tmpgrp3 scheme {" << "PAVE" << "}"  << std::endl;
              m_FileOutput << "mesh tmpgrp3" << std::endl;

            }
          else {
              m_FileOutput << "group 'tmpgrp' equals surface name \""  << szSurfTop  << "\"" << std::endl;
              m_FileOutput << "surface in tmpgrp  size {"  << szSize <<"}" << std::endl;
              m_FileOutput << "surface in tmpgrp scheme {" << "PAVE" << "}"  << std::endl;
              m_FileOutput << "mesh surface with z_coord = {Zmax}" << std::endl;
            }
        }
      // This part is for mesh top surfaces only when boundary layer surfaces are specified
      if (m_nBLAssemblyMat !=0){
          // Also look for material name in BL material list
          for (int ll=1; ll<= m_nBLAssemblyMat; ll++){
              bool duct = false;
              for (int n = 0; n < (int) m_szDuctMats.size(); n++){
                  if (strcmp(m_szDuctMats[n].c_str(), m_szBLAssmMat(ll).c_str()) == 0)
                    duct = true;
                  else
                    duct = false;
                }
              if (duct){                 //We want to use this part with pair node only for ducts and not cylinderical pins so check if this material is duct or not
                  if (m_edgeInterval != 99)
                    m_FileOutput << "curve in surf in " << m_szBLAssmMat(ll) << "_hole_surfaces interval {TOP_EDGE_INTERVAL}"<< std::endl;
                  m_FileOutput << "mesh vertex in surf in " << m_szBLAssmMat(ll) << "_hole_surfaces with z_coord = {Z_max}" << std::endl;
                  m_FileOutput << "#{corner1 = Id('node')} " << std::endl;
                  m_FileOutput << "group 'gcurves' equals curve in surface in " << m_szBLAssmMat(ll) << "_hole_surfaces'" << std::endl;
                  m_FileOutput << "#{_cntr=0} " << "\n" <<
                                  "#{_tmp_dis=0} " << "\n" <<
                                  "#{_min_dis=0}  " << "\n" <<
                                  "#{_closest_node=11} " << "\n" <<

                                  "group 'v_node' equals node in volume in surface in " <<  m_szBLAssmMat(ll) << "_hole_surfaces" << "\n" <<
                                  "group v_node remove node {corner1} " << "\n" <<
                                  "#{xc1 = Nx(corner1)} " << "\n" <<
                                  "#{yc1 = Ny(corner1)} " << "\n" <<
                                  "#{_num_nodes = NumInGrp('v_node')} " << "\n" <<
                                  "#{_min_dis = 1.e10} " << "\n" <<
                                  "#{Loop(20)} " << "\n" <<
                                  "#{_node_id = GroupMemberId('v_node', 'node', _cntr)} " << "\n" <<
                                  "#{_xni = Nx(_node_id)} " << "\n" <<
                                  "#{_yni = Ny(_node_id)} " << "\n" <<
                                  "#{_tmp_dis = (xc1 - _xni)*(xc1 -_xni) + (yc1 -_yni)*(yc1 - _yni)} " << "\n" <<
                                  "#{if(_tmp_dis < _min_dis)} " << "\n" <<
                                  "#{ _closest_node = _node_id} " << "\n" <<
                                  "# {_min_dis=_tmp_dis} " << "\n" <<
                                  "#{endif} " << "\n" <<
                                  "#{_cntr++} " << "\n" <<
                                  "#{if (_cntr >_num_nodes)} " << "\n" <<
                                  "#{break} " << "\n" <<
                                  "#{endif} " << "\n" <<
                                  "#{EndLoop} " << "\n" << std::endl;

                  // This must be used only for ducts
                  //   if (m_szBLAssmMat(ll) == duct material or it's not pin material')
                  m_FileOutput << "surf in group " << m_szBLAssmMat(ll) << "_hole_surfaces scheme hole rad_intervals "
                               << m_nBLMatIntervals(ll) << " bias " << m_dBLMatBias(ll) << " pair node {corner1} with node {_closest_node}" << std::endl;
                }
              else { // this is regular cylinder
                  m_FileOutput << "surf in group " << m_szBLAssmMat(ll) << "_hole_surfaces scheme hole rad_intervals "
                               << m_nBLMatIntervals(ll) << " bias " << m_dBLMatBias(ll) << std::endl;
                }

              m_FileOutput << "mesh surf in group " << m_szBLAssmMat(ll) << "_hole_surfaces with z_coord = {Z_max}" << std::endl;
              if(strcmp(m_szSmooth.c_str(),"on") == 0)
                m_FileOutput << "smooth surf in group " << m_szBLAssmMat(ll) << "_hole_surfaces" << std::endl;
            }
          m_FileOutput << "mesh surf in innerduct with z_coord = {Z_max}" << std::endl;
        }

      if(m_nPlanar == 0){ // volumes only
          if (m_nDuct == 1){
              m_FileOutput << "surf with z_coord > {Z_MID -.1*Z_HEIGHT}" <<
                              " and z_coord < {Z_MID + .1*Z_HEIGHT} size {AXIAL_MESH_SIZE}" << std::endl ;
              m_FileOutput << "mesh vol all" << std::endl;
            }
          else if (m_nDuct > 1){
              m_FileOutput << "### Setting Z intervals on ducts and meshing along Z " << std::endl;
              for( int p=m_nDuct; p>= 1; p--){
                  if(dMid == 0){ // z - centered
                      m_FileOutput << "surf with z_coord  > " << m_dMZAssm(p, 1) - dHeight/2.0
                                   << " and z_coord < " << m_dMZAssm(p, 2) - dHeight/2.0 << " interval " << "{BLOCK" << p << "_Z_INTERVAL}" << std::endl;
                      m_FileOutput << "mesh vol with z_coord  > " << m_dMZAssm(p, 1) - dHeight/2.0
                                   << " and z_coord < " << m_dMZAssm(p, 2) - dHeight/2.0 << std::endl;
                    }
                  else{
                      m_FileOutput << "surf with z_coord  > " << m_dMZAssm(p, 1)
                                   << " and z_coord < " << m_dMZAssm(p, 2) << " interval " << "{BLOCK" << p << "_Z_INTERVAL}" << std::endl;
                      m_FileOutput << "mesh vol with z_coord  > " << m_dMZAssm(p, 1)
                                   << " and z_coord < " << m_dMZAssm(p, 2) << std::endl;

                      m_FileOutput << "##" << std::endl;
                    }
                }
            }
        }
    }

  else if(m_szMeshType == "tet"){
      m_FileOutput << "##"<< std::endl;
      m_FileOutput << "# groupings for creating vertex groups"<< std::endl;
      for (int p=1; p<=m_nDuct; p++){
          for(int q=1;q<=m_nSides; q++){
              m_FileOutput << "group 'edge" << (m_nSides*(p-1) + q ) <<"' equals curve with name 'side_edge"
                           << (m_nSides*(p-1) + q ) << "@'" << std::endl;

              m_FileOutput << "group 'vt" <<  (m_nSides*(p-1) + q )  <<"' equals vertex with z_max == z_min in curve in edge"
                           <<  (m_nSides*(p-1) + q ) << std::endl;

            }
        }

      // creating groups for vertices on the top surface of the duct
      for (int p=1; p<=m_nDuct; p++){
          for(int q=1;q<=m_nSides; q++){

              if(q != m_nSides){
                  m_FileOutput << "group 'v" << (m_nSides*(p-1) + q ) <<"' intersect group vt" << (m_nSides*(p-1) + q )
                               << " with group vt" <<  (m_nSides*(p-1) + q + 1 )  << std::endl;
                }
              else {
                  m_FileOutput << "group 'v" << (m_nSides*(p-1) + q ) <<"' intersect group vt" << (m_nSides*(p-1) + q )
                               << " with group vt" <<  (m_nSides*(p-1) + 1 )  << std::endl;
                }
            }
        }

      // creating temp surfaces groups
      for (int p=1; p<=m_nDuct; p++){
          for(int q=1;q<=m_nSides; q++){
              m_FileOutput << "group 'st" << (m_nSides*(p-1) + q ) <<"' equals surface with z_max <> z_min in vert in v"
                           << (m_nSides*(p-1) + q ) << "'" << std::endl;
            }
        }

      // creating side curve and surface groups
      for (int p=1; p<=m_nDuct; p++){
          for(int q=1;q<=m_nSides; q++){

              m_FileOutput << "group 'c" <<  (m_nSides*(p-1) + q )  <<"' equals curve with z_max <> z_min in vert in v"
                           <<  (m_nSides*(p-1) + q ) << std::endl;

              if(q != 1){
                  m_FileOutput << "group 's" << (m_nSides*(p-1) + q ) <<"' intersect group st"  << (m_nSides*(p-1) + q )
                               << " with group st" <<  (m_nSides*(p-1) + q - 1 )  << std::endl;
                }
              else {
                  m_FileOutput << "group 's" << (m_nSides*(p-1) + q ) <<"' intersect group st" << (m_nSides*(p-1) + q )
                               << " with group st" <<  (m_nSides*(p-1) + m_nSides )  << std::endl;
                }
            }
        }

      // renaming the side surfaces for getting the split surfaces later
      for (int p=1; p<=m_nDuct; p++){
          for(int q=1;q<=m_nSides; q++){
              m_FileOutput << "surface in group s" <<  (m_nSides*(p-1) + q ) << " rename 'side_surface"
                           <<  (m_nSides*(p-1) + q ) << "'" << std::endl;
            }
        }

      // splitting the surfaces
      for (int p=1; p<=m_nDuct; p++){
          for(int q=1;q<=m_nSides; q++){

              m_FileOutput << "split surface in group s" <<  (m_nSides*(p-1) + q )  <<" direction curve in group c"
                           <<  (m_nSides*(p-1) + q ) << std::endl;
            }
        }

      // get all the split surfaces in individual groups
      for (int p=1; p<=m_nDuct; p++){
          for(int q=1;q<=m_nSides; q++){

              m_FileOutput << "group 'sname" << (m_nSides*(p-1) + q ) <<  "' equals surface with name 'side_surface"
                           <<  (m_nSides*(p-1) + q ) << "'"<< std::endl;
              m_FileOutput << "group 'svert" << (m_nSides*(p-1) + q ) <<  "' equals surface in vert in v"
                           <<  (m_nSides*(p-1) + q ) << std::endl;
              m_FileOutput << "group 'ssplit" << (m_nSides*(p-1) + q ) <<  "' intersect group sname" <<  (m_nSides*(p-1) + q )
                           << " with group svert" << (m_nSides*(p-1) + q ) << std::endl;
            }
        }

      // get all the split surfaces in individual groups
      for (int p=1; p<=m_nDuct; p++){
          for(int q=1;q<=m_nSides; q++){

              if(q != 1){
                  m_FileOutput << "group 'ssplit_" << (m_nSides*(p-1) + q ) <<"' intersect group sname"  << (m_nSides*(p-1) + q )
                               << " with group svert" <<  (m_nSides*(p-1) + q - 1 )  << std::endl;
                }
              else {
                  m_FileOutput << "group 'ssplit_" <<  (m_nSides*(p-1) + q ) <<"' intersect group sname" << (m_nSides*(p-1) + q )
                               << " with group svert" <<  (m_nSides*(p-1) + m_nSides )  << std::endl;
                }
            }
        }
      // imprint
      m_FileOutput << "#Imprint geometry" << std::endl;
      m_FileOutput << "imprint all" << std::endl;
      m_FileOutput << "#" << std::endl;
      m_FileOutput << "Merge Tolerance " << m_dMergeTol << std::endl;
      m_FileOutput << "#" << std::endl;

      // merge
      m_FileOutput << "#Merge geometry" << std::endl;
      m_FileOutput << "merge all" << std::endl;
      m_FileOutput << "#" << std::endl;

      m_FileOutput << "#Set mesh scheme and size" << std::endl;
      m_FileOutput << "volume all scheme {TET} size {TET_MESH_SIZE}" << std::endl;

      // mesh one side of each duct, such that one is flipped mesh of the other
      for (int p=1; p<=m_nDuct; p++){

          m_FileOutput << "mesh surface in group ssplit" <<  (m_nSides*(p-1) + 1) << std::endl;

          m_FileOutput << "surface in group ssplit_" << (m_nSides*(p-1) + 1) << " scheme copy source surface in group ssplit"
                       <<  (m_nSides*(p-1) + 1)
                        << " source curve in group c" <<  (m_nSides*(p-1) + 1 ) << " target curve in group c" <<  (m_nSides*(p-1) + m_nSides )
                        << " source vertex in group v" <<  (m_nSides*(p-1) + 1) << " target vertex in group v"  <<  (m_nSides*(p-1) + m_nSides )
                        << " nosmoothing" << std::endl;

          m_FileOutput << "mesh surface in group  ssplit_" << (m_nSides*(p-1) + 1) << std::endl;
        }

      // setting the copy mesh commands on the above pair of split surfaces to have all surfaces symmetrical
      for (int p=1; p<=m_nDuct; p++){
          for(int q=1;q<=m_nSides; q++){
              if(q != m_nSides){
                  m_FileOutput << "copy mesh surface in ssplit" <<  (m_nSides*(p-1) + 1)
                               << " onto surface in ssplit" << (m_nSides*(p-1) + q + 1 )
                               << " source vertex in group v" << (m_nSides*(p-1) + 1)
                               << " target vertex in group v" <<  (m_nSides*(p-1) + q + 1) << " nosmoothing" <<  std::endl;

                  m_FileOutput << "copy mesh surface in ssplit_" << (m_nSides*(p-1) + 1 )
                               << " onto surface in ssplit_" << (m_nSides*(p-1) + q +1 )
                               << " source vertex in group v" << (m_nSides*p)
                               << " target vertex in group v" <<  (m_nSides*(p-1) + q) << " nosmoothing" << std::endl;

                }
              else{
                  // do nothing
                }
            }
        }

      m_FileOutput << "# Mesh all volumes now" << std::endl;
      m_FileOutput << "mesh vol all" << std::endl;
    }

  // create and sidesets after meshing
  m_FileOutput << "#" << std::endl;
  //    }
  if(m_szSideset == "yes"){
      // top surface sidesets
      m_FileOutput << "#Creating top surface sidesets" << std::endl;
      m_FileOutput << "create group 'surfall'" << std::endl;
      for(int p=1;p<=(m_szAssmMatAlias.GetSize() - m_nBLAssemblyMat);p++){
          ++nSideset;
          szSurfTop = m_szAssmMat(p)+"_top";
          // Avoid creation if empty sideset
          m_FileOutput << "group 'tmpgrp' equals surface with name '" << szSurfTop << "' in vol in block " << m_nMaterialSetId + p -1 << std::endl;
          m_FileOutput << "sideset " << nSideset << " surface in tmpgrp " << std::endl;
          m_FileOutput << "sideset " << nSideset << " name \"" << szSurfTop << "_ss\"" << std::endl;
        }
      m_FileOutput << "#" << std::endl;
      for(int p=1;p<=m_nBLAssemblyMat;p++){
          ++nSideset;
          szSurfTop = m_szBLAssmMat(p)+"_top";
          // Avoid creation if empty sideset
          m_FileOutput << "group 'tmpgrp' equals surface with name '" << szSurfTop << std::endl;
          m_FileOutput << "sideset " << nSideset << " surface in tmpgrp " << std::endl;
          m_FileOutput << "sideset " << nSideset << " name \"" << szSurfTop << "_ss\"" << std::endl;
        }
      m_FileOutput << "#" << std::endl;
    }


  if(m_nPlanar ==0){
      if(m_szSideset == "yes"){
          // now create bot and side sideset
          m_FileOutput << "#Creating bot/side surface sidesets" << std::endl;
          for(int p=1;p<=(m_szAssmMatAlias.GetSize() - m_nBLAssemblyMat);p++){
              szSurfTop = m_szAssmMat(p)+"_bot";
              m_FileOutput << "#" << std::endl;
              ++nSideset;
              m_FileOutput << "group 'tmpgrp' equals surface with name '" << szSurfTop << "' in vol in block " << m_nMaterialSetId + p -1 << std::endl;
              m_FileOutput << "sideset " << nSideset << " surface in tmpgrp " << std::endl;
              m_FileOutput << "sideset " << nSideset << " name \"" << szSurfTop << "_ss\"" << std::endl;
            }
          for(int p=1;p<=m_nBLAssemblyMat;p++){
              ++nSideset;
              szSurfTop = m_szBLAssmMat(p)+"_bot";
              // Avoid creation if empty sideset
              m_FileOutput << "group 'tmpgrp' equals surface with name '" << szSurfTop << std::endl;
              m_FileOutput << "sideset " << nSideset << " surface in tmpgrp " << std::endl;
              m_FileOutput << "sideset " << nSideset << " name \"" << szSurfTop << "_ss\"" << std::endl;
            }
          m_FileOutput << "#" << std::endl;

          for(int p=1;p<=(m_szAssmMatAlias.GetSize() - m_nBLAssemblyMat);p++){
              szSurfSide = m_szAssmMat(p)+"_side";
              ++nSideset;
              if(m_szGeomType == "hexagonal"){
                  for (int u=1; u<=6;u++){
                      m_FileOutput << "group 'tmpgrp" << u <<"' equals surf with name '" << szSurfSide << u << "'" << std::endl;
                    }
                  m_FileOutput << "sideset " << nSideset << " surface in tmpgrp1 tmpgrp2 tmpgrp3 tmpgrp4 tmpgrp5 tmpgrp6" << std::endl;
                  m_FileOutput << "sideset " << nSideset << " name \"" << szSurfSide << "1_ss\"" << std::endl;
                  ++nSideset;
                }
              if(m_szGeomType == "hexagonal"){
                  for (int u=7; u<=12;u++){
                      m_FileOutput << "group 'tmpgrp" << u <<"' equals surf with name '" << szSurfSide << "_" << u << "'" << std::endl;
                    }
                  m_FileOutput << "sideset " << nSideset << " surface in tmpgrp7 tmpgrp8 tmpgrp9 tmpgrp10 tmpgrp11 tmpgrp12" << std::endl;
                  m_FileOutput << "sideset " << nSideset << " name \"" << szSurfSide << "2_ss\"" << std::endl;
                }
              if(m_szGeomType == "rectangular"){
                  for (int u=1; u<=4;u++){
                      m_FileOutput << "group 'tmpgrp" << u <<"' equals surf with name '" << szSurfSide << u << "'" << std::endl;
                    }
                  m_FileOutput << "sideset " << nSideset << " surface in tmpgrp1 tmpgrp2 tmpgrp3 tmpgrp4'" << std::endl;
                  m_FileOutput << "sideset " << nSideset << " name \"" << szSurfSide << "2_ss\"" << std::endl;
                  ++nSideset;
                }

              if(m_szGeomType == "rectangular"){
                  for (int u=5; u<=8;u++){
                      m_FileOutput << "group 'tmpgrp" << u <<"' equals surf with name '" << szSurfSide  << "_" << u << "'" << std::endl;
                    }
                  m_FileOutput << "sideset " << nSideset << " surface in tmpgrp5 tmpgrp6 tmpgrp7 tmpgrp8'" << std::endl;
                  m_FileOutput << "sideset " << nSideset << " name \"" << szSurfSide << "2_ss\"" << std::endl;
                }
            }


          m_FileOutput << "#" << std::endl;

          m_FileOutput << "#Creating sideset for outer most side surfaces" << std::endl;
          ++nSideset;

          m_FileOutput << "group 'tmpgrp' equals surf with name 'side_surface'" << std::endl;
          m_FileOutput << "sideset " << nSideset << " surface in tmpgrp " << std::endl;
          m_FileOutput << "sideset " << nSideset << " name \"" << "outer_side_ss\"" << std::endl;
        }
    }
  if(m_nHblock != -1){ // if more blocks are needed axially, create'em using hexes and the end
      // block creation dumps
      m_FileOutput << "#Creating blocks, Note: you might need to combine some blocks" << std::endl;
      // group creation dumps. each material has a group
      m_FileOutput << "#Creating groups" << std::endl;
      if(m_szMeshType == "hex"){
          for(int p=1;p<=(m_szAssmMatAlias.GetSize() - m_nBLAssemblyMat);p++){
              szGrp = "g_"+ m_szAssmMat(p);
              m_szAssmMat(p);
              if(m_nPlanar ==1){
                  m_FileOutput << "group \"" << szGrp << "\" add surface name \"" << m_szAssmMat(p) <<"\"" << std::endl;
                }
              else{

                  m_FileOutput << "group \"" << szGrp << "\" add body name \"" << m_szAssmMat(p) <<"\"" << std::endl;
                }
            }
          for(int p = 1; p <=  (m_szAssmMatAlias.GetSize() - m_nBLAssemblyMat);p++){
              szBlock = "b_"+ m_szAssmMat(p);
              szGrp = "g_"+ m_szAssmMat(p);
              m_FileOutput << "#{nb" << p << " =NumInGrp('" << szGrp << "')}" << std::endl;
              m_FileOutput << "#{Ifndef(nb" << p << ")}" << "\n" << "#{else}" << std::endl;
              if(m_nPlanar ==1){
                  m_FileOutput << "block " << m_nMaterialSetId + p -1 << " surface in " << szGrp  << std::endl;
                  m_FileOutput << "block " << m_nMaterialSetId + p -1 << " name \"" << szBlock <<"\""<< std::endl;
                }
              else{
                  m_FileOutput << "block " << m_nMaterialSetId + p -1 << " hex in body in " << szGrp  << std::endl;
                  m_FileOutput << "block " << m_nMaterialSetId + p -1 << " name \"" << szBlock <<"\""<< std::endl;
                }
              m_FileOutput << "#{endif}" << std::endl;
            }
          m_FileOutput << "#" << std::endl;
        }
      else{
          std::cout << "Error: Terminating journal file writing. \n Hex block (Hblock keyword) is not supported for a tet mesh." << std::endl;
          exit(1);
        }
    }

  // create super blocks
  if(m_nSuperBlocks > 0){
      for(int o = 1; o <= m_nSuperBlocks; o++){
          m_FileOutput << "block " << sb(o).m_nSuperBlockId << " vol in block ";
          for (int p = 1; p <= sb(o).m_nNumSBContents; p++){
              m_FileOutput << m_nMaterialSetId + sb(o).m_nSBContents(p) << " ";
            }
          m_FileOutput << "\n" << "block " << sb(o).m_nSuperBlockId << " name '" << sb(o).m_szSuperBlockAlias << "'" << std::endl;
          m_FileOutput << "delete block " ;
          for (int q = 1; q <= sb(o).m_nNumSBContents; q++){
              m_FileOutput << m_nMaterialSetId + sb(o).m_nSBContents(q) << " ";
            }
          m_FileOutput << "\n" << std::endl;
        }
    }



  if(m_nHblock > 0){
      // now dump the commands for making hex layers as blocks and subtracting from original
      double delta = (m_dZend - m_dZstart)/m_nHblock;
      for(int i=0; i<(m_szAssmMatAlias.GetSize() - m_nBLAssemblyMat); i++){
          m_FileOutput << "## BLOCK CREATION USING HEXES" << std::endl;
          for(int j=0; j<m_nHblock; j++){
              m_FileOutput << "group 'tmpgrp" << j+1 << "' equals hex in block " <<  m_nMaterialSetId + i
                           << " with z_coord < " << m_dZstart + (j+1)*delta << " and z_coord > "
                           << m_dZstart + j*delta << std::endl;
            }
          for(int j=0; j<m_nHblock; j++){
              m_FileOutput << "block " <<  m_nMaterialSetId+i << " group tmpgrp" << j+1 << " remove" << std::endl;
            }
          for(int j=0; j<m_nHblock; j++){
              if((m_szAssmMatAlias.GetSize() - m_nBLAssemblyMat) < 10)
                m_FileOutput << "block " << j+1 <<  m_nMaterialSetId+i << " group tmpgrp" << j+1 << std::endl;
              else
                m_FileOutput << "block " << (j+1)*10 <<  m_nMaterialSetId+i << " group tmpgrp" << j+1 << std::endl;
            }
        }
    }
  if(m_nMaterialSetId != 1)
    m_FileOutput << "renumber hex all start_id " << MAXLINES*1000 << std::endl;
  // color now
  m_FileOutput << "#Set color for different parts" << std::endl;
  if(m_nPlanar == 0){ // volumes only
      for(int p=1;p<=(m_szAssmMatAlias.GetSize() - m_nBLAssemblyMat);p++){
          szGrp = "g_"+ m_szAssmMat(p);
          if(p>20)
            nColor = 1;
          else
            nColor = p;
          m_FileOutput << "color body in " << szGrp << " " << color[nColor] << std::endl;
        }
    }
  else{ //surfaces
      // color now
      for(int p=1;p<=(m_szAssmMatAlias.GetSize() - m_nBLAssemblyMat);p++){
          szGrp = "g_"+ m_szAssmMat(p);
          if(p>20)
            nColor = 1;
          else
            nColor = p;
          m_FileOutput << "color surface in " << szGrp << " " << color[nColor] << std::endl;
        }
    }

  m_FileOutput << "delete group all" << std::endl;
  // save as .cub file dump
  m_FileOutput << "#\n#Save file" << std::endl;
  if(save_exodus){
      std::string szSave = m_szFile + ".exo";
      std::transform(szSave.begin(), szSave.end(), szSave.begin(), ::tolower);
      m_FileOutput << "export mesh '"<< szSave <<"'" << " overwrite"<<std::endl;
    }
  else{
      std::string szSave = m_szFile + ".cub";
      std::transform(szSave.begin(), szSave.end(), szSave.begin(), ::tolower);
      m_FileOutput << "save as '"<< szSave <<"'" << " overwrite"<<std::endl;
    }

  std::cout << "Schemes file created: " << m_szSchFile << std::endl;
  std::cout << "Cubit journal file created: " << m_szJouFile << std::endl;
  if(strcmp(m_szInfo.c_str(),"on") == 0)
    std::cout << "Assembly info file created: " << m_szAssmInfo << std::endl;

  m_FileOutput << "Timer Stop" << std::endl;

}
}
