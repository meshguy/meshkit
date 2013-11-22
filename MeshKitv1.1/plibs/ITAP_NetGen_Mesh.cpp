#include <meshing.hpp>

#include "ITAP_NetGen_EdgeMesh.hpp"
#include "ITAP_NetGen_SurfMesh.hpp"

#define TCL_OK 0
#define TCL_ERROR 1

  /////////////////////////////////////////////////////////////////////////////

/*
  int Mesh_Analysis(iGeom_Instance & geom, Mesh &mesh,
         int perfstepsstart, int perfstepsend )
   {
         mesh->geomtype = Mesh::GEOM_ITAP;
         mesh->SetGlobalH (mparam.maxh);
         mesh->SetMinimalH (mparam.minh);

         Array<double> maxhdom;
         maxhdom.SetSize (geom.NrSolids());
         maxhdom = mparam.maxh;

         mesh->SetMaxHDomain (maxhdom);

         Box<3> bb = geom.GetBoundingBox();
         bb.Increase (bb.Diam()/10);

         mesh->SetLocalH (bb.PMin(), bb.PMax(), 0.5);

         if (mparam.uselocalh)
         {
            mesh->SetLocalH (bb.PMin(), bb.PMax(), mparam.grading);

            int nedges = geom.emap.Extent();

            double maxedgelen = 0;
            double minedgelen = 1e99;

            // setting elements per edge

            for (int i = 1; i <= nedges; i++)
            {
               TopoDS_Edge e = TopoDS::Edge (geom.emap(i));
               if (BRep_Tool::Degenerated(e)) continue;

               GProp_GProps system;
               BRepGProp::LinearProperties(e, system);
               double len = system.Mass();

               if (len < IGNORECURVELENGTH)
               {
                  (*testout) << "ignored" << endl;
                  continue;
               }

               double localh = len/mparam.segmentsperedge;
               double s0, s1;

               // Philippose - 23/01/2009
               // Find all the parent faces of a given edge
               // and limit the mesh size of the edge based on the
               // mesh size limit of the face
               TopTools_IndexedDataMapOfShapeListOfShape edge_face_map;
               edge_face_map.Clear();

               TopExp::MapShapesAndAncestors(geom.shape, TopAbs_EDGE, TopAbs_FACE, edge_face_map);
               const TopTools_ListOfShape& parent_faces = edge_face_map.FindFromKey(e);

               TopTools_ListIteratorOfListOfShape parent_face_list;

               for(parent_face_list.Initialize(parent_faces); parent_face_list.More(); parent_face_list.Next())
               {
                  TopoDS_Face parent_face = TopoDS::Face(parent_face_list.Value());

                  int face_index = geom.fmap.FindIndex(parent_face);

                  if(face_index >= 1) localh = min(localh,geom.face_maxh[face_index - 1]);
               }

               Handle(Geom_Curve) c = BRep_Tool::Curve(e, s0, s1);

               maxedgelen = max (maxedgelen, len);
               minedgelen = min (minedgelen, len);

               // Philippose - 23/01/2009
               // Modified the calculation of maxj, because the
               // method used so far always results in maxj = 2,
               // which causes the localh to be set only at the
               // starting, mid and end of the edge.
               // Old Algorithm:
               // int maxj = 2 * (int) ceil (localh/len);
               int maxj = max((int) ceil(len/localh), 2);

               for (int j = 0; j <= maxj; j++)
               {
                  gp_Pnt pnt = c->Value (s0+double(j)/maxj*(s1-s0));
                  mesh->RestrictLocalH (Point3d(pnt.X(), pnt.Y(), pnt.Z()), localh);
               }
            }

            // setting edge curvature

            int nsections = 20;

            for (int i = 1; i <= nedges; i++)
            {
               double maxcur = 0;
               TopoDS_Edge edge = TopoDS::Edge (geom.emap(i));
               if (BRep_Tool::Degenerated(edge)) continue;
               double s0, s1;
               Handle(Geom_Curve) c = BRep_Tool::Curve(edge, s0, s1);
               BRepAdaptor_Curve brepc(edge);
               BRepLProp_CLProps prop(brepc, 2, 1e-5);

               for (int j = 1; j <= nsections; j++)
               {
                  double s = s0 + j/(double) nsections * (s1-s0);
                  prop.SetParameter (s);
                  double curvature = prop.Curvature();
                  if(curvature> maxcur) maxcur = curvature;

                  if (curvature >= 1e99)
                  continue;

                  gp_Pnt pnt = c->Value (s);

                  mesh->RestrictLocalH (Point3d(pnt.X(), pnt.Y(), pnt.Z()),
                        ComputeH (fabs(curvature)));
               }
               // (*testout) << "edge " << i << " max. curvature: " << maxcur << endl;
            }

            // setting face curvature

            int nfaces = geom.fmap.Extent();

            for (int i = 1; i <= nfaces; i++)
            {
               TopoDS_Face face = TopoDS::Face(geom.fmap(i));
               TopLoc_Location loc;
               Handle(Geom_Surface) surf = BRep_Tool::Surface (face);
               Handle(Poly_Triangulation) triangulation = BRep_Tool::Triangulation (face, loc);
               if (triangulation.IsNull()) continue;

               BRepAdaptor_Surface sf(face, Standard_True);
               BRepLProp_SLProps prop(sf, 2, 1e-5);

               int ntriangles = triangulation -> NbTriangles();
               for (int j = 1; j <= ntriangles; j++)
               {
                  gp_Pnt p[3];
                  gp_Pnt2d par[3];

                  for (int k = 1; k <=3; k++)
                  {
                     int n = triangulation->Triangles()(j)(k);
                     p[k-1] = triangulation->Nodes()(n).Transformed(loc);
                     par[k-1] = triangulation->UVNodes()(n);
                  }

                  //double maxside = 0;
                  //maxside = max (maxside, p[0].Distance(p[1]));
                  //maxside = max (maxside, p[0].Distance(p[2]));
                  //maxside = max (maxside, p[1].Distance(p[2]));
                  //cout << "\rFace " << i << " pos11 ntriangles " << ntriangles << " maxside " << maxside << flush;

                  RestrictHTriangle (par[0], par[1], par[2], &prop, *mesh, 0);
                  //cout << "\rFace " << i << " pos12 ntriangles " << ntriangles << flush;
               }
            }

            // setting close edges

            if (stlparam.resthcloseedgeenable)
            {
               int sections = 100;

               Array<Line> lines(sections*nedges);

               Box3dTree* searchtree =
               new Box3dTree (bb.PMin(), bb.PMax());

               int nlines = 0;
               for (int i = 1; i <= nedges; i++)
               {
                  TopoDS_Edge edge = TopoDS::Edge (geom.emap(i));
                  if (BRep_Tool::Degenerated(edge)) continue;

                  double s0, s1;
                  Handle(Geom_Curve) c = BRep_Tool::Curve(edge, s0, s1);
                  BRepAdaptor_Curve brepc(edge);
                  BRepLProp_CLProps prop(brepc, 1, 1e-5);
                  prop.SetParameter (s0);

                  gp_Vec d0 = prop.D1().Normalized();
                  double s_start = s0;
                  int count = 0;
                  for (int j = 1; j <= sections; j++)
                  {
                     double s = s0 + (s1-s0)*(double)j/(double)sections;
                     prop.SetParameter (s);
                     gp_Vec d1 = prop.D1().Normalized();
                     double cosalpha = fabs(d0*d1);
                     if ((j == sections) || (cosalpha < cos(10.0/180.0*M_PI)))
                     {
                        count++;
                        gp_Pnt p0 = c->Value (s_start);
                        gp_Pnt p1 = c->Value (s);
                        lines[nlines].p0 = Point<3> (p0.X(), p0.Y(), p0.Z());
                        lines[nlines].p1 = Point<3> (p1.X(), p1.Y(), p1.Z());

                        Box3d box;
                        box.SetPoint (Point3d(lines[nlines].p0));
                        box.AddPoint (Point3d(lines[nlines].p1));

                        searchtree->Insert (box.PMin(), box.PMax(), nlines+1);
                        nlines++;

                        s_start = s;
                        d0 = d1;
                     }
                  }
               }

               Array<int> linenums;

               for (int i = 0; i < nlines; i++)
               {
                  Line & line = lines[i];

                  Box3d box;
                  box.SetPoint (Point3d(line.p0));
                  box.AddPoint (Point3d(line.p1));
                  double maxhline = max (mesh->GetH(box.PMin()),
                        mesh->GetH(box.PMax()));
                  box.Increase(maxhline);

                  double mindist = 1e99;
                  linenums.SetSize(0);
                  searchtree->GetIntersecting(box.PMin(),box.PMax(),linenums);

                  for (int j = 0; j < linenums.Size(); j++)
                  {
                     int num = linenums[j]-1;
                     if (i == num) continue;
                     if ((line.p0-lines[num].p0).Length2() < 1e-15) continue;
                     if ((line.p0-lines[num].p1).Length2() < 1e-15) continue;
                     if ((line.p1-lines[num].p0).Length2() < 1e-15) continue;
                     if ((line.p1-lines[num].p1).Length2() < 1e-15) continue;
                     mindist = min (mindist, line.Dist(lines[num]));
                  }

                  mindist *= stlparam.resthcloseedgefac;

                  if (mindist < 1e-3)
                  {
                     (*testout) << "extremely small local h: " << mindist
                     << " --> setting to 1e-3" << endl;
                     (*testout) << "somewhere near " << line.p0 << " - " << line.p1 << endl;
                     mindist = 1e-3;
                  }

                  mesh->RestrictLocalHLine(line.p0, line.p1, mindist);
               }
            }

         }

         // Philippose - 09/03/2009
         // Added the capability to load the mesh size from a 
         // file also for OpenCascade Geometry
         // Note: 
         // ** If the "uselocalh" option is ticked in 
         // the "mesh options...insider" menu, the mesh 
         // size will be further modified by the topology 
         // analysis routines.
         // ** To use the mesh size file as the sole source 
         // for defining the mesh size, uncheck the "uselocalh"
         // option.
         mesh->LoadLocalMeshSize (mparam.meshsizefilename);
      }

      if (perfstepsend <= MESHCONST_ANALYSE) return TCL_OK;

}
*/

//////////////////////////////////////////////////////////////////////////////////

int ITAP_NetGen_GenerateMesh ( iGeom_Instance & geom, Mesh &mesh, 
                               int perfstepsstart, int perfstepsend)

      /////////////////////////////////////////////////////////////////////////
  
      if (perfstepsstart <= MESHCONST_MESHEDGES)
      {
         Generate_EdgeMesh(geom, mesh);
      }

      if (perfstepsend <= MESHCONST_MESHEDGES) return TCL_OK;

      /////////////////////////////////////////////////////////////////////////


      if (perfstepsstart <= MESHCONST_MESHSURFACE)
      {
         ITAPMeshSurface (geom, mesh, perfstepsend);

         MeshQuality2d (mesh);

         mesh->CalcSurfacesOfNode();
      }

      if (perfstepsend <= MESHCONST_OPTSURFACE) return TCL_OK;

      /////////////////////////////////////////////////////////////////////////

      if (perfstepsstart <= MESHCONST_MESHVOLUME)
      {
         MESHING3_RESULT res = MeshVolume (mparam, *mesh);

         if (res == MESHING3_OK) {
             RemoveIllegalElements (*mesh);
             MeshQuality3d (*mesh);
         } else 
             return TCL_ERROR;
      }

      if (perfstepsend <= MESHCONST_MESHVOLUME) return TCL_OK;

      /////////////////////////////////////////////////////////////////////////

      if (perfstepsstart <= MESHCONST_OPTVOLUME)
      {
         OptimizeVolume (mparam, *mesh);
      }

      /////////////////////////////////////////////////////////////////////////

      for (int i = 1; i <= mesh->GetNP(); i++)
      (*testout) << mesh->Point(i) << endl;

      (*testout) << endl << "NSegments: " << mesh->GetNSeg() << endl;
      for (int i = 1; i <= mesh->GetNSeg(); i++)
      (*testout) << mesh->LineSegment(i) << endl;

      return TCL_OK;
   }
}
