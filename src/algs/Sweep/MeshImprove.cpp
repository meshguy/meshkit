#include "meshkit/MeshImprove.hpp"
#include <iostream>
#include <math.h>
#include <map>

//#include "SweepWrapper.hpp"
#include "UntangleWrapper.hpp"
#include "ShapeImprover.hpp"
#include "SmartLaplacianSmoother.hpp"
#include "SweepWrapper.hpp"
#include "SmartLaplaceWrapper.hpp"
#include "ShapeImprovementWrapper.hpp"

#include "IdealWeightInverseMeanRatio.hpp"
#include "QualityAssessor.hpp"
#include "InstructionQueue.hpp"
#include "TerminationCriterion.hpp"

#include "LaplaceWrapper.hpp"
#include "SizeAdaptShapeWrapper.hpp"
#include "PaverMinEdgeLengthWrapper.hpp"
#include "DeformingDomainWrapper.hpp"
#include <MsqError.hpp>
#include <ShapeImprovementWrapper.hpp>
#include <MsqIMesh.hpp>
#include <MsqIGeom.hpp>

using namespace Mesquite; 

namespace MeshKit {

MeshImprove::MeshImprove(MKCore* core, bool isLaplacian, bool isUntangle, bool isShapeImprove, bool isSizeAdapt)
{
	mk_core = core;
	//mk = core;
	/*
	geom = geometry;
	mesh = Mesh;
	assoc = association;
	rel   = irel;
	
	int err;
	iMesh_getRootSet(mesh, &mesh_root_set, &err);
    	assert(!err);
	
	iMesh_getTagHandle(mesh, "GLOBAL_ID", &mesh_id_tag, &err, strlen("GLOBAL_ID"));
    	assert(!err);
    	
    	max = maxsize;
    	
    	int iError;
    	iMesh_save(mesh, mesh_root_set, "NotSmooth.vtk", NULL, &iError, strlen("NotSmooth.vtk"), 0);
	assert(!iError);
	*/
	IsLaplacian = isLaplacian;
	IsUntangle = isUntangle;
	IsShapeImprove = isShapeImprove;
	IsSizeAdapt = isSizeAdapt;
}

void MeshImprove::SurfMeshImprove(iBase_EntityHandle surface, iBase_EntitySetHandle surfMesh, iBase_EntityType entity_type)
{
	
	MsqError mError;
	const char* VERTEX_FIXED_TAG_NAME="MesquiteVertexFixed";
	
	iBase_TagHandle fixed_tag=0;
	iMesh::Error m_err = mk_core->imesh_instance()->getTagHandle(VERTEX_FIXED_TAG_NAME, fixed_tag);
	if (m_err)
	{
		m_err = mk_core->imesh_instance()->createTag(VERTEX_FIXED_TAG_NAME, 1, iBase_INTEGER, fixed_tag);
		IBERRCHK(m_err, "Trouble create the tag handle.");	
 	}
	
	MsqIMesh mesh_adapter(mk_core->imesh_instance()->instance(), surfMesh, entity_type, mError, &fixed_tag);
	cout << "error =" << mError << endl;
	if (mError) throw mError;
	
	
	
	//get all the vertices in surface mesh: entity_handles_out---quads     adj_entity_handles_out---vertices
	std::vector<iBase_EntityHandle> entity_handles_out, adj_entity_handles_out;
	std::vector<int> offsets_out, adj_entity_indices_out;
	
	m_err = mk_core->imesh_instance()->getAdjEntIndices(surfMesh, entity_type, iMesh_ALL_TOPOLOGIES, iBase_VERTEX, entity_handles_out, adj_entity_handles_out, adj_entity_indices_out, offsets_out);
	IBERRCHK(m_err, "Trouble get the adjacent entity indices.");

	
	cout << "number of faces is " << entity_handles_out.size() << endl;
	cout << "number of vertices is " << adj_entity_handles_out.size() << endl;
	
	//set fixed flag on all vertices
	std::vector<int> tag_data(adj_entity_handles_out.size(), 1);
	m_err = mk_core->imesh_instance()->setIntArrData(&adj_entity_handles_out[0], adj_entity_handles_out.size(), fixed_tag, &tag_data[0]);
	IBERRCHK(m_err, "Trouble set an array of int data for a list of vertices.");
		
	//clear fixed flag for vertices contained directly in set
	int count = -1;
	m_err = mk_core->imesh_instance()->getNumOfType(surfMesh, iBase_VERTEX, count);
	IBERRCHK(m_err, "Trouble get the number of vertices in the set.");

	adj_entity_handles_out.clear();
		
	cout << "Num of Vertices on the target surface is " << count << endl;
			
	m_err = mk_core->imesh_instance()->getEntities(surfMesh, iBase_VERTEX, iMesh_ALL_TOPOLOGIES, adj_entity_handles_out);
	IBERRCHK(m_err, "Trouble get the nodes from the mesh entity set.");

	
	tag_data.clear();
	tag_data.resize( adj_entity_handles_out.size(), 0 );
		
	m_err = mk_core->imesh_instance()->setIntArrData(&adj_entity_handles_out[0], adj_entity_handles_out.size(), fixed_tag, &tag_data[0]);
	IBERRCHK(m_err, "Trouble set an array of int data for mesh nodes.");
	
	//Finally smooth the mesh
	//ShapeImprovementWrapper smoother(mError);
	
	
	
	
	//SweepWrapper smoother( 1e-6,  "COORDINATES_MAP");
	if (IsUntangle){
		UntangleWrapper smoother;
		//smoother.set_cpu_time_limit(300);
		smoother.set_vertex_movement_limit_factor(0.001);	

		if (surface){
			MsqIGeom geom_adapter(mk_core->igeom_instance()->instance(), surface);
			smoother.run_instructions(&mesh_adapter, &geom_adapter, mError);
			cout << "Mesquite error in surface mesh smoothing=" << mError << endl;
		}
		else{
			smoother.run_instructions(&mesh_adapter, mError);
			cout << "Mesquite error in surface mesh smoothing=" << mError << endl;
		}
	}
	//=======================================================================================================//
	//use the smart laplace smoothing algorithm to smooth the target surface

	/*
	std::cout << "Starting the smart Laplace smoothing\n";
	IdealWeightInverseMeanRatio qa_metric;			
	QualityAssessor* qa;
	qa->add_quality_assessment( &qa_metric );
	SmartLaplacianSmoother sl_smooth;
	TerminationCriterion outer, inner;
	inner.cull_on_absolute_vertex_movement_edge_length( 0.001 );
    	sl_smooth.set_inner_termination_criterion( &inner );

	InstructionQueue q;
	MsqError err;
  	q.add_quality_assessor( qa, err ); MSQ_ERRRTN(err);
  	q.set_master_quality_improver( &sl_smooth, err ); MSQ_ERRRTN(err);
  	q.add_quality_assessor( qa, err ); MSQ_ERRRTN(err);
	if (surface){
	    MsqIGeom geom_adapter(mk_core->igeom_instance()->instance(), surface);
  	    q.run_common( &mesh_adapter, NULL, &geom_adapter, NULL, err ); MSQ_ERRRTN(err);

	}
	std::cout << "end the smart laplace smoothing\n";	
	*/
	
	
	//use the SmartLaplaceWrapper class the smooth the target surface mesh
	if (IsLaplacian){
		LaplaceWrapper sl_smooth;
		TerminationCriterion terminate;
		terminate.write_iterations( "mesquite.gpt", mError);
		if (surface){
	    		MsqIGeom geom_adapter(mk_core->igeom_instance()->instance(), surface);
            		sl_smooth.run_instructions(&mesh_adapter, &geom_adapter, mError);
           		cout << "Mesquite error in the smart Laplacian surface mesh smoothing with the geometry domain=" << mError << endl;    
		}
		else{
	    		sl_smooth.run_instructions(&mesh_adapter, mError);
            		cout << "Mesquite error in the smart Laplacian surface mesh smoothing without the geometry domain=" << mError << endl;
		}
	}
	

	//use the ShapeImprover class to smooth the target surface mesh
	if (IsShapeImprove){
		IdealWeightInverseMeanRatio extra_metric;
		ShapeImprovementWrapper smoother1;
		smoother1.quality_assessor().add_quality_assessment(&extra_metric);
		//smoother1.set_vertex_movement_limit_factor(0.01);
		//smoother1.set_cpu_time_limit(300);
		if (surface){
		    MsqIGeom geom_adapter(mk_core->igeom_instance()->instance(), surface);
		    smoother1.run_instructions(&mesh_adapter, &geom_adapter, mError);
		    cout << "Mesquite error in the ShapeImprover surface mesh smoothing=" << mError << endl;
		}
		else{
		    smoother1.run_instructions(&mesh_adapter, mError);
		    cout << "Mesquite error in the ShapeImprover surface mesh smoothing=" << mError << endl;
		}
	}
	
		
	

	/*

	dealWeightInverseMeanRatio extra_metric;
        ShapeImprovementWrapper smoother1;
        smoother1.quality_assessor().add_quality_assessment(&extra_metric);
        //smoother1.set_vertex_movement_limit_factor(0.01);
        //smoother1.set_cpu_time_limit(300);
        if (surface){
            MsqIGeom geom_adapter(mk_core->igeom_instance()->instance(), surface);
            smoother1.run_instructions(&mesh_adapter, &geom_adapter, mError);
            cout << "Mesquite error in the ShapeImprover surface mesh smoothing=" << mError << endl;
        }
        else{
            smoother1.run_instructions(&mesh_adapter, mError);
            cout << "Mesquite error in the ShapeImprover surface mesh smoothing=" << mError << endl;
        }TerminationCriterion terminate;
	terminate.write_iterations( "mesquite.gpt", mError);	
	//use the Laplace-Smoothing to smooth the target surface
        DeformingDomainWrapper smoother2;
	//smoother2.set_vertex_movement_limit_factor(0.01);
	if (surface){
	    MsqIGeom geom_adapter(mk_core->igeom_instance()->instance(), surface);
	    smoother2.run_instructions(&mesh_adapter, &geom_adapter, mError);
	    cout << "Mesquite error in LaplaceWrapper mesh smoothing=" << mError << endl;
	}
	else{
	    smoother2.run_instructions(&mesh_adapter, mError);
            cout << "Mesquite error in LaplaceWrapper mesh smoothing=" << mError << endl;
	}	
	
	
	/*
	SweepWrapper sw_smoother( 1e-6,  "COORDINATES_MAP");
	if (surface)
	{
		MsqIGeom geom_adapter(mk_core->igeom_instance()->instance(), surface);
		sw_smoother.run_instructions(&mesh_adapter, &geom_adapter, mError);
		cout << "Mesquite errorPaverMinEdgeLengthWrapper.hpp in surface mesh smoothing=" << mError << endl;
	}
	else
	{
		sw_smoother.run_instructions(&mesh_adapter, mError);
		cout << "Mesquite error in surface mesh smoothing=" << mError << endl;
	}
	*/


	//Use the Minimum Edge-Length Improvement

	if (IsSizeAdapt){
		SizeAdaptShapeWrapper Smooth_SizeAdapt(1.0e-2);
		if (surface){
		    MsqIGeom geom_adapter(mk_core->igeom_instance()->instance(), surface);
		    Smooth_SizeAdapt.run_instructions(&mesh_adapter, &geom_adapter, mError);
		    cout << "Mesquite error in the ShapeImprover surface mesh smoothing=" << mError << endl;
		}
		else{
		    Smooth_SizeAdapt.run_instructions(&mesh_adapter, mError);
		    cout << "Mesquite error in the ShapeImprover surface mesh smoothing=" << mError << endl;
		}
	}
	

}

void MeshImprove::VolumeMeshImprove(iBase_EntitySetHandle volMesh, iBase_EntityType entity_type)
{
		
	cout << "Volume smoothing is starting..." << endl;
	
	MsqError mError;
	const char* VERTEX_FIXED_TAG_NAME="MesquiteVertexFixed";
	
	iBase_TagHandle fixed_tag=0;
	
	iMesh::Error m_err = mk_core->imesh_instance()->getTagHandle(VERTEX_FIXED_TAG_NAME, fixed_tag);
	if (m_err)
	{
		m_err = mk_core->imesh_instance()->createTag(VERTEX_FIXED_TAG_NAME, 1, iBase_INTEGER, fixed_tag);
 		IBERRCHK(m_err, "Trouble create a taghandle.");
 	}
	
	MsqIMesh mesh_adapter(mk_core->imesh_instance()->instance(), volMesh, entity_type, mError, &fixed_tag);
	cout << "error =" << mError << endl;
	if (mError) throw mError;
		
	//get all the vertices in volume mesh
	int num_vtx, count;
	std::vector<iBase_EntityHandle> faces, verts;	
	std::vector<int> indices, offsets;
		
	m_err = mk_core->imesh_instance()->getAdjEntIndices(volMesh, entity_type, iMesh_ALL_TOPOLOGIES, iBase_VERTEX, faces, verts, indices, offsets);
	IBERRCHK(m_err, "Trouble get the quads and nodes on the target surface.");
	num_vtx = verts.size();
	
	cout << "number of faces is " << faces.size() << endl;
	cout << "number of vertices is " << verts.size() << endl;
	
	//set fixed flag on all vertices
	vector<int> tag_data(num_vtx, 1);
	m_err = mk_core->imesh_instance()->setIntArrData(&verts[0], verts.size(), fixed_tag, &tag_data[0]);
	IBERRCHK(m_err, "Trouble set an array of int data for nodes on the target surface.");
		
	//clear fixed flag for vertices contained directly in set
	m_err = mk_core->imesh_instance()->getNumOfType(volMesh, iBase_VERTEX, count);
	IBERRCHK(m_err, "Trouble get the number of interior nodes on the target surface.");

	//get the interior mesh nodes on the target surface
	verts.clear();
	cout << "Num of Vertices on the target surface is " << count << endl;
			
	m_err = mk_core->imesh_instance()->getEntities(volMesh, iBase_VERTEX, iMesh_ALL_TOPOLOGIES, verts);
	IBERRCHK(m_err, "Trouble get the number of interior nodes on the target surface.");
	

	tag_data.clear();
	tag_data.resize( verts.size(), 0);
	m_err = mk_core->imesh_instance()->setIntArrData(&verts[0], verts.size(), fixed_tag, &tag_data[0]);
	IBERRCHK(m_err, "Trouble set the int data for interior nodes on the target surface.");
	
	
	//using the UntangleWrapper class to smooth the mesh with the inverted elements.
	UntangleWrapper smoother1;
	
	smoother1.set_cpu_time_limit(1000);
	smoother1.set_vertex_movement_limit_factor(0.001);	
		
	smoother1.run_instructions(&mesh_adapter, mError);	

	//Finally smooth the mesh
	//Use the ShapeImprover class to smooth the target surface mesh
	ShapeImprover smoother2;
	smoother2.set_cpu_time_limit(1000);
	smoother2.set_vertex_movement_limit_factor(0.001);
	smoother2.run_instructions(&mesh_adapter, mError);
	if (mError)
	    cout << "Mesquite error in volume mesh smoothing is as follows\n" << mError << std::endl;	
	/*

	SweepWrapper smoother3( 1e-6,  "COORDINATES_MAP");
	ShapeImprovementWrapper smoother2(mError);
	
	
	smoother2.run_instructions(&mesh_adapter, mError);	
	
	if (mError)
		cout << "Mesquite error in volume mesh smoothing=" << mError << endl;
	*/
	
}

MeshImprove::~MeshImprove()
{
	cout << "It is over now in smoothing" << endl;
}

}
