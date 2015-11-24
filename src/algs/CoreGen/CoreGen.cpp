#include "meshkit/CoreGen.hpp"
#include "moab/VerdictWrapper.hpp"

#define ERROR(a) {if (iBase_SUCCESS != err) std::cerr << a << std::endl;}
#define ERRORR(a,b) {if (iBase_SUCCESS != err) {std::cerr << a << std::endl; return b;}}
namespace MeshKit
{
  // static registration of this  mesh scheme
  moab::EntityType CoreGen_tps[] = { moab::MBVERTEX,
                                     moab::MBEDGE,
                                     moab::MBTRI,
                                     moab::MBHEX,
                                     moab::MBMAXTYPE};
  const moab::EntityType* CoreGen::output_types()
  { return CoreGen_tps; }

  CoreGen::CoreGen( MKCore *mk, const MEntVector &me_vec)
    : MeshScheme( mk, me_vec),
      igeom(mk->igeom_instance()), imesh(mk->imesh_instance()),
      mb (mk->moab_instance())
  {
    err = 0;
    run_flag = 1;
    UNITCELL_DUCT = 0;
    ASSY_TYPES = 1;
    pack_type = 1;
    symm = 1;
    z_height = 1;
    z_divisions = 2;
    set_DIM = 3; // default is 3D
    PII = acos(-1.0);
    comment = "!";
    MAXCHARS = 2000;
    compute_meshtogeom = false;
    extrude_flag = false;
    mem_tflag = false;
    global_ids = true;
    merge_tol = 1.0e-4;
    do_merge = 1;
    update_sets = 0;
    merge_tag = NULL;
    nst_flag = false;
    nsb_flag = false;
    nss_flag = false;
    nssall_flag = false;
    nst_Id = 9997;
    nsb_Id = 9998;
    prob_type = "mesh";
    savefiles = "one";
    num_nsside = 0;
    linenumber = 0;
    info = "off";
    minfo = "off";
    nringsx = 0;
    nringsy = 0;

    // initialize more memory/time related variables
    ctload = 0, ctcopymove = 0, ctmerge = 0, ctextrude = 0, ctns = 0, ctgid = 0, ctsave = 0;
    tload = 0, tcopymove = 0, tmerge = 0, textrude = 0, tns = 0, tgid = 0, tsave = 0;
    ld_t = 0, ld_tload = 0, ld_tcopymove = 0, ld_tsave = 0, ld_tgid = 0, ld_tmerge = 0, ld_tns = 0;
    mem1 = 0, mem2 = 0, mem3 = 0, mem4 = 0, mem5 = 0, mem6 = 0, mem7 = 0;


  }

  CoreGen::~CoreGen()
  {}

  bool CoreGen::add_modelent(ModelEnt *model_ent)
  {
    return MeshOp::add_modelent(model_ent);
  }

  void CoreGen::setup_this()
  {
    if(rank == 0)
      logfile << "Setting-up in CoreGen meshop.." << std::endl;
    if(run_flag != 0){
        double ctload = 0;
        clock_t tload = 0;
        CClock ld_time;
        int ld_t = 0, ld_tload = 0;
        unsigned long long mem1 = 0;
        if (prob_type == "mesh") {
            if (procs == 1) {
                err = load_meshes();
                if(err !=0) {std::cout << "load meshes failed!" << std::endl; exit(1);}
              }
            else {
#ifdef USE_MPI
                err = load_meshes_parallel(rank, procs);
                if(err !=0) {std::cout << "failed to load meshes in parallel!" << std::endl; exit(1);}

                if(procs > (int) files.size()){
                    // if there are more procs than files distribute the copy/move work on each proc
                    err = distribute_mesh(rank, procs);
                    if(err !=0) {std::cout << "distribute meshes failed!" << std::endl; exit(1);}
                  }
                //Get a pcomm object
                pc = new moab::ParallelComm(mk_core()->moab_instance(), MPI::COMM_WORLD, &err);
#endif
              }

            mk_core()->moab_instance()->estimated_memory_use(0, 0, 0, &mem1);
            ld_tload = ld_time.DiffTime();
            tload = clock();
            ctload = (double) (tload - sTime)/(60*CLOCKS_PER_SEC);

            if (mem_tflag == true && procs == 1) {
                logfile << "\n" << " Clock time taken to load mesh files = " << ld_tload
                        << " seconds" << std::endl;
                logfile << " CPU time = " << ctload << " mins" << std::endl;
                logfile << " Memory used: " << mem1/1e6 << " Mb\n" << std::endl;
              }
          }

        /*********************************************/
        // load geometry files
        /*********************************************/
        else if (prob_type == "geometry" && procs == 1) {
            err = load_geometries();
            if(err !=0) {std::cout << "load geometry failed!" << std::endl; exit(1);}
          }
        else if(prob_type == "geometry" && procs > 1){
            logfile << " Parallel mode not supported for problem-type: Geometry " << std::endl;
            exit(1);
          }

        /*********************************************/
        // copy move
        /*********************************************/
        CClock ld_cm;
        err = copymove(rank, procs);
        if (prob_type == "mesh"){
            mk_core()->moab_instance()->estimated_memory_use(0, 0, 0, &mem2);
            ld_tcopymove = ld_cm.DiffTime();
            tcopymove = clock();
            ctcopymove = (double) (tcopymove - tload)/(60*CLOCKS_PER_SEC);

            if (mem_tflag == true && (strcmp(prob_type.c_str(), "mesh") == 0) && procs == 1) {
                logfile << "\n" << " Clock time taken to copy/move mesh files = " << ld_tcopymove
                        << " seconds" << std::endl;
                logfile << " CPU time = " << ctcopymove << " mins" << std::endl;
                logfile << " Memory used: " << mem2/1e6 << " Mb\n" << std::endl;
              }
          }

        for (unsigned int i = 0; i < assys.size(); i++) {
            if(prob_type =="mesh")
              cm[i]->setup_called(true);
            if(prob_type =="geometry")
              cg[i]->setup_called(true);
          }

        if (prob_type == "mesh") {
            /*********************************************/
            // merge
            /*********************************************/
            CClock ld_mm;
            if (procs == 1){
                // merge mesh now
                //std::vector<iBase_EntityHandle> ents;
                moab::Range ents;
                //int dim = imesh->getGeometricDimension();
                mb->get_entities_by_dimension(0, set_DIM, ents);
                logfile << " Merging nodes.."<< std::endl;
                moab::MergeMesh mm(mb);
                moab::ErrorCode err = mm.merge_entities(ents, merge_tol, true);
                if (moab::MB_SUCCESS != err) {
                    std::cerr << "Error in MergeMesh during merging entities" << std::endl;
                    exit(2);
                  }
              }
            else if(procs > 1){
                if (rank == 0) {
                    logfile << "Merging nodes in parallel. " << std::endl;
                  }

#ifdef USE_MPI
                //Call the resolve parallel function
                moab::ParallelMergeMesh pm(pc, merge_tol);
                err = pm.merge();
                if (err != moab::MB_SUCCESS) {
                    std::cerr << "Merge Failed" << std::endl;
                    //MPI_Abort(MPI_COMM_WORLD, 1);
                  }
#endif
              }
            mk_core()->moab_instance()->estimated_memory_use(0, 0, 0, &mem3);
            ld_tmerge = ld_mm.DiffTime();
            tmerge = clock();
            ctmerge = (double) (tmerge - tcopymove)/(60*CLOCKS_PER_SEC);

            if (mem_tflag == true && procs == 1 ) {
                logfile << "\n" << " Clock time taken to merge nodes = " << ld_tmerge
                        << " seconds" << std::endl;
                logfile << " CPU time = " << ctmerge << " mins" << std::endl;
                logfile << " Memory used: " << mem3/1e6 << " Mb\n" << std::endl;
              }
#ifdef USE_MPI
            MPI::COMM_WORLD.Barrier();
#endif
            /*********************************************/
            // extrude
            /*********************************************/
            if(procs == 1){
                if (extrude_flag == true) {
                    CClock ld_em;
                    err = extrude();
                    if(err !=0) {std::cout << "extrusion failed!" << std::endl; exit(1);}

                    mk_core()->moab_instance()->estimated_memory_use(0, 0, 0, &mem4);
                    ld_t = ld_em.DiffTime();
                    textrude = clock();
                    ctextrude = (double) (textrude - tmerge)/(60*CLOCKS_PER_SEC);

                    if (mem_tflag == true && procs == 1) {
                        logfile << "\n" << " Clock time taken to extrude = " << ld_t
                                << " seconds" << std::endl;
                        logfile << " CPU time = " << ctextrude << " mins" << std::endl;
                        logfile << " Memory used: " << mem4/1e6 << " Mb\n"
                                << std::endl;
                      }
                  }
              }
            if(extrude_flag == true)
              em->setup_called(true);
#ifdef USE_MPI
            MPI::COMM_WORLD.Barrier();
#endif
          }
      }
  }


  void CoreGen::execute_this()
  {
    if(rank == 0)
      logfile << "Executing in CoreGen meshop.." << std::endl;

    if(run_flag != 0 && prob_type != "geometry"){
        for (unsigned int i = 0; i < assys.size(); i++) {
            cm[i]->execute_called(true);
          }
        if(extrude_flag == true)
          em->execute_called(true);
        /*********************************************/
        // assign gids
        /*********************************************/
        CClock ld_gid;
        mk_core()->moab_instance()->estimated_memory_use(0, 0, 0, &mem5);
        ld_tgid = ld_gid.DiffTime();
        tgid = clock();
        ctgid = (double) (tgid-tmerge)/(60*CLOCKS_PER_SEC);

        if (mem_tflag == true && procs == 1) {
            logfile << "\n" << " Clock time taken to assign gids = " << ld_tgid
                    << " seconds" << std::endl;
            logfile << " CPU time = " << ctgid << " mins" << std::endl;
            logfile << " Memory used: " << mem5/1e6 << " Mb\n" << std::endl;
          }
        /*********************************************/
        // create neumann sets on the core model
        /*********************************************/
        if((nss_flag == true || nsb_flag == true
            || nst_flag == true) && procs == 1){
            CClock ld_ns;
            err = create_neumannset();
            if(err !=0) {std::cout << "create NeumannSet failed!" << std::endl; exit(1);}

            mk_core()->moab_instance()->estimated_memory_use(0, 0, 0, &mem6);
            ld_tns = ld_ns.DiffTime();
            tns = clock();
            ctns = (double) (tns-tgid)/(60*CLOCKS_PER_SEC);
            if (mem_tflag == true && procs == 1) {
                logfile << "\n" << " Clock time taken to create neumann sets = " << ld_tns
                        << " seconds" << std::endl;
                logfile << " CPU time = " << ctns << " mins" << std::endl;
                logfile << " Memory used: " << mem6/1e6 << " Mb\n" << std::endl;
              }
          }
      }
    if (prob_type == "mesh") {
        /*********************************************/
        // save
        /*********************************************/
        CClock ld_sv;
        if (procs == 1) {
            err = save_mesh();
            if(err !=0) {std::cout << "save mesh failed!" << std::endl; exit(1);}
          } else {
            if(savefiles != "one" && (savefiles == "multiple" || savefiles == "both")){
                err = save_mesh(rank); // uncomment to save the meshes with each proc
                if(err !=0) {std::cout << "save mesh failed!" << std::endl; exit(1);}
              }
            if(savefiles != "multiple"){
#ifdef USE_MPI
                double write_time = MPI_Wtime();
                err = save_mesh_parallel(rank, procs);
                if(err !=0) {std::cout << "save parallel mesh failed!" << std::endl; exit(1);}
                write_time = MPI_Wtime() - write_time;
                if (rank == 0)
                  logfile << "Parallel write time = " << write_time/60.0 << " mins" << std::endl;
#endif
              }
          }

        mk_core()->moab_instance()->estimated_memory_use(0, 0, 0, &mem7);
        ld_tsave = ld_sv.DiffTime();
        tsave = clock();
        ctsave = (double) (tsave - tgid)/(60*CLOCKS_PER_SEC);

        if (mem_tflag == true && procs == 1 ) {
            logfile << "\n" << " Clock time taken to save = " << ld_tsave << " seconds"
                    << std::endl;
            logfile << " CPU time = " << ctsave << " mins" << std::endl;
            logfile << " Memory used: " << mem7/1e6 << " Mb\n" << std::endl;
          }
      }
    /*********************************************/
    // geometry operations
    /*********************************************/
    else if (prob_type == "geometry") {
        err = save_geometry();
        if(err !=0) {std::cout << "save geometry failed!" << std::endl; exit(1);}
      }

    /*********************************************/
    // print memory and timing if using mpi
    /*********************************************/
    mem1/=1e6;
    mem2/=1e6;
    mem3/=1e6;
    mem5/=1e6;
    mem7/=1e6;

#ifdef USE_MPI
    unsigned long max_mem7 = 1.0;
    MPI::COMM_WORLD.Reduce( &mem7, &max_mem7, 1, MPI::UNSIGNED_LONG, MPI::MAX, 0);
#endif

#ifdef USE_MPI
    if (mem_tflag == true) {

        unsigned long max_mem1 = 1.0, max_mem2 = 1.0, max_mem3 = 1.0, max_mem5 = 1.0;

        MPI::COMM_WORLD.Reduce( &mem1, &max_mem1, 1, MPI::UNSIGNED_LONG, MPI::MAX, 0);
        MPI::COMM_WORLD.Reduce( &mem2, &max_mem2, 1, MPI::UNSIGNED_LONG, MPI::MAX, 0);
        MPI::COMM_WORLD.Reduce( &mem3, &max_mem3, 1, MPI::UNSIGNED_LONG, MPI::MAX, 0);
        MPI::COMM_WORLD.Reduce( &mem5, &max_mem5, 1, MPI::UNSIGNED_LONG, MPI::MAX, 0);

        double max_ctload = -1.0, max_ctcopymove = -1.0, max_ctgid = -1.0, max_ctsave = -1.0, max_ctmerge = -1.0;
        MPI::COMM_WORLD.Reduce( &ctload, &max_ctload, 1, MPI::DOUBLE, MPI::MAX, 0);
        MPI::COMM_WORLD.Reduce( &ctcopymove, &max_ctcopymove, 1, MPI::DOUBLE, MPI::MAX, 0);
        MPI::COMM_WORLD.Reduce( &ctmerge, &max_ctmerge, 1, MPI::DOUBLE, MPI::MAX, 0);
        MPI::COMM_WORLD.Reduce( &ctgid, &max_ctgid, 1, MPI::DOUBLE, MPI::MAX, 0);
        MPI::COMM_WORLD.Reduce( &ctsave, &max_ctsave, 1, MPI::DOUBLE, MPI::MAX, 0);

        int max_tload = -1.0, max_tcopymove = -1.0, max_tgid = -1.0, max_tsave = -1.0, max_tmerge = -1.0;
        MPI::COMM_WORLD.Reduce( &ld_tload, &max_tload, 1, MPI::INT, MPI::MAX, 0);
        MPI::COMM_WORLD.Reduce( &ld_tcopymove, &max_tcopymove, 1, MPI::INT, MPI::MAX, 0);
        MPI::COMM_WORLD.Reduce( &ld_tmerge, &max_tmerge, 1, MPI::INT, MPI::MAX, 0);
        MPI::COMM_WORLD.Reduce( &ld_tgid, &max_tgid, 1, MPI::INT, MPI::MAX, 0);
        MPI::COMM_WORLD.Reduce( &ld_tsave, &max_tsave, 1, MPI::INT, MPI::MAX, 0);

        if(rank == 0 && procs > 1){
            logfile << "\nMAXIMUM TIME TAKEN OVER ALL PROCS\nCLOCK TIME:-";
            logfile << "\n**r = " << rank<< " Time taken to load mesh files = " << max_tload
                    << " secs" << std::endl;
            logfile << "***r = : " << rank<< " Memory used: " << max_mem1 << " Mb" << std::endl;

            // copymove
            logfile << "\n**r = " << rank<< " Time taken to copy/move mesh files = " << max_tcopymove
                    << " secs" << std::endl;
            logfile << "***r = " << rank<< " Memory used: " << max_mem2 << " Mb" << std::endl;

            // merge
            logfile << "\n**r = " << rank<< " Time taken to merge nodes = " << max_tmerge
                    << " secs" << std::endl;
            logfile << "***r = " << rank<< " Memory used: " << max_mem3 << " kb" << std::endl;

            // assign gid
            logfile << "\n**r = " << rank<< " Time taken to assign gids = " << max_tgid
                    << " secs" << std::endl;
            logfile << "*** r = " << rank<< " Memory used: " << max_mem5 << " Mb" << std::endl;

            // save
            logfile << "\n**r = " << rank<< " Time taken to save = " <<  max_tsave << " secs"
                    << std::endl;
            logfile << "***r = " << rank<< " Memory used: " << max_mem7 << " Mb" << std::endl;

            // cpu times
            logfile << "\n CPU TIME:-\n" << " r = " << rank<< " Time taken to load mesh files = " << ctload
                    << " mins" << std::endl;

            logfile << " r = " << rank << " Time taken to copy/move files = " << ctcopymove
                    << " mins" << std::endl;

            logfile << " r = " << rank << " Time taken to merge = " << ctmerge
                    << " mins" << std::endl;

            logfile << " r = " << rank <<  " Time taken to assign gids = " << ctgid
                    << " mins" << std::endl;

            logfile  << " r = " << rank << " Time taken to save mesh = " << ctsave
                     << " mins" << std::endl;
          }
      }
#endif

    if (rank == 0) {
        Timer.GetDateTime(szDateTime);
        logfile << "\nXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"  << std::endl;
        logfile << "Ending at : " << szDateTime;
        logfile << "Elapsed wall clock time: " << Timer.DiffTime()
                << " seconds or " << (Timer.DiffTime()) / 60.0 << " mins\n";

        logfile << "Total CPU time used: " <<  (double) (clock() - sTime)/(CLOCKS_PER_SEC) << " seconds or " <<
                   (double) (clock() - sTime)/(60*CLOCKS_PER_SEC)
                << " mins" << std::endl;
#ifdef USE_MPI
        logfile << "Maximum memory used by a processor: " << max_mem7 <<  " Mb" << std::endl;
#endif
        if(procs == 1)
          logfile << "Maximum memory used: " << mem7 <<  " Mb" << std::endl;
        logfile << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n"  << std::endl;

      }
  }

  int CoreGen::save_mesh_parallel(const int nrank, const int numprocs)
  // -------------------------------------------------------------------------------------------
  // Function: resolve shared entitie and save mesh file in parallel (hdf5 only)
  // Input:    none
  // Output:   none
  // -------------------------------------------------------------------------------------------
  {

#ifdef USE_MPI
    // write file
    if (nrank == 0) {
        logfile << "Saving mesh file in parallel. " << std::endl;
      }

    moab::Tag mattag;
    mb->tag_get_handle( "MATERIAL_SET", 1, MB_TYPE_INTEGER, mattag );
    moab::Range matsets;
    mb->get_entities_by_type_and_tag( 0, MBENTITYSET, &mattag, 0, 1, matsets );
    pc->resolve_shared_sets( matsets, mattag );

    moab::Tag nstag;
    mb->tag_get_handle( "NEUMANN_SET", 1, MB_TYPE_INTEGER, nstag );
    moab::Range nssets;
    mb->get_entities_by_type_and_tag( 0, MBENTITYSET, &nstag, 0, 1, nssets );
    pc->resolve_shared_sets( nssets, nstag );

    //  int rval = mb->write_file(outfile.c_str() , 0,"PARALLEL=WRITE_PART;DEBUG_IO=5");
    int rval = mb->write_file(outfile.c_str() , 0,"PARALLEL=WRITE_PART");
    if(rval != moab::MB_SUCCESS) {
        std::cerr<<"Writing output file failed Code:";
        std::string foo = ""; mb->get_last_error(foo);
        std::cerr<<"File Error: "<<foo<<std::endl;
        return 1;
      }
    if (nrank == 0) {
        logfile << "Done saving mesh file: " << outfile << std::endl;
      }
#endif
    return iBase_SUCCESS;

  }

  int CoreGen::save_mesh(int nrank) {
    // export proc- nrank mesh
    std::ostringstream os;
    std::string fname;
    fname = iname;
    os << fname << nrank << ".h5m";
    fname = os.str();
    iMesh_save(imesh->instance(), root_set, fname.c_str(), NULL, &err, strlen(fname.c_str()), 0);
    ERRORR("Trouble writing output mesh.", err);
    logfile << "Saved mesh file: " << fname.c_str() << std::endl;

    return iBase_SUCCESS;
  }

  int CoreGen::save_mesh() {
    // ---------------------------------------------------------------------------
    // Function: save mesh serially
    // Input:    none
    // Output:   none
    // ---------------------------------------------------------------------------
    // export
    logfile << "Saving mesh file." << std::endl;
    iMesh_save(imesh->instance(), root_set, outfile.c_str(), NULL, &err, strlen(
                 outfile.c_str()), 0);
    ERRORR("Trouble writing output mesh.", err);
    logfile << "Saved mesh file: " << outfile.c_str() << std::endl;

    return iBase_SUCCESS;
  }


  int CoreGen::save_geometry() {
    // ---------------------------------------------------------------------------
    // Function: save geometry serially
    // Input:    none
    // Output:   none
    // ---------------------------------------------------------------------------
/*    double dTol = 1e-3;

    // getting all entities for merge and imprint
    SimpleArray<iBase_EntityHandle> entities_merge, entities_imprint;
    iGeom_getEntities(igeom->instance(), root_set, iBase_REGION,
                      ARRAY_INOUT(entities_merge), &err );
    ERRORR("Trouble writing output geometry.", err);

    // merge and imprint before save
    logfile << "Merging.." << std::endl;
    iGeom_mergeEnts(igeom->instance(), ARRAY_IN(entities_merge), dTol, &err);
    ERRORR("Trouble writing output geometry.", err);

    iGeom_getEntities( igeom->instance(), root_set, iBase_REGION, ARRAY_INOUT(entities_imprint),&err );
    ERRORR("Trouble writing output geometry.", err);

    // logfile << "Imprinting.." << std::endl;
    // iGeom_imprintEnts(igeom->instance(), ARRAY_IN(entities_imprint),&err);
    // ERRORR("Trouble writing output geometry.", err);
    // export
    logfile << "Saving geometry file: " <<  outfile << std::endl;
*/
    iGeom_save(igeom->instance(), outfile.c_str(), NULL, &err,
               strlen(outfile.c_str()), 0);
    ERRORR("Trouble writing output geometry.", err);
    logfile << "Saved geometry file: "<< outfile.c_str() <<std::endl;

    return iBase_SUCCESS;
  }

  int CoreGen::distribute_mesh(const int nrank, int numprocs)
  // -------------------------------------------------------------------------------------------
  // Function: merge the nodes within a set tolerance in the model
  // Input:    none
  // Output:   none
  // -------------------------------------------------------------------------------------------
  {
    int nback = files.size() - nassys;
    if(nrank < ((int) core_alias.size() + nback)){
        if(numprocs > (int) core_alias.size()){
            numprocs =  core_alias.size() + nback;
          }
#ifdef USE_MPI
        std::vector<int> rank_load;
        rank_load.resize(numprocs);
        int extra_procs = numprocs - files.size();
        if(numprocs >= (int) files.size() && numprocs <= (tot_assys + nback)){
            // again fill assm_meshfiles
            for(int p=0; p<nassys; p++)
              assm_meshfiles[p]=0;
            for(int p=0; p<tot_assys; p++){
                for(int q=0; q<nassys; q++){
                    if(strcmp(core_alias[p].c_str(), assm_alias[q].c_str()) ==0) {
                        assm_meshfiles[q]+=1;
                      }
                  }
              }
            //distribute
            for(int i=0; i<  (int)files.size(); i++){
                rank_load[i] = i;
              }

            int temp = 0;
            int assm_load = - 1;
            int e= 0;
            // compute the rank, mf vector for extra procs
            while(e < extra_procs){
                for(int i = 0; i < nassys; i++){
                    if (assm_meshfiles[i] > temp){
                        temp = assm_meshfiles[i];
                        assm_load = i;
                      }
                    else if (assm_load == -1){
                        logfile << "No assemblies mesh files used in core" << std::endl;
                        exit(0);
                      }
                  }
                assm_meshfiles[assm_load]-=1;
                int temp_rank = files.size()+ e;
                rank_load[temp_rank] = assm_load;
                e++;
                temp = 0;
              }
          }
        else{
            logfile << "Warning: #procs <= #assys in core, some processor will be idle" << std::endl;
          }

        std::vector<int> times_loaded(nassys);
        std::vector<std::vector<int> > meshfiles_rank (files.size());
        for(int i=0; i < (int) files.size(); i++){
            for(int j=0; j < (int) rank_load.size(); j++){
                if(rank_load[j]==i){
                    meshfiles_rank[i].push_back(j);
                    times_loaded[i]+=1;
                  }
              }
          }

        position_core.resize(numprocs);

        for(int i=0; i < (int) files.size(); i++){
            int k = 0;
            if(i < (nassys) ){
                for(int j=0; j < (int) assm_location[i].size(); j++){
                    if (k >= (int) meshfiles_rank[i].size()){
                        k = 0;
                      }
                    int p = meshfiles_rank[i][k];
                    int q = assm_location[i][j];
                    position_core[p].push_back(q);

                    ++k;
                  }
              }
            else{
                // this is background mesh set it -2, no meshfile to copy/move
                position_core[i].push_back(-2);
              }
          }

        if(nrank == 0){
            logfile << " copy/move task distribution " << std::endl;
            for(int i =0; i< numprocs; i++){
                logfile << "rank: " << i <<  " positions : ";
                for(int j=0; j< (int) position_core[i].size(); j++){
                    logfile << (int) position_core[i][j] << " ";
                  }
                logfile << "\n" << std::endl;
              }
          }
#endif
      }
    return 0;
  }

  int CoreGen::load_meshes_parallel(const int nrank, int numprocs)
  // ---------------------------------------------------------------------------
  // Function: loads all the meshes and initializes copymesh object
  // Input:    none
  // Output:   none
  // ---------------------------------------------------------------------------
  {
    int nback = files.size() - nassys;
    cm.resize(files.size());

    iMesh_getRootSet(imesh->instance(), &root_set, &err);
    ERRORR("Couldn't get the root set", err);
    if(nrank < ((int) core_alias.size() + nback)){
        if(numprocs > (int) core_alias.size()){
            numprocs =  core_alias.size() + nback;
          }

#ifdef USE_MPI
        if(numprocs > ((int) core_alias.size() + nback)){
            logfile << "Warning: #procs <= #assys in core, some processor will be idle" << std::endl;
          }

        iBase_EntitySetHandle orig_set;
        int temp_index;
        int extra_procs = numprocs - files.size();
        std::vector<int> rank_load;
        for(int i = 0; i< (int) files.size(); i++){
            temp_index = numprocs*i + nrank;
            if(temp_index >= (int) files.size()){
                if (nrank >= (int) files.size() && nrank <= (tot_assys + nback)){
                    int temp = 1;
                    int p = extra_procs;
                    // compute the rank, mf vector for extra procs
                    while(p !=0){
                        int assm_load = - 1;
                        for(int i = 0; i < nassys; i++){
                            if ((int) assm_meshfiles[i] > temp){
                                temp = assm_meshfiles[i];
                                assm_load = i;
                              }
                          }
                        if (assm_load == -1){
                            continue;
                            logfile << "Warning: #procs <= #assys in core, some processor will be idle" << std::endl;
                          }
                        assm_meshfiles[assm_load]-=1;
                        rank_load.push_back(assm_load);
                        --p;
                        temp = 1;
                      }

                    temp_index = nrank - files.size();
                    iMesh_createEntSet(imesh->instance(), 0, &orig_set, &err);
                    ERRORR("Couldn't create file set.", err);

                    if(!compute_meshtogeom){
                        // load this file
                        iMesh_load(imesh->instance(), orig_set, files[rank_load[temp_index]].c_str(), NULL, &err, strlen(files[rank_load[temp_index]].c_str()), 0);
                        ERRORR("Couldn't read mesh file.", err);
                      }
                    else{
                        load_and_compute_meshtogeom(orig_set, files[rank_load[temp_index]]);
                      }
                    logfile << "Loaded mesh file " << rank_load[temp_index] << " in processor: " << nrank << std::endl;

                    ModelEnt *me;
                    me = NULL;
                    me = new ModelEnt(mk_core(), iBase_EntitySetHandle(0), /*igeom instance*/0, (moab::EntityHandle)orig_set, 0);
                    MEntVector assm_set;
                    assm_set.push_back(me);
                    cm[i] = (CopyMesh*) mk_core()->construct_meshop("CopyMesh", assm_set);
                    cm[i]->set_name("copy_move_mesh");

                    assys.push_back(orig_set);
                    assys_index.push_back(rank_load[temp_index]);
                    break;
                  }
              }
            else{
                iMesh_createEntSet(imesh->instance(), 0, &orig_set, &err);
                ERRORR("Couldn't create file set.", err);

                if(!compute_meshtogeom){
                    // load this file
                    iMesh_load(imesh->instance(), orig_set, files[temp_index].c_str(), NULL, &err, strlen(files[temp_index].c_str()), 0);
                    ERRORR("Couldn't read mesh file.", err);
                  }
                else{
                    load_and_compute_meshtogeom(orig_set, files[temp_index]);
                  }\
                logfile << "Loaded mesh file " << temp_index << " in processor: " << nrank << std::endl;

                ModelEnt *me;
                me = NULL;
                me = new ModelEnt(mk_core(), iBase_EntitySetHandle(0), /*igeom instance*/0, (moab::EntityHandle)orig_set, 0);
                MEntVector assm_set;
                assm_set.push_back(me);
                cm[i] = (CopyMesh*) mk_core()->construct_meshop("CopyMesh", assm_set);
                cm[i]->set_name("copy_move_mesh");

                assys.push_back(orig_set);
                assys_index.push_back(temp_index);
              }
          }
#endif
      }
    return iBase_SUCCESS;
  }

  int CoreGen::load_meshes()
  // ---------------------------------------------------------------------------
  // Function: loads all the meshes and initializes copymesh and merge mesh objects
  // Input:    none
  // Output:   none
  // ---------------------------------------------------------------------------
  {
    // create cm instances for each mesh file
    cm.resize(files.size());
    iBase_EntitySetHandle orig_set;
    // loop reading all mesh files
    for (unsigned int i = 0; i < files.size(); i++) {
        iMesh_createEntSet(imesh->instance(), 0, &orig_set, &err);
        ERRORR("Couldn't create file set.", err);
        logfile << "Loading File: " << files[i].c_str() << std::endl;
        if(!compute_meshtogeom){
            iMesh_load(imesh->instance(), orig_set, files[i].c_str(), NULL, &err, strlen(files[i].c_str()), 0);
            ERRORR("Couldn't read mesh file.", err);
            mk_core()->populate_model_ents(0,0,0);
          }
        else{
            load_and_compute_meshtogeom(orig_set, files[i]);
          }
        ModelEnt *me;
        me = NULL;
        me = new ModelEnt(mk_core(), iBase_EntitySetHandle(0), /*igeom instance*/0, (moab::EntityHandle)orig_set, 0);
        MEntVector assm_set;
        //assm_set.clear();
        assm_set.push_back(me);
        cm[i] = (CopyMesh*) mk_core()->construct_meshop("CopyMesh", assm_set);
        cm[i]->set_name("copy_move_mesh");
        cm[i]->copy_sets().add_set(orig_set);

        //resize this else code will break
        //check if we've loaded the same mesh file and need to shift the Material and Neumann set start id's
        if(bsameas[i] == 0 && (int) i < nassys){
            if(all_ms_starts[i] != -1 && all_ns_starts[i] !=-1){
                if(!shift_mn_ids(orig_set, i))
                  ERRORR("Couldn't shift material and neumann set id's.", 1);
              }
          }

        assys.push_back(orig_set);
        assys_index.push_back(i);
      }
    logfile << "Loaded mesh files." << std::endl;

    return iBase_SUCCESS;
  }

  int CoreGen::load_and_compute_meshtogeom(iBase_EntitySetHandle orig_set, std::string filename)
  // ---------------------------------------------------------------------------
  // Function: loads all the meshes and initializes copymesh and merge mesh objects
  // Input:    none
  // Output:   none
  // ---------------------------------------------------------------------------
  {
    bool isCub = false;
    std::string ext;
    std::string::size_type idx;
    idx = filename.rfind('.');
    if(idx != std::string::npos)
      {
        ext = filename.substr(idx+1);
      }
    else
      {
        // No extension found
      }
    if(ext == "cub") isCub = true;

    if(isCub){
        mk_core()->igeom_instance()->deleteAll();

        iGeom_load(igeom->instance(), filename.c_str(), NULL, &err,
                   strlen(filename.c_str()), 0);
        iMesh_load(imesh->instance(), orig_set, filename.c_str(), NULL, &err, strlen(filename.c_str()), 0);
        ERRORR("Couldn't read mesh file.", err);
        // populate model eneities here so we can call irel function from MeshKit

        mk_core()->irel_pair()->inferAllRelations();
        mk_core()->populate_model_ents(0,0,0);


        // now calculate geometric and mesh volume
        VerdictWrapper vw(mb);
        // Get all the materials in this mesh file
        moab::Tag mattag, gidtag;
        moab::Range matsets;
        mb->tag_get_handle("MATERIAL_SET", 1, MB_TYPE_INTEGER, mattag);
        mb->tag_get_handle("GLOBAL_ID", 1, MB_TYPE_INTEGER, gidtag);
        mb->get_entities_by_type_and_tag((moab::EntityHandle) orig_set, MBENTITYSET, &mattag, 0, 1, matsets );

        double mtot = 0.0, volume = 0.0;
        //loop thru all the materials
        moab::Range::iterator set_it;
        for (set_it = matsets.begin(); set_it != matsets.end(); set_it++)  {
            moab::EntityHandle this_set = *set_it;

            // get the id for this set
            int material_id;
            mb->tag_get_data(mattag, &this_set, 1, &material_id);
            // add up the geometric volume for these geometric volumes
            // get the entity sets there must be one entity set for each volume
            std::vector<moab::EntityHandle> geomsets_for_gid;
            geomsets_for_gid.clear();
            mb->get_entities_by_type(this_set, MBENTITYSET, geomsets_for_gid);

            for(unsigned int volid = 0; volid < geomsets_for_gid.size(); volid++){
                // get the gid of this volume
                int my_geom_id = 0;
                mb->tag_get_data(gidtag, &geomsets_for_gid[volid], 1, &my_geom_id);
                iBase_EntityHandle ent2=NULL;
                mk_core()->irel_pair()->getSetEntRelation((iBase_EntitySetHandle) geomsets_for_gid[volid], 1, ent2);
                double myvol = 0.0;
                if(ent2 != NULL){
                    mk_core()->igeom_instance()->measure(&ent2,1,&myvol);
                    volume+=myvol;
                  }
              }
            if(volume == 0){
                logfile << "Trying to compute MESHTOGEOM, but volume obtained is zero." << std::endl;
                logfile << "Geometry engine must match file type. Aborting." << std::endl;
                exit(0);
              }

            std::vector<moab::EntityHandle> elems;
            mb->get_entities_by_dimension(this_set, 3, elems, true);
            // get all the elements in this material
            double mvol_temp = 0.0;
            for(unsigned int i=0; i<elems.size();i++){
                mvol_temp = 0.0;
                vw.quality_measure(elems[i], moab::MB_VOLUME, mvol_temp);
                mtot+=mvol_temp;
              }

            double meshtogeom = mtot/volume;
            meshtogeom_file << material_id << " "<<  meshtogeom << std::endl;
            logfile << material_id << " " << meshtogeom << std::endl;
            moab::Tag mtog_tag;
            // now set the meshtogeom tag on this set
            mb->tag_get_handle( "MESHTOGEOM", 1, MB_TYPE_DOUBLE,
                                mtog_tag, MB_TAG_SPARSE|MB_TAG_CREAT );
            mb->tag_set_data(mtog_tag, &(*set_it), 1, &meshtogeom);

            elems.clear();
            mtot = 0.0;
            volume = 0.0;
          }
        //
      }
    else{
        logfile << "Trying to load a file that doesn't have geometry: Aborting" << std::endl;
        logfile << "Try loading a .cub file." << std::endl;
        exit(0);
    }
    return iBase_SUCCESS;
  }

  int CoreGen::load_geometries()
  // ---------------------------------------------------------------------------
  // Function: loads all the meshes and initializes copymesh and merge mesh objects
  // Input:    none
  // Output:   none
  // ---------------------------------------------------------------------------
  {
    logfile << "\n--Loading geometry files." << std::endl;

    iGeom_getRootSet(igeom->instance(), &root_set, &err);
    ERRORR("Couldn't get the root set", err);

    // create cm instances for each mesh file
    cg.resize(files.size());

    iBase_EntitySetHandle orig_set, temp_set, temp_set1;

    iGeom_createEntSet(igeom->instance(), 0, &temp_set1, &err);
    ERRORR( "Problem creating entity set.",err );

    MEntVector vols_old;
    // loop reading all igeom->instance() files
    for (unsigned int i = 0; i < files.size(); i++) {
        iGeom_createEntSet(igeom->instance(), 0, &orig_set, &err);
        ERRORR( "Problem creating entity set.", err);

        iGeom_createEntSet(igeom->instance(), 0, &temp_set, &err);
        ERRORR( "Problem creating entity set.",err );

        iGeom_load(igeom->instance(), files[i].c_str(), NULL, &err,
                   strlen(files[i].c_str()), 0);
        ERRORR("Couldn't read geometry file.", err);
        mk_core()->populate_model_ents(0,-1,-1);

        iBase_EntityHandle *entities = NULL;
        int entities_ehsize = 0, entities_ehallocated = 0;

        iGeom_getEntities(igeom->instance(), root_set, iBase_REGION, &entities,
                          &entities_ehsize, &entities_ehallocated, &err);
        ERRORR( "Problem getting entities." , err);

        // add the entity
        for (int j = 0; j < entities_ehsize; j++) {
            iGeom_addEntToSet(igeom->instance(), entities[j], temp_set, &err);
            ERRORR( "Problem adding to set.", err );
          }

        iGeom_subtract(igeom->instance(), temp_set, temp_set1, &orig_set, &err);
        ERRORR( "Unable to subtract entity sets.", err );

        MEntVector vols;
        mk_core()->get_entities_by_dimension(3, vols);
        MEntVector vols1;// = vols - vols_old;
        
        std::set_difference(vols.begin(), vols.end(), vols_old.begin(), vols_old.end(), std::inserter(vols1, vols1.begin()));

        cg[i] = (CopyGeom*) mk_core()->construct_meshop("CopyGeom", vols1);
        cg[i]->set_name("copy_move_geom");

        assys.push_back(orig_set);

        assys_index.push_back(i);

        // store this set for subtraction with next entity set
        temp_set1 = temp_set;
        vols_old = vols;
      }
    logfile << "\n--Loaded geometry files.\n" << std::endl;

    return iBase_SUCCESS;
  }

  int CoreGen::shift_mn_ids(iBase_EntitySetHandle orig_set, int number)
  // ---------------------------------------------------------------------------
  // Function: loads all the meshes and initializes copymesh and merge mesh objects
  // Input:    none
  // Output:   none
  // ---------------------------------------------------------------------------
  {
    std::cout << "Swapping MS and NS ids for " << number << std::endl;
    // get all the material sets in this assembly
    moab::Tag mattag, neutag;
    mb->tag_get_handle( "MATERIAL_SET", 1, MB_TYPE_INTEGER, mattag );
    mb->tag_get_handle( "NEUMANN_SET", 1, MB_TYPE_INTEGER, neutag );

    int rval = 0;
    moab::Range matsets, neusets;

    mb->get_entities_by_type_and_tag( (moab::EntityHandle)orig_set, MBENTITYSET, &mattag, 0, 1, matsets );
    mb->get_entities_by_type_and_tag( (moab::EntityHandle)orig_set, MBENTITYSET, &neutag, 0, 1, neusets );

    int i = 0;
    moab::Range::iterator set_it;
    for (set_it = matsets.begin(); set_it != matsets.end(); set_it++)  {
        moab::EntityHandle this_set = *set_it;

        // get the id for this set
        int set_id;
        rval = mb->tag_get_data(mattag, &this_set, 1, &set_id);
        if(rval != moab::MB_SUCCESS) {
            std::cerr<<"getting tag data failed Code:";
            std::string foo = ""; mb->get_last_error(foo);
            std::cerr<<"File Error: "<<foo<<std::endl;
            return 1;
          }

        // set the new id for this set
        int new_id = all_ms_starts[number] + i;
        rval = mb->tag_set_data(mattag, &this_set, 1, &new_id);
        if(rval != moab::MB_SUCCESS) {
            std::cerr<<"getting tag data failed Code:";
            std::string foo = ""; mb->get_last_error(foo);
            std::cerr<<"File Error: "<<foo<<std::endl;
            return 1;
          }
        ++i;
      }

    int j = 0;
    for (set_it = neusets.begin(); set_it != neusets.end(); set_it++)  {
        moab::EntityHandle this_set = *set_it;

        // get the id for this set
        int set_id;
        rval = mb->tag_get_data(neutag, &this_set, 1, &set_id);
        if(rval != moab::MB_SUCCESS) {
            std::cerr<<"getting tag data failed Code:";
            std::string foo = ""; mb->get_last_error(foo);
            std::cerr<<"File Error: "<<foo<<std::endl;
            return 1;
          }

        // set the new id for this set
        int new_id = all_ns_starts[number] + j;
        rval = mb->tag_set_data(neutag, &this_set, 1, &new_id);
        if(rval != moab::MB_SUCCESS) {
            std::cerr<<"getting tag data failed Code:";
            std::string foo = ""; mb->get_last_error(foo);
            std::cerr<<"File Error: "<<foo<<std::endl;
            return 1;
          }
        ++j;
      }
    return 0;
  }

  int CoreGen::move_verts(iBase_EntitySetHandle set, const double *dx)
  // ---------------------------------------------------------------------------
  // Function: Change the coordinates for moving the assembly to its first loc.
  // Input:    none
  // Output:   none
  // ---------------------------------------------------------------------------
  {

    int verts_ents_alloc = 0, verts_ents_size = 0;
    iBase_EntityHandle *verts_ents = NULL;

    iMesh_getEntities(imesh->instance(), set, iBase_VERTEX, iMesh_ALL_TOPOLOGIES,
                      &verts_ents, &verts_ents_alloc, &verts_ents_size, &err);
    ERRORR("Failed to get any entities from original set.", iBase_FAILURE);

    double *coords = 0;
    int coords_alloc = 0, coords_size = 0;

    iMesh_getVtxArrCoords(imesh->instance(), verts_ents, verts_ents_size, iBase_INTERLEAVED,
                          &coords, &coords_alloc, &coords_size, &err);
    ERRORR("Failed to get vtx coords from set.", iBase_FAILURE);

    for (int i = 0; i < verts_ents_size; i++) {
        coords[3 * i] += dx[0];
        coords[3 * i + 1] += dx[1];
        coords[3 * i + 2] += dx[2];
      }

    iMesh_setVtxArrCoords(imesh->instance(), verts_ents, verts_ents_size, iBase_INTERLEAVED,
                          coords, coords_size, &err);
    ERRORR("Failed to set vtx coords.", iBase_FAILURE);

    return iBase_SUCCESS;
  }

  int CoreGen::move_geoms(iBase_EntitySetHandle set, const double *dx)
  // ---------------------------------------------------------------------------
  // Function: Change the coordinates for moving the assembly to its first loc.
  // Input:    none
  // Output:   none
  // ---------------------------------------------------------------------------
  {
    int entities_ehsize = 0, entities_ehallocated = 0;
    iBase_EntityHandle *entities = NULL;

    iGeom_getEntities(igeom->instance(), set, iBase_ALL_TYPES, &entities, &entities_ehsize,
                      &entities_ehallocated, &err);
    ERRORR("Failed to get entities from set.", iBase_FAILURE);

    for (int i = 0; i < entities_ehsize; i++) {
        iGeom_moveEnt(igeom->instance(), entities[i], dx[0], dx[1], dx[2], &err);
        ERRORR("Failed to move geometries.", iBase_FAILURE);
      }
    return iBase_SUCCESS;
  }

  int CoreGen::extrude() {
    // ---------------------------------------------------------------------------
    // Function: extrude 2D surface core
    // Input:    none
    // Output:   none
    // ---------------------------------------------------------------------------
    // extrude if this is a surface mesh
    if (set_DIM == 2 && extrude_flag == true) { // if surface geometry and extrude
        logfile << "Extruding surface mesh." << std::endl;
        //get entities for extrusion
        iBase_EntityHandle *ents = NULL;
        int ents_alloc = 0, ents_size;
        iMesh_getEntities(imesh->instance(), root_set, iBase_FACE, iMesh_ALL_TOPOLOGIES,
                          &ents, &ents_alloc, &ents_size, &err);
        ERRORR("Trouble getting face mesh.", err);

        // add entities for extrusion to a set
        iBase_EntitySetHandle set;
        iMesh_createEntSet(imesh->instance(), false, &set, &err);
        ERRORR("Trouble getting face mesh.", err);

        iMesh_addEntArrToSet(imesh->instance(), ents, ents_size, set, &err);
        ERRORR("Trouble getting face mesh.", err);

        ModelEnt me(mk_core(), iBase_EntitySetHandle(0), /*igeom instance*/0,
                    (moab::EntityHandle)set);
        MEntVector selection;
        selection.push_back(&me);

        // This tag needs to be set to the newly created extrude sets
        const char *tag_g1 = "GEOM_DIMENSION";
        iBase_TagHandle gtag;
        iMesh_getTagHandle(imesh->instance(), tag_g1, &gtag, &err, 14);
        ERRORR("Trouble getting geom dimension set.", err);

        // This tag needs to be set to the newly created extrude sets
        const char *tag_m1 = "MATERIAL_SET";
        iBase_TagHandle mtag;
        iMesh_getTagHandle(imesh->instance(), tag_m1, &mtag, &err, 12);
        ERRORR("Trouble getting material set.", err);

        // This tag needs to be set to the newly created extrude sets
        const char *tag_n1 = "NEUMANN_SET";
        iBase_TagHandle ntag;
        iMesh_getTagHandle(imesh->instance(), tag_n1, &ntag, &err, 11);
        ERRORR("Trouble getting neumann set.", err);


        double v[] = { 0, 0, z_height };
        int steps = z_divisions;
        em = (ExtrudeMesh*) mk_core()->construct_meshop("ExtrudeMesh", selection);
        em->set_transform(Extrude::Translate(v, steps));
        em->copy_faces(true);


        SimpleArray<iBase_EntitySetHandle> msets;
        iMesh_getEntSetsByTagsRec(imesh->instance(), root_set, &mtag, NULL, 1, 0,
                                  ARRAY_INOUT(msets), &err);
        ERRORR("Trouble getting entity set.", err);

        for (int i = 0; i < msets.size(); i++) {
            em->extrude_sets().add_set((iMesh::EntitySetHandle)msets[i]);
          }

        // some entity tag types are always copy or expand
        em->extrude_sets().add_tag("MATERIAL_SET");

        // run
        em->setup_this();
        em->execute_this();

        iMesh_destroyEntSet(imesh->instance(), set, &err);
        ERRORR("Error in destroying ent set of faces after extrusion is done.", err);

        msets.clear();
        iMesh_getEntSetsByTagsRec(imesh->instance(), root_set, &mtag, NULL, 1, 0,
                                  ARRAY_INOUT(msets), &err);
        ERRORR("Trouble getting entity set.", err);
        // now delete all the 2D material sets:
        for (int i = 0; i < msets.size(); i++) {
            int num =0;
            iMesh_getNumOfType(imesh->instance(), msets[i], iBase_REGION, &num, &err);
            ERRORR("Trouble getting num entities.", err);
            if(num == 0){ // material sets with quads
                iMesh_destroyEntSet(imesh->instance(), msets[i], &err);
                ERRORR("Trouble destroying set.", err);
              }
          }


        // Step 2: get all max. value of neumann sets, then, for newly created NS set a new value and GD =2 tag.

        iBase_EntityHandle *ents1 = NULL;
        int ents_alloc1 = 0, ents_size1 = 0;

        SimpleArray<iBase_EntitySetHandle> nsets;
        iMesh_getEntSetsByTagsRec(imesh->instance(), root_set, &ntag, NULL,
                                  1, 0, ARRAY_INOUT(nsets), &err);
        ERRORR("Trouble getting entity set.", err);

        int max_nset_value = 0;
        for (int i = 0; i < nsets.size(); i++) {
            int nvalue;
            iMesh_getEntSetIntData(imesh->instance(), nsets[i], ntag, &nvalue, &err);
            ERRORR("Trouble getting entity set.", err);
            if (nvalue > max_nset_value)
              max_nset_value = nvalue;
          }

        for (int i = 0; i < nsets.size(); i++) {

            iMesh_getEntities(imesh->instance(), nsets[i],
                              iBase_FACE, iMesh_ALL_TOPOLOGIES,
                              &ents1, &ents_alloc1, &ents_size1, &err);
            ERRORR("Trouble getting face mesh.", err);

            if(ents_size1 > 0) {
                // set GEOM_DIMENSION tag = 2 and renumber the neumann set
                const int gd = 2;
                const int nvalue = max_nset_value + i;

                iMesh_setEntSetIntData(imesh->instance(), nsets[i], ntag, nvalue, &err);
                ERRORR("Trouble getting entity set.", err);

                iMesh_setEntSetIntData(imesh->instance(),nsets[i], gtag, gd, &err);
                ERRORR("Trouble getting entity set.", err);
              }
            ents_alloc1 = 0;
            ents_size1 = 0;
            *ents1 = NULL;
          }
      }
    return iBase_SUCCESS;
  }

  int CoreGen::create_neumannset() {
    // ---------------------------------------------------------------------------
    // Function: create Neumann set on the whole core
    // Input:    none
    // Output:   none
    // ---------------------------------------------------------------------------
    if (nss_flag == true || nsb_flag == true || nst_flag == true || nssall_flag == true) {
        std::cout << "Creating NeumannSet." << std::endl;

        if (extrude_flag == true)
          set_DIM = 3;

        int err = 0, z_flag = 0, i, ents_alloc = 0, ents_size;
        double z1 = 0.0;
        iBase_TagHandle ntag1, gtag1, nametag1;
        iBase_EntityHandle *ents = NULL;
        iBase_EntitySetHandle set = NULL, set_z1 = NULL, set_z2 = NULL;
        std::vector<iBase_EntitySetHandle> set_side;

        //get entities for skinner
        if(set_DIM ==2) { // if surface geometry specified
            iMesh_getEntities(imesh->instance(), root_set,
                              iBase_FACE, iMesh_ALL_TOPOLOGIES,
                              &ents, &ents_alloc, &ents_size, &err);
          }
        else {
            iMesh_getEntities(imesh->instance(), root_set,
                              iBase_REGION, iMesh_ALL_TOPOLOGIES,
                              &ents, &ents_alloc, &ents_size, &err);
          }
        ERRORR("Trouble getting entities for specifying neumannsets via skinner.", err);

        // assign a name to the tag
        const char *ch_name1 = "NAME";
        // get tag handle
        const char *tag_neumann1 = "NEUMANN_SET";
        const char *global_id1 = "GLOBAL_ID";

        iMesh_getTagHandle(imesh->instance(), tag_neumann1, &ntag1, &err, 12);
        ERRORR("Trouble getting handle.", err);

        iMesh_getTagHandle(imesh->instance(), global_id1, &gtag1, &err, 9);
        ERRORR("Trouble getting handle.", err);

        iMesh_getTagHandle(imesh->instance(), ch_name1, &nametag1, &err, 4);
        ERRORR("Trouble getting handle.", err);

        iMesh_createEntSet(imesh->instance(),0, &set, &err); // for all other sides
        ERRORR("Trouble creating set handle.", err);

        if (set_DIM == 3) { // sets for collecting top and bottom surface
            iMesh_createEntSet(imesh->instance(),0, &set_z1, &err);
            ERRORR("Trouble creating set handle.", err);

            iMesh_createEntSet(imesh->instance(),0, &set_z2, &err);
            ERRORR("Trouble creating set handle.", err);

            set_side.resize(num_nsside);
            for(int i=0; i<num_nsside; i++){
                iMesh_createEntSet(imesh->instance(),0, &set_side[i], &err);
                ERRORR("Trouble creating set handle.", err);
              }
          }

        moab::Range tmp_elems;
        tmp_elems.insert((EntityHandle*)ents, (moab::EntityHandle*)ents + ents_size);

        // get the skin of the entities
        moab::Skinner skinner(mb);
        moab::Range skin_range;
        moab::ErrorCode result;
        moab::Range::iterator rit;

        result = skinner.find_skin(0, tmp_elems, set_DIM-1, skin_range);
        if (MB_SUCCESS != result) return result;

        for (rit = skin_range.begin(), i = 0; rit != skin_range.end(); rit++, i++) {
            if(set_DIM == 3) { // filter top and bottom
                int num_vertex=0, size_vertex =0;
                iBase_EntityHandle *vertex = NULL;
                iMesh_getEntAdj(imesh->instance(), (iBase_EntityHandle)(*rit), iBase_VERTEX, &vertex,
                                &num_vertex, &size_vertex, &err);
                ERRORR("Trouble getting number of entities after merge.", err);

                double *coords = NULL;
                int coords_alloc = 0, coords_size=0;
                iMesh_getVtxArrCoords(imesh->instance(), vertex, size_vertex, iBase_INTERLEAVED,
                                      &coords, &coords_alloc, &coords_size, &err);
                ERRORR("Trouble getting number of entities after merge.", err);
                double ztemp = coords[2];
                int flag = 0;
                for (int p=1; p<num_vertex; p++) {
                    double z1 = coords[3*p+2];
                    if( fabs(ztemp-z1) >= merge_tol) {
                        flag = 1;
                        continue;
                      }
                  }
                if(flag == 0) { // this is top or bottom surface
                    if (z_flag == 0) { // store z-coord this is the first call
                        z_flag = 1;
                        z1 = ztemp;
                      }
                    if( fabs(ztemp-z1) <= merge_tol) {
                        iMesh_addEntToSet(imesh->instance(), (iBase_EntityHandle)(*rit), set_z1, &err);
                        ERRORR("Trouble getting number of entities after merge.", err);
                      }
                    else {
                        iMesh_addEntToSet(imesh->instance(), (iBase_EntityHandle)(*rit), set_z2, &err);
                        ERRORR("Trouble getting number of entities after merge.", err);
                      }
                  }
                else if (flag == 1) { // add the faces that are not top or bottom surface
                    // filter the sidesets based on their x and y coords

                    for(int k=0; k<num_nsside; k++){
                        if ( fabs((coords[0])*nsx[k] + (coords[1])*nsy[k] + nsc[k]) <= merge_tol
                             && fabs((coords[3])*nsx[k] + (coords[4])*nsy[k] + nsc[k]) <= merge_tol
                             && fabs((coords[6])*nsx[k] + (coords[7])*nsy[k] + nsc[k]) <= merge_tol) {
                            iMesh_addEntToSet(imesh->instance(), (iBase_EntityHandle)(*rit), (iBase_EntitySetHandle) set_side[k], &err);
                            ERRORR("Trouble getting number of entities after merge.", err);
                            continue;
                          }
                        else{ // outside the specified
                            iMesh_addEntToSet(imesh->instance(), (iBase_EntityHandle)(*rit), set, &err);
                            ERRORR("Trouble getting number of entities after merge.", err);
                            continue;
                          }
                      }
                    if(num_nsside == 0){
                        iMesh_addEntToSet(imesh->instance(), (iBase_EntityHandle)(*rit), set, &err);
                        ERRORR("Trouble getting number of entities after merge.", err);
                      }

                  }
              }
            else if(set_DIM == 2) { // edges add all for sideset
                iMesh_addEntToSet(imesh->instance(), (iBase_EntityHandle)(*rit), set, &err);
                ERRORR("Trouble getting number of entities after merge.", err);
              }
          }

        if (set_DIM == 3) {
            if (nst_flag == true || nsb_flag == true) {

                iMesh_setEntSetIntData( imesh->instance(), set_z1, ntag1, nst_Id, &err);
                ERRORR("Trouble getting handle.", err);

                iMesh_setEntSetIntData( imesh->instance(), set_z1, gtag1, nst_Id, &err);
                ERRORR("Trouble getting handle.", err);

                std::string name1 = "core_top_ss";
                iMesh_setEntSetData( imesh->instance(), set_z1, nametag1, name1.c_str(), 11, &err);
                ERRORR("Trouble getting handle.", err);

                iMesh_setEntSetIntData( imesh->instance(), set_z2, ntag1, nsb_Id, &err);
                ERRORR("Trouble getting handle.", err);

                iMesh_setEntSetIntData( imesh->instance(), set_z2, gtag1, nsb_Id, &err);
                ERRORR("Trouble getting handle.", err);

                std::string name2 = "core_bottom_ss";
                iMesh_setEntSetData( imesh->instance(), set_z2, nametag1, name2.c_str(), 14, &err);
                ERRORR("Trouble getting handle.", err);

                for(int j=0; j<num_nsside; j++){
                    iMesh_setEntSetIntData( imesh->instance(), set_side[j], ntag1, nss_Id[j], &err);
                    ERRORR("Trouble getting handle.", err);

                    iMesh_setEntSetIntData( imesh->instance(), set_side[j], gtag1, nss_Id[j], &err);
                    ERRORR("Trouble getting handle.", err);

                    std::stringstream ss;
                    ss << j;
                    std::string name3 = "side" + ss.str();
                    iMesh_setEntSetData( imesh->instance(), set_side[j], nametag1, name3.c_str(), 8, &err);
                    ERRORR("Trouble getting handle.", err);
                  }
              }
          }
        // same for both 2D and 3D models
        if (nssall_flag == true) {
            iMesh_setEntSetIntData( imesh->instance(), set, ntag1, nssall_Id, &err);
            ERRORR("Trouble getting handle.", err);

            iMesh_setEntSetIntData( imesh->instance(), set, gtag1, nssall_Id, &err);
            ERRORR("Trouble getting handle.", err);

            std::string name = "all_sides";
            iMesh_setEntSetData( imesh->instance(), set, nametag1, name.c_str(), 9, &err);
            ERRORR("Trouble getting handle.", err);
          }
      }
    return iBase_SUCCESS;
  }

  void CoreGen::IOErrorHandler (ErrorStates ECode) const
  // ---------------------------------------------------------------------------
  // Function: displays error messages related to input file parsing data
  // Input:    error code
  // Output:   none
  // ---------------------------------------------------------------------------
  {
    std::cerr << '\n';
    if (ECode == INVALIDINPUT) // invalid input
      std::cerr << "Invalid input.";
    else if (ECode == ENEGATIVE) // invalid input
      std::cerr << "Unexpected negative value.";
    else
      std::cerr << "Unknown error ...?";

    std::cerr << '\n' << "Error reading input file, line : " << linenumber;
    std::cerr << std::endl;
    exit (1);
  }

  int CoreGen::write_minfofile()
  // ---------------------------------------------------------------------------
  // Function: write the spreadsheet mesh info file based on inputs read from mesh & input file
  // Input:    none
  // Output:   none
  // ---------------------------------------------------------------------------
  {
    logfile << "Writing mesh info file indicating elements and pin number" << std::endl;

    moab::Tag ntag;
    mb->tag_get_handle(NAME_TAG_NAME, NAME_TAG_SIZE, MB_TYPE_OPAQUE, ntag);
    moab::Tag mattag;
    mb->tag_get_handle( "MATERIAL_SET", 1, MB_TYPE_INTEGER, mattag );
    int rval = 0;
    moab::Range matsets;
    std::vector <EntityHandle> set_ents;

    mb->get_entities_by_type_and_tag( 0, MBENTITYSET, &mattag, 0, 1, matsets );

    moab::Range::iterator set_it;
    for (set_it = matsets.begin(); set_it != matsets.end(); set_it++)  {
        moab::EntityHandle this_set = *set_it;

        // get the id for this set
        int set_id;
        rval = mb->tag_get_data(mattag, &this_set, 1, &set_id);
        if(rval != moab::MB_SUCCESS) {
            std::cerr<<"getting tag data failed Code:";
            std::string foo = ""; mb->get_last_error(foo);
            std::cerr<<"File Error: "<<foo<<std::endl;
            return 1;
          }
        char name[NAME_TAG_SIZE];
        rval = mb->tag_get_data(ntag, &this_set, 1, &name);
        if(rval != moab::MB_SUCCESS) {
            std::cerr<<"getting tag data failed Code:";
            std::string foo = ""; mb->get_last_error(foo);
            std::cerr<<"File Error: "<<foo<<std::endl;
            return 1;
          }
        // check for the special block _xp created in AssyGen stage
        // now print out elements for each pin on the mesh info file
        if(name[0]=='_' && name[1]=='x' && name[2] == 'p'){

            // get the entities in the set, recursively
            rval = mb->get_entities_by_dimension(this_set, 3, set_ents, true);

            logfile << "Block: " << set_id << " has "
                    << set_ents.size() << " entities. Name = " << name;

            // loop thro' all the elements in all the sets
            for (int i=0;i<int(set_ents.size());i++){

                std::vector<EntityHandle> conn;
                EntityHandle handle = set_ents[i];

                // get the connectivity of this element
                rval = mb->get_connectivity(&handle, 1, conn);
                double coords[3];
                double x_sum = 0.0, y_sum = 0.0, z_sum = 0.0;
                for (int j = 0; j<int(conn.size()); ++j){
                    rval = mb->get_coords(&conn[j], 1, coords);
                    x_sum+=coords[0];
                    y_sum+=coords[1];
                    z_sum+=coords[2];
                  }
                int p = 3;
                while(name[p]!='\0'){
                    minfo_file << name[p];
                    ++p;
                  }
                minfo_file << " \t" << x_sum/conn.size() << " \t" << y_sum/conn.size() << " \t" <<  z_sum/conn.size() <<  std::endl;
              }
            logfile << ". Deleting block " << set_id << std::endl;
            mb->delete_entities(&this_set, 1);
            set_ents.clear();
          }

      }
    return 0;
  }


  int CoreGen::prepareIO(int argc, char *argv[], int nrank, int nprocs, std::string  TestDir)
  // -----------------------------------------------------------------------------------
  // Function: Obtains file names and opens input/output files and then read/write them
  // Input:    command line arguments
  // Output:   none
  // -----------------------------------------------------------------------------------
  {
    // set rank and total number of processors
    rank = nrank;
    procs = nprocs;
    sTime = clock();
    testdir = TestDir;
    if (argc > 1) {
        if (argv[1][0] == '-' && argv[1][1] == 'm') {
            // set to zero, when run_flag = 1, program runs and does copy, move, merge, extrude, assign gids, save and close
            run_flag = 0;
          }
      }

    bool bDone = false;
    do {
        if (2 == argc) {
            if (argv[1][0] == '-' && nrank == 0) {
                if (argv[1][1] == 'h') {
                    logfile << "Usage: coregen [-t -m -h] <coregen input file>"
                            << std::endl;
                    logfile << "        -t print timing and memory usage info in each step"
                            << std::endl;
                    logfile << "        -m create makefile only" << std::endl;
                    logfile << "        -h print help" << std::endl;
                    logfile << "\nInstruction on writing coregen input file can be found at: "
                            << std::endl;
                    logfile << "        http://press3.mcs.anl.gov/sigma/meshkit/rgg/coregen-input-file-keyword-definitions/"
                            << std::endl;
                    exit(0);
                  }
              }

            iname = argv[1];
            ifile = iname + ".inp";
            outfile = iname + ".h5m";
            mfile = iname + ".makefile";
            meshtogeomfile = "meshtogeom." + iname;
            infofile = iname + "_info.csv";
            minfofile = iname + "_mesh_info.csv";
            logfilename = iname + ".log";
          } else if (3 == argc) {
            int i = 1;// will loop through arguments, and process them
            for (i = 1; i < argc - 1; i++) {
                if (argv[i][0] == '-') {
                    switch (argv[i][1]) {
                      case 'm': {
                          if (nrank == 0) {
                              logfile << "Creating Make/Info file Only" << std::endl;
                            }
                          // only makefile creation specified
                          iname = argv[2];
                          ifile = iname + ".inp";
                          outfile = iname + ".h5m";
                          mfile = iname + ".makefile";
                          meshtogeomfile = "meshtogeom." + iname;
                          infofile = iname + "_info.csv";
                          minfofile = iname + "_mesh_info.csv";
                          logfilename = iname + ".log";
                          break;
                        }
                      case 't': {
                          mem_tflag = true;
                          iname = argv[2];
                          ifile = iname + ".inp";
                          outfile = iname + ".h5m";
                          mfile = iname + ".makefile";
                          meshtogeomfile = "meshtogeom." + iname;
                          infofile = iname + "_info.csv";
                          minfofile = iname + "_mesh_info.csv";
                          logfilename = iname + ".log";
                          break;
                        }
                      case 'h': {
                          if (nrank == 0) {
                              logfile << "Usage: coregen [-t -m -h] <coregen input file>"
                                      << std::endl;
                              logfile << "        -t print timing and memory usage info in each step"
                                      << std::endl;
                              logfile << "        -m create makefile only"
                                      << std::endl;
                              logfile << "        -h print help" << std::endl;
                              logfile << "\nInstruction on writing coregen input file can also be found at: "
                                      << std::endl;
                              logfile << "        http://press3.mcs.anl.gov/sigma/meshkit/rgg/coregen-input-file-keyword-definitions/"
                                      << std::endl;
                              exit(0);
                              break;
                            }
                        }
                      }
                  }
              }
          } else { //default case
            if (nrank == 0) {
                logfile << "Usage: " << argv[0]
                        << " <input file> WITHOUT EXTENSION" << std::endl;
              }
            iname = TestDir + "/" + COREGEN_DEFAULT_TEST_FILE;
            ifile = iname + ".inp";
            std::string temp = CTEST_FILE_NAME;
            outfile = temp + ".h5m";
            mfile = temp + ".makefile";
            meshtogeomfile = "meshtogeom." + temp;
            infofile = temp + "_info.csv";
            minfofile = temp + "_mesh_info.csv";
            logfilename = temp + ".log";
          }

        // open the file
        file_input.open(ifile.c_str(), std::ios::in);
        logfile.coss.open(logfilename.c_str(), std::ios::out);

        /*********************************************/
        // Print banner on logfile and std out
        /*********************************************/
        if (rank == 0) {
            banner();
            Timer.GetDateTime(szDateTime);
            logfile << "\nStarting out at : " << szDateTime << "\n";
          }

        if (!file_input) {
            if (nrank == 0) {
                logfile << "Default case input file located here: <MeshKit/data>" << std::endl;
                logfile << "Usage: coregen [-t -m -h] <coregen input file>"
                        << std::endl;
                logfile << "        -t print timing and memory usage info in each step"
                        << std::endl;
                logfile << "        -m create makefile only" << std::endl;
                logfile << "        -h print help" << std::endl;
                logfile << "\nInstruction on writing coregen input file can be found at: "
                        << std::endl;
                logfile << "        http://press3.mcs.anl.gov/sigma/meshkit/rgg/coregen-input-file-keyword-definitions/"
                        << std::endl;
              }
            file_input.clear();
            exit(1);
          } else
          bDone = true; // file opened successfully
      } while (!bDone);

    // open Makefile-rgg
    do {
        make_file.open(mfile.c_str(), std::ios::out);
        if (!make_file) {
            if (nrank == 0) {
                logfile << "Unable to open makefile for writing" << std::endl;
              }
            make_file.clear();
          } else
          bDone = true; // file opened successfully
      } while (!bDone);

    if (nrank == 0) {
        logfile << "\nEntered input file name: " << ifile << std::endl;
      }

    // now call the functions to read and write
    err = read_inputs_phase1(argc, argv);
    ERRORR("Failed to read inputs in phase1.", 1);

    err = read_inputs_phase2(argc, argv);
    ERRORR("Failed to read inputs in phase2.", 1);

    // open meshtogeomfile
    if(compute_meshtogeom){
        meshtogeom_file.coss.open(meshtogeomfile.c_str(), std::ios::out);
        logfile << "Created meshtogeom file: " << meshtogeomfile << std::endl;
    }
    // open info file
    if(strcmp(info.c_str(),"on") == 0 && nrank == 0){
        do {
            info_file.open(infofile.c_str(), std::ios::out);
            if (!info_file) {
                if (nrank == 0) {
                    logfile << "Unable to open makefile for writing" << std::endl;
                  }
                info_file.clear();
              } else
              bDone = true; // file opened successfully
            logfile << "Created core info file: " << infofile << std::endl;
          } while (!bDone);

        info_file << "assm index"  << " \t" << "assm number" << " \t" << "dX" << " \t" << "dY" << " \t" << "dZ"  << " \t" << "rank" << std::endl;
      }

    // open mesh info file
    if(strcmp(minfo.c_str(),"on") == 0 && nrank == 0){
        do {
            minfo_file.open(minfofile.c_str(), std::ios::out);
            if (!info_file) {
                if (nrank == 0) {
                    logfile << "Unable to open minfofile for writing" << std::endl;
                  }
                minfo_file.clear();
              } else
              bDone = true; // file opened successfully
            logfile << "Created mesh details info file: " << minfofile << std::endl;
          } while (!bDone);
        minfo_file << "pin_number"  << " \t" << "x_centroid" << " \t" << "y_centroid" << " \t" << "z_centroid" << std::endl;
      }
    if (nrank == 0) {
        err = write_makefile();
        ERRORR("Failed to write a makefile.", 1);
      }
    return 0;
  }


  int CoreGen::read_inputs_phase1(int argc, char *argv[]) {
    // ---------------------------------------------------------------------------
    // Function: Reads the dimension and symmetry of the problem
    // Input:    none
    // Output:   none
    // ---------------------------------------------------------------------------
    CParser parse;
    for (;;) {
        if (!parse.ReadNextLine(file_input, linenumber, input_string, MAXCHARS,
                                comment))
          ERRORR("Reading input file failed",1);
        //    logfile << input_string << std::endl;
        if (input_string.substr(0, 11) == "problemtype") {
            std::istringstream formatString(input_string);
            formatString >> card >> prob_type;
            if(((strcmp (prob_type.c_str(), "geometry") != 0)
                && (strcmp (prob_type.c_str(), "mesh") != 0)) || formatString.fail())
              IOErrorHandler (INVALIDINPUT);
            if ((strcmp(prob_type.c_str(), "geometry") == 0)) {
                prob_type = "geometry";
              }
          }
        if (input_string.substr(0, 8) == "geometry" && input_string.substr(0, 12) != "geometrytype") {
            std::istringstream formatString(input_string);
            formatString >> card >> geometry;
            if(((strcmp (geometry.c_str(), "volume") != 0)
                && (strcmp (geometry.c_str(), "surface") != 0)) || formatString.fail())
              IOErrorHandler (INVALIDINPUT);
            if ((strcmp(geometry.c_str(), "surface") == 0)) {
                set_DIM = 2;
              }
          }

        // igeom->instance() engine
        if (input_string.substr(0, 10) == "geomengine") {
            std::istringstream formatString(input_string);
            formatString >> card >> geom_engine;
            if(((strcmp (geom_engine.c_str(), "acis") != 0)
                && (strcmp (geom_engine.c_str(), "occ") != 0)) || formatString.fail())
              IOErrorHandler (INVALIDINPUT);
          }

        // symmetry
        if (input_string.substr(0, 8) == "symmetry") {
            std::istringstream formatString(input_string);
            formatString >> card >> symm;
            if((symm !=1 && symm !=6 && symm !=12) || formatString.fail())
              IOErrorHandler (INVALIDINPUT);
          }

        // merge tolerance
        if (input_string.substr(0, 14) == "mergetolerance") {
            std::istringstream formatString(input_string);
            formatString >> card >> merge_tol;
            if(merge_tol < 0 || formatString.fail())
              IOErrorHandler (ENEGATIVE);
          }

        // save onefile for each proc (multiple) flag
        if (input_string.substr(0, 12) == "saveparallel") {
            std::istringstream formatString(input_string);
            formatString >> card >> savefiles;
            if(formatString.fail())
              IOErrorHandler (INVALIDINPUT);
          }
        // info flag
        if (input_string.substr(0, 4) == "info") {
            std::istringstream formatString(input_string);
            formatString >> card >> info;
            if(formatString.fail())
              IOErrorHandler (INVALIDINPUT);
          }
        // info flag
        if (input_string.substr(0, 8) == "meshinfo") {
            std::istringstream formatString(input_string);
            formatString >> card >> minfo;
            if(formatString.fail())
              IOErrorHandler (INVALIDINPUT);
          }
        // neumannset card
        if (input_string.substr(0, 10) == "neumannset") {
            std::istringstream formatString(input_string);
            std::string nsLoc = "", temp1, temp2, temp3;
            double x, y, c;
            int nsId = 0;
            formatString >> card >> nsLoc >> nsId;
            if(nsId < 0 || formatString.fail())
              IOErrorHandler (INVALIDINPUT);
            if ((strcmp(nsLoc.c_str(), "top") == 0)) {
                nst_flag = true;
                nst_Id = nsId;
              } else if ((strcmp(nsLoc.c_str(), "bot") == 0)) {
                nsb_flag = true;
                nsb_Id = nsId;
              } else if ((strcmp(nsLoc.c_str(), "sideall") == 0)) {
                nssall_flag = true;
                nssall_Id = nsId;
              } else if ((strcmp(nsLoc.c_str(), "side") == 0)) {
                nss_Id.push_back(nsId);
                // we are reading the equation of a straight line ax + by + c = 0
                formatString >> temp1 >> x >> temp2 >> y >> temp3 >> c;
                if(formatString.fail())
                  IOErrorHandler (INVALIDINPUT);
                nsx.push_back(x);
                nsy.push_back(y);
                nsc.push_back(c);

                ++num_nsside;
                nss_flag = true;
              } else {
                std::cout << "Invalid Neumann set specification" << std::endl;
              }
          }

        // breaking condition
        if (input_string.substr(0, 3) == "end") {
            std::istringstream formatstring(input_string);
            break;
          }
      }
    return iBase_SUCCESS;
  }

  int CoreGen::read_inputs_phase2(int argc, char *argv[])
  // ---------------------------------------------------------------------------
  // Function: read all the inputs
  // Input:    command line arguments
  // Output:   none
  // ---------------------------------------------------------------------------
  {
    //Rewind the input file
    file_input.clear(std::ios_base::goodbit);
    file_input.seekg(0L, std::ios::beg);
    linenumber = 0;

    CParser parse;
    for (;;) {
        if (!parse.ReadNextLine(file_input, linenumber, input_string, MAXCHARS,
                                comment))
          ERRORR("Reading input file failed",1);

        if (input_string.substr(0, 12) == "geometrytype" ) {
            std::istringstream formatString(input_string);
            formatString >> card >> geom_type;
            if(formatString.fail())
              IOErrorHandler (INVALIDINPUT);

            if(geom_type.substr(0,3) == "hex"){ // for all hex type assemblies read the pitch and mesh files names
                bool reading_assemblies = false;
                do{
                    if (!parse.ReadNextLine(file_input, linenumber, input_string,
                                            MAXCHARS, comment))
                      ERRORR("Reading input file failed",1);
                    std::istringstream formatString(input_string);
                    if (input_string.substr(0, 10) == "assemblies"){
                        formatString >> card >> nassys >> pitch;
                        if(nassys < 0 || formatString.fail())
                          IOErrorHandler (INVALIDINPUT);

                        // reading file and alias names
                        if(!parse_assembly_names(parse, argc, argv))
                          ERRORR("error parsing names of assemblies",1);
                        reading_assemblies = true;
                      }
                  } while (reading_assemblies == false) ;
              }
            // read lattice info for all assemblies, also assemblies for rect assemblies.
            if (geom_type == "hexvertex" && symm == 6) {

                // reading lattice
                bool reading_lattice = false;
                do{
                    if (!parse.ReadNextLine(file_input, linenumber, input_string,
                                            MAXCHARS, comment))
                      ERRORR("Reading input file failed",1);
                    std::istringstream formatString1(input_string);
                    if (input_string.substr(0,7) == "lattice"){
                        reading_lattice = true;
                        formatString1 >> card >> nrings;
                        if(nrings < 0 || formatString1.fail())
                          IOErrorHandler (INVALIDINPUT);
                        if (nrings % 2 == 0)
                          tot_assys = (nrings * (nrings)) / 2;
                        else
                          tot_assys = ((nrings * (nrings - 1)) / 2)
                              + (nrings + 1) / 2;
                      }
                  } while (reading_lattice == false) ;

                // now reading the arrangement
                if (!parse.ReadNextLine(file_input, linenumber, input_string,
                                        MAXCHARS, comment))
                  ERRORR("Reading input file failed",1);
                std::istringstream formatString2(input_string);
                for (int i = 1; i <= tot_assys; i++) {
                    formatString2 >> temp_alias;
                    if(formatString2.fail())
                      IOErrorHandler (INVALIDINPUT);
                    core_alias.push_back(temp_alias);
                  }
              }

            else if (geom_type == "rectangular" && symm == 1) {

                    bool reading_assemblies = false;
                    do{
                        if (!parse.ReadNextLine(file_input, linenumber, input_string,
                                                MAXCHARS, comment))
                          ERRORR("Reading input file failed",1);
                        std::istringstream formatString(input_string);
                        if (input_string.substr(0, 10) == "assemblies"){
                            formatString >> card >> nassys >> pitchx >> pitchy;
                            if(nassys < 0 || pitchx < 0 || pitchy< 0 || formatString.fail())
                              IOErrorHandler (INVALIDINPUT);

                            // reading file and alias names
                            if(!parse_assembly_names(parse, argc, argv))
                              ERRORR("error parsing names of assemblies",1);
                            reading_assemblies = true;
                          }
                      } while (reading_assemblies == false) ;

                    // reading lattice
                    bool reading_lattice = false;
                    do{
                        if (!parse.ReadNextLine(file_input, linenumber, input_string,
                                                MAXCHARS, comment))
                          ERRORR("Reading input file failed",1);
                        std::istringstream formatString1(input_string);
                        if (input_string.substr(0,7) == "lattice"){
                            reading_lattice = true;
                            formatString1 >> card >> nringsx >> nringsy;
                            if(nringsx <= 0 || nringsy <= 0 || formatString1.fail())
                              IOErrorHandler (INVALIDINPUT);
                            tot_assys = nringsx * nringsy;
                          }
                      } while (reading_lattice == false) ;

                // now reading the arrangement
                if (!parse.ReadNextLine(file_input, linenumber, input_string,
                                        MAXCHARS, comment))
                  ERRORR("Reading input file failed",1);
                std::istringstream formatString2(input_string);
                for (int i = 1; i <= tot_assys; i++) {
                    formatString2 >> temp_alias;
                    if(formatString2.fail())
                      IOErrorHandler (INVALIDINPUT);
                    core_alias.push_back(temp_alias);
                  }
              }

            else if (geom_type == "hexflat" && symm == 6) {

                // reading lattice
                bool reading_lattice = false;
                do{
                    if (!parse.ReadNextLine(file_input, linenumber, input_string,
                                            MAXCHARS, comment))
                      ERRORR("Reading input file failed",1);
                    std::istringstream formatString1(input_string);
                    if (input_string.substr(0,7) == "lattice"){
                        reading_lattice = true;
                        formatString1 >> card >> nrings;
                        if(nrings < 0 || formatString1.fail())
                          IOErrorHandler (INVALIDINPUT);
                        tot_assys = (nrings * (nrings + 1)) / 2;
                      }
                  } while (reading_lattice == false) ;


                // now reading the arrangement
                if (!parse.ReadNextLine(file_input, linenumber, input_string,
                                        MAXCHARS, comment))
                  ERRORR("Reading input file failed",1);
                std::istringstream formatString2(input_string);
                for (int i = 1; i <= tot_assys; i++) {
                    formatString2 >> temp_alias;
                    if(formatString2.fail())
                      IOErrorHandler (INVALIDINPUT);
                    core_alias.push_back(temp_alias);
                  }
              }

            else if (geom_type == "hexflat" && symm == 1) {
                // reading lattice
                bool reading_lattice = false;
                do{
                    if (!parse.ReadNextLine(file_input, linenumber, input_string,
                                            MAXCHARS, comment))
                      ERRORR("Reading input file failed",1);
                    std::istringstream formatString1(input_string);
                    if (input_string.substr(0,7) == "lattice"){
                        reading_lattice = true;
                        formatString1 >> card >> nrings;
                        if(nrings < 0 || formatString1.fail())
                          IOErrorHandler (INVALIDINPUT);
                        tot_assys = 3 * (nrings * (nrings - 1)) + 1;
                      }
                  } while (reading_lattice == false) ;

                // now reading the arrangement
                if (!parse.ReadNextLine(file_input, linenumber, input_string,
                                        MAXCHARS, comment))
                  ERRORR("Reading input file failed",1);
                std::istringstream formatString2(input_string);
                for (int i = 1; i <= tot_assys; i++) {
                    formatString2 >> temp_alias;
                    if(formatString2.fail())
                      IOErrorHandler (INVALIDINPUT);
                    core_alias.push_back(temp_alias);
                  }
              }

            else if (geom_type == "hexflat" && symm == 12) {
                // reading lattice
                bool reading_lattice = false;
                do{
                    if (!parse.ReadNextLine(file_input, linenumber, input_string,
                                            MAXCHARS, comment))
                      ERRORR("Reading input file failed",1);
                    std::istringstream formatString1(input_string);
                    if (input_string.substr(0,7) == "lattice"){
                        reading_lattice = true;
                        formatString1 >> card >> nrings;
                        if(nrings < 0 || formatString.fail())
                          IOErrorHandler (INVALIDINPUT);
                        if (nrings % 2 == 0)
                          tot_assys = (nrings * (nrings + 2)) / 4;
                        else
                          tot_assys = ((nrings + 1) * (nrings + 1)) / 4;
                      }
                  } while (reading_lattice == false) ;

                // now reading the arrangement
                if (!parse.ReadNextLine(file_input, linenumber, input_string,
                                        MAXCHARS, comment))
                  ERRORR("Reading input file failed",1);
                std::istringstream formatString2(input_string);
                for (int i = 1; i <= tot_assys; i++) {
                    formatString2 >> temp_alias;
                    if(formatString2.fail())
                      IOErrorHandler (INVALIDINPUT);
                    core_alias.push_back(temp_alias);
                  }
              }

            else {
                ERRORR("Invalid geometry type",1);
              }
          }
        // background mesh
        if (input_string.substr(0, 10) == "background") {
            std::istringstream formatString(input_string);
            formatString >> card >> back_meshfile;
            if(formatString.fail())
              IOErrorHandler (INVALIDINPUT);

            all_meshfiles.push_back(back_meshfile);

            if (iname == COREGEN_DEFAULT_TEST_FILE){
                back_meshfile = testdir + back_meshfile;
              }
            files.push_back(back_meshfile);
            back_mesh = true;
          }
        // z-height and z-divisions
        if (input_string.substr(0, 7) == "extrude") {
            extrude_flag = true;
            std::istringstream formatString(input_string);
            formatString >> card >> z_height >> z_divisions;
            if(z_divisions < 0 || formatString.fail())
              IOErrorHandler (INVALIDINPUT);
          }
        // if keyword MESHTOGEOM is specified set tags on material sets and write a meshtogeom.txt file with #material id and #meshtogeom ratio
        if (input_string.substr(0, 10) == "meshtogeom") {
            compute_meshtogeom = true;
          }
        // OutputFileName
        if (input_string.substr(0, 14) == "outputfilename") {
            std::istringstream formatString(input_string);
            formatString >> card >> outfile;
            if(formatString.fail())
              IOErrorHandler (INVALIDINPUT);
          }

        // breaking condition
        if (input_string.substr(0, 3) == "end") {
            std::istringstream formatstring(input_string);
            break;
          }
      }
    // set some variables
    assm_meshfiles.resize(nassys);
    load_per_assm.resize(nassys);
    size_mf.resize(nassys);
    times_loaded.resize(nassys);
    assm_location.resize(nassys);
    bsameas.resize(nassys);

    for(int i = 0; i < tot_assys; i++){
        for (int j = 0; j < nassys; j++){
            if (strcmp(core_alias[i].c_str(), assm_alias[j].c_str()) == 0) {
                assm_meshfiles[j]+=1;
                assm_location[j].push_back(i);
                break;
              }
          }
      }
    return iBase_SUCCESS;
  }

  int CoreGen::parse_assembly_names(CParser parse, int argc, char *argv[] )
  // ---------------------------------------------------------------------------
  // Function: Reads all the assemblies from CoreGen input file
  // Input:    none
  // Output:   none
  // ---------------------------------------------------------------------------
  {

    // reading file and alias names
    for (int i = 1; i <= nassys; i++) {
        if (!parse.ReadNextLine(file_input, linenumber,
                                input_string, MAXCHARS, comment, false))
          ERRORR("Reading input file failed",1);
        std::istringstream formatString(input_string);
        formatString >> meshfile >> mf_alias >> same_as >> reloading_mf >> ms_startid >> ns_startid;
        // we don't check for formatting since same_as and parameters after it may not be present.
        // variable gets populated correctly in the file

        // if meshfile variable is a path then only convert the filename to lower case
        unsigned pos = 0;
        pos = meshfile.find_last_of("/\\");
        if (pos > 0 && pos < meshfile.length()){
            std::string filename = "", path="";
            path = meshfile.substr(0,pos);
            filename = meshfile.substr(pos+1);
            std::transform(filename.begin(), filename.end(), filename.begin(), ::tolower);
            meshfile = path+"/"+ filename;
          }
        else{
            std::transform(meshfile.begin(), meshfile.end(), meshfile.begin(), ::tolower);
          }
        // also convert the alias and reloading_mf name to lower case, since we've changed the actual reading.
        std::transform(mf_alias.begin(), mf_alias.end(), mf_alias.begin(), ::tolower);
        std::transform(reloading_mf.begin(), reloading_mf.end(), reloading_mf.begin(), ::tolower);

        if (same_as == "same_as")
          bsameas.push_back(0);
        else
          bsameas.push_back(1);
        if(bsameas[i-1] == 1){
            all_meshfiles.push_back(meshfile);
            if (argc == 1){
                meshfile = testdir + "/" + meshfile;
              }
            files.push_back(meshfile);
            assm_alias.push_back(mf_alias);
            all_ms_starts.push_back(-1);
            all_ns_starts.push_back(-1);
          }
        else{
            all_meshfiles.push_back(reloading_mf);
            if (iname == COREGEN_DEFAULT_TEST_FILE){
                meshfile = testdir + reloading_mf;
              }
            files.push_back(reloading_mf);
            assm_alias.push_back(mf_alias);
            all_ms_starts.push_back(ms_startid);
            all_ns_starts.push_back(ns_startid);
          }
        //  bsameas = false;
      }
    return 0;
  }

  int CoreGen::find_assm(const int i, int &assm_index)
  // ---------------------------------------------------------------------------
  // Function: find the assembly index (0 to n) for n assemblies for core alias i
  // Input:    none
  // Output:   none
  // ---------------------------------------------------------------------------
  {
    int flag = 0;
    for (int j = 0; j < nassys; j++)
      if (strcmp(core_alias[i].c_str(), assm_alias[j].c_str()) == 0) {
          assm_index = j;
          flag = 1;
          break;
        }
    if (flag == 0)//nothing found return -1 or no assembly
      assm_index = -1;
    return iBase_SUCCESS;
  }

  int CoreGen::write_makefile()
  // ---------------------------------------------------------------------------
  // Function: write the makefile based on inputs read from input file
  // Input:    none
  // Output:   none
  // ---------------------------------------------------------------------------
  {
    std::string name;
    std::vector<std::string> f_no_ext, f_sat, f_inp, f_jou, f_injou;
    make_file << "##" << std::endl;
    make_file
        << "## This makefile is automatically generated by coregen program"
        << std::endl;
    make_file << "##" << std::endl;
    make_file << "## Check your coregen, assygen and cubit location"
              << std::endl;
    make_file << "##" << std::endl;
    make_file << "\nCUBIT = cubit\n" << std::endl;
    make_file << "COREGEN = ../../coregen\n" << std::endl;
    make_file << "ASSYGEN = ../../assygen\n" << std::endl;

    // remove the ./ if run from the current working TestDirectory
    make_file << "MESH_FILES = ";
    std::string filename;
    for (unsigned int i = 0; i < files.size(); i++) {
        if (files[i][0] == '.' && files[i][1] == '/') {
            filename = files[i].substr(2, files[i].length());
          } else if (files[i][0] != '.' || files[i][1] != '/') {
            int loc1 = files[i].find_last_of(".");
            int loc2 = files[i].find_last_of("/");
            filename = files[i].substr(loc2 + 1, loc1);
          } else {
            filename = files[i];
          }
        mk_files.push_back(filename);
        make_file << all_meshfiles[i] << "  ";
      }

    // get file names without extension
    for (unsigned int i = 0; i < mk_files.size(); i++) {
        int loc = mk_files[i].find_first_of(".");
        f_no_ext.push_back(mk_files[i].substr(0, loc));
      }

    make_file << "\n\nGEOM_FILES = ";
    for (unsigned int i = 0; i < mk_files.size(); i++) {
        if (geom_engine == "occ")
          name = f_no_ext[i] + ".stp";
        else
          name = f_no_ext[i] + ".sat";
        f_sat.push_back(name);
        make_file << name << "  ";
        name = "";
      }

    make_file << "\n\nJOU_FILES = ";
    for (unsigned int i = 0; i < mk_files.size(); i++) {
        name = f_no_ext[i] + ".jou";
        f_jou.push_back(name);
        make_file << name << "  ";
        name = "";
      }

    make_file << "\n\nINJOU_FILES = ";
    for (unsigned int i = 0; i < mk_files.size(); i++) {
        name = f_no_ext[i] + ".template.jou";
        f_injou.push_back(name);
        make_file << name << "  ";
        name = "";
      }

    make_file << "\n\nASSYGEN_FILES = ";
    for (unsigned int i = 0; i < mk_files.size(); i++) {
        name = f_no_ext[i] + ".inp";
        f_inp.push_back(name);
        make_file << name << "  ";
        name = "";
      }

    make_file << "\n\n" << outfile << " : ${MESH_FILES} " << ifile << std::endl;
    make_file << "\t" << "${COREGEN} " << iname << std::endl;
    for (unsigned int i = 0; i < mk_files.size(); i++) {
        make_file << all_meshfiles[i] << " : " << f_sat[i] << "  " << f_jou[i]
                     << "  " << f_injou[i] << std::endl;
        make_file << "\t" << "${CUBIT} -batch " << f_jou[i] << "\n"
                  << std::endl;

        make_file << f_sat[i] << " " << f_jou[i] << " " << f_injou[i] << " : "
                  << f_inp[i] << std::endl;
        make_file << "\t" << "${ASSYGEN} " << f_no_ext[i] << "\n" << std::endl;
      }

    make_file.close();
    logfile << "Created makefile: " << mfile << std::endl;
    return 0;
  }

  void CoreGen::banner()
  // ---------------------------------------------------------------------------
  // Function: display the program banner
  // Input:    none
  // Output:   none
  // ---------------------------------------------------------------------------
  {
    if(rank == 0 ){
        logfile << '\n';
        logfile
            << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
            << '\n';
        logfile
            << "Program to Assemble Nuclear Reactor Assembly Meshes and Form a Core     "
            << '\n';
        logfile << "\t\t\tArgonne National Laboratory" << '\n';
        logfile << "\t\t\t        2015         " << '\n';
        logfile
            << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
            << '\n';
      }
  }

} // namespace MeshKit
