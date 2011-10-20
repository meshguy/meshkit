/*********************************************
 June,10
 Reactor Mesh Assembler
 Argonne National Laboratory

 Driver program to generate core mesh by
 copy/move/merge'ing assembly meshes
 as specified in the input file
*********************************************/
#include "crgen.hpp"
#include "clock.hpp"

int main(int argc, char *argv[]) {

  int rank = 0, nprocs = 1;

  //Initialize MPI
#ifdef USE_MPI
  MPI::Init(argc, argv);
  nprocs = MPI::COMM_WORLD.Get_size();
  rank = MPI::COMM_WORLD.Get_rank();
#endif

  // start program timer and declare timing variables
  CClock Timer;
  std::string szDateTime;
  clock_t sTime = clock();
  double ctload = 0, ctcopymove = 0, ctmerge = 0, ctextrude = 0, ctns = 0, ctgid = 0, ctsave = 0;
  clock_t tload = 0, tcopymove = 0, tmerge = 0, textrude = 0, tns = 0, tgid = 0, tsave = 0;

  CCrgen TheCore;
  int err = 0;
  int run_flag = 1;

  // more memory/time related variables
  int ld_t = 0, ld_tload = 0, ld_tcopymove = 0, ld_tsave = 0, ld_tgid = 0, ld_tmerge = 0, ld_tns = 0;
  unsigned long mem1 = 0, mem2 = 0, mem3 = 0, mem4 = 0, mem5 = 0, mem6 = 0, mem7 = 0;

  /*********************************************/
  // Print banner on standard output
  /*********************************************/
  if (rank == 0) {
    err = TheCore.banner();
    ERRORR("Failed in creating banner", 1);

    Timer.GetDateTime(szDateTime);
    std::cout << "\nStarting out at : " << szDateTime << "\n";
  }

  /***********************************************************/
  // read inputs from input file; do this on all processors
  /***********************************************************/
  err = TheCore.prepareIO(argc, argv, rank, nprocs);
  ERRORR("Failed in preparing i/o.", 1);

  if (argc > 1) {
    if (argv[1][0] == '-' && argv[1][1] == 'm') {
      // when run_flag = 1, program runs and does copy, move, merge, extrude, assign gids, save and close
      run_flag = 0;
      // when run_flag = 1, program only generates a makefile
    }
  }

  if (run_flag == 1 ) {   
    /*********************************************/
    // load mesh files
    /*********************************************/
    CClock ld_time;
    if (TheCore.prob_type == "mesh") {
      if (nprocs == 1) {
  	err = TheCore.load_meshes();
  	ERRORR("Failed to load meshes.", 1);
      } 
      else {
#ifdef USE_MPI	
	err = TheCore.load_meshes_parallel(rank, nprocs);
	ERRORR("Failed to load meshes.", 1);

	if(nprocs > (int) TheCore.files.size()){
	  // if there are more procs than files distribute the copy/move work on each proc
	  err = TheCore.distribute_mesh(rank, nprocs);
	  ERRORR("Failed to load meshes.", 1);
	}
  	//Get a pcomm object
	TheCore.pc = new moab::ParallelComm(TheCore.mbImpl());
#endif
      }

      TheCore.mbImpl()->estimated_memory_use(0, 0, 0, &mem1);
      ld_tload = ld_time.DiffTime();
      tload = clock();
      ctload = (double) (tload - sTime)/(60*CLOCKS_PER_SEC);

      if (TheCore.mem_tflag == true && nprocs == 1) {
        std::cout << "\n" << " Clock time taken to load mesh files = " << ld_tload
		  << " seconds" << std::endl;
	std::cout << " CPU time = " << ctload << " mins" << std::endl;
        std::cout << " Memory used: " << mem1/1e6 << " Mb\n" << std::endl;
      }
    }
    /*********************************************/
    // load geometry files
    /*********************************************/
    else if (TheCore.prob_type == "geometry" && nprocs == 1) {
      err = TheCore.load_geometries();
      ERRORR("Failed to load geometries", 1);
    }
    else if(TheCore.prob_type == "geometry" && nprocs > 1){
      std::cout << " Parallel mode not supported for problem-type: Geometry " << std::endl;
      exit(1);
    }

#ifdef USE_MPI
    MPI::COMM_WORLD.Barrier();
#endif
    /*********************************************/
    // copy move
    /*********************************************/
    CClock ld_cm;
    if(TheCore.prob_type == "mesh"){ 
      err = TheCore.copymove_parallel(rank, nprocs);
      ERRORR("Failed in copy move routine.", 1);

      TheCore.mbImpl()->estimated_memory_use(0, 0, 0, &mem2);
      ld_tcopymove = ld_cm.DiffTime();
      tcopymove = clock();
      ctcopymove = (double) (tcopymove - tload)/(60*CLOCKS_PER_SEC);

      if (TheCore.mem_tflag == true && (strcmp(TheCore.prob_type.c_str(), "mesh") == 0) && nprocs == 1) {
	std::cout << "\n" << " Clock time taken to copy/move mesh files = " << ld_tcopymove
		  << " seconds" << std::endl;
	std::cout << " CPU time = " << ctcopymove << " mins" << std::endl;
	std::cout << " Memory used: " << mem2/1e6 << " Mb\n" << std::endl;
      }
    }
    else{
      err = TheCore.copy_move();
    }

#ifdef USE_MPI
    MPI::COMM_WORLD.Barrier();
#endif
    if (TheCore.prob_type == "mesh") {
      /*********************************************/
      // merge
      /*********************************************/
      CClock ld_mm;
      if (nprocs == 1) {
	std::cout << "Merging.." << std::endl;
    	err = TheCore.merge_nodes();
    	ERRORR("Failed to merge nodes.", 1);
      } else {
	err = TheCore.merge_nodes_parallel(rank, nprocs);
    	ERRORR("Failed to merge nodes in parallel.", 1);
      }

      TheCore.mbImpl()->estimated_memory_use(0, 0, 0, &mem3);
      ld_tmerge = ld_mm.DiffTime();
      tmerge = clock();
      ctmerge = (double) (tmerge - tcopymove)/(60*CLOCKS_PER_SEC);

      if (TheCore.mem_tflag == true && nprocs == 1 ) {
    	std::cout << "\n" << " Clock time taken to merge nodes = " << ld_tmerge
    		  << " seconds" << std::endl;
	std::cout << " CPU time = " << ctmerge << " mins" << std::endl;
    	std::cout << " Memory used: " << mem3/1e6 << " Mb\n" << std::endl;
      }
#ifdef USE_MPI
      MPI::COMM_WORLD.Barrier();
#endif
      /*********************************************/
      // extrude
      /*********************************************/
      if(nprocs == 1){
	if (TheCore.extrude_flag == true) {

	  // assign global ids after copy/move step
	  if (rank == 0)
	    err = TheCore.assign_gids();

	  CClock ld_em;
	  err = TheCore.extrude();
	  ERRORR("Failed to extrude.", 1);

	  TheCore.mbImpl()->estimated_memory_use(0, 0, 0, &mem4);
	  ld_t = ld_em.DiffTime();
	  textrude = clock();
	  ctextrude = (double) (textrude - tmerge)/(60*CLOCKS_PER_SEC);

	  if (TheCore.mem_tflag == true && nprocs == 1) {
	    std::cout << "\n" << " Clock time taken to extrude = " << ld_t
		      << " seconds" << std::endl;
	    std::cout << " CPU time = " << ctextrude << " mins" << std::endl;
	    std::cout << " Memory used: " << mem4/1e6 << " Mb\n"
		      << std::endl;
	  }
	}
      }
#ifdef USE_MPI
      MPI::COMM_WORLD.Barrier();
#endif
      /*********************************************/
      // assign gids
      /*********************************************/
      CClock ld_gid;
      if (nprocs == 1) {
    	err = TheCore.assign_gids();
    	ERRORR("Failed to assign global ids.", 1);
      }
      else{
	// err = TheCore.assign_gids_parallel(rank, nprocs);
    	// ERRORR("Failed to assign global ids.", 1);
      }

      TheCore.mbImpl()->estimated_memory_use(0, 0, 0, &mem5);
      ld_tgid = ld_gid.DiffTime();
      tgid = clock();
      ctgid = (double) (tgid-tmerge)/(60*CLOCKS_PER_SEC);

      if (TheCore.mem_tflag == true && nprocs == 1) {
      	std::cout << "\n" << " Clock time taken to assign gids = " << ld_tgid
      		  << " seconds" << std::endl;
      	std::cout << " CPU time = " << ctgid << " mins" << std::endl;
      	std::cout << " Memory used: " << mem5/1e6 << " Mb\n" << std::endl;
      }
      /*********************************************/
      // create neumann sets on the core model
      /*********************************************/
      if((TheCore.nss_flag == true || TheCore.nsb_flag == true 
	  || TheCore.nst_flag == true) && nprocs == 1){
	CClock ld_ns;
	err = TheCore.create_neumannset();
	ERRORR("Failed to create neumann set.", 1);

	TheCore.mbImpl()->estimated_memory_use(0, 0, 0, &mem6);
	ld_tns = ld_ns.DiffTime();
	tns = clock();
	ctns = (double) (tns-tgid)/(60*CLOCKS_PER_SEC);
	if (TheCore.mem_tflag == true && nprocs == 1) {
	  std::cout << "\n" << " Clock time taken to create neumann sets = " << ld_tns
		    << " seconds" << std::endl;
	  std::cout << " CPU time = " << ctns << " mins" << std::endl;
	  std::cout << " Memory used: " << mem6/1e6 << " Mb\n" << std::endl;
	}
      }
#ifdef USE_MPI
      MPI::COMM_WORLD.Barrier();
#endif
      /*********************************************/
      // save
      /*********************************************/
      CClock ld_sv;
      if (nprocs == 1) {
      	err = TheCore.save_mesh();
      	ERRORR("Failed to save o/p file.", 1);
      } else {
	if(TheCore.savefiles != "one" && (TheCore.savefiles == "multiple" || TheCore.savefiles == "both")){
	  err = TheCore.save_mesh(rank); // uncomment to save the meshes with each proc
	  ERRORR("Failed to save o/p file.", 1);
	}
	if(TheCore.savefiles != "multiple"){
#ifdef USE_MPI
	  double write_time = MPI_Wtime();
	  err = TheCore.save_mesh_parallel(rank, nprocs);
	  ERRORR("Failed to save o/p file.", 1);
	  write_time = MPI_Wtime() - write_time;
	  if (rank == 0)
	    std::cout << "Parallel write time = " << write_time/60.0 << " mins" << std::endl;
#endif	 
	}
      }

      TheCore.mbImpl()->estimated_memory_use(0, 0, 0, &mem7);
      ld_tsave = ld_sv.DiffTime();
      tsave = clock();
      ctsave = (double) (tsave - tgid)/(60*CLOCKS_PER_SEC);

      if (TheCore.mem_tflag == true && nprocs == 1 ) {
	std::cout << "\n" << " Clock time taken to save = " << ld_tsave << " seconds"
		  << std::endl;
	std::cout << " CPU time = " << ctsave << " mins" << std::endl;
	std::cout << " Memory used: " << mem7/1e6 << " Mb\n" << std::endl;
      }
    }
    /*********************************************/
    // geometry operations
    /*********************************************/
    else if (TheCore.prob_type == "geometry") {
      err = TheCore.save_geometry();
      ERRORR("Failed to save o/p file.", 1);
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
    if (TheCore.mem_tflag == true) {

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

      if(rank == 0 && nprocs > 1){
	std::cout << "\nMAXIMUM TIME TAKEN OVER ALL PROCS\nCLOCK TIME:-";
	std::cout << "\n**r = " << rank<< " Time taken to load mesh files = " << max_tload
		  << " secs" << std::endl;
	std::cout << "***r = : " << rank<< " Memory used: " << max_mem1 << " Mb" << std::endl;
    
	// copymove
	std::cout << "\n**r = " << rank<< " Time taken to copy/move mesh files = " << max_tcopymove
		  << " secs" << std::endl;
	std::cout << "***r = " << rank<< " Memory used: " << max_mem2 << " Mb" << std::endl;
    
	// merge
	std::cout << "\n**r = " << rank<< " Time taken to merge nodes = " << max_tmerge
		  << " secs" << std::endl;
	std::cout << "***r = " << rank<< " Memory used: " << max_mem3 << " kb" << std::endl;

	// assign gid
	std::cout << "\n**r = " << rank<< " Time taken to assign gids = " << max_tgid
		  << " secs" << std::endl;
	std::cout << "*** r = " << rank<< " Memory used: " << max_mem5 << " Mb" << std::endl;
    
	// save
	std::cout << "\n**r = " << rank<< " Time taken to save = " <<  max_tsave << " secs"
		  << std::endl;
	std::cout << "***r = " << rank<< " Memory used: " << max_mem7 << " Mb" << std::endl;

	// cpu times
	std::cout << "\n CPU TIME:-\n" << " r = " << rank<< " Time taken to load mesh files = " << ctload
		  << " mins" << std::endl;

	std::cout << " r = " << rank << " Time taken to copy/move files = " << ctcopymove
		  << " mins" << std::endl;

	std::cout << " r = " << rank << " Time taken to merge = " << ctmerge
		  << " mins" << std::endl;

	std::cout << " r = " << rank <<  " Time taken to assign gids = " << ctgid
		  << " mins" << std::endl;

	std::cout  << " r = " << rank << " Time taken to save mesh = " << ctsave
		   << " mins" << std::endl;
      }
    }
#endif

    if (rank == 0) {
      Timer.GetDateTime(szDateTime);
      std::cout << "\nXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"  << std::endl;
      std::cout << "Ending at : " << szDateTime;
      std::cout << "Elapsed wall clock time: " << Timer.DiffTime()
    		<< " seconds or " << (Timer.DiffTime()) / 60.0 << " mins\n";

      std::cout << "Total CPU time used: " <<  (double) (clock() - sTime)/(CLOCKS_PER_SEC) << " seconds or " << 
	(double) (clock() - sTime)/(60*CLOCKS_PER_SEC) 
		<< " mins" << std::endl;
#ifdef USE_MPI
      std::cout << "Maximum memory used by a processor: " << max_mem7 <<  " Mb" << std::endl;
#endif
      if(nprocs == 1)
	std::cout << "Maximum memory used: " << mem7 <<  " Mb" << std::endl;
      std::cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n"  << std::endl;

    }
    /*********************************************/
    // close
    /*********************************************/
    if (run_flag == 1 && nprocs == 1) {
      err = TheCore.close();
      ERRORR("Failed to dellocate.", 1);
    } else {
      err = TheCore.close_parallel(rank, nprocs);
      ERRORR("Failed to dellocate.", 1);
    }
  }
  else{
    // makefile already generated
  }

#ifdef USE_MPI
  MPI::COMM_WORLD.Barrier();
  MPI::Finalize();
#endif

  return 0;
}
