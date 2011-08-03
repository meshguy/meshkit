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
  //Initialize MPI
  int rank = 0, nprocs = 1;

#ifdef USE_MPI
  MPI::Init(argc, argv);
  nprocs = MPI::COMM_WORLD.Get_size();
  rank = MPI::COMM_WORLD.Get_rank();
#endif

  // start program timer
  CClock Timer;
  std::string szDateTime;
  clock_t sTime = clock();

  CCrgen TheCore;
  int err = 0;
  int run_flag = 1;

  // memory related variables
  double ld_t = 0;
  unsigned long mem1, mem2, mem3, mem4, mem5, mem6, mem7;

  if (rank == 0) {
    err = TheCore.banner();
    ERRORR("Failed in creating banner", 1);

    Timer.GetDateTime(szDateTime);
    std::cout << "\nStarting out at : " << szDateTime << "\n";
  }
  /*********************************************/
  // read inputs and create makefile, do this on all processors
  /*********************************************/
  err = TheCore.prepareIO(argc, argv, rank, nprocs);
  ERRORR("Failed in preparing i/o.", 1);

  if (argc > 1) {
    if (argv[1][0] == '-' && argv[1][1] == 'm') {
      run_flag = 0;
    }
  }

  // copy, move, merge, extrude, assign gids, save and close
  if (run_flag == 1 ) {
    /*********************************************/
    // load mesh files
    /*********************************************/
    CClock ld_time;
    if (TheCore.prob_type == "mesh") {
      if (nprocs == 1) {
  	err = TheCore.load_meshes();
  	ERRORR("Failed to load meshes.", 1);
	//  	TheCore.pc = new moab::ParallelComm(TheCore.mbImpl());	
      } else {
#ifdef USE_MPI	
	err = TheCore.load_meshes_parallel(rank, nprocs);
	ERRORR("Failed to load meshes.", 1);

	MPI::COMM_WORLD.Barrier();

	if(nprocs > (int) TheCore.files.size()){
	  // if there are more procs than files distribute the copy/move work on each proc
	  err = TheCore.distribute_mesh(rank, nprocs);
	  ERRORR("Failed to load meshes.", 1);
	}
  	//Get a pcomm object
	TheCore.pc = new moab::ParallelComm(TheCore.mbImpl());
#endif
      }
      if (TheCore.mem_tflag == true) {
        TheCore.mbImpl()->estimated_memory_use(0, 0, 0, &mem1);

        ld_t = ld_time.DiffTime();
        std::cout << "\n**from rank: " << rank<< " Time taken to load mesh files = " << ld_t
		  << " seconds" << std::endl;
        std::cout << "***from rank: " << rank<< " Memory used: " << mem1 << " kb\n" << std::endl;
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

    /*********************************************/
    // copy move
    /*********************************************/
    CClock ld_cm;
    if(TheCore.prob_type == "mesh"){ 
      if (nprocs <= (int) TheCore.files.size()) {
	err = TheCore.copy_move_parallel(rank, nprocs);
	ERRORR("Failed in copy move routine.", 1);
      } else {
	err = TheCore.copymove_parallel(rank, nprocs);
	ERRORR("Failed in copy move routine.", 1);
      }
    }
    else{
      err = TheCore.copy_move();
    }
    if (TheCore.mem_tflag == true && (strcmp(TheCore.prob_type.c_str(), "mesh") == 0)) {
      TheCore.mbImpl()->estimated_memory_use(0, 0, 0, &mem2);

      ld_t = ld_cm.DiffTime();
      std::cout << "\n**from rank: " << rank<< " Time taken to copy/move mesh files = " << ld_t
    		<< " seconds" << std::endl;
      std::cout << "***from rank: " << rank<< " Memory used: " << mem2 << " kb\n" << std::endl;
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
    	err = TheCore.merge_nodes();
    	ERRORR("Failed to merge nodes.", 1);
      } else {
	err = TheCore.merge_nodes_parallel(rank, nprocs);
    	ERRORR("Failed to merge nodes in parallel.", 1);
      }

      if (TheCore.mem_tflag == true) {
    	TheCore.mbImpl()->estimated_memory_use(0, 0, 0, &mem3);

    	ld_t = ld_mm.DiffTime();
    	std::cout << "\n**from rank: " << rank<< " Time taken to merge nodes = " << ld_t
    		  << " seconds" << std::endl;
    	std::cout << "***from rank: " << rank<< " Memory used: " << mem3 << " kb\n" << std::endl;
      }

      /*********************************************/
      // extrude
      /*********************************************/
      if (TheCore.extrude_flag == true) {

    	// assign global ids after copy/move step
    	if (rank == 0)
    	  err = TheCore.assign_gids();

    	CClock ld_em;
    	err = TheCore.extrude();
    	ERRORR("Failed to extrude.", 1);

    	if (TheCore.mem_tflag == true) {
    	  TheCore.mbImpl()->estimated_memory_use(0, 0, 0, &mem4);

    	  ld_t = ld_em.DiffTime();
    	  std::cout << "\n**from rank: " << rank<< " Time taken to extrude = " << ld_t
    		    << " seconds" << std::endl;
    	  std::cout << "***from rank: " << rank<< " Memory used: " << mem4 << " kb\n"
    		    << std::endl;
    	}
      }
      /*********************************************/
      // assign gids
      /*********************************************/
#ifdef USE_MPI   
      MPI::COMM_WORLD.Barrier();
#endif
      CClock ld_gid;
      if (nprocs == 1) {
    	err = TheCore.assign_gids();
    	ERRORR("Failed to assign global ids.", 1);
      }
      else{
	err = TheCore.assign_gids_parallel(rank, nprocs);
    	ERRORR("Failed to assign global ids.", 1);
      }

      if (TheCore.mem_tflag == true) {
    	TheCore.mbImpl()->estimated_memory_use(0, 0, 0, &mem5);

    	ld_t = ld_gid.DiffTime();
    	std::cout << "\n**from rank: " << rank<< " Time taken to assign gids = " << ld_t
    		  << " seconds" << std::endl;
    	std::cout << "***from rank: " << rank<< " Memory used: " << mem5 << " kb\n" << std::endl;
      }
      /*********************************************/
      // create neumann sets on the core model
      /*********************************************/
      CClock ld_ns;
      err = TheCore.create_neumannset();
      ERRORR("Failed to create neumann set.", 1);

      if (TheCore.mem_tflag == true) {
    	TheCore.mbImpl()->estimated_memory_use(0, 0, 0, &mem6);

    	ld_t = ld_ns.DiffTime();
    	std::cout << "\n**from rank: " << rank<< " Time taken to create neumann sets = " << ld_t
    		  << " seconds" << std::endl;
    	std::cout << "***from rank: " << rank<< " Memory used: " << mem6 << " kb\n" << std::endl;
      }
      /*********************************************/
      // save
      /*********************************************/
      CClock ld_sv;
      if (nprocs == 1) {
      	err = TheCore.save_mesh();
      	ERRORR("Failed to save o/p file.", 1);
      } else {
	//	err = TheCore.save_mesh(rank); // uncomment to save the meshes with each proc
	//      	ERRORR("Failed to save o/p file.", 1);
#ifdef USE_MPI
	double write_time = MPI_Wtime();
	err = TheCore.save_mesh_parallel(rank, nprocs);
	ERRORR("Failed to save o/p file.", 1);
	write_time = MPI_Wtime() - write_time;
	if (rank == 0){
	  std::cout << "Parallel write time = " << write_time << " seconds" << std::endl;
	}
#endif
      }
      if (TheCore.mem_tflag == true) {
	TheCore.mbImpl()->estimated_memory_use(0, 0, 0, &mem7);

	ld_t = ld_sv.DiffTime();
	std::cout << "\n**from rank: " << rank<< " Time taken to save = " << ld_t << " seconds"
		  << std::endl;
	std::cout << "***from rank: " << rank<< " Memory used: " << mem7 << " kb\n" << std::endl;
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
    // close
    /*********************************************/
    if (run_flag == 1 && nprocs == 1) {
      err = TheCore.close();
      ERRORR("Failed to dellocate.", 1);
    } else {
      err = TheCore.close_parallel(rank, nprocs);
      ERRORR("Failed to dellocate.", 1);
    }
  
    // compute the elapsed time
#ifdef USE_MPI   
    MPI::COMM_WORLD.Barrier();
#endif
    if (rank == 0) {
      Timer.GetDateTime(szDateTime);
      std::cout << "Ending at : " << szDateTime;
      std::cout << "Elapsed wall clock time: " << Timer.DiffTime()
    		<< " seconds or " << (Timer.DiffTime()) / 60.0 << " mins\n";

      std::cout << "Total CPU time used: " << (double) (clock() - sTime)/CLOCKS_PER_SEC \
		<< " seconds" << std::endl;
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
