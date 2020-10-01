#include "battery.h"


#define DIMS 3
using namespace dealii;


int main(int argc, char **argv){
  try{
		Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv,1);
		deallog.depth_console (0);
					
		battery<DIMS> problem(argv[1]);	
    problem.run ();
  }
  catch (std::exception &exc){
    std::cerr << std::endl << std::endl
	      << "----------------------------------------------------"
	      << std::endl;
    std::cerr << "Exception on processing: " << std::endl
	      << exc.what() << std::endl
	      << "Aborting!" << std::endl
	      << "----------------------------------------------------"
	      << std::endl;

    return 1;
  }
  catch (...){
    std::cerr << std::endl << std::endl
	      << "----------------------------------------------------"
	      << std::endl;
    std::cerr << "Unknown exception!" << std::endl
	      << "Aborting!" << std::endl
	      << "----------------------------------------------------"
	      << std::endl;
    return 1;
  }

  return 0;
}
