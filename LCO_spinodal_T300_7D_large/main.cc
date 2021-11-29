#include "CahnHilliard.h"
#include <chrono>

#define DIMS 2
using namespace dealii;


int main(int argc, char **argv){
  try{
    Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv,1);
    deallog.depth_console (0);
    CahnHilliard<DIMS> problem;	

    std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now(); 
    problem.run ();
    std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();

    std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);

    std::cout << "It took me " << time_span.count() << " seconds.";
    std::cout << std::endl;
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
