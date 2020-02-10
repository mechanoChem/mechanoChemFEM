#include "L2projection.h"

#define DIMS 2
using namespace dealii;


int main(int argc, char **argv){
  try{
		Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv,1);
		deallog.depth_console (0);

		//main fields 
		std::vector<std::vector<std::string> > primary_variables_add(4);		
	  primary_variables_add[0].push_back("C1"); primary_variables_add[0].push_back("component_is_scalar");
	  primary_variables_add[1].push_back("mu1"); primary_variables_add[1].push_back("component_is_scalar");
	  primary_variables_add[2].push_back("C2"); primary_variables_add[2].push_back("component_is_scalar");
	  primary_variables_add[3].push_back("mu2"); primary_variables_add[3].push_back("component_is_scalar");
			
		int number_domain=1;
		int mech_degree=1;
		
		std::vector<std::vector<int> > FE_support_add(number_domain);// store order of polynomial basis functions, 0 means FE_Nothing	
		for(unsigned int i=0;i<4;i++) FE_support_add[0].push_back(1);	
		L2projection<DIMS> problem(primary_variables_add,FE_support_add, params);	
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
