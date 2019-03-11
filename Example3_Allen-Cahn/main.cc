#include "AllenCahn.h"


#define DIMS 2
using namespace dealii;


int main(int argc, char **argv){
  try{
		Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv,1);
		deallog.depth_console (0);

		//main fields 
		std::vector<std::vector<std::string> > primary_variables(1);		
	  primary_variables[0].push_back("mu"); primary_variables[0].push_back("component_is_scalar");
		
		int number_domain=1;
		int diff_degree=1;
		
		std::vector<std::vector<int> > FE_support(number_domain);// store order of polynomial basis functions, 0 means FE_Nothing	
		FE_support[0].push_back(diff_degree);
		
		ParameterHandler params;
		params.enter_subsection("Parameters");	
		params.declare_entry("omega","0",Patterns::Double() );
		params.declare_entry("c_alpha","0",Patterns::Double() );
		params.declare_entry("c_beta","0",Patterns::Double() );
		params.declare_entry("M","0",Patterns::Double() );
		params.declare_entry("kappa","0",Patterns::Double() );
		params.leave_subsection();		
		
			
		AllenCahn<DIMS> problem(primary_variables,FE_support, params);	
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
