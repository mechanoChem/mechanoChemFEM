#include "diffusion_reaction.h"


#define DIMS 2
using namespace dealii;


int main(int argc, char **argv){
  try{
		Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv,1);
		deallog.depth_console (0);

		//main fields 
		std::vector<std::vector<std::string> > primary_variables(2);		
	  primary_variables[0].push_back("c1"); primary_variables[0].push_back("component_is_scalar");
	  primary_variables[1].push_back("c2"); primary_variables[1].push_back("component_is_scalar");
		
		int number_domain=1;
		int diff_degree=1;
		
		std::vector<std::vector<int> > FE_support(number_domain);// store order of polynomial basis functions, 0 means FE_Nothing	
		FE_support[0].push_back(diff_degree);
		FE_support[0].push_back(diff_degree);
		
		ParameterHandler params;
		params.enter_subsection("Concentration");	
		params.declare_entry("D_1","0",Patterns::Double() );
		params.declare_entry("D_2","0",Patterns::Double() );
		params.declare_entry("R_10","0",Patterns::Double() );
		params.declare_entry("R_11","0",Patterns::Double() );
		params.declare_entry("R_12","0",Patterns::Double() );
		params.declare_entry("R_13","0",Patterns::Double() );
		params.declare_entry("R_20","0",Patterns::Double() );
		params.declare_entry("R_21","0",Patterns::Double() );
		params.declare_entry("R_22","0",Patterns::Double() );
		params.declare_entry("R_23","0",Patterns::Double() );
		params.declare_entry("jn","0",Patterns::Double() );
		params.leave_subsection();	
		
			
		diffusion_reaction<2> problem(primary_variables,FE_support, params);	
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
