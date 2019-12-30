#include "CahnHilliard.h"


#define DIMS 2
using namespace dealii;


int main(int argc, char **argv){
  try{
		Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv,1);
		deallog.depth_console (0);

		//main fields 
		std::vector<std::vector<std::string> > primary_variables(4);		
	  primary_variables[0].push_back("c1"); primary_variables[0].push_back("component_is_scalar");
	  primary_variables[1].push_back("mu1"); primary_variables[1].push_back("component_is_scalar");
	  primary_variables[2].push_back("c2"); primary_variables[2].push_back("component_is_scalar");
	  primary_variables[3].push_back("mu2"); primary_variables[3].push_back("component_is_scalar");
		
		int number_domain=1;
		int diff_degree=1;
		
		std::vector<std::vector<int> > FE_support(number_domain);// store order of polynomial basis functions, 0 means FE_Nothing	
		for (unsigned int i=0;i<4;i++) FE_support[0].push_back(diff_degree);
		
		ParameterHandler params;
		params.enter_subsection("Problem");
		params.declare_entry("output_w_theta","true",Patterns::Bool());
		params.leave_subsection();
		params.enter_subsection("Concentration");
		params.declare_entry("c1_ini","0",Patterns::Double() );
		params.declare_entry("c2_ini","0",Patterns::Double() );
		params.declare_entry("mobility_1","0",Patterns::Double() );
		params.declare_entry("mobility_2","0",Patterns::Double() );
		params.declare_entry("kappa_1","0",Patterns::Double() );
		params.declare_entry("kappa_2","0",Patterns::Double() );
		
		params.declare_entry("d","0",Patterns::Double() );
		params.declare_entry("s","0",Patterns::Double() );
		
		params.leave_subsection();		
		
			
		CahnHilliard<DIMS> problem(primary_variables,FE_support, params);	
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
