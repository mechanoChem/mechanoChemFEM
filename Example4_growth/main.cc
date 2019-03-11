#include "growth.h"


#define DIMS 3
using namespace dealii;


int main(int argc, char **argv){
  try{
		Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv,1);
		deallog.depth_console (0);

		//main fields 
		std::vector<std::vector<std::string> > primary_variables(2);		
	  primary_variables[0].push_back("c1"); primary_variables[0].push_back("component_is_scalar");
	  primary_variables[1].push_back("u"); primary_variables[1].push_back("component_is_vector");
		
		int number_domain=2;
		int basis_order=1;
		
		std::vector<std::vector<int> > FE_support(number_domain);// store order of polynomial basis functions, 0 means FE_Nothing	
		FE_support[0].push_back(basis_order);
		FE_support[0].push_back(0);
		
		FE_support[1].push_back(basis_order);
		FE_support[1].push_back(basis_order);
		
		ParameterHandler params;
		params.enter_subsection("parameters");
		params.declare_entry("youngsModulus","0",Patterns::Double() );
		params.declare_entry("poissonRatio","0",Patterns::Double() );
		params.declare_entry("c_ini","0",Patterns::Double() );
		params.declare_entry("M","0",Patterns::Double() );
		params.leave_subsection();	
		
			
		growth<DIMS> problem(primary_variables,FE_support, params);	
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
