#include "initBoundValProbs.h"


#define DIMS 3
using namespace dealii;


int main(int argc, char **argv){
  try{
		Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv,1);
		deallog.depth_console (0);

		//main fields 
		std::vector<std::vector<std::string> > primary_variables(3);		
		primary_variables[0].push_back("u"); primary_variables[0].push_back("component_is_vector");
	        primary_variables[1].push_back("c1"); primary_variables[1].push_back("component_is_scalar");
	        primary_variables[2].push_back("c2"); primary_variables[2].push_back("component_is_scalar");
		// primary_variables[1].push_back("advection"); primary_variables[2].push_back("component_is_scalar");
		
		int number_domain=3;
		int mech_degree=1;
		int diff_degree=1;
		std::vector<std::vector<int> > FE_support(number_domain);// store order of polynomial basis functions, 0 means FE_Nothing	
		for(unsigned int i=0;i<2;i++){
			FE_support[i].push_back(mech_degree);
			FE_support[i].push_back(diff_degree);
			FE_support[i].push_back(diff_degree);
			//FE_support[i].push_back(diff_degree);
		}
		FE_support[2].push_back(0);
		FE_support[2].push_back(0);
		FE_support[2].push_back(0);
		//FE_support[2].push_back(0);

		
		std::vector<std::vector<std::string> > variables_add(3);	
		variables_add[0].push_back("advection"); variables_add[0].push_back("component_is_vector");
		variables_add[1].push_back("c1"); variables_add[1].push_back("component_is_scalar");
		variables_add[2].push_back("c2"); variables_add[2].push_back("component_is_scalar");
		std::vector<std::vector<int> > FE_support_add(number_domain);// store order of polynomial basis functions, 0 means FE_Nothing	
		for(unsigned int i=0;i<2;i++){
			FE_support_add[i].push_back(mech_degree);
			FE_support_add[i].push_back(diff_degree);
			FE_support_add[i].push_back(diff_degree);
		}
		FE_support_add[2].push_back(0);
		FE_support_add[2].push_back(0);
		FE_support_add[2].push_back(0);


		ParameterHandler params;
    initBoundValProbs<DIMS> problem(primary_variables,FE_support,variables_add, FE_support_add, params);
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
