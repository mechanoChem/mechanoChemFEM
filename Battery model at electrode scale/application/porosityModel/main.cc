#include "initBoundValProbs.h"

#define DIMS 3
using namespace dealii;


int main(int argc, char **argv){
  try{
		Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv,1);
		deallog.depth_console (0);
		
		ParameterHandler params;
		initBoundValProbs<DIMS> problem(params);
		ElectricChemo<Sacado::Fad::DFad<double>,DIMS> _electricChemoFormula(params);
		problem.electricChemoFormula=&_electricChemoFormula;
		problem.declare_parameters();
		params.read_input ("parameters.prm");
		
		params.enter_subsection("Problem");	
		int separator_fe=params.get_integer("separator_fe");
		int electrode_fe=params.get_integer("electrode_fe");
		params.leave_subsection();
		
		//main fields 
		std::vector<std::vector<std::string> > primary_variables(6);//u; c_li_plus;phi_e;c_li;phi_s;T	
		primary_variables[0].push_back("u"); primary_variables[0].push_back("component_is_vector");
		primary_variables[1].push_back("C_li_plus"); primary_variables[1].push_back("component_is_scalar");
		primary_variables[2].push_back("phi_e"); primary_variables[2].push_back("component_is_scalar");
		primary_variables[3].push_back("C_li"); primary_variables[3].push_back("component_is_scalar");
		primary_variables[4].push_back("phi_s"); primary_variables[4].push_back("component_is_scalar");
		primary_variables[5].push_back("T"); primary_variables[5].push_back("component_is_scalar");
		
		//active material, separator, currentCollector, pure solid
		int number_domain=2;
		std::vector<std::vector<int> > FE_support(number_domain);// store order of polynomial basis functions, 0 means FE_Nothing	
		//electrode domain
		FE_support[separator_fe].push_back(1);//u
		FE_support[separator_fe].push_back(1);//C_li_plus
		FE_support[separator_fe].push_back(1);//phi_e
		FE_support[separator_fe].push_back(0);//C_li
		FE_support[separator_fe].push_back(0);//phi_s
		FE_support[separator_fe].push_back(1);//T
	
		//separator
		FE_support[electrode_fe].push_back(1);//u
		FE_support[electrode_fe].push_back(1);//C_li_plus
		FE_support[electrode_fe].push_back(1);//phi_e
		FE_support[electrode_fe].push_back(1);//C_li
		FE_support[electrode_fe].push_back(1);//phi_s
		FE_support[electrode_fe].push_back(1);//T

    problem.run (primary_variables, FE_support);
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
