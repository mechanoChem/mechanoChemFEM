#include "precipitate_PF.h"


#define DIMS 2
using namespace dealii;


int main(int argc, char **argv){
  try{
    Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv,1);
    deallog.depth_console (0);

    /**************************************
     * Define unknown fields
     *************************************/
    std::vector<std::vector<std::string> > primary_variables(6);
    primary_variables[0].push_back("c"); primary_variables[0].push_back("component_is_scalar");
    primary_variables[1].push_back("mu"); primary_variables[1].push_back("component_is_scalar");
    primary_variables[2].push_back("eta1"); primary_variables[2].push_back("component_is_scalar");
    primary_variables[3].push_back("eta2"); primary_variables[3].push_back("component_is_scalar");
    primary_variables[4].push_back("eta3"); primary_variables[4].push_back("component_is_scalar");
    primary_variables[5].push_back("u"); primary_variables[5].push_back("component_is_vector");
		
    int number_domain=1;
    int basis_order=1;
		
    std::vector<std::vector<int> > FE_support(number_domain);// store order of polynomial basis functions, 0 means FE_Nothing	
    FE_support[0].push_back(basis_order);
    FE_support[0].push_back(basis_order);
    FE_support[0].push_back(basis_order);
    FE_support[0].push_back(basis_order);
    FE_support[0].push_back(basis_order);
    FE_support[0].push_back(basis_order);

    /**************************************
     * Define parameters
     *************************************/
    ParameterHandler params;
    params.enter_subsection("parameters");

    //Chemistry
    params.declare_entry("mobility","0",
			 Patterns::Double() );
    params.declare_entry("kinetic_coeff","0",
			 Patterns::Double() );
    params.declare_entry("sqrt_kappa1","0",
			 Patterns::Double() );
    params.declare_entry("sqrt_kappa2","0",
			 Patterns::Double() );
    params.declare_entry("c_avg","0",
			 Patterns::Double() );
    params.leave_subsection();	


    //Mesh refinement
    params.enter_subsection("mesh_refinement");
    params.declare_entry("n_init_global_refine","0",
			 Patterns::Integer() );
    params.declare_entry("n_init_local_refine","0",
			 Patterns::Integer() );
    params.declare_entry("max_refine_level","0",
			 Patterns::Integer() );
    params.leave_subsection();

    /**************************************
     * Run the simulation
     *************************************/
    precipitate_PF<DIMS> problem(primary_variables,FE_support, params);	
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
