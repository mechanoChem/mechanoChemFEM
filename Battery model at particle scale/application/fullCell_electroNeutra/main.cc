#include "initBoundValProbs.h"


#define DIMS 2
using namespace dealii;


int main(int argc, char **argv){
  try{
		Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv,1);
		deallog.depth_console (0);

		//main fields 
		std::vector<std::vector<std::string> > primary_variables(9);//u;v;p;_li_plus;phi_e;c_li;phi_s;T;u_mesh		
		primary_variables[0].push_back("u"); primary_variables[0].push_back("component_is_vector");
		primary_variables[1].push_back("v"); primary_variables[1].push_back("component_is_vector");
		primary_variables[2].push_back("p"); primary_variables[2].push_back("component_is_scalar");
		primary_variables[3].push_back("C_li_plus"); primary_variables[3].push_back("component_is_scalar");
		primary_variables[4].push_back("phi_e"); primary_variables[4].push_back("component_is_scalar");
		primary_variables[5].push_back("C_li"); primary_variables[5].push_back("component_is_scalar");
		primary_variables[6].push_back("phi_s"); primary_variables[6].push_back("component_is_scalar");
		primary_variables[7].push_back("T"); primary_variables[7].push_back("component_is_scalar");
		primary_variables[8].push_back("u_mesh"); primary_variables[8].push_back("component_is_vector");
		
		//active material, electrolyte, currentCollector, pure solid
		int number_domain=4;
		std::vector<std::vector<int> > FE_support(number_domain);// store order of polynomial basis functions, 0 means FE_Nothing	
		//active material domain
		FE_support[0].push_back(1);//u
		FE_support[0].push_back(0);//v
		FE_support[0].push_back(0);//p
		FE_support[0].push_back(0);//C_li_plus
		FE_support[0].push_back(0);//phi_e
		FE_support[0].push_back(1);//C_li
		FE_support[0].push_back(1);//phi_s
		FE_support[0].push_back(1);//T
		FE_support[0].push_back(0);//u_mesh
		//electrolyte
		FE_support[1].push_back(0);//u
		FE_support[1].push_back(2);//v
		FE_support[1].push_back(1);//p
		FE_support[1].push_back(1);//C_li_plus
		FE_support[1].push_back(1);//phi_e
		FE_support[1].push_back(0);//C_li
		FE_support[1].push_back(0);//phi_s
		FE_support[1].push_back(1);//T
		FE_support[1].push_back(1);//u_mesh
		//currentCollector and binder (currentCollector_id=binder_id)
		FE_support[2].push_back(1);//u
		FE_support[2].push_back(0);//v
		FE_support[2].push_back(0);//p
		FE_support[2].push_back(0);//C_li_plus
		FE_support[2].push_back(0);//phi_e
		FE_support[2].push_back(0);//C_li
		FE_support[2].push_back(1);//phi_s
		FE_support[2].push_back(1);//T
		FE_support[2].push_back(0);//u_mesh
		//pure solid
		FE_support[3].push_back(1);//u
		FE_support[3].push_back(0);//v
		FE_support[3].push_back(0);//p
		FE_support[3].push_back(0);//C_li_plus
		FE_support[3].push_back(0);//phi_e
		FE_support[3].push_back(0);//C_li
		FE_support[3].push_back(0);//phi_s
		FE_support[3].push_back(1);//T
		FE_support[3].push_back(0);//u_mesh

		
		ParameterHandler params;
		
    initBoundValProbs<DIMS> problem(primary_variables,FE_support,params);
		
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
