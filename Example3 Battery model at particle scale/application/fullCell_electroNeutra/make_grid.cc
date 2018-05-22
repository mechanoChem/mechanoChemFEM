#include "initBoundValProbs.h"

template <int dim>
void initBoundValProbs<dim>::make_grid()
{	
	params->enter_subsection("Problem");
	std::string mesh_directory=params->get("mesh");
	this->pcout << "reading external  mesh:"<<mesh_directory<<std::endl;
	params->leave_subsection();	
  GridIn<dim> gridin;
  gridin.attach_triangulation(this->triangulation);
  std::ifstream f(mesh_directory);
  gridin.read_abaqus(f);
	//gridin.read_msh(f);
	
}


template class initBoundValProbs<1>;
template class initBoundValProbs<2>;
template class initBoundValProbs<3>;
