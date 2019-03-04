#include "initBoundValProbs.h"

template <int dim>
void initBoundValProbs<dim>::make_grid()
{	
	
	params->enter_subsection("Problem");
	std::string meshType=params->get("meshType");
	int initial_global_refine=params->get_integer("initial_global_refine");
	if(std::strcmp(meshType.c_str(),"external")==0 ){
		this->pcout << "reading external  mesh\n";
		
		std::string mesh_directory=params->get("mesh");
  	GridIn<dim> gridin; 
		gridin.attach_triangulation(this->triangulation);
  	std::ifstream f(mesh_directory);
  	//gridin.read_abaqus(f);
		gridin.read_msh(f);
	}
	else if(std::strcmp(meshType.c_str(),"half_hyperShell")==0 ){
		double x=params->get_double("Shell_center_X");
		double y=params->get_double("Shell_center_Y");
		double z=params->get_double("Shell_center_Z");
		double inner_radius=params->get_double("inner_radius");
		double outer_radius=params->get_double("outer_radius");
		const Point<dim> center(x,y,z);
		GridGenerator::half_hyper_shell(this->triangulation, center, inner_radius, outer_radius, 0, false);
    //this->triangulation.set_all_manifold_ids_on_boundary(1);
    static SphericalManifold<dim> manifold(center);
    static HyperShellBoundary<dim> boundary(center);
    this->triangulation.set_manifold (0,manifold);
		this->triangulation.set_all_manifold_ids(0);
		//this->triangulation.set_manifold (1, boundary);	
	}
	this->triangulation.refine_global (initial_global_refine);
	params->leave_subsection();	
	
}


template class initBoundValProbs<1>;
template class initBoundValProbs<2>;
template class initBoundValProbs<3>;
