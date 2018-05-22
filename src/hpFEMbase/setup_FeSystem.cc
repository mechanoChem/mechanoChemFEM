#include"../../include/hpFEM.h"

template <int dim>
void hpFEM<dim>::setup_FeSystem(std::vector<std::shared_ptr<FESystem<dim>> >& fe_system, hp::FECollection<dim>& fe_collection, hp::QCollection<dim>& q_collection, std::vector<unsigned int >& primary_variables_dof, 
									 std::vector<std::vector<std::string> >& primary_variables, std::vector<std::vector<int> >& FE_support, const QGauss<dim>& volume_quadrature)
{
	
	std::vector<unsigned int> fe_multiplicities;
	int tem_dof=0;
	for(unsigned int i=0;i<primary_variables.size();i++){
		if(std::strcmp(primary_variables[i][1].c_str(),"component_is_scalar")==0) {
			fe_multiplicities.push_back (1);
			primary_variables_dof.push_back(tem_dof);
			tem_dof++;
		}
		else if(std::strcmp(primary_variables[i][1].c_str(),"component_is_vector")==0){
			fe_multiplicities.push_back (dim);
			primary_variables_dof.push_back(tem_dof);
			tem_dof+=dim;
		}	
		else{
			printf ("component type not supported \n"); exit(1);
		}
	}
	
	for(unsigned int i=0;i<FE_support.size();i++){
		std::vector<const FiniteElement<dim>*> tem_fe_component;
		for(unsigned int j=0;j<primary_variables.size();j++){
			if(FE_support[i][j]>0){
				tem_fe_component.push_back (new FE_Q<dim>(FE_support[i][j]));
			}
			else if(FE_support[i][j]==0){
				tem_fe_component.push_back (new FE_Nothing<dim>());
			}
			else{
				printf ("negative order of polynomial basis functions \n"); exit(1);
			}
		}
			
		std::shared_ptr<FESystem<dim>> tem_fe(new FESystem<dim>(tem_fe_component, fe_multiplicities));//T	
		fe_system.push_back(tem_fe);
		
		fe_collection.push_back(*fe_system[i]);
		q_collection.push_back (volume_quadrature);
		
		for (unsigned int i=0; i<tem_fe_component.size(); ++i) {
			delete tem_fe_component[i];
		}	
	}				 
}

template <int dim>
int hpFEM<dim>::totalDOF(std::vector<std::vector<std::string> >& primary_variables)
{
	int tem_dof=0;
	for(unsigned int i=0;i<primary_variables.size();i++){
		if(std::strcmp(primary_variables[i][1].c_str(),"component_is_scalar")==0) {
			tem_dof++;
		}
		else if(std::strcmp(primary_variables[i][1].c_str(),"component_is_vector")==0){
			tem_dof+=dim;
		}	
		else{
			printf ("component type not supported \n"); exit(1);
		}
	}
	return tem_dof;
}

template <int dim>
void hpFEM<dim>::set_active_fe_indices (std::vector<std::vector<int> >& FE_support, hp::DoFHandler<dim>& local_dof_handler, int domain)
{
  typename hp::DoFHandler<dim>::active_cell_iterator cell = local_dof_handler.begin_active(), endc=local_dof_handler.end();
  for (;cell!=endc; ++cell){
		for(unsigned int i=0;i<FE_support.size();i++){
			if (cell->material_id()==i+domain){
				cell->set_active_fe_index(i);
    	}
		}
	}
}

template class hpFEM<1>;
template class hpFEM<2>;
template class hpFEM<3>;