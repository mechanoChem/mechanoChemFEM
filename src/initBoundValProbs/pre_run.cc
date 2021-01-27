/*
zhenlin wang 2019
*/

#include"../../include/mechanoChemFEM.h"
#include <cstdlib>

template <int dim>
void mechanoChemFEM<dim>::define_primary_fields()
{
	//define_primary_fields from parameters file.
	std::vector<std::string> primary_variables_s;
	std::vector<int> FE_support_list_v;
	
	if(this->use_ParameterHandler){
		params_mechanoChemFEM->enter_subsection("Problem");
		std::string primary_variables_list=params_mechanoChemFEM->get("primary_variables_list");
		std::string FE_support_list=params_mechanoChemFEM->get("FE_support_list");
		primary_variables_s=Utilities::split_string_list(primary_variables_list);
		std::vector<std::string> FE_support_list_s=Utilities::split_string_list(FE_support_list);
		FE_support_list_v=Utilities::string_to_int(FE_support_list_s);
		
		params_mechanoChemFEM->leave_subsection();		
	}
	if(this->use_ParameterJson){
		primary_variables_s=(*params_mechanoChemFEM_json)["Problem"]["primary_variables_list"].get<std::vector<std::string>>();
		FE_support_list_v=(*params_mechanoChemFEM_json)["Problem"]["FE_support_list"].get<std::vector<int>>();
	}
	
	int num_variables=primary_variables_s.size()/2;
	int num_domain=FE_support_list_v.size()/num_variables;
	primary_variables.resize(num_variables);
	FE_support.resize(num_domain);
	pcout<<"=========== Defining Primary Fields ==========="<<std::endl;
	pcout<<"primary variables are:"<<std::endl;
	for(unsigned int i=0;i< num_variables;i++){
		pcout<<primary_variables_s[i*2]<<" "<<primary_variables_s[i*2+1]<<std::endl;
		primary_variables[i].push_back(primary_variables_s[i*2]);
		primary_variables[i].push_back(primary_variables_s[i*2+1]);
	}
	pcout<<std::endl;
	for(unsigned int i=0;i< num_domain;i++){
		pcout<<"for domain number id "<<i<<" the basis order is: "<<std::endl;
		for(unsigned int j=0;j< num_variables;j++){
			pcout<<FE_support_list_v[i*num_variables+j]<<"  ";
			FE_support[i].push_back(FE_support_list_v[i*num_variables+j]);
		}
		pcout<<std::endl;
	}
	pcout<<std::endl;
	pcout<<"=========== Primary Fields Defined ==========="<<std::endl;
	if(primary_variables_s.size()%2!=0 or FE_support_list_v.size()%num_variables!=0){pcout<<"variable and domains are not consistent!!"<<std::endl; exit(1);}
	
}

template<int dim>
void mechanoChemFEM<dim>::init_ibvp()
{
	FEMdata_out.set_output_name(primary_variables);	
	this->setup_FeSystem(fe_system, fe_collection, q_collection, primary_variables_dof,primary_variables,FE_support,*volume_quadrature);
	make_grid();
	refine_grid();
	setMultDomain();
  mark_boundary();
	setup_linear_system();
	apply_boundary_condition();
  pcout << "   Number of active cells:       " << hpFEM<dim>::triangulation.n_active_cells() << std::endl;
  pcout << "   Number of degrees of freedom: " << hpFEM<dim>::dof_handler.n_dofs() << std::endl; 
	if(!resuming_from_snapshot) {
		apply_initial_condition();		
	  std::string output_path = output_directory+"output-0.vtk";
	  FEMdata_out.write_vtk(solution, output_path);
	  if(save_snapshot){
			std::string snapshot_path = snapshot_directory+"snapshot-"+std::to_string(0)+".dat";
	  	FEMdata_out.create_vector_snapshot(solution, snapshot_path);
		}
	}
	else {
	  pcout<<"resuming from snapshot"<<std::endl;
	  FEMdata_out.resume_vector_from_snapshot(solution,snapfile);
	  std::string output_path = output_directory+"output-resume.vtk";
    FEMdata_out.write_vtk(solution, output_path);
	  solution_prev=solution;
	}
}

template class mechanoChemFEM<1>;
template class mechanoChemFEM<2>;
template class mechanoChemFEM<3>;
