/*
zhenlin wang 2019
*/

#include"../../include/mechanoChemFEM.h"
#include <cstdlib>

template <int dim>
void mechanoChemFEM<dim>::define_primary_fields()
{
	//define_primary_fields from parameters file.
	params_mechanoChemFEM->enter_subsection("Problem");
	std::string primary_variables_list=params_mechanoChemFEM->get("primary_variables_list");
	std::string FE_support_list=params_mechanoChemFEM->get("FE_support_list");
	params_mechanoChemFEM->leave_subsection();		
	
	std::vector<std::string> primary_variables_s=Utilities::split_string_list(primary_variables_list);
	std::vector<std::string> FE_support_list_s=Utilities::split_string_list(FE_support_list);
	std::vector<int> FE_support_list_v=Utilities::string_to_int(FE_support_list_s);
	
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
			FE_support[i].push_back(FE_support_list_v[j]);
		}
	}
	pcout<<std::endl;
	pcout<<"=========== Primary Fields Defined ==========="<<std::endl;
	if(primary_variables_s.size()%2!=0 or FE_support_list_v.size()%2!=0){pcout<<"variable and domains are not consistent!!"<<std::endl; exit(1);}
	
}

template <int dim>
void mechanoChemFEM<dim>::init_ibvp()
{
	FEMdata_out.set_output_name(primary_variables);	
	this->setup_FeSystem(fe_system, fe_collection, q_collection, primary_variables_dof,primary_variables,FE_support,*volume_quadrature);
	make_grid();
	refine_grid();
	setMultDomain();
  mark_boundary();
	setup_linear_system();
}

template class mechanoChemFEM<1>;
template class mechanoChemFEM<2>;
template class mechanoChemFEM<3>;
