#include"../../include/supplementary/computedField.h"

template <int dim>
computedField<dim>::computedField (){}
	
template <int dim>
computedField<dim>::~computedField (){}

template <int dim>
void computedField<dim>::setupComputedField(std::vector<std::vector<std::string> > _primary_variables)
{
	primary_variables=_primary_variables;
}



template <int dim>
std::vector<std::string> computedField<dim>::get_names() const
{
  std::vector<std::string> solution_names;
	for(unsigned int i=0;i<primary_variables.size();i++){
		if(std::strcmp(primary_variables[i][1].c_str(),"component_is_scalar")==0) {
			solution_names.push_back(primary_variables[i][0]);
		}
		else if(std::strcmp(primary_variables[i][1].c_str(),"component_is_vector")==0){
		  for (unsigned int j=0; j<dim; ++j){
		    solution_names.push_back(primary_variables[i][0]);
			}
		}
		else{ std::cout<<"primary_variables component type does not support \n"; exit(1);}
	}
  return solution_names;
}



template <int dim>
std::vector<DataComponentInterpretation::DataComponentInterpretation>
computedField<dim>::get_data_component_interpretation () const
{
  std::vector<DataComponentInterpretation::DataComponentInterpretation> interpretation;
	for(unsigned int i=0;i<primary_variables.size();i++){
		if(std::strcmp(primary_variables[i][1].c_str(),"component_is_scalar")==0) {
			interpretation.push_back(dealii::DataComponentInterpretation::component_is_scalar);
		}
		else if(std::strcmp(primary_variables[i][1].c_str(),"component_is_vector")==0){
		  for (unsigned int j=0; j<dim; ++j){
		    interpretation.push_back(dealii::DataComponentInterpretation::component_is_part_of_vector);
			}
		}
		else{ std::cout<<"primary_variables component type does not support \n"; exit(1);}
	}
  return interpretation;
}

template <int dim>
UpdateFlags
computedField<dim>::get_needed_update_flags() const
{
  return update_values | update_gradients | update_q_points;
}

template class computedField<1>;
template class computedField<2>;
template class computedField<3>;
