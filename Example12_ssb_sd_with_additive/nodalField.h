/**
*author zhenlin wang, 2018
*/
#ifndef nodalField_h
#define nodalField_h

#include <supplementary/computedField.h>
#include "battery_package/include/Battery_fields.h"

template <int dim>
class nodalField : public computedField<dim>
{
public:
	nlohmann::json* params_json;
	Battery_fields<dim>* battery_fields; 
	void setup(Battery_fields<dim>& _battery_fields, nlohmann::json& _params){battery_fields=&_battery_fields; params_json=&_params;}
	void evaluate_vector_field(const DataPostprocessorInputs::Vector< dim > &input_data, std::vector< Vector< double >> &computed_quantities) const;
								 							 
	std::vector<unsigned int > primary_variables_dof;
};

#endif
