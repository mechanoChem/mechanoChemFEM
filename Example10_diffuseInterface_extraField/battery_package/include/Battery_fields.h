/*
zhenlin wang 2020
*module transportation
*/
#ifndef Battery_fields_h
#define Battery_fields_h

#include "mechanoChemFEM.h"
template <int dim>
class Battery_fields
{
	public:
		Battery_fields();
		int num_primary_fields;
		struct Fields {
		  dealii::Table<1,double>  value_conv;
		  dealii::Table<1,Sacado::Fad::DFad<double> > value;
			dealii::Table<2,Sacado::Fad::DFad<double> > value_grad;
		};
		std::map<std::string,int> active_fields_index;
		std::vector<struct Fields> quad_fields;
		void declare_parameters(nlohmann::json& _params);
		virtual void set_up_active_fields(std::vector<std::vector<std::string> > primary_variables);
		void update_fields(const typename hp::DoFHandler<dim>::active_cell_iterator &cell, const FEValues<dim>& fe_values, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<1, double >& ULocalConv);
		
		nlohmann::json* params_json;
		const typename hp::DoFHandler<dim>::active_cell_iterator *current_cell;
		int pos_electrode_domain_id=0, separator_domain_id=1, neg_electrode_domain_id=2;
		int current_domain_id=0;
		ConditionalOStream pcout;
};
#endif