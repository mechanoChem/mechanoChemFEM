/*
zhenlin wang 2020
*module transportation
*/
#ifndef Battery_fields_h
#define Battery_fields_h

#include "mechanoChemFEM.h"
class Battery_fields
{
	public:
		Battery_fields();
		Battery_fields(int _dim);
		int num_primary_fields;
		struct Fields {
		  dealii::Table<1,double>  value_conv;
		  dealii::Table<1,Sacado::Fad::DFad<double> > value;
			dealii::Table<2,Sacado::Fad::DFad<double> > value_grad;
		};
		std::map<std::string,int> active_fields_index;
		
		std::vector<struct Fields> quad_fields;
		void declare_parameters_batteryFields();
		virtual void set_up_active_fields(std::vector<std::vector<std::string> > primary_variables, int dim);
		ParameterHandler* params_battery_fileds;
		int dim;
		
		ConditionalOStream pcout;
};
#endif