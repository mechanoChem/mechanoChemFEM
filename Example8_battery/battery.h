/*
zhenlin wang 2019
*coupled diffusion reaction
*/
#include "mechanoChemFEM.h"
#include "battery_package/include/battery_components.h"

template <int dim>
class battery: public mechanoChemFEM<dim>
{
	public:
		battery();
		//this is a overloaded function 
		void apply_boundary_condition();
		void get_residual(const typename hp::DoFHandler<dim>::active_cell_iterator &cell, const FEValues<dim>& fe_values, Table<1, Sacado::Fad::DFad<double> >& R, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<1, double >& ULocalConv);
		ParameterHandler* params;		
		ConstraintMatrix* constraints;
		
		Battery_fields battery_fields;
		Lithium<dim> lithium;
		Displacement<dim> displacement;
		
		
};
template <int dim>
battery<dim>::battery()
{
	//pass the pointer to "constraints" in that defined in mechanoChemFEM
	constraints=this->constraints_mechanoChemFEM;
	//This let you use one params to get all parameters pre-defined in the mechanoChemFEM
	params=this->params_mechanoChemFEM;
	
	lithium.declare_parameters(*params);
	displacement.declare_parameters(*params);
	//Declear the parameters before load it
	this->load_parameters("../parameters.prm");		
	this->define_primary_fields();
	

	//set active battery fields
	battery_fields.set_up_active_fields(this->primary_variables,dim);	
	if(battery_fields.active_fields_index["Lithium"]>-1) lithium.set_up_fields(battery_fields, this->ResidualEq, battery_fields.active_fields_index["Lithium"]);
	if(battery_fields.active_fields_index["Displacement"]>-1) displacement.set_up_fields(battery_fields,this->ResidualEq, battery_fields.active_fields_index["Displacement"]);

	this->init_ibvp();
}



template <int dim>
void battery<dim>::get_residual(const typename hp::DoFHandler<dim>::active_cell_iterator &cell, const FEValues<dim>& fe_values, Table<1, Sacado::Fad::DFad<double> >& R, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<1, double >& ULocalConv)
{
	if(battery_fields.active_fields_index["Lithium"]>-1) lithium.r_get_residual(cell,fe_values,R,ULocal, ULocalConv);
	if(battery_fields.active_fields_index["Displacement"]>-1) displacement.r_get_residual(cell,fe_values,R,ULocal, ULocalConv);
	
}


//set Dirichlet BC
template <int dim>
void battery<dim>::apply_boundary_condition()
{
	constraints->clear ();
	DoFTools::make_hanging_node_constraints (this->dof_handler, *constraints);
	
	int totalDOF=this->totalDOF(this->primary_variables);
  std::vector<bool> All_component (totalDOF, true);	
	VectorTools:: interpolate_boundary_values (this->dof_handler, dim, ZeroFunction<dim> (totalDOF),*constraints, All_component);
	
	constraints->close ();
}

template <int dim>
void InitialConditions<dim>::vector_value (const Point<dim>   &p, Vector<double>   &values) const{
	Assert (values.size() == totalDOF, ExcDimensionMismatch (values.size(), totalDOF));
	values=0;
	for(unsigned int i=0;i<primary_variables.size();i++){
		if(std::strcmp(primary_variables[i][0].c_str(),"Lithium")==0){
			values(primary_variables_dof[i])= 1.0;
			if(p[2]==0) values(primary_variables_dof[i])= 1.5;
		}
		//if(std::strcmp(primary_variables[i][0].c_str(),"Displacement")==0) values(primary_variables_dof[i])= 1.0;
	}
}

template class InitialConditions<1>;
template class InitialConditions<2>;
template class InitialConditions<3>;