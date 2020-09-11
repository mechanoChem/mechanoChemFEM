/*
zhenlin wang 2019
*coupled diffusion reaction
*/
#include "mechanoChemFEM.h"
template <int dim>
class restart: public mechanoChemFEM<dim>
{
	public:
		restart();
		//this is a overloaded function 
		void get_residual(const typename hp::DoFHandler<dim>::active_cell_iterator &cell, const FEValues<dim>& fe_values, Table<1, Sacado::Fad::DFad<double> >& R, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<1, double >& ULocalConv);
		ParameterHandler* params;	
		void save_mesh();	
};
template <int dim>
restart<dim>::restart()
{
	//This let you use one params to get all parameters pre-defined in the mechanoChemFEM
	params=this->params_mechanoChemFEM;	
	//Declear the parameters before load it
	this->load_parameters("../parameters.prm");
	
	//define main fields from parameter file.
	this->define_primary_fields();
	//Set up the ibvp.
	this->init_ibvp();
	save_mesh();
}

template <int dim>
void restart<dim>::save_mesh()
{
	 std::string triangulation_name = "mesh_restart";
	 this->triangulation.save(triangulation_name.c_str());
	 
}

template <int dim>
void restart<dim>::get_residual(const typename hp::DoFHandler<dim>::active_cell_iterator &cell, const FEValues<dim>& fe_values, Table<1, Sacado::Fad::DFad<double> >& R, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<1, double >& ULocalConv)
{
	//evaluate primary fields
	unsigned int n_q_points= fe_values.n_quadrature_points;
	dealii::Table<1,Sacado::Fad::DFad<double> > c_1(n_q_points), target(n_q_points);
	
	evaluateScalarFunction<Sacado::Fad::DFad<double>,dim>(fe_values, 0, ULocal, c_1);	
	for (unsigned int q=0;q<n_q_points;q++) target[q]=c_1[q]-1;
	this->ResidualEq.residualForEqualityEq(fe_values, 0, R,target);
}



template <int dim>
void InitialConditions<dim>::vector_value (const Point<dim>   &p, Vector<double>   &values) const{
  Assert (values.size() == 2, ExcDimensionMismatch (values.size(), 2));
  values(1) = 0;    
}
template class InitialConditions<1>;
template class InitialConditions<2>;
template class InitialConditions<3>;
