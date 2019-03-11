/*
zhenlin wang 2019
*CahnHilliard
*/
#include "initBoundValProbs.h"
template <int dim>
class CahnHilliard: public initBoundValProbs<dim>
{
	public:
		CahnHilliard(std::vector<std::vector<std::string> > _primary_variables, std::vector<std::vector<int> > _FE_support, ParameterHandler& _params);
		//this is a overloaded function 
		void get_residual(const typename hp::DoFHandler<dim>::active_cell_iterator &cell, const FEValues<dim>& fe_values, Table<1, Sacado::Fad::DFad<double> >& R, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<1, double >& ULocalConv);
		ParameterHandler* params;		
};
template <int dim>
CahnHilliard<dim>::CahnHilliard(std::vector<std::vector<std::string> > _primary_variables, std::vector<std::vector<int> > _FE_support, ParameterHandler& _params)
	:initBoundValProbs<dim>(_primary_variables, _FE_support, _params),params(&_params){}

template <int dim>
void CahnHilliard<dim>::get_residual(const typename hp::DoFHandler<dim>::active_cell_iterator &cell, const FEValues<dim>& fe_values, Table<1, Sacado::Fad::DFad<double> >& R, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<1, double >& ULocalConv)
{
	//evaluate primary fields
	params->enter_subsection("Concentration");
	double M=params->get_double("M");
	double omega=params->get_double("omega");
	double c_alpha=params->get_double("c_alpha");
	double c_beta=params->get_double("c_beta");
	double kappa=params->get_double("kappa");
	params->leave_subsection();
	unsigned int n_q_points= fe_values.n_quadrature_points;
	int c_dof=0, mu_dof=1;
		
	dealii::Table<1,double>  c_1_conv(n_q_points);
	dealii::Table<1,Sacado::Fad::DFad<double> > c_1(n_q_points), mu(n_q_points);
	dealii::Table<2,Sacado::Fad::DFad<double> >  c_1_grad(n_q_points, dim), mu_grad(n_q_points, dim);
	
	evaluateScalarFunction<double,dim>(fe_values, c_dof, ULocalConv, c_1_conv);
	evaluateScalarFunction<Sacado::Fad::DFad<double>,dim>(fe_values, c_dof, ULocal, c_1);	
	evaluateScalarFunctionGradient<Sacado::Fad::DFad<double>,dim>(fe_values, c_dof, ULocal, c_1_grad);
	
	evaluateScalarFunction<Sacado::Fad::DFad<double>,dim>(fe_values, mu_dof, ULocal, mu);	
	evaluateScalarFunctionGradient<Sacado::Fad::DFad<double>,dim>(fe_values, mu_dof, ULocal, mu_grad);
	
	//evaluate diffusion and reaction term
	dealii::Table<1,Sacado::Fad::DFad<double> > rhs_mu(n_q_points);
	dealii::Table<2,Sacado::Fad::DFad<double> > j_c_1(n_q_points, dim), kappa_c_1_grad(n_q_points, dim);
	
	j_c_1=table_scaling<Sacado::Fad::DFad<double>, dim>(mu_grad,-M);//-D_1*c_1_grad
	kappa_c_1_grad=table_scaling<Sacado::Fad::DFad<double>, dim>(c_1_grad,kappa);
	
	for(unsigned int q=0; q<n_q_points;q++) rhs_mu[q]=2*omega*(c_1[q]-c_alpha)*(c_1[q]-c_beta)*(2*c_1[q]-c_alpha-c_beta)-mu[q];
	
	//call residual functions
	this->ResidualEq.residualForDiffusionEq(fe_values, c_dof, R, c_1, c_1_conv, j_c_1);
	this->ResidualEq.residualForPoissonEq(fe_values, mu_dof, R, kappa_c_1_grad, rhs_mu);
	
}

template <int dim>
void InitialConditions<dim>::vector_value (const Point<dim>   &p, Vector<double>   &values) const{
  Assert (values.size() == 2, ExcDimensionMismatch (values.size(), 2));
  values(1) = 0;    
 values(0)= 0.5 + 0.04*(static_cast <double> (rand())/(static_cast <double>(RAND_MAX))-0.5);
}
template class InitialConditions<1>;
template class InitialConditions<2>;
template class InitialConditions<3>;