/*
zhenlin wang 2019
*AllenCahn
*/
#include "initBoundValProbs.h"
template <int dim>
class AllenCahn: public initBoundValProbs<dim>
{
	public:
		AllenCahn(std::vector<std::vector<std::string> > _primary_variables, std::vector<std::vector<int> > _FE_support, ParameterHandler& _params);
		//this is a overloaded function 
		void get_residual(const typename hp::DoFHandler<dim>::active_cell_iterator &cell, const FEValues<dim>& fe_values, Table<1, Sacado::Fad::DFad<double> >& R, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<1, double >& ULocalConv);
		ParameterHandler* params;		
};
template <int dim>
AllenCahn<dim>::AllenCahn(std::vector<std::vector<std::string> > _primary_variables, std::vector<std::vector<int> > _FE_support, ParameterHandler& _params)
	:initBoundValProbs<dim>(_primary_variables, _FE_support, _params),params(&_params){}

template <int dim>
void AllenCahn<dim>::get_residual(const typename hp::DoFHandler<dim>::active_cell_iterator &cell, const FEValues<dim>& fe_values, Table<1, Sacado::Fad::DFad<double> >& R, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<1, double >& ULocalConv)
{
	//evaluate primary fields
	params->enter_subsection("Parameters");
	double M=params->get_double("M");
	double omega=params->get_double("omega");
	double c_alpha=params->get_double("c_alpha");
	double c_beta=params->get_double("c_beta");
	double kappa=params->get_double("kappa");
	params->leave_subsection();
	unsigned int n_q_points= fe_values.n_quadrature_points;
		
	dealii::Table<1,double>  mu_conv(n_q_points);
	dealii::Table<1,Sacado::Fad::DFad<double> >  mu(n_q_points);
	dealii::Table<2,Sacado::Fad::DFad<double> >  mu_grad(n_q_points, dim);
	
	evaluateScalarFunction<double,dim>(fe_values, 0, ULocalConv, mu_conv);
	evaluateScalarFunction<Sacado::Fad::DFad<double>,dim>(fe_values, 0, ULocal, mu);	
	evaluateScalarFunctionGradient<Sacado::Fad::DFad<double>,dim>(fe_values, 0, ULocal, mu_grad);
	
	//evaluate diffusion and reaction term
	dealii::Table<1,Sacado::Fad::DFad<double> > rhs_mu(n_q_points);
	dealii::Table<2,Sacado::Fad::DFad<double> > j_mu(n_q_points, dim);
	
	j_mu=table_scaling<Sacado::Fad::DFad<double>, dim>(mu_grad,-kappa*M);
	
	for(unsigned int q=0; q<n_q_points;q++){
		 rhs_mu[q]=-M*(2*omega*(mu[q]-c_alpha)*(mu[q]-c_beta)*(2*mu[q]-c_alpha-c_beta));
	 }
	
	//call residual functions
	this->ResidualEq.residualForDiff_ReacEq(fe_values, 0, R, mu, mu_conv, j_mu, rhs_mu);
	
}
template <int dim>
double InitialConditions<dim>::value(const Point<dim>   &p, const unsigned int 	component) const{
  return 0.5 + 0.2*(static_cast <double> (rand())/(static_cast <double>(RAND_MAX))-0.5); 
}

template class InitialConditions<1>;
template class InitialConditions<2>;
template class InitialConditions<3>;