/*
zhenlin wang 2019
*CahnHilliard
*/
#include "mechanoChemFEM.h"
template <int dim>
class CahnHilliard: public mechanoChemFEM<dim>
{
	public:
		CahnHilliard(std::vector<std::vector<std::string> > _primary_variables, std::vector<std::vector<int> > _FE_support, ParameterHandler& _params);
		//this is a overloaded function 
		void get_residual(const typename hp::DoFHandler<dim>::active_cell_iterator &cell, const FEValues<dim>& fe_values, Table<1, Sacado::Fad::DFad<double> >& R, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<1, double >& ULocalConv);
		void solve_ibvp();
		int iter_count=0;
		ParameterHandler* params;		
};
template <int dim>
CahnHilliard<dim>::CahnHilliard(std::vector<std::vector<std::string> > _primary_variables, std::vector<std::vector<int> > _FE_support, ParameterHandler& _params)
	:mechanoChemFEM<dim>(_primary_variables, _FE_support, _params),params(&_params){
		this->pcout<<"CahnHilliard initiated"<<std::endl;
	}

template <int dim>
void CahnHilliard<dim>::solve_ibvp()
{		
	int reduce_time=0;
	while(reduce_time<10) {
		bool converge_flag=this->nonlinearSolve(this->solution);
		//reset iter_count to 0
		if (converge_flag) break;
		else{
			iter_count=0;
			this->pcout<<"not converge, reduce dt by half"<<std::endl;
			this->current_dt *= 0.5;
		}		
	}
	iter_count++;
	if(iter_count>5) {
		iter_count=0;
		this->current_dt *= 2;//double dt
	}
	this->solution_prev=this->solution;
}

template <int dim>
void CahnHilliard<dim>::get_residual(const typename hp::DoFHandler<dim>::active_cell_iterator &cell, const FEValues<dim>& fe_values, Table<1, Sacado::Fad::DFad<double> >& R, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<1, double >& ULocalConv)
{
	//evaluate primary fields
	params->enter_subsection("Concentration");
	double M1=params->get_double("mobility_1");
	double M2=params->get_double("mobility_2");
	double k1=params->get_double("kappa_1");
	double k2=params->get_double("kappa_2");
	
	double d=params->get_double("d");
	double s=params->get_double("s");

	params->leave_subsection();
	unsigned int n_q_points= fe_values.n_quadrature_points;
	int c1_dof=0, mu1_dof=1,c2_dof=2, mu2_dof=3;
  //define fields
	dealii::Table<1,double>  c1_conv(n_q_points);
	dealii::Table<1,Sacado::Fad::DFad<double> > c1(n_q_points), mu1(n_q_points);
	dealii::Table<2,Sacado::Fad::DFad<double> >  c1_grad(n_q_points, dim), mu1_grad(n_q_points, dim);
	
	dealii::Table<1,double>  c2_conv(n_q_points);
	dealii::Table<1,Sacado::Fad::DFad<double> > c2(n_q_points), mu2(n_q_points);
	dealii::Table<2,Sacado::Fad::DFad<double> >  c2_grad(n_q_points, dim), mu2_grad(n_q_points, dim);
  //evaluate fields
	evaluateScalarFunction<double,dim>(fe_values, c1_dof, ULocalConv, c1_conv);
	evaluateScalarFunction<Sacado::Fad::DFad<double>,dim>(fe_values, c1_dof, ULocal, c1);	
	evaluateScalarFunctionGradient<Sacado::Fad::DFad<double>,dim>(fe_values, c1_dof, ULocal, c1_grad);
	evaluateScalarFunction<double,dim>(fe_values, c2_dof, ULocalConv, c2_conv);
	evaluateScalarFunction<Sacado::Fad::DFad<double>,dim>(fe_values, c2_dof, ULocal, c2);	
	evaluateScalarFunctionGradient<Sacado::Fad::DFad<double>,dim>(fe_values, c2_dof, ULocal, c2_grad);
	
	evaluateScalarFunction<Sacado::Fad::DFad<double>,dim>(fe_values, mu1_dof, ULocal, mu1);	
	evaluateScalarFunctionGradient<Sacado::Fad::DFad<double>,dim>(fe_values, mu1_dof, ULocal, mu1_grad);
	evaluateScalarFunction<Sacado::Fad::DFad<double>,dim>(fe_values, mu2_dof, ULocal, mu2);	
	evaluateScalarFunctionGradient<Sacado::Fad::DFad<double>,dim>(fe_values, mu2_dof, ULocal, mu2_grad);
	
	
	//evaluate diffusion and reaction term
	dealii::Table<1,Sacado::Fad::DFad<double> > rhs_mu1(n_q_points);
	dealii::Table<2,Sacado::Fad::DFad<double> > j_c1(n_q_points, dim), kappa_c1_grad(n_q_points, dim);
	dealii::Table<1,Sacado::Fad::DFad<double> > rhs_mu2(n_q_points);
	dealii::Table<2,Sacado::Fad::DFad<double> > j_c2(n_q_points, dim), kappa_c2_grad(n_q_points, dim);
	
	j_c1=table_scaling<Sacado::Fad::DFad<double>, dim>(mu1_grad,-M1);//-D_1*c_1_grad
	j_c2=table_scaling<Sacado::Fad::DFad<double>, dim>(mu2_grad,-M2);//-D_1*c_1_grad
	kappa_c1_grad=table_scaling<Sacado::Fad::DFad<double>, dim>(c1_grad,k1);
	kappa_c2_grad=table_scaling<Sacado::Fad::DFad<double>, dim>(c2_grad,k2);
	
	for(unsigned int q=0; q<n_q_points;q++){
		 Sacado::Fad::DFad<double> F_c1=6*d/std::pow(s,4)*(c1[q]*c1[q]+c2[q]*c2[q])*c1[q]-6*d/std::pow(s,3)*c1[q]*c2[q]-3*d/std::pow(s,2)*c1[q];
		 Sacado::Fad::DFad<double> F_c2=6*d/std::pow(s,4)*(c1[q]*c1[q]+c2[q]*c2[q])*c2[q]+3*d/std::pow(s,3)*c2[q]*c2[q]-3*d/std::pow(s,2)*c2[q];
		
		 rhs_mu1[q]=F_c1-mu1[q];
		 rhs_mu2[q]=F_c2-mu2[q];
	}
	
	//call residual functions
	this->ResidualEq.residualForDiffusionEq(fe_values, c1_dof, R, c1, c1_conv, j_c1);
	this->ResidualEq.residualForPoissonEq(fe_values, mu1_dof, R, kappa_c1_grad, rhs_mu1);
	this->ResidualEq.residualForDiffusionEq(fe_values, c2_dof, R, c2, c2_conv, j_c2);
	this->ResidualEq.residualForPoissonEq(fe_values, mu2_dof, R, kappa_c2_grad, rhs_mu2);
	
}

template <int dim>
void InitialConditions<dim>::vector_value (const Point<dim>   &p, Vector<double>   &values) const{
	params->enter_subsection("Concentration");
 	values(0)= params->get_double("c1_ini") + 0.01*(static_cast <double> (rand())/(static_cast <double>(RAND_MAX))-0.5);
  values(1) = 0;    
	values(2)= params->get_double("c2_ini") + 0.01*(static_cast <double> (rand())/(static_cast <double>(RAND_MAX))-0.5);
 	values(3) = 0;    
	params->leave_subsection();
}
template class InitialConditions<1>;
template class InitialConditions<2>;
template class InitialConditions<3>;