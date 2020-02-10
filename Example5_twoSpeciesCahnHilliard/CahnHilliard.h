/*
zhenlin wang 2019
*CahnHilliard
*/
#include "mechanoChemFEM.h"
#include "supplementary/computedField.h"
//#include <deal.II/lac/affine_constraints.h>

template <int dim>
class nodalField : public computedField<dim>
{
public:
	nodalField(dealii::ParameterHandler& _params);
	
	dealii::ParameterHandler* params;
	void compute_derived_quantities_vector(const std::vector<Vector<double> > &uh,
					       const std::vector<std::vector<Tensor<1,dim> > > &duh,
					       const std::vector<std::vector<Tensor<2,dim> > > &dduh,
					       const std::vector<Point<dim> >                  &normals,
					       const std::vector<Point<dim> >                  &evaluation_points,
					       std::vector<Vector<double> >                    &computed_quantities) const;								 							 
};

template <int dim>
nodalField<dim>::nodalField(dealii::ParameterHandler& _params):params(&_params)
{}

template <int dim>
void nodalField<dim>::compute_derived_quantities_vector(const std::vector<Vector<double> > &uh,
					       const std::vector<std::vector<Tensor<1,dim> > > &duh,
					       const std::vector<std::vector<Tensor<2,dim> > > &dduh,
					       const std::vector<Point<dim> >                  &normals,
					       const std::vector<Point<dim> >                  &evaluation_points,
					       std::vector<Vector<double> >                    &computed_quantities) const
{
	
	const unsigned int n_q_points = uh.size();
	for(unsigned int q=0; q<n_q_points;q++){
		double c1=uh[q][0], c2=uh[q][2];
		if (c2+0.866*c1 > 0 and c1 >= 0) computed_quantities[q][0]=1;
		else if (c2-0.866*c1>= 0 and c1 < 0) computed_quantities[q][0]=0;
		else computed_quantities[q][0]=-1;
	}
}

template <int dim>
class CahnHilliard: public mechanoChemFEM<dim>
{
	public:
		CahnHilliard();
		//this is a overloaded function 
		void ini_updateLinearSystem();
		void get_residual(const typename hp::DoFHandler<dim>::active_cell_iterator &cell, const FEValues<dim>& fe_values, Table<1, Sacado::Fad::DFad<double> >& R, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<1, double >& ULocalConv);
		void solve_ibvp();
		void save_features();
		int iter_count=0;
		ParameterHandler* params;		
		nodalField<dim> computedNodalField;
		std::vector<double> local_features;
		std::vector<double> features;
		bool output_w_theta;
		
		//AffineConstraints<double> null_constraints;
		

};
template <int dim>
CahnHilliard<dim>::CahnHilliard():computedNodalField(*params)
{
		this->pcout<<"CahnHilliard initiated"<<std::endl;
		//This let you use one params to get all parameters pre-defined in the mechanoChemFEM
		params=this->params_mechanoChemFEM;
		
		params->enter_subsection("Problem");
		params->declare_entry("output_w_theta","true",Patterns::Bool());
		params->leave_subsection();
		params->enter_subsection("Concentration");
		params->declare_entry("c1_ini","0",Patterns::Double() );
		params->declare_entry("c2_ini","0",Patterns::Double() );
		params->declare_entry("mobility_1","0",Patterns::Double() );
		params->declare_entry("mobility_2","0",Patterns::Double() );
		params->declare_entry("kappa_1","0",Patterns::Double() );
		params->declare_entry("kappa_2","0",Patterns::Double() );
		
		params->declare_entry("d","0",Patterns::Double() );
		params->declare_entry("s","0",Patterns::Double() );
		
		params->leave_subsection();		
		
		//Declear the parameters before load it
		this->load_parameters("../parameters.prm");
		//initiate computed fields.
		std::vector<std::vector<std::string> > computed_primary_variables={ {"theta", "component_is_scalar"}};
		computedNodalField.setupComputedField(computed_primary_variables);
		params->enter_subsection("Problem");
		output_w_theta=params->get_bool("output_w_theta");
		params->leave_subsection();
		local_features.resize(4,0.0);
		features.resize(4,0.0);
		
		//define main fields from parameter file.
		this->define_primary_fields();
		//Set up the ibvp.
		this->init_ibvp();
	}


template <int dim>
void CahnHilliard<dim>::ini_updateLinearSystem()
{
	for (unsigned int i=0;i<4;i++)
	{
		local_features[i]=0;
		features[i]=0;
	}
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
	MPI_Reduce(&local_features[0], &features[0], 4, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);	
	
	save_features();
	if(output_w_theta){
		Vector<double> localized_U(this->solution_prev);
		this->FEMdata_out.clear_data_vectors();
		Vector<float> material_id(this->triangulation.n_active_cells()); 
		this->FEMdata_out.data_out.add_data_vector(localized_U, computedNodalField);
		std::string output_path = this->output_directory+"output-"+std::to_string(this->current_increment+this->off_output_index)+".vtk";
		this->FEMdata_out.write_vtk(this->solution_prev,output_path);
		this->save_output=false;
	}


	
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
	
	//calculate features
	dealii::Table<1,Sacado::Fad::DFad<double> > theta1(n_q_points);
	dealii::Table<1,Sacado::Fad::DFad<double> > theta2(n_q_points);
	dealii::Table<1,Sacado::Fad::DFad<double> > theta3(n_q_points);
	dealii::Table<1,Sacado::Fad::DFad<double> > gamma(n_q_points);
	for(unsigned int q=0; q<n_q_points;q++){
		if (c2[q]+0.866*c1[q] > 0 and c1[q] >= 0) theta1[q]=1;
		else if (c2[q]-0.866*c1[q] >= 0 and c1[q] < 0) theta2[q]=1;
		else theta3[q]=1;
		gamma[q]=0.5*(c1_grad[q][0]*c1_grad[q][0]+c1_grad[q][1]*c1_grad[q][1]);
	}
	local_features[0]+=this->ResidualEq.volumeIntegration(fe_values, theta1);
	local_features[1]+=this->ResidualEq.volumeIntegration(fe_values, theta2);
	local_features[2]+=this->ResidualEq.volumeIntegration(fe_values, theta3);
	local_features[3]+=this->ResidualEq.volumeIntegration(fe_values, gamma);
}


template <int dim>
void CahnHilliard<dim>::save_features(){
	if (this->this_mpi_process == 0 ){
	  std::ofstream myfile;
	  myfile.open ( (this->output_directory+"Features.dat").c_str(),std::ios_base::app);
		if(!myfile.is_open()) {std::cout<<"file failed to open!"; exit(1);}
		myfile <<this->current_increment<<" ";
		for(unsigned int i=0;i<4;i++) {
		  myfile <<features[i]<<" ";
		}
		myfile <<this->current_time<<" ";
		myfile <<"\n";
	}
	
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
