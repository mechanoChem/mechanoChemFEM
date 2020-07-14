/*
zhenlin wang 2019
*coupled diffusion reaction
*/
#include "mechanoChemFEM.h"
template <int dim>
class diffusion_reaction: public mechanoChemFEM<dim>
{
	public:
		diffusion_reaction();
		//this is a overloaded function 
		void get_residual(const typename hp::DoFHandler<dim>::active_cell_iterator &cell, const FEValues<dim>& fe_values, Table<1, Sacado::Fad::DFad<double> >& R, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<1, double >& ULocalConv);
		ParameterHandler* params;		
};
template <int dim>
diffusion_reaction<dim>::diffusion_reaction()
{
	//This let you use one params to get all parameters pre-defined in the mechanoChemFEM
	params=this->params_mechanoChemFEM;
	
	params->enter_subsection("Concentration");	
	params->declare_entry("D_1","0",Patterns::Double() );
	params->declare_entry("D_2","0",Patterns::Double() );
	params->declare_entry("R_10","0",Patterns::Double() );
	params->declare_entry("R_11","0",Patterns::Double() );
	params->declare_entry("R_12","0",Patterns::Double() );
	params->declare_entry("R_13","0",Patterns::Double() );
	params->declare_entry("R_20","0",Patterns::Double() );
	params->declare_entry("R_21","0",Patterns::Double() );
	params->declare_entry("R_22","0",Patterns::Double() );
	params->declare_entry("R_23","0",Patterns::Double() );
	params->declare_entry("jn","0",Patterns::Double() );
	params->leave_subsection();
	
	//Declear the parameters before load it
	this->load_parameters("parameters.prm");
	
	//define main fields from parameter file.
	this->define_primary_fields();
	//Set up the ibvp.
	this->init_ibvp();
}


template <int dim>
void diffusion_reaction<dim>::get_residual(const typename hp::DoFHandler<dim>::active_cell_iterator &cell, const FEValues<dim>& fe_values, Table<1, Sacado::Fad::DFad<double> >& R, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<1, double >& ULocalConv)
{
	//evaluate primary fields
	params->enter_subsection("Concentration");
	double D_1=params->get_double("D_1");
	double D_2=params->get_double("D_2");
	double R_10=params->get_double("R_10");
	double R_11=params->get_double("R_11");
	double R_12=params->get_double("R_12");
	double R_13=params->get_double("R_13");
	double R_20=params->get_double("R_20");
	double R_21=params->get_double("R_21");
	double R_22=params->get_double("R_22");
	double R_23=params->get_double("R_23");
	double jn=params->get_double("jn");
	params->leave_subsection();	
	unsigned int n_q_points= fe_values.n_quadrature_points;
	int c_1_dof=0, c_2_dof=1;
		
	dealii::Table<1,double>  c_1_conv(n_q_points), c_2_conv(n_q_points);
	dealii::Table<1,Sacado::Fad::DFad<double> > c_1(n_q_points), c_2(n_q_points);
	dealii::Table<2,Sacado::Fad::DFad<double> >  c_1_grad(n_q_points, dim), c_2_grad(n_q_points, dim);
	
	evaluateScalarFunction<double,dim>(fe_values, c_1_dof, ULocalConv, c_1_conv);
	evaluateScalarFunction<Sacado::Fad::DFad<double>,dim>(fe_values, c_1_dof, ULocal, c_1);	
	evaluateScalarFunctionGradient<Sacado::Fad::DFad<double>,dim>(fe_values, c_1_dof, ULocal, c_1_grad);
	
	evaluateScalarFunction<double,dim>(fe_values, c_2_dof, ULocalConv, c_2_conv);
	evaluateScalarFunction<Sacado::Fad::DFad<double>,dim>(fe_values, c_2_dof, ULocal, c_2);	
	evaluateScalarFunctionGradient<Sacado::Fad::DFad<double>,dim>(fe_values, c_2_dof, ULocal, c_2_grad);

	
	//evaluate diffusion and reaction term
	dealii::Table<1,Sacado::Fad::DFad<double> > reaction_1(n_q_points), reaction_2(n_q_points);
	dealii::Table<2,Sacado::Fad::DFad<double> > j_c_1(n_q_points, dim),j_c_2(n_q_points, dim);
	
	j_c_1=table_scaling<Sacado::Fad::DFad<double>, dim>(c_1_grad,-D_1);//-D_1*c_1_grad
	j_c_2=table_scaling<Sacado::Fad::DFad<double>, dim>(c_2_grad,-D_2);//-D_2*c_2_grad
	
	for(unsigned int q=0; q<n_q_points;q++){
		reaction_1[q]=R_10+R_11*c_1[q]+R_12*c_2[q]+R_13*c_1[q]*c_1[q]*c_2[q];
		reaction_2[q]=R_20+R_21*c_1[q]+R_22*c_2[q]+R_23*c_1[q]*c_1[q]*c_2[q];
	}
	//call residual functions
	this->ResidualEq.residualForDiff_ReacEq(fe_values, c_1_dof, R, c_1, c_1_conv, j_c_1, reaction_1);
	this->ResidualEq.residualForDiff_ReacEq(fe_values, c_2_dof, R, c_2, c_2_conv, j_c_2, reaction_2);
	
	//BC
	for (unsigned int faceID=0; faceID<2*dim; faceID++){
		if(cell->face(faceID)->boundary_id()==dim*2 ){
		  FEFaceValues<dim> fe_face_values(fe_values.get_fe(), *(this->common_face_quadrature), update_values | update_quadrature_points | update_JxW_values);
			fe_face_values.reinit(cell,faceID);
			this->ResidualEq.residualForNeummanBC(fe_values, fe_face_values, c_1_dof, R, jn);
		}
	}
	
}



template <int dim>
void InitialConditions<dim>::vector_value (const Point<dim>   &p, Vector<double>   &values) const{
  Assert (values.size() == 2, ExcDimensionMismatch (values.size(), 2));
  values(1) = 0;    
  values(0)= 0.5 + 0.1*static_cast <double> (rand())/(static_cast <double>(RAND_MAX/2.0))/2;
}
template class InitialConditions<1>;
template class InitialConditions<2>;
template class InitialConditions<3>;
