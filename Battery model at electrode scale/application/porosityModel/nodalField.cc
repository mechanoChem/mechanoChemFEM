#include "nodalField.h"

template <int dim>
nodalField<dim>::nodalField(dealii::ParameterHandler& _params):params(&_params){}

template <int dim>
nodalField<dim>::~nodalField(){}

template <int dim>
void nodalField<dim>::compute_derived_quantities_vector(const std::vector<Vector<double> > &uh,
					       const std::vector<std::vector<Tensor<1,dim> > > &duh,
					       const std::vector<std::vector<Tensor<2,dim> > > &dduh,
					       const std::vector<Point<dim> >                  &normals,
					       const std::vector<Point<dim> >                  &evaluation_points,
					       std::vector<Vector<double> >                    &computed_quantities) const
{
	ElectricChemo<double,dim> electricChemoFormula;
 	electricChemoFormula.params=params;
 	electricChemoFormula.setParametersFromHandler();

	params->enter_subsection("Geometry");
	double electrode_Y1=params->get_double("electrode_Y1");
	double electrode_Y2=params->get_double("electrode_Y2");
  params->leave_subsection();
	
	int u_dof=primary_variables_dof[0];
	int c_li_plus_dof=primary_variables_dof[1];
	int phi_e_dof=primary_variables_dof[2];
	int c_li_dof=primary_variables_dof[3];
	int phi_s_dof=primary_variables_dof[4];
	int T_dof=primary_variables_dof[5];
	
	const unsigned int dof_per_node = uh.size();
	for (unsigned int q=0; q<dof_per_node; ++q){
		computed_quantities[q][0]=0;
		
		double c_li_plus, phi_e, c_li, phi_s, T;
		c_li_plus=uh[q][c_li_plus_dof];
		phi_e=uh[q][phi_e_dof];
		c_li=uh[q][c_li_dof];
		phi_s=uh[q][phi_s_dof];
		T=uh[q][T_dof];
		int domain;
		
		if(evaluation_points[q][1]<=electrode_Y1) domain=-1;
		else if(evaluation_points[q][1]>=electrode_Y2)domain=1;
		else domain=0;
		computed_quantities[q][0]=electricChemoFormula.formula_jn(T, c_li, c_li_plus, phi_s, phi_e, domain);		
	}	
}

template class nodalField<1>;
template class nodalField<2>;
template class nodalField<3>;