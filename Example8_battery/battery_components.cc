/*
zhenlin wang 2019
*coupled diffusion reaction
*/
#include "battery_package/include/battery_components.h"
template <int dim>
void Lithium<dim>::set_diffusion_term(dealii::Table<2,Sacado::Fad::DFad<double> >& diffu){
	params->enter_subsection("parameters_Lithium");
	double M=params->get_double("lithium_diffusivity");
	params->leave_subsection();
	diffu=table_scaling<Sacado::Fad::DFad<double>, dim>(this->Solution->quad_fields[this->primiary_dof].value_grad,-M);
}

template class Lithium<1>;
template class Lithium<2>;
template class Lithium<3>;

template <int dim>
void Displacement<dim>::set_stress(dealii::Table<3,Sacado::Fad::DFad<double> >& F, dealii::Table<3,Sacado::Fad::DFad<double> >& P){
	params->enter_subsection("parameters_Mechanics");
	this->youngsModulus=params->get_double("youngsModulus");
	this->poissonRatio=params->get_double("poissonRatio");
	this->ResidualEq->setLameParametersByYoungsModulusPoissonRatio(this->youngsModulus, this->poissonRatio);
	params->leave_subsection();
	
	unsigned int n_q_points= P.size(0);
	int lithium_index=this->Solution->active_fields_index["Lithium"];
		
	
	if (lithium_index==-1){
		this->ResidualEq->evaluateNeoHookeanStress(P, F);
	}
	else {
		dealii::Table<3,Sacado::Fad::DFad<double> > Fe(n_q_points,dim,dim);
		params->enter_subsection("parameters_Mechanics");
		double C_0=params->get_double("lithium_0");
		params->leave_subsection();
		for(unsigned int q=0; q<n_q_points;q++){			
			for (unsigned int i=0; i<dim; ++i){
				for (unsigned int j=0; j<dim; ++j){
		 			Fe[q][i][j]=F[q][i][j]/std::pow((this->Solution->quad_fields[lithium_index].value[q]/C_0), 1.0/3.0); //Isotropic growth
				}
			}
		}
		this->ResidualEq->evaluateNeoHookeanStress(P, Fe);
	}
}

template class Displacement<1>;
template class Displacement<2>;
template class Displacement<3>;