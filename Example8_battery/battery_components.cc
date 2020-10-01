/*
zhenlin wang 2019
*coupled diffusion reaction
*/


/*
*Lithium
*/
#include "battery_package/include/battery_components.h"
template <int dim>
void Lithium<dim>::set_diffusion_term(dealii::Table<2,Sacado::Fad::DFad<double> >& diffu){
	double M=(*params_json)["Lithium"]["lithium_diffusivity"];
	M=1;
	int phaesField_index=this->Solution->active_fields_index["Lithium_phaseField"];
	if(phaesField_index==-1) diffu=table_scaling<2,Sacado::Fad::DFad<double>,double >(this->Solution->quad_fields[this->primiary_dof].value_grad,-M);
	else{
		diffu=table_scaling<2,Sacado::Fad::DFad<double>,double >(this->Solution->quad_fields[phaesField_index].value_grad,-M);
	}
}

template class Lithium<1>;
template class Lithium<2>;
template class Lithium<3>;


/*
*Lithium phaseField
*/
template <int dim>
void Lithium_phaseField<dim>::set_field_and_source_term(dealii::Table<2,Sacado::Fad::DFad<double> >& field, dealii::Table<1,Sacado::Fad::DFad<double> >& source)
{
	int Lithium_index=this->Solution->active_fields_index["Lithium_phaseField"];
	if (Lithium_index==-1){std::cout<<"Lithium_phaseField is defined, but Lithium is not defined"<<std::endl; exit(-1);}
	double kappa=(*params_json)["Lithium_phaseField"]["kappa"];
	double C_alpha=(*params_json)["Lithium_phaseField"]["c_alpha"];
	double C_beta=(*params_json)["Lithium_phaseField"]["c_beta"];
	double omega=(*params_json)["Lithium_phaseField"]["omega"];
	
	field=table_scaling<2,Sacado::Fad::DFad<double>,double >(this->Solution->quad_fields[Lithium_index].value_grad,-kappa);
	unsigned int n_q_points= source.size(0);
	dealii::Table<1,Sacado::Fad::DFad<double> > C_li(this->Solution->quad_fields[Lithium_index].value);
	for (unsigned int q=0;q<n_q_points;q++){
		source[q]=2*omega*(C_li[q]-C_alpha)*(C_li[q]-C_beta)*(2*C_li[q]-C_alpha-C_beta)-this->Solution->quad_fields[this->primiary_dof].value[q];
	}
}

template class Lithium_phaseField<1>;
template class Lithium_phaseField<2>;
template class Lithium_phaseField<3>;



/*
*Displacement
*/
template <int dim>
void Displacement<dim>::set_stress(dealii::Table<3,Sacado::Fad::DFad<double> >& F, dealii::Table<3,Sacado::Fad::DFad<double> >& P){
	this->youngsModulus=(*params_json)["Mechanics"]["youngs_modulus"];
	this->poissonRatio=(*params_json)["Mechanics"]["poisson_ratio"];
	this->ResidualEq->setLameParametersByYoungsModulusPoissonRatio(this->youngsModulus, this->poissonRatio);
	
	unsigned int n_q_points= P.size(0);
	int lithium_index=this->Solution->active_fields_index["Lithium"];
		
	
	if (lithium_index==-1){
		this->ResidualEq->evaluateNeoHookeanStress(P, F);
	}
	else {
		dealii::Table<3,Sacado::Fad::DFad<double> > Fe(n_q_points,dim,dim);
		double C_a=(*params_json)["Mechanics"]["lithium_a"];
		double C_b=(*params_json)["Mechanics"]["lithium_b"];
		dealii::Table<2,Sacado::Fad::DFad<double>> Feiga(dim,dim);
		dealii::Table<2,Sacado::Fad::DFad<double>> Feigba(dim,dim);
		Feiga[0][0]=(*params_json)["Mechanics"]["Feiga_11"];
		Feigba[0][0]=(*params_json)["Mechanics"]["Feigb_11"];
		Feigba[0][0]-=(*params_json)["Mechanics"]["Feiga_11"].get<double>();
		if(dim>=2){
			Feiga[1][1]=(*params_json)["Mechanics"]["Feiga_22"];
			Feigba[1][1]=(*params_json)["Mechanics"]["Feigb_22"];
			Feigba[1][1]-=(*params_json)["Mechanics"]["Feiga_22"].get<double>();
		}
		if(dim==3){
			Feiga[2][2]=(*params_json)["Mechanics"]["Feiga_33"];
			Feigba[2][2]=(*params_json)["Mechanics"]["Feigb_33"];
			Feigba[2][2]-=(*params_json)["Mechanics"]["Feiga_33"].get<double>();
		}
		
		for(unsigned int q=0; q<n_q_points;q++){
			Sacado::Fad::DFad<double> C_q=this->Solution->quad_fields[lithium_index].value[q];
			
					
			dealii::Table<2,Sacado::Fad::DFad<double> > Feig(dim,dim);
			dealii::Table<2,Sacado::Fad::DFad<double>> invFeig(dim,dim);
			Feig=table_scaling<2,Sacado::Fad::DFad<double>,Sacado::Fad::DFad<double> > (Feigba, (C_q-C_a)/(C_b-C_a) );   
			Feig=table_add<2,Sacado::Fad::DFad<double>,Sacado::Fad::DFad<double> > (Feig, Feiga);
			getInverse<Sacado::Fad::DFad<double>,dim> (Feig,invFeig);
			for (unsigned int i=0; i<dim; ++i){
				for (unsigned int j=0; j<dim; ++j){
					for (unsigned int k=0; k<dim; ++k){
		 				Fe[q][i][j]+=F[q][i][k]*invFeig[k][j];
					}
				}
			}
		}
		this->ResidualEq->evaluateNeoHookeanStress(P, Fe);
	}
}

template class Displacement<1>;
template class Displacement<2>;
template class Displacement<3>;