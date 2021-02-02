/*
zhenlin wang 2019
*displacement
*/
#include "../../include/battery_components.h"

template <int dim>
void Displacement<dim>::declare_parameters(nlohmann::json& _params){
	params_json=&_params;
}

template <int dim>
void Displacement<dim>::set_stress(dealii::Table<3,Sacado::Fad::DFad<double> >& F, dealii::Table<3,Sacado::Fad::DFad<double> >& P){
	int mat_id = (* this->battery_fields->current_cell)->material_id();
	this->youngsModulus=(*params_json)["Mechanics"]["youngs_modulus_particle"];
	if (mat_id == 2) this->youngsModulus=(*params_json)["Mechanics"]["youngs_modulus_electrolyte"];
	this->poissonRatio=(*params_json)["Mechanics"]["poisson_ratio"];
	this->ResidualEq->setLameParametersByYoungsModulusPoissonRatio(this->youngsModulus, this->poissonRatio);
	
	unsigned int n_q_points= P.size(0);
	int lithium_index=this->battery_fields->active_fields_index["Lithium"];
		
	
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
		double eps_0=1.0e-5;
		int interface_index=this->battery_fields->active_fields_index["Diffuse_interface"];
		if(this->battery_fields->quad_fields[interface_index].value[0]>=0.5 ){
			for(unsigned int q=0; q<n_q_points;q++){
				Sacado::Fad::DFad<double> C_q=this->battery_fields->quad_fields[lithium_index].value[q];
			
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
			else{this->ResidualEq->evaluateNeoHookeanStress(P, F);}
		
	}
}

template class Displacement<1>;
template class Displacement<2>;
template class Displacement<3>;