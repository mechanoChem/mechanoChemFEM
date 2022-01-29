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
	double separator_line=(*params_json)["ElectroChemo"]["separator_line"];
	int orientation=(*params_json)["ElectroChemo"]["orientation"];
	this->youngsModulus=(*params_json)["Mechanics"]["youngs_modulus_particle"];
	if (mat_id == 1) this->youngsModulus=(*params_json)["Mechanics"]["youngs_modulus_electrolyte"];
	this->poissonRatio=(*params_json)["Mechanics"]["poisson_ratio"];
	this->ResidualEq->setLameParametersByYoungsModulusPoissonRatio(this->youngsModulus, this->poissonRatio);
	
	Point<dim> center=this->battery_fields->cell_center;

//electrolyte_id=0, active_particle_id=1, interface_id=2, li_metal_id=3, additive_id=4, li_metal_interface_id=5, additive_interface_id=6
//
	int electrolyte_id=(*params_json)["Mechanics"]["electrolyte_id"];
	int active_particle_id=(*params_json)["Mechanics"]["active_particle_id"];
	int li_metal_id=(*params_json)["Mechanics"]["li_metal_id"];
	
	unsigned int n_q_points= P.size(0);
	int lithium_index=this->battery_fields->active_fields_index["Lithium"];
		
	if (lithium_index==-1 or mat_id== electrolyte_id){
		this->ResidualEq->evaluateNeoHookeanStress(P, F);
	}
	else {
		dealii::Table<3,Sacado::Fad::DFad<double> > Fe(n_q_points,dim,dim);
		double C_a=(*params_json)["Mechanics"]["lithium_a"];
		double C_b=(*params_json)["Mechanics"]["lithium_b"];
		dealii::Table<2,Sacado::Fad::DFad<double>> Feiga(dim,dim);
		dealii::Table<2,Sacado::Fad::DFad<double>> Feigba(dim,dim);

    //---------------------------------------------------------------------
    // zhenlin's old implementation, difficult to generate cracks as the 1st time step requires a equilibrium process
		//Feiga[0][0]=(*params_json)["Mechanics"]["Feiga_11"]; // = 1
		//Feigba[0][0]=(*params_json)["Mechanics"]["Feigb_11"]; // 0.9
		//Feigba[0][0]-=(*params_json)["Mechanics"]["Feiga_11"].get<double>(); // -0.1
		//if(dim>=2){
			//Feiga[1][1]=(*params_json)["Mechanics"]["Feiga_22"]; // = 1
			//Feigba[1][1]=(*params_json)["Mechanics"]["Feigb_22"]; // 0.9
			//Feigba[1][1]-=(*params_json)["Mechanics"]["Feiga_22"].get<double>(); // -0.1
		//}
		//if(dim==3){
			//Feiga[2][2]=(*params_json)["Mechanics"]["Feiga_33"]; // = 1
			//Feigba[2][2]=(*params_json)["Mechanics"]["Feigb_33"]; // 0.9
			//Feigba[2][2]-=(*params_json)["Mechanics"]["Feiga_33"].get<double>(); // -0.1
		//}
		//double eps_0=1.0e-5;

		//if(mat_id==0){
			//for(unsigned int q=0; q<n_q_points;q++){
				//Sacado::Fad::DFad<double> C_q=this->battery_fields->quad_fields[lithium_index].value[q];
				////std::cout<<"tem"<<(C_q.val()-C_a)/(C_b-C_a) <<std::endl;
				//dealii::Table<2,Sacado::Fad::DFad<double> > Feig(dim,dim);
				//dealii::Table<2,Sacado::Fad::DFad<double>> invFeig(dim,dim);
				//Feig=table_scaling<2,Sacado::Fad::DFad<double>,Sacado::Fad::DFad<double> > (Feigba, (C_q-C_a)/(C_b-C_a) );  // scale * (-0.1) 
				//Feig=table_add<2,Sacado::Fad::DFad<double>,Sacado::Fad::DFad<double> > (Feig, Feiga); // + 2nd order identity
				//getInverse<Sacado::Fad::DFad<double>,dim> (Feig,invFeig);
					//for (unsigned int i=0; i<dim; ++i){
						//for (unsigned int j=0; j<dim; ++j){
							//for (unsigned int k=0; k<dim; ++k){
								 //Fe[q][i][j]+=F[q][i][k]*invFeig[k][j];
							//}
						//}
					//}
				//}
				//this->ResidualEq->evaluateNeoHookeanStress(P, Fe);				
			//}
			//else{this->ResidualEq->evaluateNeoHookeanStress(P, F);}
    //------------------------------------------------------------


    double swell_ratio = 0;
    double C_li_100_neg=(*params_json)["ElectroChemo"]["c_li_100_neg"];
    double C_li_100_pos=(*params_json)["ElectroChemo"]["c_li_100_pos"];
    double C_li_max_neg=(*params_json)["ElectroChemo"]["c_li_max_neg"];
    double C_li_max_pos=(*params_json)["ElectroChemo"]["c_li_max_pos"];

    if (mat_id == li_metal_id)
    {
        double swell_ratio_anode = (*params_json)["Mechanics"]["SwellRatioAnode"];
        swell_ratio = swell_ratio_anode;
    }
    else if (mat_id == active_particle_id)
    {
        double swell_ratio_cathode = (*params_json)["Mechanics"]["SwellRatioCathode"];
        swell_ratio = swell_ratio_cathode;
    }

    double C_0 = 0.0;
	C_0=C_li_100_neg * C_li_max_pos;
	Point<dim> center=(* this->battery_fields->current_cell)->center();
	if (center[orientation]>separator_line){
		C_0=C_li_100_pos * C_li_max_pos;
	}

    if(mat_id==active_particle_id or mat_id == li_metal_id){
      for(unsigned int q=0; q<n_q_points;q++){
        Sacado::Fad::DFad<double> C_q=this->battery_fields->quad_fields[lithium_index].value[q];
        //std::cout << "C_q" << C_q << std::endl;
        //std::cout<<"tem"<<(C_q.val()-C_a)/(C_b-C_a) <<std::endl;
        dealii::Table<2,Sacado::Fad::DFad<double> > Feig(dim,dim);
        dealii::Table<2,Sacado::Fad::DFad<double>> invFeig(dim,dim);

        Feig[0][0]=(C_q - C_0) * swell_ratio + 1.0; 
        if(dim>=2){
          Feig[1][1]=(C_q - C_0) * swell_ratio + 1.0; 
        }
        if(dim==3){
          Feig[2][2]=(C_q - C_0) * swell_ratio + 1.0; 
        }
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
