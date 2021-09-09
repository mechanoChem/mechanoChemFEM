/*
zhenlin wang 2019
*module transportation
*/
#include "../include/Mechanics.h"

template <int dim>
Mechanics<dim>::Mechanics(){}

template <int dim>
Mechanics<dim>::Mechanics(Battery_fields<dim>& _fields, Residual<Sacado::Fad::DFad<double>,dim>& _ResidualEq)
{
	battery_fields=&_fields;
	ResidualEq=&_ResidualEq;
}

template <int dim>
Mechanics<dim>::Mechanics(Battery_fields<dim>& _fields,Residual<Sacado::Fad::DFad<double>,dim>& _ResidualEq, int _primiary_dof)
{
	battery_fields=&_fields;
	ResidualEq=&_ResidualEq;
	primiary_dof=_primiary_dof;
}

template <int dim>
void Mechanics<dim>::set_up_fields(Battery_fields<dim>& _battery_fields, Residual<Sacado::Fad::DFad<double>,dim>& _ResidualEq, int _primiary_dof)
{
	battery_fields=&_battery_fields;
	ResidualEq=&_ResidualEq;
	primiary_dof=_primiary_dof;
}

template <int dim>
void Mechanics<dim>::set_up_fields(Battery_fields<dim>& _battery_fields, ElectricChemo<dim, Sacado::Fad::DFad<double>>& _electricChemoFormula, Residual<Sacado::Fad::DFad<double>,dim>& _ResidualEq, int _primiary_dof)
{
	battery_fields=&_battery_fields;
	electricChemoFormula=&_electricChemoFormula;
	ResidualEq=&_ResidualEq;
	primiary_dof=_primiary_dof;
}

template <int dim>
void Mechanics<dim>::set_primiary_dof(int _primiary_dof)
{
	primiary_dof=_primiary_dof;
}

template <int dim>
void Mechanics<dim>::set_cell_id(int _cell_id)
{
	cell_id=_cell_id;
}

template <int dim>
void Mechanics<dim>::r_get_residual(const FEValues<dim>& fe_values, Table<1, Sacado::Fad::DFad<double> >& R, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<1, double >& ULocalConv, std::vector<std::vector<double>> &pressure)
{
//evaluate primary fields	
	ResidualEq->setLameParametersByYoungsModulusPoissonRatio(youngsModulus, poissonRatio);
	unsigned int n_q_points= fe_values.n_quadrature_points;
	dealii::Table<3, Sacado::Fad::DFad<double> >Fe(n_q_points,dim,dim);
	deformationMap<Sacado::Fad::DFad<double>, dim> defMap(n_q_points); 
	getDeformationMap<Sacado::Fad::DFad<double>, dim>(fe_values, primiary_dof, ULocal, defMap);
	dealii::Table<3, Sacado::Fad::DFad<double> > P(n_q_points,dim,dim);
	set_stress(defMap.F, P);
  
  // compute pressure
	dealii::Table<2, Sacado::Fad::DFad<double> > sigma(dim,dim);

  pressure[cell_id].resize(n_q_points);
  for (unsigned int q = 0; q < n_q_points; ++q){
		for(unsigned int i=0;i<dim;i++){
			for (unsigned int j=0; j<dim; ++j){
        sigma[i][j] = 0.0;
				for (unsigned int k=0; k<dim; ++k){
					sigma[i][j]+=P[q][i][k]*defMap.F[q][j][k]/defMap.detF[q];
				} // k
			} // j
		} //i
    double _pressure = 0.0;
		for(unsigned int i=0;i<dim;i++){
      _pressure += 1.0/3.0 * sigma[i][i].val();
    }
    pressure[cell_id][q] = _pressure;
    //std::cout << "pressure: " << _pressure << " q " << q << std::endl;
	} // q

	  
	//chemo
	ResidualEq->residualForMechanics(fe_values, primiary_dof, R, P);
	
}

template <int dim>
void Mechanics<dim>::set_stress(dealii::Table<3,Sacado::Fad::DFad<double> >& F, dealii::Table<3,Sacado::Fad::DFad<double> >& P)
{
	ResidualEq->evaluateNeoHookeanStress(P, F);
}


template class Mechanics<1>;
template class Mechanics<2>;
template class Mechanics<3>;
