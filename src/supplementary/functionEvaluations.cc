/**
 * @defgroup EvaluationFunctions Evaluation Functions
 */
#include"../../include/supplementary/functionEvaluations.h"


/**
 * @ingroup EvaluationFunctions
 */
template <class T, int dim>
void evaluateScalarFunction(const FEValues<dim>& fe_values, unsigned int DOF, Table<1, T>& ULocal, Table<1, T>& U)
{
	unsigned int dofs_per_cell= fe_values.dofs_per_cell;
	unsigned int n_q_points= fe_values.n_quadrature_points;

	//Loop over quadrature points
	for (unsigned int q=0; q<n_q_points; ++q){
		U[q]=0.0; //U
		for (unsigned int k=0; k<dofs_per_cell; ++k){
			if (fe_values.get_fe().system_to_component_index(k).first==DOF){
				U[q]+=ULocal[k]*fe_values.shape_value(k, q); //U
			}
		}
	}
}
template void evaluateScalarFunction<Sacado::Fad::DFad<double>, 1>(const FEValues<1>& fe_values, unsigned int DOF, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<1, Sacado::Fad::DFad<double>>& U);
template void evaluateScalarFunction<Sacado::Fad::DFad<double>, 2>(const FEValues<2>& fe_values, unsigned int DOF, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<1, Sacado::Fad::DFad<double>>& U);
template void evaluateScalarFunction<Sacado::Fad::DFad<double>, 3>(const FEValues<3>& fe_values, unsigned int DOF, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<1, Sacado::Fad::DFad<double>>& U);
template void evaluateScalarFunction<double, 1>(const FEValues<1>& fe_values, unsigned int DOF, Table<1, double>& ULocal, Table<1, double>& U);
template void evaluateScalarFunction<double, 2>(const FEValues<2>& fe_values, unsigned int DOF, Table<1, double>& ULocal, Table<1, double>& U);
template void evaluateScalarFunction<double, 3>(const FEValues<3>& fe_values, unsigned int DOF, Table<1, double>& ULocal, Table<1, double>& U);

/**
*@ingroup EvaluationFunctions
*evaluate value at surface
*/
template <class T, int dim>
void evaluateScalarFunction(const FEValues<dim>& fe_values, const FEFaceValues<dim>& fe_face_values, unsigned int DOF, Table<1, T>& ULocal, Table<1, T>& U)
{
	unsigned int dofs_per_cell= fe_values.dofs_per_cell;
	unsigned int n_q_points= fe_face_values.n_quadrature_points;

	//Loop over quadrature points
	for (unsigned int q=0; q<n_q_points; ++q){
		U[q]=0.0; //U
		for (unsigned int k=0; k<dofs_per_cell; ++k){
			if (fe_values.get_fe().system_to_component_index(k).first==DOF){
				U[q]+=ULocal[k]*fe_face_values.shape_value(k, q); //U
			}
		}
	}
}
template void evaluateScalarFunction<Sacado::Fad::DFad<double>, 1>(const FEValues<1>& fe_values, const FEFaceValues<1>& fe_face_values, unsigned int DOF, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<1, Sacado::Fad::DFad<double>>& U);
template void evaluateScalarFunction<Sacado::Fad::DFad<double>, 2>(const FEValues<2>& fe_values, const FEFaceValues<2>& fe_face_values, unsigned int DOF, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<1, Sacado::Fad::DFad<double>>& U);
template void evaluateScalarFunction<Sacado::Fad::DFad<double>, 3>(const FEValues<3>& fe_values, const FEFaceValues<3>& fe_face_values, unsigned int DOF, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<1, Sacado::Fad::DFad<double>>& U);
template void evaluateScalarFunction<double, 1>(const FEValues<1>& fe_values, const FEFaceValues<1>& fe_face_values, unsigned int DOF, Table<1, double>& ULocal, Table<1, double>& U);
template void evaluateScalarFunction<double, 2>(const FEValues<2>& fe_values, const FEFaceValues<2>& fe_face_values, unsigned int DOF, Table<1, double>& ULocal, Table<1, double>& U);
template void evaluateScalarFunction<double, 3>(const FEValues<3>& fe_values, const FEFaceValues<3>& fe_face_values, unsigned int DOF, Table<1, double>& ULocal, Table<1, double>& U);


/**
*@ingroup EvaluationFunctions
*evaluate scalar function gradient
*/
template <class T, int dim>
void evaluateScalarFunctionGradient(const FEValues<dim>& fe_values, unsigned int DOF, Table<1, T>& ULocal, Table<2, T>& gradU)
{
	unsigned int dofs_per_cell= fe_values.dofs_per_cell;
	unsigned int n_q_points= fe_values.n_quadrature_points;
	//Loop over quadrature points
	for (unsigned int q=0; q<n_q_points; ++q){
		for (unsigned int i=0;i<dim;i++) {gradU[q][i]=0;}
		for (unsigned int k=0; k<dofs_per_cell; ++k){
			if (fe_values.get_fe().system_to_component_index(k).first==DOF){
				for (unsigned int i=0; i<dim; ++i){
					gradU[q][i]+=ULocal[k]*fe_values.shape_grad(k, q)[i]; //gradU
				}
			}	
		}	
	}
}
template void evaluateScalarFunctionGradient<Sacado::Fad::DFad<double>, 1>(const FEValues<1>& fe_values, unsigned int DOF, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<2, Sacado::Fad::DFad<double>>& gradU);
template void evaluateScalarFunctionGradient<Sacado::Fad::DFad<double>, 2>(const FEValues<2>& fe_values, unsigned int DOF, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<2, Sacado::Fad::DFad<double>>& gradU);
template void evaluateScalarFunctionGradient<Sacado::Fad::DFad<double>, 3>(const FEValues<3>& fe_values, unsigned int DOF, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<2, Sacado::Fad::DFad<double>>& gradU);
template void evaluateScalarFunctionGradient<double, 1>(const FEValues<1>& fe_values, unsigned int DOF, Table<1, double>& ULocal, Table<2, double>& gradU);
template void evaluateScalarFunctionGradient<double, 2>(const FEValues<2>& fe_values, unsigned int DOF, Table<1, double>& ULocal, Table<2, double>& gradU);
template void evaluateScalarFunctionGradient<double, 3>(const FEValues<3>& fe_values, unsigned int DOF, Table<1, double>& ULocal, Table<2, double>& gradU);

/**
*@ingroup EvaluationFunctions
*evaluate scalar function gradient
*/
template <class T, int dim>
void evaluateScalarFunctionGradient(const FEValues<dim>& fe_values, unsigned int DOF, Table<1, T>& ULocal, Table<2, T>& gradU, deformationMap<T, dim>& defMap)
{
	unsigned int dofs_per_cell= fe_values.dofs_per_cell;
	unsigned int n_q_points= fe_values.n_quadrature_points;
	Table<1, T> refGradU(dim);
	//Loop over quadrature points
	for (unsigned int q=0; q<n_q_points; ++q){
		for (unsigned int i=0;i<dim;i++) {gradU[q][i]=0;refGradU[i]=0;}
		for (unsigned int k=0; k<dofs_per_cell; ++k){
			if (fe_values.get_fe().system_to_component_index(k).first==DOF){
				for (unsigned int i=0; i<dim; ++i){
					refGradU[i]+=ULocal[k]*fe_values.shape_grad(k, q)[i]; //gradU
				}
			}
			
		}
		//Transform gradient to current configuration. gradW=(F^-T)*GradW
		for (unsigned int i=0; i<dim; ++i){
			for (unsigned int j=0; j<dim; ++j){
				gradU[q][i]+=defMap.invF[q][j][i]*refGradU[j];
			}
		}		
	}
}
template void evaluateScalarFunctionGradient<Sacado::Fad::DFad<double>, 1>(const FEValues<1>& fe_values, unsigned int DOF, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<2, Sacado::Fad::DFad<double>>& gradU, deformationMap<Sacado::Fad::DFad<double>, 1>& defMap);
template void evaluateScalarFunctionGradient<Sacado::Fad::DFad<double>, 2>(const FEValues<2>& fe_values, unsigned int DOF, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<2, Sacado::Fad::DFad<double>>& gradU, deformationMap<Sacado::Fad::DFad<double>, 2>& defMap);
template void evaluateScalarFunctionGradient<Sacado::Fad::DFad<double>, 3>(const FEValues<3>& fe_values, unsigned int DOF, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<2, Sacado::Fad::DFad<double>>& gradU, deformationMap<Sacado::Fad::DFad<double>, 3>& defMap);
template void evaluateScalarFunctionGradient<double, 1>(const FEValues<1>& fe_values, unsigned int DOF, Table<1, double>& ULocal, Table<2, double>& gradU, deformationMap<double, 1>& defMap);
template void evaluateScalarFunctionGradient<double, 2>(const FEValues<2>& fe_values, unsigned int DOF, Table<1, double>& ULocal, Table<2, double>& gradU, deformationMap<double, 2>& defMap);
template void evaluateScalarFunctionGradient<double, 3>(const FEValues<3>& fe_values, unsigned int DOF, Table<1, double>& ULocal, Table<2, double>& gradU, deformationMap<double, 3>& defMap);

/**
*@ingroup EvaluationFunctions
*evaluate scalar function gradient at surface
*/
template <class T, int dim>
void evaluateScalarFunctionGradient(const FEValues<dim>& fe_values, const FEFaceValues<dim>& fe_face_values, unsigned int DOF, Table<1, T>& ULocal, Table<2, T>& gradU)
{
	unsigned int dofs_per_cell= fe_values.dofs_per_cell;
	unsigned int n_q_points= fe_face_values.n_quadrature_points;
	//Loop over quadrature points
	for (unsigned int q=0; q<n_q_points; ++q){
		for (unsigned int i=0;i<dim;i++) {gradU[q][i]=0;}
		for (unsigned int k=0; k<dofs_per_cell; ++k){
			if (fe_values.get_fe().system_to_component_index(k).first==DOF){
				for (unsigned int i=0; i<dim; ++i){
					gradU[q][i]+=ULocal[k]*fe_face_values.shape_grad(k, q)[i]; //gradU
				}
			}
		}		
	}
}
template void evaluateScalarFunctionGradient<Sacado::Fad::DFad<double>, 1>(const FEValues<1>& fe_values, const FEFaceValues<1>& fe_face_values, unsigned int DOF, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<2, Sacado::Fad::DFad<double>>& gradU);
template void evaluateScalarFunctionGradient<Sacado::Fad::DFad<double>, 2>(const FEValues<2>& fe_values, const FEFaceValues<2>& fe_face_values, unsigned int DOF, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<2, Sacado::Fad::DFad<double>>& gradU);
template void evaluateScalarFunctionGradient<Sacado::Fad::DFad<double>, 3>(const FEValues<3>& fe_values, const FEFaceValues<3>& fe_face_values, unsigned int DOF, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<2, Sacado::Fad::DFad<double>>& gradU);
template void evaluateScalarFunctionGradient<double, 1>(const FEValues<1>& fe_values, const FEFaceValues<1>& fe_face_values, unsigned int DOF, Table<1, double>& ULocal, Table<2, double>& gradU);
template void evaluateScalarFunctionGradient<double, 2>(const FEValues<2>& fe_values, const FEFaceValues<2>& fe_face_values, unsigned int DOF, Table<1, double>& ULocal, Table<2, double>& gradU);
template void evaluateScalarFunctionGradient<double, 3>(const FEValues<3>& fe_values, const FEFaceValues<3>& fe_face_values, unsigned int DOF, Table<1, double>& ULocal, Table<2, double>& gradU);

/**
*@ingroup EvaluationFunctions
*evaluate scalar function gradient
*/
template <class T, int dim>
void evaluateScalarFunctionGradient(const FEValues<dim>& fe_values, const FEFaceValues<dim>& fe_face_values, unsigned int DOF, Table<1, T>& ULocal, Table<2, T>& gradU, deformationMap<T, dim>& defMap)
{
	unsigned int dofs_per_cell= fe_values.dofs_per_cell;
	unsigned int n_q_points= fe_face_values.n_quadrature_points;
	Table<1, T> refGradU(dim);
	//Loop over quadrature points
	for (unsigned int q=0; q<n_q_points; ++q){
		for (unsigned int i=0;i<dim;i++) {gradU[q][i]=0;refGradU[i]=0;}
		for (unsigned int k=0; k<dofs_per_cell; ++k){
			if (fe_values.get_fe().system_to_component_index(k).first==DOF){
				for (unsigned int i=0; i<dim; ++i){
					refGradU[i]+=ULocal[k]*fe_face_values.shape_grad(k, q)[i]; //gradU
				}
			}
		}
		//Transform gradient to current configuration. gradW=(F^-T)*GradW
		for (unsigned int i=0; i<dim; ++i){
			for (unsigned int j=0; j<dim; ++j){
				gradU[q][i]+=defMap.invF[q][j][i]*refGradU[j];
			}
		}
	}
}
template void evaluateScalarFunctionGradient<Sacado::Fad::DFad<double>, 1>(const FEValues<1>& fe_values, const FEFaceValues<1>& fe_face_values, unsigned int DOF, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<2, Sacado::Fad::DFad<double>>& gradU, deformationMap<Sacado::Fad::DFad<double>, 1>& defMap);
template void evaluateScalarFunctionGradient<Sacado::Fad::DFad<double>, 2>(const FEValues<2>& fe_values, const FEFaceValues<2>& fe_face_values, unsigned int DOF, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<2, Sacado::Fad::DFad<double>>& gradU, deformationMap<Sacado::Fad::DFad<double>, 2>& defMap);
template void evaluateScalarFunctionGradient<Sacado::Fad::DFad<double>, 3>(const FEValues<3>& fe_values, const FEFaceValues<3>& fe_face_values, unsigned int DOF, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<2, Sacado::Fad::DFad<double>>& gradU, deformationMap<Sacado::Fad::DFad<double>, 3>& defMap);
template void evaluateScalarFunctionGradient<double, 1>(const FEValues<1>& fe_values, const FEFaceValues<1>& fe_face_values, unsigned int DOF, Table<1, double>& ULocal, Table<2, double>& gradU, deformationMap<double, 1>& defMap);
template void evaluateScalarFunctionGradient<double, 2>(const FEValues<2>& fe_values, const FEFaceValues<2>& fe_face_values, unsigned int DOF, Table<1, double>& ULocal, Table<2, double>& gradU, deformationMap<double, 2>& defMap);
template void evaluateScalarFunctionGradient<double, 3>(const FEValues<3>& fe_values, const FEFaceValues<3>& fe_face_values, unsigned int DOF, Table<1, double>& ULocal, Table<2, double>& gradU, deformationMap<double, 3>& defMap);




/**
*@ingroup EvaluationFunctions
*evaluate vector function
*/
template <class T, int dim>
void evaluateVectorFunction(const FEValues<dim>& fe_values, unsigned int DOF, Table<1, T>& ULocal, Table<2, T>& U){
	unsigned int dofs_per_cell= fe_values.dofs_per_cell;
	unsigned int n_q_points= fe_values.n_quadrature_points;
	//Loop over quadrature points
	for (unsigned int q=0; q<n_q_points; ++q){
		for (unsigned int k=0; k<dim; ++k) U[q][k]=0; 
		for (unsigned int k=0; k<dofs_per_cell; ++k){
			unsigned int ck = fe_values.get_fe().system_to_component_index(k).first - DOF;
			if (ck>=0 && ck<dim){
				U[q][ck]+=ULocal[k]*fe_values.shape_value(k, q); //U
			}
		}
	}
}
template void evaluateVectorFunction<Sacado::Fad::DFad<double>, 1>(const FEValues<1>& fe_values, unsigned int DOF, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<2, Sacado::Fad::DFad<double>>& U);
template void evaluateVectorFunction<Sacado::Fad::DFad<double>, 2>(const FEValues<2>& fe_values, unsigned int DOF, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<2, Sacado::Fad::DFad<double>>& U);
template void evaluateVectorFunction<Sacado::Fad::DFad<double>, 3>(const FEValues<3>& fe_values, unsigned int DOF, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<2, Sacado::Fad::DFad<double>>& U);
template void evaluateVectorFunction<double, 1>(const FEValues<1>& fe_values, unsigned int DOF, Table<1, double>& ULocal, Table<2, double>& U);
template void evaluateVectorFunction<double, 2>(const FEValues<2>& fe_values, unsigned int DOF, Table<1, double>& ULocal, Table<2, double>& U);
template void evaluateVectorFunction<double, 3>(const FEValues<3>& fe_values, unsigned int DOF, Table<1, double>& ULocal, Table<2, double>& U);


/**
*@ingroup EvaluationFunctions
*evaluateVectorFunctionGradient at volume quadrature points
*/
template <class T, int dim>
void evaluateVectorFunctionGradient(const FEValues<dim>& fe_values, unsigned int DOF, Table<1, T>& ULocal, Table<3, T>& gradU){
	unsigned int dofs_per_cell= fe_values.dofs_per_cell;
	unsigned int n_q_points= fe_values.n_quadrature_points;
	for (unsigned int q=0; q<n_q_points; ++q){
		for (unsigned int i=0; i<dim; ++i){
			for(unsigned int j=0;j<dim;++j){
					gradU[q][i][j]=0.0; //gradU
			}
		}
		//gradU.fill(0.0);
		//Loop over quadrature points
		for (unsigned int k=0; k<dofs_per_cell; ++k){
			unsigned int ck = fe_values.get_fe().system_to_component_index(k).first - DOF;
			if (ck>=0 && ck<dim){
				for (unsigned int i=0; i<dim; ++i){
					gradU[q][ck][i]+=ULocal[k]*fe_values.shape_grad(k, q)[i]; //gradU
				}
			}
		}
	}
}
template void evaluateVectorFunctionGradient<Sacado::Fad::DFad<double>, 1>(const FEValues<1>& fe_values, unsigned int DOF, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<3, Sacado::Fad::DFad<double>>& gradU);
template void evaluateVectorFunctionGradient<Sacado::Fad::DFad<double>, 2>(const FEValues<2>& fe_values, unsigned int DOF, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<3, Sacado::Fad::DFad<double>>& gradU);
template void evaluateVectorFunctionGradient<Sacado::Fad::DFad<double>, 3>(const FEValues<3>& fe_values, unsigned int DOF, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<3, Sacado::Fad::DFad<double>>& gradU);
template void evaluateVectorFunctionGradient<double, 1>(const FEValues<1>& fe_values, unsigned int DOF, Table<1, double>& ULocal, Table<3, double>& gradU);
template void evaluateVectorFunctionGradient<double, 2>(const FEValues<2>& fe_values, unsigned int DOF, Table<1, double>& ULocal, Table<3, double>& gradU);
template void evaluateVectorFunctionGradient<double, 3>(const FEValues<3>& fe_values, unsigned int DOF, Table<1, double>& ULocal, Table<3, double>& gradU);




/**
*@ingroup EvaluationFunctions
*evaluateVectorFunctionGradient at volume quadrature points at current configuration
*/
template <class T, int dim>
void evaluateVectorFunctionGradient(const FEValues<dim>& fe_values, unsigned int DOF, Table<1, T>& ULocal, Table<3, T>& gradU, deformationMap<T, dim>& defMap)
{
	unsigned int dofs_per_cell= fe_values.dofs_per_cell;
	unsigned int n_q_points= fe_values.n_quadrature_points;
	
	Table<2, T> refGradU(dim,dim);
	
	for (unsigned int q=0; q<n_q_points; ++q){
		for (unsigned int i=0; i<dim; ++i){
			for(unsigned int j=0;j<dim;++j){
					gradU[q][i][j]=0.0; //gradU
					refGradU[i][j]=0.0;
			}
		}
		for (unsigned int k=0; k<dofs_per_cell; ++k){
			unsigned int ck = fe_values.get_fe().system_to_component_index(k).first - DOF;
			if (ck>=0 && ck<dim){
				for (unsigned int i=0; i<dim; ++i){
					refGradU[ck][i]+=ULocal[k]*fe_values.shape_grad(k, q)[i]; //gradU
				}
			}
		}
		
		for(unsigned int k=0;k<dim;k++){
			for (unsigned int i=0; i<dim; ++i){
				for (unsigned int j=0; j<dim; ++j){
					gradU[q][k][i]+=defMap.invF[q][j][i]*refGradU[k][j];
				}
			}
		}
				
	}
	
	
}
template void evaluateVectorFunctionGradient<Sacado::Fad::DFad<double>, 1>(const FEValues<1>& fe_values, unsigned int DOF, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<3, Sacado::Fad::DFad<double>>& gradU,deformationMap<Sacado::Fad::DFad<double>, 1>& defMap);
template void evaluateVectorFunctionGradient<Sacado::Fad::DFad<double>, 2>(const FEValues<2>& fe_values, unsigned int DOF, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<3, Sacado::Fad::DFad<double>>& gradU,deformationMap<Sacado::Fad::DFad<double>, 2>& defMap);
template void evaluateVectorFunctionGradient<Sacado::Fad::DFad<double>, 3>(const FEValues<3>& fe_values, unsigned int DOF, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<3, Sacado::Fad::DFad<double>>& gradU,deformationMap<Sacado::Fad::DFad<double>, 3>& defMap);
template void evaluateVectorFunctionGradient<double, 1>(const FEValues<1>& fe_values, unsigned int DOF, Table<1, double>& ULocal, Table<3, double>& gradU,deformationMap<double, 1>& defMap);
template void evaluateVectorFunctionGradient<double, 2>(const FEValues<2>& fe_values, unsigned int DOF, Table<1, double>& ULocal, Table<3, double>& gradU,deformationMap<double, 2>& defMap);
template void evaluateVectorFunctionGradient<double, 3>(const FEValues<3>& fe_values, unsigned int DOF, Table<1, double>& ULocal, Table<3, double>& gradU,deformationMap<double, 3>& defMap);




/**
*@ingroup EvaluationFunctions
*evaluateVectorFunctionGradient at surface
*/
template <class T, int dim>
void evaluateVectorFunctionGradient(const FEValues<dim>& fe_values, const FEFaceValues<dim>& fe_face_values, unsigned int DOF, Table<1, T>& ULocal, Table<3, T>& gradU)
{
	unsigned int dofs_per_cell= fe_values.dofs_per_cell;
	unsigned int n_q_points= fe_face_values.n_quadrature_points;
	for (unsigned int q=0; q<n_q_points; ++q){
		for (unsigned int i=0; i<dim; ++i){
			for(unsigned int j=0;j<dim;++j){
					gradU[q][i][j]=0.0; //gradU
			}
		}
		//gradU.fill(0.0);
		//Loop over quadrature points
		for (unsigned int k=0; k<dofs_per_cell; ++k){
			unsigned int ck = fe_values.get_fe().system_to_component_index(k).first - DOF;
			if (ck>=0 && ck<dim){
				for (unsigned int i=0; i<dim; ++i){
					gradU[q][ck][i]+=ULocal[k]*fe_face_values.shape_grad(k, q)[i]; //gradU
				}
			}
		}
	}
}
template void evaluateVectorFunctionGradient<Sacado::Fad::DFad<double>, 1>(const FEValues<1>& fe_values, const FEFaceValues<1>& fe_face_values, unsigned int DOF, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<3, Sacado::Fad::DFad<double>>& gradU);
template void evaluateVectorFunctionGradient<Sacado::Fad::DFad<double>, 2>(const FEValues<2>& fe_values, const FEFaceValues<2>& fe_face_values, unsigned int DOF, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<3, Sacado::Fad::DFad<double>>& gradU);
template void evaluateVectorFunctionGradient<Sacado::Fad::DFad<double>, 3>(const FEValues<3>& fe_values, const FEFaceValues<3>& fe_face_values, unsigned int DOF, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<3, Sacado::Fad::DFad<double>>& gradU);
template void evaluateVectorFunctionGradient<double, 1>(const FEValues<1>& fe_values, const FEFaceValues<1>& fe_face_values, unsigned int DOF, Table<1, double>& ULocal, Table<3, double>& gradU);
template void evaluateVectorFunctionGradient<double, 2>(const FEValues<2>& fe_values, const FEFaceValues<2>& fe_face_values, unsigned int DOF, Table<1, double>& ULocal, Table<3, double>& gradU);
template void evaluateVectorFunctionGradient<double, 3>(const FEValues<3>& fe_values, const FEFaceValues<3>& fe_face_values, unsigned int DOF, Table<1, double>& ULocal, Table<3, double>& gradU);


/**
*@ingroup EvaluationFunctions
*evaluateVectorFunctionGradient at surface at current configuration
*/
template <class T, int dim>
void evaluateVectorFunctionGradient(const FEValues<dim>& fe_values, const FEFaceValues<dim>& fe_face_values, unsigned int DOF, Table<1, T>& ULocal, Table<3, T>& gradU, deformationMap<T, dim>& defMap)
{
	unsigned int dofs_per_cell= fe_values.dofs_per_cell;
	unsigned int n_q_points= fe_face_values.n_quadrature_points;
	
	Table<2, T> refGradU(dim,dim);
	
	for (unsigned int q=0; q<n_q_points; ++q){
		for (unsigned int i=0; i<dim; ++i){
			for(unsigned int j=0;j<dim;++j){
					gradU[q][i][j]=0.0; //gradU
					refGradU[i][j]=0.0;
			}
		}
	//gradU.fill(0.0);
	//Loop over quadrature points
		for (unsigned int k=0; k<dofs_per_cell; ++k){
			unsigned int ck = fe_values.get_fe().system_to_component_index(k).first - DOF;
			if (ck>=0 && ck<dim){
				for (unsigned int i=0; i<dim; ++i){
					refGradU[ck][i]+=ULocal[k]*fe_face_values.shape_grad(k, q)[i]; //gradU
				}
			}
		}
		
		for(unsigned int k=0;k<dim;k++){
			for (unsigned int i=0; i<dim; ++i){
				for (unsigned int j=0; j<dim; ++j){
					gradU[q][k][i]+=defMap.invF[q][j][i]*refGradU[k][j];
				}
			}
		}
	}
}
template void evaluateVectorFunctionGradient<Sacado::Fad::DFad<double>, 1>(const FEValues<1>& fe_values, const FEFaceValues<1>& fe_face_values, unsigned int DOF, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<3, Sacado::Fad::DFad<double>>& gradU,deformationMap<Sacado::Fad::DFad<double>, 1>& defMap);
template void evaluateVectorFunctionGradient<Sacado::Fad::DFad<double>, 2>(const FEValues<2>& fe_values, const FEFaceValues<2>& fe_face_values, unsigned int DOF, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<3, Sacado::Fad::DFad<double>>& gradU,deformationMap<Sacado::Fad::DFad<double>, 2>& defMap);
template void evaluateVectorFunctionGradient<Sacado::Fad::DFad<double>, 3>(const FEValues<3>& fe_values, const FEFaceValues<3>& fe_face_values, unsigned int DOF, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<3, Sacado::Fad::DFad<double>>& gradU,deformationMap<Sacado::Fad::DFad<double>, 3>& defMap);
template void evaluateVectorFunctionGradient<double, 1>(const FEValues<1>& fe_values, const FEFaceValues<1>& fe_face_values, unsigned int DOF, Table<1, double>& ULocal, Table<3, double>& gradU,deformationMap<double, 1>& defMap);
template void evaluateVectorFunctionGradient<double, 2>(const FEValues<2>& fe_values, const FEFaceValues<2>& fe_face_values, unsigned int DOF, Table<1, double>& ULocal, Table<3, double>& gradU,deformationMap<double, 2>& defMap);
template void evaluateVectorFunctionGradient<double, 3>(const FEValues<3>& fe_values, const FEFaceValues<3>& fe_face_values, unsigned int DOF, Table<1, double>& ULocal, Table<3, double>& gradU,deformationMap<double, 3>& defMap);



/**
*@ingroup EvaluationFunctions
*get volume deformatonMap at volume quadrature point
*/
template <class T, int dim>
void getDeformationMap(const FEValues<dim>& fe_values, unsigned int DOF, Table<1, T>& ULocal, deformationMap<T, dim>& defMap,unsigned int iteration)
{
  unsigned int n_q_points= fe_values.n_quadrature_points;
	if(n_q_points!=defMap.n_q_points) {printf("number of quadrature points mismatch!"); exit(-1);}
  //evaluate dx/dX
  Table<3, T> gradU(n_q_points, dim, dim);
  evaluateVectorFunctionGradient<T, dim>(fe_values, DOF, ULocal, gradU);

  //Loop over quadrature points
  for (unsigned int q=0; q<n_q_points; ++q){
    Table<2, T > Fq(dim, dim), invFq(dim, dim); T detFq;
    for (unsigned int i=0; i<dim; ++i){
      for (unsigned int j=0; j<dim; ++j){
				defMap.F[q][i][j] = Fq[i][j] = (i==j) + gradU[q][i][j]; //F (as double value)
      }
    }
    getInverse<T, dim>(Fq, invFq, detFq); //get inverse(F)
    defMap.detF[q] = detFq;
    for (unsigned int i=0; i<dim; ++i){
      for (unsigned int j=0; j<dim; ++j){
				defMap.invF[q][i][j] = invFq[i][j];
      }
    }
  }
}
template void getDeformationMap<Sacado::Fad::DFad<double>, 1>(const FEValues<1>& fe_values, unsigned int DOF, Table<1, Sacado::Fad::DFad<double>>& ULocal, deformationMap<Sacado::Fad::DFad<double>, 1>& defMap, unsigned int iteration);
template void getDeformationMap<Sacado::Fad::DFad<double>, 2>(const FEValues<2>& fe_values, unsigned int DOF, Table<1, Sacado::Fad::DFad<double>>& ULocal, deformationMap<Sacado::Fad::DFad<double>, 2>& defMap, unsigned int iteration);
template void getDeformationMap<Sacado::Fad::DFad<double>, 3>(const FEValues<3>& fe_values, unsigned int DOF, Table<1, Sacado::Fad::DFad<double>>& ULocal, deformationMap<Sacado::Fad::DFad<double>, 3>& defMap, unsigned int iteration);
template void getDeformationMap<double, 1>(const FEValues<1>& fe_values, unsigned int DOF, Table<1, double>& ULocal, deformationMap<double, 1>& defMap, unsigned int iteration);
template void getDeformationMap<double, 2>(const FEValues<2>& fe_values, unsigned int DOF, Table<1, double>& ULocal, deformationMap<double, 2>& defMap, unsigned int iteration);
template void getDeformationMap<double, 3>(const FEValues<3>& fe_values, unsigned int DOF, Table<1, double>& ULocal, deformationMap<double, 3>& defMap, unsigned int iteration);


/**
*@ingroup EvaluationFunctions
*get volume deformatonMap at surface quadrature point
*/
template <class T, int dim>
void getDeformationMap(const FEValues<dim>& fe_values, const FEFaceValues<dim>& fe_face_values, unsigned int DOF, Table<1, T>& ULocal, deformationMap<T, dim>& defMap,unsigned int iteration)
{
  unsigned int n_q_points= fe_face_values.n_quadrature_points;
	if(n_q_points!=defMap.n_q_points) {printf("number of quadrature points mismatch!"); exit(-1);}
  //evaluate dx/dX
  Table<3, T> gradU(n_q_points, dim, dim);
  evaluateVectorFunctionGradient<T, dim>(fe_values, fe_face_values, DOF, ULocal, gradU);

  //Loop over quadrature points
  for (unsigned int q=0; q<n_q_points; ++q){
    Table<2, T > Fq(dim, dim), invFq(dim, dim); T detFq;
    for (unsigned int i=0; i<dim; ++i){
      for (unsigned int j=0; j<dim; ++j){
				defMap.F[q][i][j] = Fq[i][j] = (i==j) + gradU[q][i][j]; //F (as double value)
      }
    }
    getInverse<T, dim>(Fq, invFq, detFq); //get inverse(F)
    defMap.detF[q] = detFq;
    for (unsigned int i=0; i<dim; ++i){
      for (unsigned int j=0; j<dim; ++j){
				defMap.invF[q][i][j] = invFq[i][j];
      }
    }	
	}
}
template void getDeformationMap<Sacado::Fad::DFad<double>, 1>(const FEValues<1>& fe_values, const FEFaceValues<1>& fe_face_values, unsigned int DOF, Table<1, Sacado::Fad::DFad<double>>& ULocal, deformationMap<Sacado::Fad::DFad<double>, 1>& defMap, unsigned int iteration);
template void getDeformationMap<Sacado::Fad::DFad<double>, 2>(const FEValues<2>& fe_values, const FEFaceValues<2>& fe_face_values, unsigned int DOF, Table<1, Sacado::Fad::DFad<double>>& ULocal, deformationMap<Sacado::Fad::DFad<double>, 2>& defMap, unsigned int iteration);
template void getDeformationMap<Sacado::Fad::DFad<double>, 3>(const FEValues<3>& fe_values, const FEFaceValues<3>& fe_face_values, unsigned int DOF, Table<1, Sacado::Fad::DFad<double>>& ULocal, deformationMap<Sacado::Fad::DFad<double>, 3>& defMap, unsigned int iteration);
template void getDeformationMap<double, 1>(const FEValues<1>& fe_values, const FEFaceValues<1>& fe_face_values, unsigned int DOF, Table<1, double>& ULocal, deformationMap<double, 1>& defMap, unsigned int iteration);
template void getDeformationMap<double, 2>(const FEValues<2>& fe_values, const FEFaceValues<2>& fe_face_values, unsigned int DOF, Table<1, double>& ULocal, deformationMap<double, 2>& defMap, unsigned int iteration);
template void getDeformationMap<double, 3>(const FEValues<3>& fe_values, const FEFaceValues<3>& fe_face_values, unsigned int DOF, Table<1, double>& ULocal, deformationMap<double, 3>& defMap, unsigned int iteration);


/**
*@ingroup EvaluationFunctions
*/
template <class T, int dim>
  void getDeformationMapTranspose(deformationMap<T, dim>& defMap_new, deformationMap<T, dim>& defMap_old){
 //Loop over quadrature points
	unsigned int n_q_points=defMap_old.n_q_points;
  for (unsigned int q=0; q<n_q_points; ++q){
		defMap_new.detF[q]=defMap_old.detF[q];
		for(unsigned int i=0;i<dim;i++){
			for (unsigned int j=0; j<dim; ++j){
				defMap_new.F[q][i][j]=defMap_old.F[q][j][i];
				defMap_new.invF[q][i][j]=defMap_old.invF[q][i][j];
			}
		}
  }
}

template void getDeformationMapTranspose<Sacado::Fad::DFad<double>, 1>(deformationMap<Sacado::Fad::DFad<double>, 1>& defMap_new, deformationMap<Sacado::Fad::DFad<double>, 1>& defMap_old);
template void getDeformationMapTranspose<Sacado::Fad::DFad<double>, 2>(deformationMap<Sacado::Fad::DFad<double>, 2>& defMap_new, deformationMap<Sacado::Fad::DFad<double>, 2>& defMap_old);
template void getDeformationMapTranspose<Sacado::Fad::DFad<double>, 3>(deformationMap<Sacado::Fad::DFad<double>, 3>& defMap_new, deformationMap<Sacado::Fad::DFad<double>, 3>& defMap_old);
template void getDeformationMapTranspose<double, 1>(deformationMap<double, 1>& defMap_new, deformationMap<double, 1>& defMap_old);
template void getDeformationMapTranspose<double, 2>(deformationMap<double, 2>& defMap_new, deformationMap<double, 2>& defMap_old);
template void getDeformationMapTranspose<double, 3>(deformationMap<double, 3>& defMap_new, deformationMap<double, 3>& defMap_old);




/**
*@ingroup EvaluationFunctions
*/
template <class T, int dim>
void evaluateSpatialgradients(Table<2, T>& gradU, Table<2, T> gradU_spat, deformationMap<T, dim>& defMap)
{
	unsigned int n_q_points= gradU.size(0);
  for (unsigned int q = 0; q < n_q_points; ++q){
		for (unsigned int i=0; i<dim; ++i){
			gradU_spat[q][i]=0;
			for (unsigned int j=0; j<dim; ++j){
				gradU_spat[q][i]+=defMap.invF[q][j][i]*gradU[q][j];
			}
		}
	}
}
template void evaluateSpatialgradients<Sacado::Fad::DFad<double>, 1>(Table<2, Sacado::Fad::DFad<double>>& gradU, Table<2, Sacado::Fad::DFad<double>> gradU_spat, deformationMap<Sacado::Fad::DFad<double>, 1>& defMap);
template void evaluateSpatialgradients<Sacado::Fad::DFad<double>, 2>(Table<2, Sacado::Fad::DFad<double>>& gradU, Table<2, Sacado::Fad::DFad<double>> gradU_spat, deformationMap<Sacado::Fad::DFad<double>, 2>& defMap);
template void evaluateSpatialgradients<Sacado::Fad::DFad<double>, 3>(Table<2, Sacado::Fad::DFad<double>>& gradU, Table<2, Sacado::Fad::DFad<double>> gradU_spat, deformationMap<Sacado::Fad::DFad<double>, 3>& defMap);
template void evaluateSpatialgradients<double, 1>(Table<2, double >& gradU, Table<2, double> gradU_spat, deformationMap<double , 1>& defMap);
template void evaluateSpatialgradients<double, 2>(Table<2, double >& gradU, Table<2, double> gradU_spat, deformationMap<double , 2>& defMap);
template void evaluateSpatialgradients<double, 3>(Table<2, double >& gradU, Table<2, double> gradU_spat, deformationMap<double , 3>& defMap);



/**
*@ingroup EvaluationFunctions
*/
template <class T, int dim>
void evaluateSpatialgradients(Table<3, T>& gradU, Table<3, T> gradU_spat, deformationMap<T, dim>& defMap)
{
	unsigned int n_q_points= gradU.size(0);
  for (unsigned int q = 0; q < n_q_points; ++q){
		for(unsigned int k=0;k<dim;k++){
			for (unsigned int i=0; i<dim; ++i){
				gradU_spat[q][k][i]=0;
				for (unsigned int j=0; j<dim; ++j){
					gradU_spat[q][k][i]+=defMap.invF[q][j][i]*gradU[q][k][j];
				}
			}
		}
	}
}
template void evaluateSpatialgradients<Sacado::Fad::DFad<double>, 1>(Table<3, Sacado::Fad::DFad<double>>& gradU, Table<3, Sacado::Fad::DFad<double>> gradU_spat, deformationMap<Sacado::Fad::DFad<double>, 1>& defMap);
template void evaluateSpatialgradients<Sacado::Fad::DFad<double>, 2>(Table<3, Sacado::Fad::DFad<double>>& gradU, Table<3, Sacado::Fad::DFad<double>> gradU_spat, deformationMap<Sacado::Fad::DFad<double>, 2>& defMap);
template void evaluateSpatialgradients<Sacado::Fad::DFad<double>, 3>(Table<3, Sacado::Fad::DFad<double>>& gradU, Table<3, Sacado::Fad::DFad<double>> gradU_spat, deformationMap<Sacado::Fad::DFad<double>, 3>& defMap);
template void evaluateSpatialgradients<double, 1>(Table<3, double >& gradU, Table<3, double> gradU_spat, deformationMap<double , 1>& defMap);
template void evaluateSpatialgradients<double, 2>(Table<3, double >& gradU, Table<3, double> gradU_spat, deformationMap<double , 2>& defMap);
template void evaluateSpatialgradients<double, 3>(Table<3, double >& gradU, Table<3, double> gradU_spat, deformationMap<double , 3>& defMap);