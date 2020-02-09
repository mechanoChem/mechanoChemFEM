#ifndef Residual_h
#define Residual_h
#include <Sacado.hpp>
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/table.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/lac/vector.h>
#include "supplementary/dataStruct.h"
#include "supplementary/supplementaryFunctions.h"
#include "supplementary/functionEvaluations.h"

template <class T, int dim>
class Residual
{
public:
  Residual ();
  ~Residual();
		
	double currentTime, dt;
	//unsigned int iteration;
	//material properity
	double lambda, mu, viscosity, density;	
	std::string constitutiveModel;
	
	
	void setLameParametersByYoungsModulusPoissonRatio(double youngsModulus, double poissonRatio);
	/**
	*Saint-Venant Kirchhoff constitutive model in 3D
	*/
	double SVK3D(unsigned int i, unsigned int j, unsigned int k, unsigned int l);

	/**
	*Saint-Venant Kirchhoff constitutive model in 2D
	*/
	double SVK2D(unsigned int i, unsigned int j);
		

	/**
	*evaluate strain tensor with infinitesimal_strain_indicator from residualForMechanics
	*/
	void evaluateStrain(dealii::Table<3, T >& Fe, dealii::Table<3, T >& E, deformationMap<T, dim>& defMap, bool infinitesimal_strain_indicator=false);

	/**
	*evaluate stress tensor using Saint-Venant Kirchhoff constitutive model
	*/
	void evaluateSaint_Venant_KirchhoffStress(dealii::Table<3, T >& P, dealii::Table<3, T > &Fe, dealii::Table<3, T >& E);
	
	/**
	*evaluate stress tensor using neoHookean constitutive model
	*/
	void evaluateNeoHookeanStress(dealii::Table<3, T >& P, dealii::Table<3, T > &Fe);
	
	/**
	*assemble residual for mechanics governing equation, \nabla P=0, with infinitesimal_strain_indicator (false as default)
	*/
	void residualForMechanics(const FEValues<dim>& fe_values, unsigned int DOF, Table<1, T >& R,dealii::Table<3, T > P);

	/**
	*assemble residual for Neumman boundary condtion  \bP\cdot\bn=gradn
	*/
	void residualForNeummanBC(const FEValues<dim>& fe_values, const FEFaceValues<dim>& fe_face_values, unsigned int DOF, dealii::Table<1, T >& R, dealii::Table<2,T >& gradn);


	/**
	*assemble residual for Stokes equation with fixed mesh of fluid doamin
	*/
	void residualForStokesEq(const FEValues<dim>& fe_values, unsigned int DOF, Table<1, T >& R, Table<3, T>& gradV, dealii::Table<1, T >& pressure);
	
	/**
	*assemble residual for Stokes equation with moving mesh of fluid doamin stored in deformationMap
	*/
	void residualForStokesEq(const FEValues<dim>& fe_values, unsigned int DOF, Table<1, T >& R, deformationMap<T, dim>& defMap, Table<3, T>& gradV, dealii::Table<1, T >& pressure);
	
	/**
	*assemble residual for Navier_Stokes equation with moving mesh of fluid doamin stored in deformationMap
	*/
	void residualForNavier_StokesEq(const FEValues<dim>& fe_values, unsigned int DOF, Table<1, T >& R, deformationMap<T, dim>& defMap, Table<2, T>& V, dealii::Table<2,double>& V_conv, Table<3, T>& gradV, dealii::Table<1, T >& pressure);
	
	/**
	*evaluate Stokes Stress
	*/
	void evaluateStokesStress(dealii::Table<3, T > P_stoke, const FEValues<dim>& fe_values, Table<3, T>& gradV,dealii::Table<1, T >& pressure);
	
	/**
	*assemble residual for Continuity equation with fixed mesh of fluid doamin
	*/
	void residualForContinuityEq(const FEValues<dim>& fe_values, unsigned int DOF, Table<1, T >& R, Table<1, T >&divVelocity);

	/**
	*assemble residual for Continuity equation with movint mesh of fluid doamin stored in deformationMap
	*/
	void residualForContinuityEq(const FEValues<dim>& fe_values, unsigned int DOF, Table<1, T >& R, deformationMap<T, dim>& defMap, Table<1, T >&divVelocity);

	/**
	*assemble residual for diffusion equation at reference configuration:
	* general fick's law of diffusion
	* \partial c/\partial t + \nabla flux=0;
	*/
  void residualForDiffusionEq(const FEValues<dim>& fe_values, unsigned int DOF, dealii::Table<1, T >& R, dealii::Table<1, T >& c, dealii::Table<1,double>& c_conv, dealii::Table<2,T >& flux);
 
	/**
	*same as above but assemble residual for diffusion equation at current configuration by deformationMap
	*/
	
  void residualForDiffusionEq(const FEValues<dim>& fe_values, unsigned int DOF, dealii::Table<1, T >& R,deformationMap<T, dim>& defMap, dealii::Table<1, T >& c, dealii::Table<1,double>& c_conv, dealii::Table<2,T >& flux);
	
	/**
	*assemble residual for diffusion-reaction equation at reference configuration:
	* general diffusion-reaction equation
	* \partial c/\partial t + \nabla flux=reaction;
	*/
	void residualForDiff_ReacEq(const FEValues<dim>& fe_values, unsigned int DOF, dealii::Table<1, T >& R, dealii::Table<1, T >& c, dealii::Table<1,double>& c_conv, dealii::Table<2,T >& flux, dealii::Table<1, T >& reaction);
	
	/**
	*same as above but assemble residual for diffusion-reaction equation at current configuration by deformationMap
	*/
	void residualForDiff_ReacEq(const FEValues<dim>& fe_values, unsigned int DOF, dealii::Table<1, T >& R, deformationMap<T, dim>& defMap, dealii::Table<1, T >& c, dealii::Table<1,double>& c_conv, dealii::Table<2,T >& flux, dealii::Table<1, T >& reaction);
	/**
	*assemble residual for Poisson equation at reference configuration:
	* \nabla^2 phi=rhs;
	*/
  void residualForPoissonEq(const FEValues<dim>& fe_values, unsigned int DOF, dealii::Table<1,T >& R, dealii::Table<2, T >& phi_grad, dealii::Table<1, T >& rhs);
  
	/**
	*same as before but assemble residual for Poisson equation at current configuration by deformationMap
	*/
	void residualForPoissonEq(const FEValues<dim>& fe_values, unsigned int DOF, dealii::Table<1,T >& R, deformationMap<T, dim>& defMap, dealii::Table<2, T >& phi_grad, dealii::Table<1, T >& rhs);
  /*
	*Poisson equation for vector variables
	*/
  void residualForPoissonEq(const FEValues<dim>& fe_values, unsigned int DOF, dealii::Table<1,T >& R, dealii::Table<3, T >& phi_grad, dealii::Table<2, T >& rhs);
  void residualForPoissonEq(const FEValues<dim>& fe_values, unsigned int DOF, dealii::Table<1,T >& R, deformationMap<T, dim>& defMap, dealii::Table<3, T >& phi_grad, dealii::Table<2, T >& rhs);
	
  /*
	*gneral weak form
	*/
  void residualForGeneralWeakForm(const FEValues<dim>& fe_values, unsigned int DOF, dealii::Table<1,T >& R, dealii::Table<1, T >& scalar_term, dealii::Table<2, T >& vector_term);
  void residualForGeneralWeakForm(const FEValues<dim>& fe_values, unsigned int DOF, dealii::Table<1,T >& R, deformationMap<T, dim>& defMap, dealii::Table<1, T >& scalar_term, dealii::Table<2, T >& vector_term);
	
	/**
	*applying Neumman boundary condition \gradc *n=gradn
	*/
	void residualForNeummanBC(const FEValues<dim>& fe_values, const FEFaceValues<dim>& fe_face_values, unsigned int DOF, dealii::Table<1, T >& R, dealii::Table<1,T >& gradn);
	/**
	*applying Neumman boundary condition \grad c *n=gradn at current configuration
	*/
	void residualForNeummanBC(const FEValues<dim>& fe_values, const FEFaceValues<dim>& fe_face_values, unsigned int DOF, dealii::Table<1, T >& R, deformationMap<T, dim>& defMap, dealii::Table<1,T >& gradn);
	
	/**
	*applying Neumman boundary condition \grad c *n=gradn for constant gradn
	*/
	void residualForNeummanBC(const FEValues<dim>& fe_values, const FEFaceValues<dim>& fe_face_values, unsigned int DOF, dealii::Table<1, T >& R, double gradn);
	
	
	/**
	*assemble constant equation for scalar varibable value=0
	*/
	void residualForEqualityEq(const FEValues<dim>& fe_values, unsigned int DOF, dealii::Table<1, T >& R, dealii::Table<1,T >& value);
	/**
	*assemble constant equation for vector variable 
	*/
	void residualForEqualityEq(const FEValues<dim>& fe_values, unsigned int DOF, dealii::Table<1, T >& R, dealii::Table<2,T >& value);
	
	
	/**
	*scale the residual by scallingFactor
	*/
	void scalling(const FEValues<dim>& fe_values, unsigned int DOF, Table<1, Sacado::Fad::DFad<double> >& R, double scallingFactor);
	
	/**
	*doing integration of scalar variable over volume of the cell and return the value
	*/
	double volumeIntegration(const FEValues<dim>& fe_values, double value_quad);
	
	/**
	*same as above but at current configuration
	*/
	double volumeIntegration(const FEValues<dim>& fe_values, double value_quad, deformationMap<T, dim>& defMap);
	
	
	/**
	*doing integration of vector variable over volume of the cell and return the value
	*/
	double volumeIntegration(const FEValues<dim>& fe_values, Table<1, Sacado::Fad::DFad<double> >& value_quad);
	/**
	*same as above but at current configuration
	*/
	double volumeIntegration(const FEValues<dim>& fe_values, Table<1, Sacado::Fad::DFad<double> >& value_quad, deformationMap<T, dim>& defMap);
	
	/**
	*same as above, but with doulbe type input
	*/
	double volumeIntegration(const FEValues<dim>& fe_values, Table<1, double >& value_quad);
	/**
	*same as above but at current configuration
	*/
	double volumeIntegration(const FEValues<dim>& fe_values, Table<1, double >& value_quad, deformationMap<T, dim>& defMap);
	
	/**
	*doing integration of scalar variable over surface of the cell and return the value
	*/
	double surfaceIntegration(const FEFaceValues<dim>& fe_face_values, double value_quad);
	/**
	*same as above, but at current configuration using Nanson's formula
	*/
	double surfaceIntegration(const FEFaceValues<dim>& fe_face_values, double value_quad, deformationMap<T, dim>& defMap_face);
	
	/**
	*doing integration of vector variable over surface of the cell and return the value
	*/
	double surfaceIntegration(const FEFaceValues<dim>& fe_face_values, dealii::Table<1, Sacado::Fad::DFad<double> >& value_quad);
	/**
	*same as above, but at current configuration using Nanson's formula
	*/
	double surfaceIntegration(const FEFaceValues<dim>& fe_face_values, dealii::Table<1, Sacado::Fad::DFad<double> >& value_quad,deformationMap<T, dim>& defMap_face);
	
	/**
	*same as above, but with doulbe type input
	*/
	double surfaceIntegration(const FEFaceValues<dim>& fe_face_values, dealii::Table<1, double >& value_quad);
	/**
	*same as above, but at current configuration using Nanson's formula
	*/
	double surfaceIntegration(const FEFaceValues<dim>& fe_face_values, dealii::Table<1, double >& value_quad, deformationMap<T, dim>& defMap_face);
	
	
};

#endif