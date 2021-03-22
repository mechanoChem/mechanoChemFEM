#include "mechanoChemFEM.h"
template <int dim>
class plasticity: public mechanoChemFEM<dim>
{
 public:
  plasticity();
  //this is a overloaded function 
  void apply_boundary_condition();
  void get_residual(const typename hp::DoFHandler<dim>::active_cell_iterator &cell, const FEValues<dim>& fe_values, Table<1, Sacado::Fad::DFad<double> >& R, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<1, double >& ULocalConv);
  ParameterHandler* params;		
  ConstraintMatrix* constraints;
};

template <int dim>
plasticity<dim>::plasticity()
{
  //pass the pointer to "constraints" in that defined in mechanoChemFEM
  constraints=this->constraints_mechanoChemFEM;
  //This let you use one params to get all parameters pre-defined in the mechanoChemFEM
  params=this->params_mechanoChemFEM;
  params->enter_subsection("parameters");
  params->declare_entry("youngsModulus","0",Patterns::Double() );
  params->declare_entry("poissonRatio","0",Patterns::Double() );
  params->leave_subsection();	
	
  //Declear the parameters before load it
  this->load_parameters("../parameters.prm");
	
  //define main fields from parameter file.
  this->define_primary_fields();

  //Set up the ibvp.
  this->init_ibvp();
}

template <int dim>
void plasticity<dim>::get_residual(const typename hp::DoFHandler<dim>::active_cell_iterator &cell, const FEValues<dim>& fe_values, Table<1, Sacado::Fad::DFad<double> >& R, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<1, double >& ULocalConv)
{
  //evaluate primary fields
  params->enter_subsection("parameters");
  double youngsModulus=params->get_double("youngsModulus");
  double poissonRatio=params->get_double("poissonRatio");	
  params->leave_subsection();
	
  unsigned int n_q_points= fe_values.n_quadrature_points;
  int u_dof=0;
	  
  //mechanics
  deformationMap<Sacado::Fad::DFad<double>, dim> defMap(n_q_points); 
  getDeformationMap<Sacado::Fad::DFad<double>, dim>(fe_values, u_dof, ULocal, defMap);	
  dealii::Table<3, Sacado::Fad::DFad<double> > P_FpinvT(n_q_points,dim,dim), Cpinv_conv(n_q_points,dim,dim), Cpinv(n_q_points,dim,dim);
  dealii::Table<1, Sacado::Fad::DFad<double> > alpha_conv(n_q_points), alpha(n_q_points);
  for (int q=0; q<n_q_points; ++q){
    alpha[q] = 0.;
    alpha_conv[q] = 0.;
    for (int i=0; i<dim; ++i){
      for (int j=0; j<dim; ++j){
	Cpinv_conv[q][i][j] = (i==j);
	Cpinv[q][i][j] = (i==j);
      }
    }
  }

  //call residual functions
  this->ResidualEq.setLameParametersByYoungsModulusPoissonRatio(youngsModulus, poissonRatio);
  this->ResidualEq.evaluateQuadLogStress(P_FpinvT, defMap.F, Cpinv_conv, alpha_conv, Cpinv, alpha);
  //this->ResidualEq.evaluateNeoHookeanStress(P_FpinvT, defMap.F);
  this->ResidualEq.residualForMechanics(fe_values, u_dof, R, P_FpinvT);	
	
  //BC
  //*
  for (unsigned int faceID=0; faceID<2*dim; faceID++){
    if(cell->face(faceID)->boundary_id()==4 ){
      FEFaceValues<dim> fe_face_values(fe_values.get_fe(), *(this->common_face_quadrature), update_values | update_quadrature_points | update_JxW_values);
      fe_face_values.reinit(cell,faceID);
      this->ResidualEq.residualForNeummanBC(fe_values, fe_face_values, 1, R, -1.e6);
      //this->ResidualEq.residualForNeummanBC(fe_values, fe_face_values, 0, R, -1.e6);
    }
  }
  // */
	
}

//set Dirichlet BC
template <int dim>
void plasticity<dim>::apply_boundary_condition()
{
  this->pcout<<"setup_constraints"<<std::endl;
  constraints->clear ();
  DoFTools::make_hanging_node_constraints (this->dof_handler, *constraints);

  std::vector<bool> fixed (dim, true); 
  std::vector<bool> move (dim, false); move[0]=true; 
  //apply constraints on boundary
  VectorTools:: interpolate_boundary_values (this->dof_handler, 1, ZeroFunction<dim> (dim),*constraints, fixed);
  //VectorTools:: interpolate_boundary_values (this->dof_handler, 4, ConstantFunction<dim> (0.00001,dim),*constraints, move);
	
  constraints->close ();
}

template class InitialConditions<1>;
template class InitialConditions<2>;
template class InitialConditions<3>;
