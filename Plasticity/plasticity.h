#include "mechanoChemFEM.h"
#include <vector>
template <int dim>
class plasticity: public mechanoChemFEM<dim>
{
 public:
  plasticity();
  //this is a overloaded function 
  void apply_boundary_condition();
  void get_residual(const typename hp::DoFHandler<dim>::active_cell_iterator &cell, const FEValues<dim>& fe_values, Table<1, Sacado::Fad::DFad<double> >& R, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<1, double >& ULocalConv);
  void solve_ibvp();
  ParameterHandler* params;		
  ConstraintMatrix* constraints;
  int cell_label;
  std::vector<dealii::Table<3, Sacado::Fad::DFad<double> > > Cpinv;
  std::vector<dealii::Table<3, double > > Cpinv_conv;
  std::vector<dealii::Table<1, Sacado::Fad::DFad<double> > > alpha;
  std::vector<dealii::Table<1, double > > alpha_conv;
  bool plasFlag = false;
  int iterCount = 0;
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

  //Initialize plasticity history variables
  int n_q_points = this->volume_quadrature->size();
  dealii::Table<3, Sacado::Fad::DFad<double> > tmp_Cpinv(n_q_points,dim,dim);
  dealii::Table<3, double > tmp_Cpinv_conv(n_q_points,dim,dim);
  dealii::Table<1, Sacado::Fad::DFad<double> > tmp_alpha(n_q_points);
  dealii::Table<1, double > tmp_alpha_conv(n_q_points);
  for (int q=0; q<n_q_points; ++q){
    tmp_alpha_conv[q] = 0.;
    tmp_alpha[q] = 0.;
    for (int i=0; i<dim; ++i){
      for (int j=0; j<dim; ++j){
	tmp_Cpinv_conv[q][i][j] = (i==j);
	tmp_Cpinv[q][i][j] = (i==j);
      }
    }
  }	
  int cell_label = 0;
  typename hp::DoFHandler<dim>::active_cell_iterator cell = this->dof_handler.begin_active(), endc=this->dof_handler.end();
  for (;cell!=endc; ++cell){
    if (cell->is_locally_owned()){
      Cpinv.push_back(tmp_Cpinv);
      Cpinv_conv.push_back(tmp_Cpinv_conv);
      alpha.push_back(tmp_alpha);
      alpha_conv.push_back(tmp_alpha_conv);
      cell->set_user_index(cell_label);
      cell_label++;
    }
  }
  
}

template <int dim>
void plasticity<dim>::get_residual(const typename hp::DoFHandler<dim>::active_cell_iterator &cell, const FEValues<dim>& fe_values, Table<1, Sacado::Fad::DFad<double> >& R, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<1, double >& ULocalConv)
{
  
  //evaluate primary fields
  params->enter_subsection("parameters");
  double youngsModulus=params->get_double("youngsModulus");
  double poissonRatio=params->get_double("poissonRatio");
  double tau_y = 500.e6;
  params->leave_subsection();
	
  unsigned int n_q_points= fe_values.n_quadrature_points;
  int u_dof=0;
  int alpha_dof=dim;
	  
  //mechanics
  deformationMap<Sacado::Fad::DFad<double>, dim> defMap(n_q_points); 
  getDeformationMap<Sacado::Fad::DFad<double>, dim>(fe_values, u_dof, ULocal, defMap);
  for (int i=0; i<ULocal.size(0); ++i){
    if(std::isnan(ULocal[i].val())){
      std::cout << "nan in Ulocal\n" << std::endl;
      for (int k=0; k<ULocal.size(0); ++k){
	std::cout << ULocal[k].val() << " ";
      }
      std::cout << std::endl;
      break;
    }
  }
  /*
  for (int q=0; q<n_q_points; ++q){
    for (int i=0; i<dim; ++i){
      for (int j=0; j<dim; ++j){
	if(std::isnan(defMap.F[q][i][j].val())){
	  std::cout << "nan in f\n";
	  for (int k=0; k<ULocal.size(0); ++k){
	    std::cout << ULocal[k] << " " << std::endl;
	  }
	  std::cout << std::endl;
	}
      }
    }
  }
  */
  dealii::Table<3, Sacado::Fad::DFad<double> > P_FpinvT(n_q_points,dim,dim);
  dealii::Table<1, Sacado::Fad::DFad<double> > proj_alpha(n_q_points), local_alpha(n_q_points);
  evaluateScalarFunction<Sacado::Fad::DFad<double>,dim>(fe_values, alpha_dof, ULocal, local_alpha);

  //call residual function
  int cell_label = cell->user_index();
  this->ResidualEq.setLameParametersByYoungsModulusPoissonRatio(youngsModulus, poissonRatio);
  this->ResidualEq.evaluateQuadLogStress(P_FpinvT, defMap.F, Cpinv_conv[cell_label], alpha_conv[cell_label], Cpinv[cell_label], alpha[cell_label], tau_y, 0);
  /*  
  if (this->current_increment == 1){
    this->ResidualEq.evaluateQuadLogStress(P_FpinvT, defMap.F, Cpinv_conv[cell_label], alpha_conv[cell_label], Cpinv[cell_label], alpha[cell_label], tau_y, this->currentIteration);
  }
  else{
    this->ResidualEq.evaluateQuadLogStress(P_FpinvT, defMap.F, Cpinv_conv[cell_label], alpha_conv[cell_label], Cpinv[cell_label], alpha[cell_label], tau_y, 2);
  }
  */
  for (int q=0; q<n_q_points; ++q){
    if (alpha[cell_label][q] > 0 && !plasFlag){
      this->pcout << "Onset of plasticity\n\n";
      plasFlag = true;
    }
    proj_alpha[q] = local_alpha[q] - alpha[cell_label][q];
  }
  //this->ResidualEq.evaluateNeoHookeanStress(P_FpinvT, defMap.F);
  this->ResidualEq.residualForMechanics(fe_values, u_dof, R, P_FpinvT);	
  this->ResidualEq.residualForEqualityEq(fe_values, alpha_dof, R, proj_alpha);	
	
  //BC
  /*
  for (unsigned int faceID=0; faceID<2*dim; faceID++){
    if(cell->face(faceID)->boundary_id()==4 ){
      FEFaceValues<dim> fe_face_values(fe_values.get_fe(), *(this->common_face_quadrature), update_values | update_quadrature_points | update_JxW_values);
      fe_face_values.reinit(cell,faceID);
      //this->ResidualEq.residualForNeummanBC(fe_values, fe_face_values, 1, R, -1.e6);
      this->ResidualEq.residualForNeummanBC(fe_values, fe_face_values, 0, R, -1.e8*this->current_increment);
    }
  }
  // */
	
}


template <int dim>
void plasticity<dim>::solve_ibvp(){

  int reduce_time=0;
  int converge_iter;
  bool converged;
  while(reduce_time<8) {
  
    plasFlag = false;
  
    // apply nonhomogeneous Dirichlet boundary conditions
    // NOTE: this is only half of the process
    std::vector<types::global_dof_index> local_face_dof_indices (this->fe_system[0]->dofs_per_face);  
    typename hp::DoFHandler<dim>::active_cell_iterator cell = this->dof_handler.begin_active(), endc=this->dof_handler.end();
    for (;cell!=endc; ++cell){
      for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f){
	if (cell->face(f)->boundary_id()==4){
	  cell->face(f)->get_dof_indices (local_face_dof_indices,0);
	  for (unsigned int i=0; i<local_face_dof_indices.size(); ++i){
	    const unsigned int ck = this->fe_system[0]->face_system_to_component_index(i).first;
	    if(ck==0 and this->solution.in_local_range(local_face_dof_indices[i])){
	      this->solution(local_face_dof_indices[i]) = 0.001*(this->current_time+this->current_dt)/this->total_time; //add given value for x_displacement dof
	    }
	  }
	}
      }
    }
    this->solution.compress(VectorOperation::insert); // Since it is a distributed object
    
    converged = this->nonlinearSolve(this->solution);

    if (!converged) {
      iterCount = 0;
      this->pcout << "not converged; reduce loading by half\n";
      this->current_dt *= 0.5;
      this->solution=this->solution_prev;
      reduce_time++;
    }
    else
      break;
  }
  if (this->currentIteration < 4){
    iterCount++;
  }
  else
    iterCount = 0;
  if (iterCount > 3){
    iterCount = 0;
    this->current_dt *= 2.;
  }

  /*
  std::cout << std::endl;
  for (int i=0; i<this->solution.size(); ++i){
    std::cout << this->solution[i] << " ";
  }
  std::cout << std::endl;
  */
  
  //update solution
  this->solution_prev=this->solution;
  //*
  //update plasticity history variables
  for (int cell=0; cell<alpha.size(); ++cell){
    for (int q=0; q<alpha[cell].size(0); ++q){
      alpha_conv[cell][q] = alpha[cell][q].val();
      for (int i=0; i<dim; ++i){
	for (int j=0; j<dim; ++j){
	  Cpinv_conv[cell][q][i][j] = Cpinv[cell][q][i][j].val();
	}
      }
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
  //apply constraints on boundary
  //VectorTools:: interpolate_boundary_values (this->dof_handler, 1, ZeroFunction<dim> (dim),*constraints, fixed);
  
  std::vector<types::global_dof_index> local_face_dof_indices (this->fe_system[0]->dofs_per_face);  
  typename hp::DoFHandler<dim>::active_cell_iterator cell = this->dof_handler.begin_active(), endc=this->dof_handler.end();
  for (;cell!=endc; ++cell){
    for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f){
      if (cell->face(f)->boundary_id()==4){
	cell->face(f)->get_dof_indices (local_face_dof_indices,0);
	for (unsigned int i=0; i<local_face_dof_indices.size(); ++i){
	  const unsigned int ck = this->fe_system[0]->face_system_to_component_index(i).first;
	  if(ck==0 || ck==1){
	    constraints->add_line (local_face_dof_indices[i]); //add constraint line for x_displacement dof
	  }
	}
      }
      else if (cell->face(f)->boundary_id()==1){
	cell->face(f)->get_dof_indices (local_face_dof_indices,0);
	for (unsigned int i=0; i<local_face_dof_indices.size(); ++i){
	  constraints->add_line (local_face_dof_indices[i]); //add fixed (zero) constraint for all dofs
	}
      }
    }
  }
  constraints->close ();
}

template class InitialConditions<1>;
template class InitialConditions<2>;
template class InitialConditions<3>;
