/*
zhenlin wang 2019
*/
#include "mechanoChemFEM.h"
template <int dim>
class L2projection: public mechanoChemFEM<dim>
{
	public:
		L2projection(std::vector<std::vector<std::string> > _primary_variables, std::vector<std::vector<int> > _FE_support,std::vector<std::vector<std::string> > _variables_add, std::vector<std::vector<int> > _FE_support_add, ParameterHandler& _params);
		~L2projection();
		//these are overloaded functions 
		void updateLinearSystem();
		void apply_boundary_condition();
		ParameterHandler* params;			
		/*
		--------------------------------------------------------------------------------------------
		* FEM for additional data 
		--------------------------------------------------------------------------------------------
		*/
		void setup_addtionalSystem();		
		std::vector<std::shared_ptr<FESystem<dim>> > fe_system_add;		
    hp::FECollection<dim> fe_collection_add;
    hp::QCollection<dim>  q_collection_add;
		hp::DoFHandler<dim>*   dof_handler_add;
		
		std::vector<std::vector<std::string> > variables_add;
		std::vector<std::vector<int> > FE_support_add;
		
		PETScWrappers::MPI::Vector additional_data;
};

template <int dim>
L2projection<dim>::L2projection(std::vector<std::vector<std::string> > _primary_variables, std::vector<std::vector<int> > _FE_support,std::vector<std::vector<std::string> > _variables_add, std::vector<std::vector<int> > _FE_support_add, ParameterHandler& _params)
	:mechanoChemFEM<dim>(_primary_variables, _FE_support, _params),variables_add(_variables_add), FE_support_add(_FE_support_add),params(&_params){
  dof_handler_add=new hp::DoFHandler<dim>(this->triangulation);
	this->pcout<<"L2 projection initialized"<<std::endl;
}
template <int dim>
L2projection<dim>::~L2projection(){delete dof_handler_add; }

template <int dim>
void L2projection<dim>::apply_boundary_condition()
{
	setup_addtionalSystem();
	
	hpFEM<dim>::constraints.clear ();
  DoFTools::make_hanging_node_constraints (hpFEM<dim>::dof_handler, hpFEM<dim>::constraints);
  hpFEM<dim>::constraints.close ();  
}



template <int dim>
void L2projection<dim>::setup_addtionalSystem()
{
	//adjoint Fe_system with same volume_quadrature	
  this->pcout<<"setup_addtionalSystem"<<std::endl;

	std::vector<unsigned int > variables_dof_tem;
	this->setup_FeSystem(fe_system_add,fe_collection_add, q_collection_add, variables_dof_tem,variables_add,FE_support_add,*(this->volume_quadrature) );
	dof_handler_add->distribute_dofs (fe_collection_add);
	DoFRenumbering::component_wise (*dof_handler_add);
	const types::global_dof_index n_local_dofs = DoFTools::count_dofs_with_subdomain_association(*dof_handler_add, this->this_mpi_process);
	const types::global_dof_index n_total_dofs=dof_handler_add->n_dofs();
										
	additional_data.reinit (this->mpi_communicator,n_total_dofs,n_local_dofs); 	
	this->pcout<<"setup_addtionalSystem done!"<<std::endl;
}

template <int dim>
void L2projection<dim>::updateLinearSystem()
{
  this->pcout<<"updateLinearSystem"<<std::endl;	
	params->enter_subsection("Problem");
	std::string additional_data_snapshot=params->get("additional_data_snapshot");
	params->leave_subsection();	
	FEMdata<dim,PETScWrappers::MPI::Vector> FEMdata_out_add(*dof_handler_add);
	FEMdata_out_add.set_output_name(variables_add);
	FEMdata_out_add.resume_vector_from_snapshot(additional_data,additional_data_snapshot+"snapshot-"+std::to_string(this->current_increment)+".dat");
	
	//initialize 
	this->reinitLinearSystem();
	
  hp::FEValues<dim> hp_fe_values (this->fe_collection, this->q_collection, update_values | update_quadrature_points  | update_JxW_values | update_gradients);  
  hp::FEValues<dim> hp_fe_values_add (fe_collection_add, q_collection_add, update_values | update_quadrature_points  | update_JxW_values | update_gradients);  
 
  FullMatrix<double> local_matrix;
  Vector<double>            local_rhs;
	local_matrix = 0; local_rhs = 0; 
  std::vector<types::global_dof_index> local_dof_indices;
	std::vector<types::global_dof_index> local_dof_indices_add;
  //loop over cells
  PETScWrappers::Vector localized_U(this->solution);
	 PETScWrappers::Vector localized_Un(this->solution_prev);
	PETScWrappers::Vector localized_U_add(additional_data);
	
	this->ResidualEq.dt=this->current_dt;		
  typename hp::DoFHandler<dim>::active_cell_iterator cell = hpFEM<dim>::dof_handler.begin_active(), endc=hpFEM<dim>::dof_handler.end();
	typename hp::DoFHandler<dim>::active_cell_iterator cell_add = dof_handler_add->begin_active();
	for (;cell!=endc; ++cell, ++cell_add){
		if (cell->subdomain_id() == this->this_mpi_process){
			hp_fe_values.reinit (cell);
			hp_fe_values_add.reinit (cell_add);
			
			const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
    	const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values();
			
			const unsigned int dofs_per_cell_add = cell_add->get_fe().dofs_per_cell;
    	const FEValues<dim> &fe_values_add = hp_fe_values_add.get_present_fe_values();
			
			//n_q_points should be same for both dof_handler!
    	unsigned int n_q_points= fe_values.n_quadrature_points;
    	local_matrix.reinit (dofs_per_cell, dofs_per_cell);
    	local_rhs.reinit (dofs_per_cell);
    	local_dof_indices.resize (dofs_per_cell);
			cell->get_dof_indices (local_dof_indices);
			
    	local_dof_indices_add.resize (dofs_per_cell_add);
			cell_add->get_dof_indices (local_dof_indices_add);
			
    	//AD variables
    	Table<1, Sacado::Fad::DFad<double> > ULocal(dofs_per_cell);
			Table<1, double > ULocalConv(dofs_per_cell);
			
			Table<1, double > ULocal_add(dofs_per_cell_add);
			
	  	for (unsigned int i=0; i<dofs_per_cell; ++i){
				if (std::abs(localized_U(local_dof_indices[i]))<1.0e-16) ULocal[i]=0.0;
				else{ULocal[i]=localized_U(local_dof_indices[i]);}
				ULocal[i].diff (i, dofs_per_cell);
				ULocalConv[i]= localized_Un(local_dof_indices[i]);
	  	}
			
	  	for (unsigned int i=0; i<dofs_per_cell_add; ++i){
				if (std::abs(localized_U_add(local_dof_indices_add[i]))<1.0e-16) ULocal_add[i]=0.0;
				else{ULocal_add[i]=localized_U_add(local_dof_indices_add[i]);}
	  	}
		
    	Table<1, Sacado::Fad::DFad<double> > R(dofs_per_cell); 
			for(unsigned int i=0;i<dofs_per_cell;i++) R[i]=0.0;
			
			dealii::Table<1,double>  c1(n_q_points),c2(n_q_points);
			evaluateScalarFunction<double,dim>(fe_values_add, 0, ULocal_add, c1);
			evaluateScalarFunction<double,dim>(fe_values_add, 2, ULocal_add, c2);
			
			dealii::Table<1,Sacado::Fad::DFad<double>>  theta(n_q_points),norm(n_q_points);
			evaluateScalarFunction<Sacado::Fad::DFad<double>,dim>(fe_values, 0, ULocal, theta);
			
			
			for(unsigned int q=0; q<n_q_points;q++){
				if (c2[q]+0.866*c1[q] > 0 and c1[q] >= 0) norm[q]=1;
				else if (c2[q]-0.866*c1[q]>= 0 and c1[q] < 0) norm[q]=0;
				else norm[q]=-1;
				norm[q]=norm[q]-theta[q];
			}
			
			///this->ResidualEq.residualForGeneralWeakForm(fe_values, 0, R,norm,vector_term);
		  for (unsigned int i=0; i<dofs_per_cell; ++i) {
		    for (unsigned int q=0; q<n_q_points; ++q){
					R[i] +=  fe_values.shape_value(i, q)*norm[q]*fe_values.JxW(q);
				}
			}
    	//Residual(R) and Jacobian(R')		
    	for (unsigned int i=0; i<dofs_per_cell; ++i) {
      	for (unsigned int j=0; j<dofs_per_cell; ++j){
					// R' by AD
					local_matrix(i,j)= R[i].dx(j);
      	}//g
      	//R
      	local_rhs(i) = -R[i].val(); 
   	  }
			this->distribute_local_to_global(local_matrix, local_rhs, local_dof_indices);
		}
	}
	this->LinearSystemCompressAdd();
	this->pcout<<"updateLinearSystem done"<<std::endl;
}




