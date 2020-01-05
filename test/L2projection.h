/*
zhenlin wang 2019
*/
#include "mechanoChemFEM.h"
template <int dim>
class L2projection: public mechanoChemFEM<dim>
{
	public:
		L2projection(std::vector<std::vector<std::string> > _primary_variables, std::vector<std::vector<int> > _FE_support, ParameterHandler& _params);
		L2projection(std::vector<std::vector<std::string> > _primary_variables, std::vector<std::vector<int> > _FE_support,std::vector<std::vector<std::string> > _variables_add, std::vector<std::vector<int> > _FE_support_add, ParameterHandler& _params);
		~L2projection();
		//these are overloaded functions 
		void apply_boundary_condition();
		void get_residual(const typename hp::DoFHandler<dim>::active_cell_iterator &cell, const FEValues<dim>& fe_values, Table<1, Sacado::Fad::DFad<double> >& R, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<1, double >& ULocalConv);
		ParameterHandler* params;		
		
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
L2projection<dim>::L2projection(std::vector<std::vector<std::string> > _primary_variables, std::vector<std::vector<int> > _FE_support, ParameterHandler& _params)
	:mechanoChemFEM<dim>(_primary_variables, _FE_support, _params),params(&_params){

	this->pcout<<"L2 projection initialized"<<std::endl;
}

template <int dim>
L2projection<dim>::L2projection(std::vector<std::vector<std::string> > _primary_variables, std::vector<std::vector<int> > _FE_support,std::vector<std::vector<std::string> > _variables_add, std::vector<std::vector<int> > _FE_support_add, ParameterHandler& _params)
	:mechanoChemFEM<dim>(_primary_variables, _FE_support, _params),variables_add(_variables_add), FE_support_add(_FE_support_add),params(&_params){
  dof_handler_add=new hp::DoFHandler<dim>(this->triangulation);
	this->pcout<<"L2 projection initialized"<<std::endl;
}

template <int dim>
L2projection<dim>::~L2projection()
{
	delete dof_handler_add;
}

template <int dim>
void L2projection<dim>::apply_boundary_condition()
{
	//setup_addtionalSystem();

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
void L2projection<dim>::get_residual(const typename hp::DoFHandler<dim>::active_cell_iterator &cell, const FEValues<dim>& fe_values, Table<1, Sacado::Fad::DFad<double> >& R, Table<1, Sacado::Fad::DFad<double>>& ULocal, Table<1, double >& ULocalConv)
{
	unsigned int n_q_points= fe_values.n_quadrature_points;
	const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
	dealii::Table<1,Sacado::Fad::DFad<double>>  theta(n_q_points),norm(n_q_points);
	evaluateScalarFunction<Sacado::Fad::DFad<double>,dim>(fe_values, 0, ULocal, theta);
	for(unsigned int q=0; q<n_q_points;q++) norm[q]=1-theta[q];
	  
	for (unsigned int i=0; i<dofs_per_cell; ++i) {
		for (unsigned int q=0; q<n_q_points; ++q){
			R[i] +=  fe_values.shape_value(i, q)*norm[q]*fe_values.JxW(q);
		}
	}
}