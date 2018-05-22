#include "../../include/FEMdata.h"

template  <int dim, class vectorType>
void FEMdata<dim, vectorType>::create_vector_snapshot(vectorType& U, std::string out_snap_file, unsigned int precision, const bool scientific, const bool across) const
{		
	dealii::Vector<double> Un(U);
	int this_mpi_process;
	MPI_Comm_rank(mpi_communicator, &this_mpi_process);
	if(this_mpi_process == 0){
		std::ofstream myfile;
  	myfile.open (out_snap_file.c_str());
		if(!myfile.is_open()) {std::cout<<"file failed to open!"; exit(1);}
		Un.print(myfile,precision,scientific,across);
	 	myfile.close();
	}
}


template  <int dim, class vectorType>
void FEMdata<dim, vectorType>::resume_vector_from_snapshot(vectorType& Un, std::string in_snap_file)
{
	int this_mpi_process;
	MPI_Comm_rank(mpi_communicator, &this_mpi_process);
	std::vector<double> Un_snapshot;
	if(this_mpi_process == 0){
	
  	std::ifstream ifile(in_snap_file.c_str(), std::ios::in);
 	
  	if (!ifile.is_open()) {
      std::cout << "There was a problem opening the input file!\n";
      exit(1);
  	}
		double num = 0.0;
  	while (ifile >> num) {
      Un_snapshot.push_back(num);		
 	 	}
		if(!(Un_snapshot.size()==dof_handler->n_dofs()))  {
      std::cout << "Un_snapshot.size="<<Un_snapshot.size()<<" is not equal to dof_handler->n_dofs()="<<dof_handler->n_dofs()<<", processor #="<<this_mpi_process<<std::endl;
      exit(1);
  	}
		
	}
	else {Un_snapshot.resize(dof_handler->n_dofs());}
	MPI_Bcast(&Un_snapshot[0], dof_handler->n_dofs(), MPI_DOUBLE,0, mpi_communicator);

	//assign value to Un parallelly 
	Un=0;
  typename dealii::hp::DoFHandler<dim>::active_cell_iterator cell =dof_handler->begin_active(), endc=dof_handler->end();
  for (;cell!=endc; ++cell){	
		if (cell->subdomain_id() ==this_mpi_process){
    	const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
			std::vector<unsigned int> local_dof_indices (dofs_per_cell);
			cell->get_dof_indices (local_dof_indices);
			
			for (unsigned int i=0; i<dofs_per_cell; ++i) {
				Un(local_dof_indices[i])=Un_snapshot[local_dof_indices[i]];
			}	
		}
	}    
	Un.compress(dealii::VectorOperation::insert);
	
}


template class FEMdata<1, dealii::Vector<double> >;
template class FEMdata<2, dealii::Vector<double> >;
template class FEMdata<3, dealii::Vector<double> >;

template class FEMdata<1, dealii::PETScWrappers::MPI::Vector>;
template class FEMdata<2, dealii::PETScWrappers::MPI::Vector>;
template class FEMdata<3, dealii::PETScWrappers::MPI::Vector>;
