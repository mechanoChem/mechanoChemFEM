#ifndef FEMdata_h
#define FEMdata_h
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string.h>
#include <ctime>
#include <cmath>
#include <vector>
#include <map>
#include <mpi.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/solution_transfer.h>
#include <deal.II/hp/dof_handler.h>

#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/petsc_parallel_vector.h>

template <int dim, class vectorType>
class FEMdata
{
public:
	FEMdata();
	FEMdata(dealii::hp::DoFHandler<dim>& _dof_handler);
	~FEMdata();
	/**
	*set nodal_solution_names and nodal_data_component_interpretation in deal.ii convention
	*/
	void set_output_name(std::vector<std::vector<std::string> > primary_variables);
	
	/*
	* wrap for dealii::DataOut::clear_data_vectors
	*/
	void clear_data_vectors();
	/**
	*write vtk using solution vector corresponding to this dof_handler into path.
	*path is the full directory with file name e.g. '/home/output/output.vtk'. 
	*/
	void write_vtk(vectorType& solution, std::string path);
	
	/**
	*create snapshot of vector
	*basically do printing of the vector
	*use it with resume_vector_from_snapshot to do restart
	*/
	void create_vector_snapshot(vectorType& U, std::string out_snap_file, unsigned int precision=16, const bool scientific = true, const bool across = true) const;

	/**
	*resume snapshot of vector stored in 'in_snap_file' to vector 'Un' according to this dof_handler.
	*only works when two vector initialized based on same mesh, otherwise dof sequence and partition would not be recoverd
	*/
	void resume_vector_from_snapshot(vectorType& Un, std::string in_snap_file);
	
	/**
	*deal.ii dof_handler
	*/
	dealii::hp::DoFHandler<dim>* dof_handler;
	dealii::DataOut<dim,dealii::hp::DoFHandler<dim> >data_out; 
	
	std::vector<std::string> nodal_solution_names;
	std::vector<dealii::DataComponentInterpretation::DataComponentInterpretation> nodal_data_component_interpretation;
	
	MPI_Comm mpi_communicator;
};

#endif
	