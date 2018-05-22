#include "solve_test_PETSc.h"


template <int dim>
void solve_test_PETSc<dim>::updateLinearSystem()
{
	this->reinitLinearSystem();
	
	for(unsigned int i=0;i<n_local_dofs;i++)
	{
		int d=this_mpi_process*n_local_dofs+i;
		this->system_rhs(d)=-(solution(d)+solution(d)*solution(d));
		this->system_matrix.set(d,d,1+2*solution(d));
		std::cout<<d<<" solution="<<solution(d)<<std::endl;
	}
	
	
	this->system_rhs.compress(VectorOperation::insert);
	this->system_matrix.compress(VectorOperation::insert);
		
	this->system_matrix.print(std::cout);
}

template class solve_test_PETSc<1>;
template class solve_test_PETSc<2>;
template class solve_test_PETSc<3>;