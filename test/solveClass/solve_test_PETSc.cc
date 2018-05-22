#include "solve_test_PETSc.h"


template <int dim>
solve_test_PETSc<dim>::solve_test_PETSc()
{}


template <int dim>
solve_test_PETSc<dim>::~solve_test_PETSc (){}

template class solve_test_PETSc<1>;
template class solve_test_PETSc<2>;
template class solve_test_PETSc<3>;
