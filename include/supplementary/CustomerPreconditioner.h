#ifndef CustomerPreconditioner_H_
#define CustomerPreconditioner_H_
#include <deal.II/lac/petsc_parallel_sparse_matrix.h>
#include <deal.II/lac/petsc_parallel_vector.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>

using namespace dealii;

template<class matrixType, class vectorType>
class CustomerPreconditioner : public Subscriptor
{
	public:
  	CustomerPreconditioner(const matrixType &A);
  	void vmult (vectorType &dst, vectorType &src) const;
	private:
  const SmartPointer<const matrixType> system_matrix;

};
#endif	