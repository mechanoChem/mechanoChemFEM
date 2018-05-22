/**
*author zhenlin wang, 2018
*/
#ifndef nodalField_h
#define nodalField_h

#include <deal.II/base/parameter_handler.h>
#include <supplementary/computedField.h>
#include "electricChemo.h"

template <int dim>
class nodalField : public computedField<dim>
{
public:
	nodalField(dealii::ParameterHandler& _params);
	~nodalField();
	
	dealii::ParameterHandler* params;
	void compute_derived_quantities_vector(const std::vector<Vector<double> > &uh,
					       const std::vector<std::vector<Tensor<1,dim> > > &duh,
					       const std::vector<std::vector<Tensor<2,dim> > > &dduh,
					       const std::vector<Point<dim> >                  &normals,
					       const std::vector<Point<dim> >                  &evaluation_points,
					       std::vector<Vector<double> >                    &computed_quantities) const;
								 
	
								 //std::vector<std::string> get_names () const;
								 //std::vector<DataComponentInterpretation::DataComponentInterpretation> get_data_component_interpretation () const;
								 //virtual UpdateFlags get_needed_update_flags () const;
								 //void setupComputedField(std::vector<std::vector<std::string> > _primary_variables);
								 							 
	std::vector<unsigned int > primary_variables_dof;
	//std::vector<std::vector<std::string> > primary_variables;
};

#endif