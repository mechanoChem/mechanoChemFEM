/**
*author zhenlin wang
*/
#ifndef computedField_h
#define computedField_h

#include <deal.II/base/function.h>
#include <deal.II/lac/vector.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/data_postprocessor.h>
using namespace dealii;

template <int dim>
class computedField : public DataPostprocessor<dim>
{
public:
	computedField ();
	~computedField ();

	/**
	*setup names and variable types:scalar, vector 
	*/
	void setupComputedField(std::vector<std::vector<std::string> > _primary_variables);
	
 std::vector<std::string> get_names () const;
 std::vector<DataComponentInterpretation::DataComponentInterpretation> get_data_component_interpretation () const;
 virtual UpdateFlags get_needed_update_flags () const;
 
 std::vector<std::vector<std::string> > primary_variables;

};

#endif