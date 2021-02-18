#ifndef SDdata_h
#define SDdata_h
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>

// create a data structure to store the strong discontinuity related information
//
template <int dim>
struct SDdata {
  // iso-value/surface related information
  int cell_id;
  int lnode_1;
  int lnode_2;
  bool is_interface_element = false;
  FullMatrix<double> Kcc, Kxic, Kcxi, Kxixi;
  Vector<double> rlocal, xi_old;
  
  Point<dim> edge1_node1, edge1_node2, edge1_node ;
  Point<dim> edge2_node1, edge2_node2, edge2_node ;
  double edge1_local_s = 0.0, edge2_local_s = 0.0;

  // plus side: value >= iso-value
  std::vector<int> lnode_plus;
  // minus side: value < iso-value
  std::vector<int> lnode_minus;

  // 1-d shape function information
  FullMatrix<double> shape_value_1d;
  Vector<double> jxw_1d;

  // Reaction rate flow into the plus side (active particle) is positive. Thus it is flow out relative to the electrolyte. The reaction rate would be negative
  // If particles flow output of the plus side (active particle), the value will be negative, and the value would be positive for the electrolyte.
  double reaction_rate = 0.0;

  // reaction rate in the opposite_flux_dof will have a negative sign.
  int opposite_flux_dof = -1;

  // iso-value/surface related information
};
#endif
