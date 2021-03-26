#ifndef SDdata_h
#define SDdata_h
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/base/table.h>


// create a data structure to store the strong discontinuity related information
//
using namespace dealii;
template <int dim>
struct SDdata {
  // iso-value/surface related information
  int cell_id;
  int lnode_1;
  int lnode_2;
  bool is_interface_element = false;

  FullMatrix<double> Kcc, Kxic, Kcxi, Kxixi_inv;
  Vector<double> rlocal, xi_old, xi_conv;

  FullMatrix<double> Kcc_c_e, Kxic_c_e, Kcxi_c_e, Kxixi_inv_c_e;
  Vector<double> rlocal_c_e, xi_old_c_e, xi_conv_c_e;

  FullMatrix<double> Kcc_phi_s, Kxic_phi_s, Kcxi_phi_s, Kxixi_inv_phi_s;
  Vector<double> rlocal_phi_s, xi_old_phi_s, xi_conv_phi_s;

  FullMatrix<double> Kcc_phi_e, Kxic_phi_e, Kcxi_phi_e, Kxixi_inv_phi_e;
  Vector<double> rlocal_phi_e, xi_old_phi_e, xi_conv_phi_e;

  Point<dim, double> edge1_node1, edge1_node2, edge1_node;
  Point<dim, double> edge2_node1, edge2_node2, edge2_node;
  Point<dim, double> one_plus_node;
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
  Sacado::Fad::DFad<double> reaction_rate_li = 0.0;
  Sacado::Fad::DFad<double> reaction_rate_potential = 0.0;

  // reaction rate in the opposite_flux_dof will have a negative sign.
  int opposite_flux_dof_li = -1;
  int opposite_flux_dof_potential = -1;

  // iso-value/surface related information

  // update of xi
  std::vector<double> crk_n{1.0, 0.0};  // this should be auto-determined by the iso-surfaces.
  double area_elem = 0.0;
  double interface_length = 0.0;
  double computed_area = 0.0;

  Vector<double> ULocal_k;  // create a larger vector to store previous step values
  Vector<double> C_Li_plus_old;  // enhanced c Li plus at the interface
  Vector<double> C_Li_plus_new;  // 
};
#endif
