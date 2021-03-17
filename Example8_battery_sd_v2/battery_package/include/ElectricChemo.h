/*
zhenlin wang 2020
*ElectricChemo
*/

#ifndef ElectricChemo_h
#define ElectricChemo_h

#include "mechanoChemFEM.h"
#include "Battery_fields.h"
template <int dim, class T>
class ElectricChemo
{
public:
  ElectricChemo ();
	Battery_fields<dim>* battery_fields; 
	void declare_parameters(nlohmann::json& _params);
	void init(Battery_fields<dim>& _battery_fields);	
	/*
	*expression, formula and equations
	*/
  T formula_jn(T Temp, T c_li, T c_li_plus, T phi_s, T phi_e, int domainflag);
	T formula_j0(T c_li, T c_li_plus, int domainflag);
	T formula_Usc(T x, int domainflag);
	
	
	T formula_dUdt(T UnitC);
	
	dealii::Table<1,T > D_li_plus(int type=1);
	dealii::Table<1,T > sigma_e(int type=1);

	dealii::Table<1,double > sigma_e_interface(dealii::Table<1,double> &C_li_plus_q, int type=1);
	dealii::Table<1,double > D_li_plus_interface(dealii::Table<1,double> &C_li_plus_q, int type=1);
	
	

	nlohmann::json* params_ElectricChemo_json;
	
	double F, Rr;
	double k_neg, k_pos, alpha_neg, alpha_pos, c_max_neg, c_max_pos;
	
};

#endif
