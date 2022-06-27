#include "battery.h"
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>
#include <iostream>
#include <string.h>
#include <vector>
#include <fstream>

using namespace dealii;
using namespace std;

// add other types ...

template <class T>
void write_std_vector(ostream& fo, std::vector<T> &std_vector0)
{
  fo.write(reinterpret_cast<char*>(&std_vector0[0]), std_vector0.size()* sizeof(T));
}

template <class T>
void read_std_vector(istream& fi, std::vector<T> &std_vector0)
{
  fi.read(reinterpret_cast<char*>(&std_vector0[0]), std_vector0.size()*sizeof(T));
}

// not working 
//template <class T>
//void write_std_vector_vector(ostream& fo, std::vector< std::vector<T> > &std_vector0)
//{
  //fo.write(reinterpret_cast<char*>(&std_vector0[0][0]), std_vector0.size()*std_vector0[0].size()*sizeof(T));
//}

//template <class T>
//void read_std_vector_vector(istream& fi, std::vector< std::vector<T> > &std_vector0)
//{
  //fi.read(reinterpret_cast<char*>(&std_vector0[0][0]), std_vector0.size()*std_vector0[0].size()*sizeof(T));
//}

void write_vector(ostream& fo, Vector<double> &vector0)
{
  fo.write(reinterpret_cast<char*>(&vector0(0)), vector0.size()* sizeof(double));
}

void read_vector(istream& fi, Vector<double> &vector0)
{
  fi.read(reinterpret_cast<char*>(&vector0(0)), vector0.size()*sizeof(double));
}

void write_full_matrix(ostream& fo, FullMatrix<double> &full_matrix0)
{
  fo.write(reinterpret_cast<char*>(&full_matrix0(0,0)), full_matrix0.size()[0] * full_matrix0.size()[1] * sizeof(double));
}

void read_full_matrix(istream& fi, FullMatrix<double> &full_matrix0)
{
  fi.read(reinterpret_cast<char*>(&full_matrix0(0,0)), full_matrix0.size()[0] * full_matrix0.size()[1] *sizeof(double));
}

void write_bool(ostream& fo, bool &bool0)
{
  fo.write(reinterpret_cast<char*>(&bool0), sizeof(bool));
}

void read_bool(istream& fi, bool &bool0)
{
  fi.read(reinterpret_cast<char*>(&bool0), sizeof(bool));
}

void write_int(ostream& fo, int &int0)
{
  fo.write(reinterpret_cast<char*>(&int0), sizeof(int));
}

void read_int(istream& fi, int &int0)
{
  fi.read(reinterpret_cast<char*>(&int0), sizeof(int));
}

void write_double(ostream& fo, double &double0)
{
  fo.write(reinterpret_cast<char*>(&double0), sizeof(double));
}

void read_double(istream& fi, double &double0)
{
  fi.read(reinterpret_cast<char*>(&double0), sizeof(double));
}

//void write_sacado_double(ostream& fo, Sacado::Fad::DFad<double> &sacado_double0)
//{
  //fo.write(reinterpret_cast<char*>(&sacado_double0), sizeof(Sacado::Fad::DFad<double>));
//}

//void read_sacado_double(istream& fi, Sacado::Fad::DFad<double> &sacado_double0)
//{
  //fi.read(reinterpret_cast<char*>(&sacado_double0), sizeof(Sacado::Fad::DFad<double>));
//}

template <int dim>
void battery<dim>::save_sd_data()
{
    std::ofstream fout(this->snapshot_directory+"restart-"+std::to_string(this->current_increment +this->off_output_index)+".dat" + std::to_string(this->this_mpi_process), std::ios::out | std::ios::binary);
    for (unsigned int count = 0; count < cell_SDdata.size(); count++)
    {
      if (cell_SDdata[count].is_interface_element)
      {
        //std::cout << " mpi " << this->this_mpi_process << " cell_id " << cell_SDdata[count].cell_id << std::endl;
         write_int(fout, cell_SDdata[count].cell_id);
         write_int(fout, cell_SDdata[count].lnode_1);
         write_int(fout, cell_SDdata[count].lnode_2);

         write_bool(fout, cell_SDdata[count].is_interface_element);
         write_bool(fout, cell_SDdata[count].is_fractured);

         write_double(fout, cell_SDdata[count].max_Tn);
         write_double(fout, cell_SDdata[count].max_Tm);
         write_double(fout, cell_SDdata[count].Tn_old);
         write_double(fout, cell_SDdata[count].Tn_new);

         write_full_matrix(fout, cell_SDdata[count].Kuu_wd);
         write_full_matrix(fout, cell_SDdata[count].Kxiu_wd);
         write_full_matrix(fout, cell_SDdata[count].Kuxi_wd);
         write_full_matrix(fout, cell_SDdata[count].Kxixi_inv_u_wd);

         write_vector(fout, cell_SDdata[count].rlocal_u_wd);
         write_vector(fout, cell_SDdata[count].xi_old_u_wd);
         write_vector(fout, cell_SDdata[count].xi_conv_u_wd);

         write_full_matrix(fout, cell_SDdata[count].Kuu_sd);
         write_full_matrix(fout, cell_SDdata[count].Kxiu_sd);
         write_full_matrix(fout, cell_SDdata[count].Kuxi_sd);
         write_full_matrix(fout, cell_SDdata[count].Kxixi_inv_u_sd);

         write_vector(fout, cell_SDdata[count].rlocal_u_sd);
         write_vector(fout, cell_SDdata[count].xi_old_u_sd);
         write_vector(fout, cell_SDdata[count].xi_conv_u_sd);

         write_full_matrix(fout, cell_SDdata[count].Kcc);
         write_full_matrix(fout, cell_SDdata[count].Kxic);
         write_full_matrix(fout, cell_SDdata[count].Kcxi);
         write_full_matrix(fout, cell_SDdata[count].Kxixi_inv);

         write_vector(fout, cell_SDdata[count].rlocal);
         write_vector(fout, cell_SDdata[count].xi_old);
         write_vector(fout, cell_SDdata[count].xi_conv);

         write_full_matrix(fout, cell_SDdata[count].Kcc_c_e);
         write_full_matrix(fout, cell_SDdata[count].Kxic_c_e);
         write_full_matrix(fout, cell_SDdata[count].Kcxi_c_e);
         write_full_matrix(fout, cell_SDdata[count].Kxixi_inv_c_e);

         write_vector(fout, cell_SDdata[count].rlocal_c_e);
         write_vector(fout, cell_SDdata[count].xi_old_c_e);
         write_vector(fout, cell_SDdata[count].xi_conv_c_e);

         write_full_matrix(fout, cell_SDdata[count].Kcc_phi_s);
         write_full_matrix(fout, cell_SDdata[count].Kxic_phi_s);
         write_full_matrix(fout, cell_SDdata[count].Kcxi_phi_s);
         write_full_matrix(fout, cell_SDdata[count].Kxixi_inv_phi_s);

         write_vector(fout, cell_SDdata[count].rlocal_phi_s);
         write_vector(fout, cell_SDdata[count].xi_old_phi_s);
         write_vector(fout, cell_SDdata[count].xi_conv_phi_s);

         write_full_matrix(fout, cell_SDdata[count].Kcc_phi_e);
         write_full_matrix(fout, cell_SDdata[count].Kxic_phi_e);
         write_full_matrix(fout, cell_SDdata[count].Kcxi_phi_e);
         write_full_matrix(fout, cell_SDdata[count].Kxixi_inv_phi_e);

         write_vector(fout, cell_SDdata[count].rlocal_phi_e);
         write_vector(fout, cell_SDdata[count].xi_old_phi_e);
         write_vector(fout, cell_SDdata[count].xi_conv_phi_e);

         write_double(fout, cell_SDdata[count].edge1_local_s);
         write_double(fout, cell_SDdata[count].edge2_local_s);

         write_std_vector<int>(fout, cell_SDdata[count].lnode_plus);
         write_std_vector<int>(fout, cell_SDdata[count].lnode_minus);

         write_full_matrix(fout, cell_SDdata[count].shape_value_1d);

         write_vector(fout, cell_SDdata[count].jxw_1d);

         // sacado cannot be properly loaded
         //write_sacado_double(fout, cell_SDdata[count].reaction_rate_li);
         //write_sacado_double(fout, cell_SDdata[count].reaction_rate_potential);

         write_int(fout, cell_SDdata[count].opposite_flux_dof_li);
         write_int(fout, cell_SDdata[count].opposite_flux_dof_potential);

         write_std_vector<double>(fout, cell_SDdata[count].crk_n);

         write_double(fout, cell_SDdata[count].area_elem);
         write_double(fout, cell_SDdata[count].interface_length);
         write_double(fout, cell_SDdata[count].computed_area);

         write_vector(fout, cell_SDdata[count].ULocal_k);
         write_vector(fout, cell_SDdata[count].C_Li_plus_old);
         write_vector(fout, cell_SDdata[count].C_Li_plus_new);

         write_int(fout, cell_SDdata[count].crack_id);
         write_double(fout, cell_SDdata[count].T_n);
      }
    }
    for (unsigned int count = 0;count < pressure.size(); count++)
    {
      write_std_vector<double>(fout, pressure[count]);
    }
    for (unsigned int count = 0;count < pressure_old.size(); count++)
    {
      write_std_vector<double>(fout, pressure_old[count]);
    }

    fout.close();
}

template <int dim>
void battery<dim>::load_sd_data()
{
  std::string sdfile=(*params_json)["Problem"]["sd_data_file"];
  sdfile = sdfile + std::to_string(this->this_mpi_process);
  // Note (WARNING): it is assumed that the code is not changed while restarting. So, all the SDdata info has been properly initialized. 
  // If a segmentation error is encountered here, most likely there is a mismatch of the size of some vector or matrix.
    std::fstream fin;
    fin.open(sdfile, std::ios::binary | std::ios::in);
    if (fin.fail())
    {
	    std::cout<<"\n FILE DOES NOT EXIST \n";
    }
    
    if(fin.good())
    {
    
      for (unsigned int count = 0;count < cell_SDdata.size(); count++)
      {
        if (cell_SDdata[count].is_interface_element)
        {
          // ------- debugging purpose ------
          //int new_cell_id = -1;
          //double new_max_Tn = -1.0;
          //bool new_is_interface_element = false;
          //FullMatrix<double> new_Kuu_sd;
          //Vector<double> new_rlocal_u_sd;
          //new_Kuu_sd.reinit(8,8);
          //new_Kuu_sd(0,0)=1.0;
          //new_rlocal_u_sd.reinit(3);
          //new_rlocal_u_sd=99.0;
          //std::vector<double> new_crk_n{0.0, 0.0};
          //Sacado::Fad::DFad<double> new_reaction_rate_li = 2.0;
          //read_int(fin, new_cell_id);
          //read_double(fin, new_max_Tn);
          //read_bool(fin, new_is_interface_element);
          //read_full_matrix(fin, new_Kuu_sd);
          //read_vector(fin, new_rlocal_u_sd);
          //read_std_vector(fin, new_crk_n);
          //read_sacado_double(fin, new_reaction_rate_li);

          //std::cout << count 
            //<< " cell_id: " << cell_SDdata[count].cell_id 
            //<< " new_cell_id " << new_cell_id
            //<< " max_Tn: " << cell_SDdata[count].max_Tn 
            //<< " new_max_Tn " << new_max_Tn
            //<< " is_interface_element: " << cell_SDdata[count].is_interface_element 
            //<< " new_is_interface_element " << new_is_interface_element
            //<< " Kuu_sd(0,0): " << cell_SDdata[count].Kuu_sd(0,0) 
            //<< " new_Kuu_sd(0,0) " << new_Kuu_sd(0,0)
            //<< " rlocal_u_sd(0): " << cell_SDdata[count].rlocal_u_sd(0) 
            //<< " new_rlocal_u_sd(0) " << new_rlocal_u_sd(0)
            //<< " crk_n(0): " << cell_SDdata[count].crk_n[0]
            //<< " new_crk_n(0) " << new_crk_n[0]
            //<< " reaction_rate_li: " << cell_SDdata[count].reaction_rate_li
            //<< " new_reaction_rate_li " << new_reaction_rate_li
            ////<< " is_interface: " << cell_SDdata[count].is_interface_element 
            ////<< " ULocal_k: " << cell_SDdata[count].ULocal_k.size()
            ////<< " crk_n: " << cell_SDdata[count].crk_n[0] << cell_SDdata[count].crk_n[1] 
            ////<< " sizeof(SDdata) " << sizeof(SDdata<dim>) << " another size: " << sizeof(cell_SDdata[count])
            //<< std::endl;
          //--------------------
         read_int(fin, cell_SDdata[count].cell_id);
         read_int(fin, cell_SDdata[count].lnode_1);
         read_int(fin, cell_SDdata[count].lnode_2);

         read_bool(fin, cell_SDdata[count].is_interface_element);
         read_bool(fin, cell_SDdata[count].is_fractured);

         read_double(fin, cell_SDdata[count].max_Tn);
         read_double(fin, cell_SDdata[count].max_Tm);
         read_double(fin, cell_SDdata[count].Tn_old);
         read_double(fin, cell_SDdata[count].Tn_new);

         read_full_matrix(fin, cell_SDdata[count].Kuu_wd);
         read_full_matrix(fin, cell_SDdata[count].Kxiu_wd);
         read_full_matrix(fin, cell_SDdata[count].Kuxi_wd);
         read_full_matrix(fin, cell_SDdata[count].Kxixi_inv_u_wd);

         read_vector(fin, cell_SDdata[count].rlocal_u_wd);
         read_vector(fin, cell_SDdata[count].xi_old_u_wd);
         read_vector(fin, cell_SDdata[count].xi_conv_u_wd);

         read_full_matrix(fin, cell_SDdata[count].Kuu_sd);
         read_full_matrix(fin, cell_SDdata[count].Kxiu_sd);
         read_full_matrix(fin, cell_SDdata[count].Kuxi_sd);
         read_full_matrix(fin, cell_SDdata[count].Kxixi_inv_u_sd);

         read_vector(fin, cell_SDdata[count].rlocal_u_sd);
         read_vector(fin, cell_SDdata[count].xi_old_u_sd);
         read_vector(fin, cell_SDdata[count].xi_conv_u_sd);

         read_full_matrix(fin, cell_SDdata[count].Kcc);
         read_full_matrix(fin, cell_SDdata[count].Kxic);
         read_full_matrix(fin, cell_SDdata[count].Kcxi);
         read_full_matrix(fin, cell_SDdata[count].Kxixi_inv);

         read_vector(fin, cell_SDdata[count].rlocal);
         read_vector(fin, cell_SDdata[count].xi_old);
         read_vector(fin, cell_SDdata[count].xi_conv);

         read_full_matrix(fin, cell_SDdata[count].Kcc_c_e);
         read_full_matrix(fin, cell_SDdata[count].Kxic_c_e);
         read_full_matrix(fin, cell_SDdata[count].Kcxi_c_e);
         read_full_matrix(fin, cell_SDdata[count].Kxixi_inv_c_e);

         read_vector(fin, cell_SDdata[count].rlocal_c_e);
         read_vector(fin, cell_SDdata[count].xi_old_c_e);
         read_vector(fin, cell_SDdata[count].xi_conv_c_e);

         read_full_matrix(fin, cell_SDdata[count].Kcc_phi_s);
         read_full_matrix(fin, cell_SDdata[count].Kxic_phi_s);
         read_full_matrix(fin, cell_SDdata[count].Kcxi_phi_s);
         read_full_matrix(fin, cell_SDdata[count].Kxixi_inv_phi_s);

         read_vector(fin, cell_SDdata[count].rlocal_phi_s);
         read_vector(fin, cell_SDdata[count].xi_old_phi_s);
         read_vector(fin, cell_SDdata[count].xi_conv_phi_s);

         read_full_matrix(fin, cell_SDdata[count].Kcc_phi_e);
         read_full_matrix(fin, cell_SDdata[count].Kxic_phi_e);
         read_full_matrix(fin, cell_SDdata[count].Kcxi_phi_e);
         read_full_matrix(fin, cell_SDdata[count].Kxixi_inv_phi_e);

         read_vector(fin, cell_SDdata[count].rlocal_phi_e);
         read_vector(fin, cell_SDdata[count].xi_old_phi_e);
         read_vector(fin, cell_SDdata[count].xi_conv_phi_e);

         read_double(fin, cell_SDdata[count].edge1_local_s);
         read_double(fin, cell_SDdata[count].edge2_local_s);

         read_std_vector<int>(fin, cell_SDdata[count].lnode_plus);
         read_std_vector<int>(fin, cell_SDdata[count].lnode_minus);

         read_full_matrix(fin, cell_SDdata[count].shape_value_1d);

         read_vector(fin, cell_SDdata[count].jxw_1d);

         // sacado cannot be properly loaded
         //read_sacado_double(fin, cell_SDdata[count].reaction_rate_li);
         //read_sacado_double(fin, cell_SDdata[count].reaction_rate_potential);

         read_int(fin, cell_SDdata[count].opposite_flux_dof_li);
         read_int(fin, cell_SDdata[count].opposite_flux_dof_potential);

         read_std_vector<double>(fin, cell_SDdata[count].crk_n);

         read_double(fin, cell_SDdata[count].area_elem);
         read_double(fin, cell_SDdata[count].interface_length);
         read_double(fin, cell_SDdata[count].computed_area);

         read_vector(fin, cell_SDdata[count].ULocal_k);
         read_vector(fin, cell_SDdata[count].C_Li_plus_old);
         read_vector(fin, cell_SDdata[count].C_Li_plus_new);

         read_int(fin, cell_SDdata[count].crack_id);
         read_double(fin, cell_SDdata[count].T_n);
        }
      }
      for (unsigned int count = 0;count < pressure.size(); count++)
      {
        read_std_vector<double>(fin, pressure[count]);
      }
      for (unsigned int count = 0;count < pressure_old.size(); count++)
      {
        read_std_vector<double>(fin, pressure_old[count]);
    }
      this->pcout << "finish the loading of SDdata!"  << std::endl;
      //std::cout << "finish the loading of SDdata!" << " pressure " << pressure[0][0] << " pressure_old " << pressure_old[0][0] << std::endl;
    }
}

template class battery<1>;
template class battery<2>;
template class battery<3>;
