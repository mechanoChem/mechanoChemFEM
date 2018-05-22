/**
 * @defgroup Supplementary Supplementary functions
 */
#include"../../include/supplementary/supplementaryFunctions.h"

/**
 * @ingroup Supplementary
 */
template <int dim>
void print_mesh_info(const Triangulation<dim> &tria)
{
  std::cout << "Mesh info:" << std::endl
             << " dimension: " << dim << std::endl
             << " no. of cells: " << tria.n_active_cells() << std::endl;
   {
     std::map<unsigned int, unsigned int> boundary_count;
     typename Triangulation<dim>::active_cell_iterator
     cell = tria.begin_active(),
     endc = tria.end();
     for (; cell!=endc; ++cell)
       {
         for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
           {
             if (cell->face(face)->at_boundary())
               boundary_count[cell->face(face)->boundary_id()]++;
           }
       }
     std::cout << " boundary indicators: ";
     for (std::map<unsigned int, unsigned int>::iterator it=boundary_count.begin();
          it!=boundary_count.end();
          ++it)
       {
         std::cout << it->first << "(" << it->second << " times) ";
       }
     std::cout << std::endl;
   }
}

template void print_mesh_info<1>(const Triangulation<1> &tria);
template void print_mesh_info<2>(const Triangulation<2> &tria);
template void print_mesh_info<3>(const Triangulation<3> &tria);

/**
 * @ingroup Supplementary
 */
template <int dim>
void output_mesh(const Triangulation<dim> &tria, std::string path)
{	
	 GridOut grid_out;
	 grid_out.write_vtk (tria, path.c_str());	
}

/**
 * @ingroup Supplementary
 */
void move_file (const std::string &old_name, const std::string &new_name)
{
    int error = system (("mv " + old_name + " " + new_name).c_str());
    
    // If the above call failed, e.g. because there is no command-line
    // available, try with internal functions.
		/* deal.ii does not have fexists!!
		
    if (error != 0)
    {
        if (Utilities::fexists(new_name))
        {
            error = remove(new_name.c_str());
            AssertThrow (error == 0, ExcMessage(std::string ("Unable to remove file: "
                                                             + new_name
                                                             + ", although it seems to exist. "
                                                             + "The error code is "
                                                             + Utilities::to_string(error) + ".")));
        }
        
        error = rename(old_name.c_str(),new_name.c_str());
        AssertThrow (error == 0, ExcMessage(std::string ("Unable to rename files: ")
                                            +
                                            old_name + " -> " + new_name
                                            + ". The error code is "
                                            + Utilities::to_string(error) + "."));
    }
		*/
}