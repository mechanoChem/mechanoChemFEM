#include"initBoundValProbs.h"

template <int dim>
void initBoundValProbs<dim>::setMultDomain()
{
	params->enter_subsection("Geometry");
	double currentCollector_Y1=params->get_double("currentCollector_Y1");
	double currentCollector_Y2=params->get_double("currentCollector_Y2");
	double particle_R=params->get_double("particle_R");
	params->leave_subsection();
	//for(unsigned int i=0;i<num_particle*2;i++) interval[i]=0.1*i;
	bool particle;
	typename Triangulation<dim>::active_cell_iterator cell = hpFEM<dim>::dof_handler.begin_active(), endc = hpFEM<dim>::dof_handler.end();
	for (;cell!=endc; ++cell){
	  typename Triangulation<dim>::active_cell_iterator cell = hpFEM<dim>::dof_handler.begin_active(), endc = hpFEM<dim>::dof_handler.end();
	  for (;cell!=endc; ++cell){
			bool particle=false;
			bool binder=false;
	    const Point<dim> cell_center = cell->center();//face n=y-;
			if(cell_center[1]<=currentCollector_Y1 or cell_center[1]>=currentCollector_Y2) cell->set_material_id (current_collector_id);
			else if(cell_center[1]<=electrode_Y1 or cell_center[1]>=electrode_Y2){
				for(unsigned int i=0;i<5;i++){
					const Point<dim> origin1(0,currentCollector_Y1+i*15);
					const Point<dim> origin10(0,currentCollector_Y1+i*15+7.5+1000);
					const Point<dim> origin2(0,electrode_Y2+i*15);
					const Point<dim> origin20(0,electrode_Y2+i*15+7.5+1000);
					if(cell_center.distance(origin1)<=particle_R or cell_center.distance(origin10)<=particle_R){
						if( abs(cell_center[1]-origin1[1]-7.5)<=2 or abs(cell_center[1]-origin1[1]+7.5)<=2  ){binder=true;}
						else {particle=true;}
						break;
					}
					else if(cell_center.distance(origin2)<=particle_R or cell_center.distance(origin20)<=particle_R){
						if( abs(cell_center[1]-origin2[1]-7.5)<2  or abs(cell_center[1]-origin2[1]+7.5)<2) {binder=true;}
						else {particle=true;}
						break;
					}
					
				}
				if(binder) cell->set_material_id (current_collector_id);
				else if(particle) cell->set_material_id (active_material_id);
				else cell->set_material_id (electrolyte_id);
			}
			else{
				if(cell_center[0]<=-(13-6) or cell_center[0]>=13-6){
					if(cell_center[1]>=electrode_Y1+3.5 and cell_center[1]<=electrode_Y1+3.5+4 ) cell->set_material_id (solid_id);
					else if(cell_center[1]>=electrode_Y1+3.5+1*2+4 and cell_center[1]<=electrode_Y1+3.5+1*2+4*2 ) cell->set_material_id (solid_id);
					else if(cell_center[1]>=electrode_Y1+3.5+2*2+4*2 and cell_center[1]<=electrode_Y1+3.5+2*2+4*3 ) cell->set_material_id (solid_id);
					else cell->set_material_id (electrolyte_id);
				}
				else {
					cell->set_material_id (electrolyte_id);
					for(unsigned int i=0;i<3;i++){
						const Point<dim> origin1(-(13-6),electrode_Y1+5.5+2*i+4*i);
						const Point<dim> origin2(13-6,electrode_Y1+5.5+2*i+4*i);
						if(cell_center.distance(origin1)<=2) {cell->set_material_id (solid_id); break;}
						else if(cell_center.distance(origin2)<=2) {cell->set_material_id (solid_id); break;}
					}
				}
			}	
	  }
	}
}

template class initBoundValProbs<1>;
template class initBoundValProbs<2>;
template class initBoundValProbs<3>;
