/*
zhenlin wang 2019
*/

#include"../../include/mechanoChemFEM.h"
template <int dim>
void mechanoChemFEM<dim>::declare_parameters_mechanoChemFEM()
{	
	params_mechanoChemFEM->enter_subsection("Problem");
	params_mechanoChemFEM->declare_entry("print_parameter","true",Patterns::Bool());
	params_mechanoChemFEM->declare_entry("primary_variables_list","u , component_is_scalar",Patterns::FileName() );
	params_mechanoChemFEM->declare_entry("FE_support_list","1",Patterns::FileName() );
	
	params_mechanoChemFEM->declare_entry("dt","0",Patterns::Double() );
	params_mechanoChemFEM->declare_entry("totalTime","0",Patterns::Double() );
	params_mechanoChemFEM->declare_entry("current_increment","0",Patterns::Integer());
	params_mechanoChemFEM->declare_entry("current_time","0",Patterns::Double());
	params_mechanoChemFEM->declare_entry("resuming_from_snapshot","false",Patterns::Bool());
	params_mechanoChemFEM->declare_entry("save_output","true",Patterns::Bool());
	params_mechanoChemFEM->declare_entry("save_snapshot","false",Patterns::Bool());
	
	params_mechanoChemFEM->declare_entry("off_output_index","0",Patterns::Integer() );
	params_mechanoChemFEM->declare_entry("skip_output","1",Patterns::Integer() );
	
	params_mechanoChemFEM->declare_entry("mesh","1",Patterns::FileName() );
	params_mechanoChemFEM->declare_entry("snapshot_file","1",Patterns::DirectoryName() );
	params_mechanoChemFEM->declare_entry("output_directory","1",Patterns::DirectoryName() );
	params_mechanoChemFEM->declare_entry("snapshot_directory","1",Patterns::DirectoryName() );
	
	//FEM
	params_mechanoChemFEM->declare_entry("volume_quadrature","3",Patterns::Integer());
	params_mechanoChemFEM->declare_entry("face_quadrature","2",Patterns::Integer() );
	params_mechanoChemFEM->leave_subsection();	
	
	params_mechanoChemFEM->enter_subsection("Geometry");
	params_mechanoChemFEM->declare_entry("x_min","0",Patterns::Double() );
	params_mechanoChemFEM->declare_entry("y_min","0",Patterns::Double() );
	params_mechanoChemFEM->declare_entry("z_min","0",Patterns::Double() );
	params_mechanoChemFEM->declare_entry("x_max","0",Patterns::Double() );
	params_mechanoChemFEM->declare_entry("y_max","0",Patterns::Double() );
	params_mechanoChemFEM->declare_entry("z_max","0",Patterns::Double() );
	
	params_mechanoChemFEM->declare_entry("num_elem_x","0",Patterns::Integer() );
	params_mechanoChemFEM->declare_entry("num_elem_y","0",Patterns::Integer() );
	params_mechanoChemFEM->declare_entry("num_elem_z","0",Patterns::Integer() );
	params_mechanoChemFEM->leave_subsection();		
	//Json file
	(*params_mechanoChemFEM_json)["Problem"]["print_parameter"]=true;
	(*params_mechanoChemFEM_json)["Problem"]["primary_variables_list"]={"u","component_is_scalar"};
	(*params_mechanoChemFEM_json)["Problem"]["FE_support_list"]={1};
	(*params_mechanoChemFEM_json)["Problem"]["dt"]=0;
	(*params_mechanoChemFEM_json)["Problem"]["totalTime"]=0;
	(*params_mechanoChemFEM_json)["Problem"]["current_increment"]=0;	
	
	(*params_mechanoChemFEM_json)["Problem"]["current_time"]=0;
	(*params_mechanoChemFEM_json)["Problem"]["resuming_from_snapshot"]=false;
	(*params_mechanoChemFEM_json)["Problem"]["save_output"]=true;
	(*params_mechanoChemFEM_json)["Problem"]["save_snapshot"]=false;
	(*params_mechanoChemFEM_json)["Problem"]["off_output_index"]=0;
	(*params_mechanoChemFEM_json)["Problem"]["skip_output"]=1;
	(*params_mechanoChemFEM_json)["Problem"]["mesh"]="emPty";
	(*params_mechanoChemFEM_json)["Problem"]["snapshot_file"]="emPty";
	(*params_mechanoChemFEM_json)["Problem"]["output_directory"]="emPty";
	(*params_mechanoChemFEM_json)["Problem"]["snapshot_directory"]="emPty";
		
	(*params_mechanoChemFEM_json)["Problem"]["volume_quadrature"]=3;
	(*params_mechanoChemFEM_json)["Problem"]["face_quadrature"]=2;
	
	(*params_mechanoChemFEM_json)["Geometry"]["x_min"]=0;
	(*params_mechanoChemFEM_json)["Geometry"]["y_min"]=0;
	(*params_mechanoChemFEM_json)["Geometry"]["z_min"]=0;
	(*params_mechanoChemFEM_json)["Geometry"]["x_max"]=0;
	(*params_mechanoChemFEM_json)["Geometry"]["y_max"]=0;
	(*params_mechanoChemFEM_json)["Geometry"]["z_max"]=0;
	
	(*params_mechanoChemFEM_json)["Geometry"]["num_elem_x"]=0;
	(*params_mechanoChemFEM_json)["Geometry"]["num_elem_y"]=0;
	(*params_mechanoChemFEM_json)["Geometry"]["num_elem_z"]=0;
	
}

template <int dim>
void mechanoChemFEM<dim>::load_parameters(std::string parametersfile, std::string paramepterFile_type)
{	
	if(std::strcmp(paramepterFile_type.c_str(),"auto")==0 ){
	  paramepterFile_type= parametersfile.substr(parametersfile.find_last_of(".")+1 );
		
		if(std::strcmp(paramepterFile_type.c_str(),"prm")==0 ){
			paramepterFile_type="dealii";
			pcout<<"Detect parameter file as deal.ii type"<<std::endl;
		}
		else if(std::strcmp(paramepterFile_type.c_str(),"json")==0 ){
			paramepterFile_type="json";
			pcout<<"Detect parameter file as json type"<<std::endl;
		}
	}
	
	if(std::strcmp(paramepterFile_type.c_str(),"dealii")==0 ){
		this->use_ParameterHandler=true;
		params_mechanoChemFEM->parse_input (parametersfile);
		params_mechanoChemFEM->enter_subsection("Problem");
		bool printParameter=params_mechanoChemFEM->get_bool("print_parameter");
		output_directory=params_mechanoChemFEM->get("output_directory");
		snapshot_directory=params_mechanoChemFEM->get("snapshot_directory");
		skip_output=params_mechanoChemFEM->get_integer("skip_output");
		snapfile=params_mechanoChemFEM->get("snapshot_file");
		current_dt=params_mechanoChemFEM->get_double("dt");
		total_time=params_mechanoChemFEM->get_double("totalTime");
		current_increment=params_mechanoChemFEM->get_integer("current_increment");
		current_time=params_mechanoChemFEM->get_double("current_time");
		resuming_from_snapshot=params_mechanoChemFEM->get_bool("resuming_from_snapshot");
		save_snapshot=params_mechanoChemFEM->get_bool("save_snapshot");
		save_output=params_mechanoChemFEM->get_bool("save_output");
		off_output_index=params_mechanoChemFEM->get_integer("off_output_index");

		volume_quadrature= new const QGauss<dim>(params_mechanoChemFEM->get_integer("volume_quadrature"));
		common_face_quadrature= new const QGauss<dim-1>(params_mechanoChemFEM->get_integer("face_quadrature"));
		params_mechanoChemFEM->leave_subsection();	
	
		if(printParameter) {
			if(this_mpi_process == 0) params_mechanoChemFEM->print_parameters (std::cout, ParameterHandler::Text);
		}
	}
	else if(std::strcmp(paramepterFile_type.c_str(),"json")==0 ){
		this->use_ParameterJson=true;
		std::fstream file(parametersfile);
	  if (!file){
			std::cout<<"No parameters.json found"<<std::endl;
			exit(-1);
	  }
		
	  if(file){
	    nlohmann::json j_tmp;
	    file >> j_tmp;
			j_tmp=j_tmp.flatten();
			
			(*params_mechanoChemFEM_json)=params_mechanoChemFEM_json->flatten();
	    params_mechanoChemFEM_json->update(j_tmp); //Update values from parameters file
			(*params_mechanoChemFEM_json)=params_mechanoChemFEM_json->unflatten();
	    file.close();
		}
				
		bool printParameter=(*params_mechanoChemFEM_json)["Problem"]["print_parameter"];
		output_directory=(*params_mechanoChemFEM_json)["Problem"]["output_directory"];
		snapshot_directory=(*params_mechanoChemFEM_json)["Problem"]["snapshot_directory"];
		skip_output=(*params_mechanoChemFEM_json)["Problem"]["skip_output"].get<int>();
		snapfile=(*params_mechanoChemFEM_json)["Problem"]["snapshot_file"];
		current_dt=(*params_mechanoChemFEM_json)["Problem"]["dt"];
		total_time=(*params_mechanoChemFEM_json)["Problem"]["totalTime"];
		current_increment=(*params_mechanoChemFEM_json)["Problem"]["current_increment"].get<int>();
		current_time=(*params_mechanoChemFEM_json)["Problem"]["current_time"];
		resuming_from_snapshot=(*params_mechanoChemFEM_json)["Problem"]["resuming_from_snapshot"];
		save_snapshot=(*params_mechanoChemFEM_json)["Problem"]["save_snapshot"];
		save_output=(*params_mechanoChemFEM_json)["Problem"]["save_output"];
		off_output_index=(*params_mechanoChemFEM_json)["Problem"]["off_output_index"].get<int>();

		volume_quadrature= new const QGauss<dim>((*params_mechanoChemFEM_json)["Problem"]["volume_quadrature"].get<int>());
		common_face_quadrature= new const QGauss<dim-1>((*params_mechanoChemFEM_json)["Problem"]["face_quadrature"].get<int>());
		
		if(printParameter) {
			if(this_mpi_process == 0) 	std::cout << std::setw(4) << (*params_mechanoChemFEM_json) << '\n';
		}
	}
	else if (std::strcmp(paramepterFile_type.c_str(),"test")==0 ){
	}
	else{
		pcout<<"Wrong parameter file"<<std::endl;
		exit(-1);
	}
  const int dir_err1 = system(("mkdir -p " + output_directory).c_str());
  const int dir_err2 = system(("mkdir -p " + snapshot_directory).c_str());
  if (dir_err1 == -1 or dir_err2 == -1)
  {
    printf("Error creating directory!\n");
    exit(1);
  }
}

template <int dim>
void mechanoChemFEM<dim>::set_parameter(std::vector<std::string> names,double val){
  std::string one_name;
  for (auto i : names){
    one_name += "/" + i;
  }
  std::cout<<"double"<<std::endl;
  nlohmann::json tmp_val = val;
  nlohmann::json tmp_old = params_mechanoChemFEM_json->flatten();
	
	std::cout<<tmp_val.type_name()<<std::endl;
	
  //Throw error if parameter name does not exist
  if ( tmp_old[one_name].type() == nlohmann::json::value_t::null ){
    tmp_old[one_name] = val;
    params_mechanoChemFEM_json->update(tmp_old.unflatten());
    //throw std::invalid_argument( "Given parameter: '" + one_name + "' does not exist.");
  }
  if ( tmp_old[one_name].type_name() == tmp_val.type_name() ){
    // New parameter type matches old parameter type
    tmp_old[one_name] = val;
    params_mechanoChemFEM_json->update(tmp_old.unflatten());
  }
	else if(tmp_old[one_name].type_name()=="number (integer)" and fmod(val,1) == 0 ){
		set_parameter(names,int(val));
	}
  else{
    std::string new_type_name(tmp_val.type_name());
    std::string true_type_name(tmp_old[one_name].type_name());
    
    throw std::invalid_argument( "Incorrect value type: '" + new_type_name + "' for given parameter: '" + one_name + "'. Should be of type: '" + true_type_name + "'");
  }
}

template <int dim>
void mechanoChemFEM<dim>::set_parameter(std::vector<std::string> names,int val){  
  std::string one_name;
  for (auto i : names){
    one_name += "/" + i;
  }
  
  nlohmann::json tmp_val = val;
  nlohmann::json tmp_old = params_mechanoChemFEM_json->flatten();

  //Throw error if parameter name does not exist
  if ( tmp_old[one_name].type() == nlohmann::json::value_t::null ){
    tmp_old[one_name] = val;
    params_mechanoChemFEM_json->update(tmp_old.unflatten());
    //throw std::invalid_argument( "Given parameter: '" + one_name + "' does not exist.");
  }
  if ( tmp_old[one_name].type_name() == tmp_val.type_name() ){
    // New parameter type matches old parameter type
    tmp_old[one_name] = val;
    params_mechanoChemFEM_json->update(tmp_old.unflatten());
  }
  else{
    std::string new_type_name(tmp_val.type_name());
    std::string true_type_name(tmp_old[one_name].type_name());
    
    throw std::invalid_argument( "Incorrect value type: '" + new_type_name + "' for given parameter: '" + one_name + "'. Should be of type: '" + true_type_name + "'");
  }
}

template <int dim>
void mechanoChemFEM<dim>::set_parameter(std::vector<std::string> names,bool val){
  
  std::string one_name;
  for (auto i : names){
    one_name += "/" + i;
  }
  
  nlohmann::json tmp_val = val;
  nlohmann::json tmp_old = params_mechanoChemFEM_json->flatten();

  //Throw error if parameter name does not exist
  if ( tmp_old[one_name].type() == nlohmann::json::value_t::null ){
    tmp_old[one_name] = val;
    params_mechanoChemFEM_json->update(tmp_old.unflatten());
    //throw std::invalid_argument( "Given parameter: '" + one_name + "' does not exist.");
  }
  if ( tmp_old[one_name].type_name() == tmp_val.type_name() ){
    // New parameter type matches old parameter type
    tmp_old[one_name] = val;
    params_mechanoChemFEM_json->update(tmp_old.unflatten());
  }
  else{
    std::string new_type_name(tmp_val.type_name());
    std::string true_type_name(tmp_old[one_name].type_name());
    
    throw std::invalid_argument( "Incorrect value type: '" + new_type_name + "' for given parameter: '" + one_name + "'. Should be of type: '" + true_type_name + "'");
  }
}

template <int dim>
void mechanoChemFEM<dim>::set_parameter(std::vector<std::string> names,std::string val){
  
  std::string one_name;
  for (auto i : names){
    one_name += "/" + i;
  }
  
  nlohmann::json tmp_val = val;
  nlohmann::json tmp_old = params_mechanoChemFEM_json->flatten();

  //Throw error if parameter name does not exist
  if ( tmp_old[one_name].type() == nlohmann::json::value_t::null ){
    tmp_old[one_name] = val;
    params_mechanoChemFEM_json->update(tmp_old.unflatten());
    //throw std::invalid_argument( "Given parameter: '" + one_name + "' does not exist.");
  }
  if ( tmp_old[one_name].type_name() == tmp_val.type_name() ){
    // New parameter type matches old parameter type
    tmp_old[one_name] = val;
    params_mechanoChemFEM_json->update(tmp_old.unflatten());
  }
  else{
    std::string new_type_name(tmp_val.type_name());
    std::string true_type_name(tmp_old[one_name].type_name());
    
    throw std::invalid_argument( "Incorrect value type: '" + new_type_name + "' for given parameter: '" + one_name + "'. Should be of type: '" + true_type_name + "'");
  }
}

template class mechanoChemFEM<1>;
template class mechanoChemFEM<2>;
template class mechanoChemFEM<3>;
