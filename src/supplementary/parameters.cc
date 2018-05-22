#include "../../include/supplementary/parameters.h"
/*
*class for storage parameters
*type could be: double, int, bool, string, dealii::Point
*/
template<int dim>
void parametersClass<dim>::setDouble(std::string param, double value, bool print){
  pDouble[param.c_str()]=value;
  if(print) printf("%s:%12.6e\n", param.c_str(), value);  
}
template<int dim>
void parametersClass<dim>::setInt(std::string param, int value, bool print){
  pInt[param.c_str()]=value;
  if(print) printf("%s:%u\n", param.c_str(), value);
}
template<int dim>
void parametersClass<dim>::setUnsignedInt(std::string param, unsigned int value, bool print){
  pUnsignedInt[param.c_str()]=value;
  if(print) printf("%s:%u\n", param.c_str(), value);
}
template<int dim>
void parametersClass<dim>::setBool(std::string param, bool value, bool print){
  pBool[param.c_str()]=value;
  if(print) printf("%s:%u\n", param.c_str(), value);
}
template<int dim>
void parametersClass<dim>::setString(std::string param, std::string value, bool print){
  pString[param.c_str()]=value;
  if(print) printf("%s:%s\n", param.c_str(), value.c_str());
}

template<int dim>
void parametersClass<dim>::setPoint(std::string param, double value[dim], bool print){
  for(unsigned int i=0;i<dim;i++){
    pPoint[param.c_str()][i]=value[i];
    if(print) printf("%s[%u]=%%12.6e\n",param.c_str(),i, value[i]);
  }
}
template<int dim>
double parametersClass<dim>::getDouble(std::string param){
  if (pDouble.count(param.c_str())==0){printf("unknown parameter '%s' requested\n", param.c_str()); exit(-1);}
  return pDouble[param.c_str()];
}
template<int dim>
int parametersClass<dim>::getInt(std::string param){
  if (pInt.count(param.c_str())==0){printf("unknown parameter '%s' requested\n", param.c_str()); exit(-1);}
  return pInt[param.c_str()];
}
template<int dim>
unsigned int parametersClass<dim>::getUnsignedInt(std::string param){
  if (pUnsignedInt.count(param.c_str())==0){printf("unknown parameter '%s' requested\n", param.c_str()); exit(-1);}
  return pUnsignedInt[param.c_str()];
}
template<int dim>
bool parametersClass<dim>::getBool(std::string param){
  if (pBool.count(param.c_str())==0){return false;}
  return pBool[param.c_str()];
}
template<int dim>
std::string parametersClass<dim>::getString(std::string param){
  if (pString.count(param.c_str())==0){printf("unknown parameter '%s' requested\n", param.c_str()); exit(-1);}
  return pString[param.c_str()];
}

template<int dim>
dealii::Point<dim> parametersClass<dim>::getPoint(std::string param){
  if (pPoint.count(param.c_str())==0){printf("unknown parameter '%s' requested\n", param.c_str()); exit(-1);}
  return pPoint[param.c_str()];
}

template class parametersClass<1>;
template class parametersClass<2>;
template class parametersClass<3>;
/*
void parametersClass::readInParameters(std::string fileName){
  //read parameters from file
  printf("\nreading parameter file: %s\n", fileName.c_str());
  std::ifstream inputFile(fileName.c_str());
  std::string line;
  while (std::getline(inputFile, line)){
    std::istringstream ss(line);
    std::string key; ss >> key;
    if (pInt.find(key)!=pInt.end()) {
      int temp;
      ss >> temp; setInt(key, temp);
    }
    else if (pDouble.find(key)!=pDouble.end()) {
      double temp;
      ss >> temp; setDouble(key, temp);
    }
    else if (pBool.find(key)!=pBool.end()) {
      bool temp;
      ss >> temp; setBool(key, temp);
    }
    else if (pString.find(key)!=pString.end()) {
      std::string temp;
      ss >> temp; setString(key, temp);
    }
    else {
      std::cout << "unknown parameter '" << key << "' in parameter file\n"; exit(-1);
    }
  }
  printf("\n");
}
*/
