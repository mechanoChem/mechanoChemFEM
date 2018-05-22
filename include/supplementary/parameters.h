/**
* Provide parameters storage,
*/

#ifndef PARAMETERS_H_
#define PARAMETERS_H_
#include <iostream>
#include <sstream>
#include <string.h>
#include <map>
#include <deal.II/base/point.h>
template<int dim>
class parametersClass{
 public:
	 /**
	 *
	 **/
  void setDouble(std::string param, double value, bool print=false);
  void setInt(std::string param, int value, bool print=false);
	void setUnsignedInt(std::string param, unsigned int value, bool print=false);
  void setBool(std::string param, bool value, bool print=false);
  void setString(std::string param, std::string value, bool print=false);
	
	void setPoint(std::string param, double value[dim], bool print=false);
	
  double getDouble(std::string param);
  int getInt(std::string param);
	unsigned int getUnsignedInt(std::string param);
  bool getBool(std::string param);
  std::string getString(std::string param);
	dealii::Point<dim> getPoint(std::string param);
  //void readInParameters(std::string fileName);
 private:
  std::map<std::string, int> pInt;
	std::map<std::string, unsigned int> pUnsignedInt;
  std::map<std::string, double> pDouble;
  std::map<std::string, bool> pBool;
  std::map<std::string, std::string> pString;
	std::map<std::string, dealii::Point<dim> > pPoint;
};

#endif