/**
*author zhenlin wang
*/
#include <deal.II/base/parameter_handler.h>

#ifndef mechanoChemPrimitive_h
#define mechanoChemPrimitive_h

/**
*Virtual class, with virtual functions.Initially designed to communicate with the python wrapper. and for params_ref
*/
class mechanoChemPrimitive
{
	
  public:
		/**
		* abstract class.
		*/		
    mechanoChemPrimitive(){};
    ~mechanoChemPrimitive(){};
		
		dealii::ParameterHandler params_ref;

};

#endif
