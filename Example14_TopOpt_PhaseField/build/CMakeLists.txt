##
#  CMake script
##

# Set the name of the project and target:
SET(TARGET "main")
#SET_PROPERTY(GLOBAL PROPERTY TARGET_SUPPORTS_SHARED_LIBS TRUE)

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)

#SET(ENV{DEAL_II_DIR} "/Applications/deal.II.app/Contents/Resources/brew/")
FIND_PACKAGE(deal.II 8.0 QUIET
  HINTS ${deal.II_DIR} ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
  )
IF(NOT ${deal.II_FOUND})
  MESSAGE(FATAL_ERROR "\n"
    "*** Could not locate deal.II. ***\n\n"
    "You may want to either pass a flag -DDEAL_II_DIR=/path/to/deal.II to cmake\n"
    "or set an environment variable \"DEAL_II_DIR\" that contains this path."
    )
ENDIF()


#path on flux
#SET(dealiiBase_inc_DIR "/nfs/turbo/compphys/dealiiBase/include")
#SET(dealiiBase_lib "/nfs/turbo/compphys/dealiiBase/lib/dealMultiPhysics.a")

#SET(src_DIR "/home/wzhenlin/workspace/dataDrivenIdentifyGovenEq/FEMcode/src")
#SET(inc_DIR "/home/wzhenlin/workspace/dataDrivenIdentifyGovenEq/FEMcode/include")

#path on my mac
SET(mechanoChemFEM_INC $ENV{mechanoChemFEM_DIR}/v0.7/include)
SET(mechanoChemFEM_LIB $ENV{mechanoChemFEM_DIR}/v0.7/lib/mechanoChemFEM.a)

SET (sacado_INC $ENV{TRILINOS_DIR}/include) 
#SET (sacado_LIB $ENV{TRILINOS_DIR}/lib/libsacado.so) 
SET (sacado_LIB $ENV{TRILINOS_DIR}/lib/libkokkoscore.so $ENV{TRILINOS_DIR}/lib/libsacado.so)
INCLUDE_DIRECTORIES(${sacado_INC}) 


FILE (GLOB TARGET_SRC "../*.cc")

DEAL_II_INITIALIZE_CACHED_VARIABLES()
PROJECT(${TARGET})
DEAL_II_INVOKE_AUTOPILOT()

#export mechanoChemFEM_DIR=/Users/wzhenlin/numerics/mechanochem
INCLUDE_DIRECTORIES(${mechanoChemFEM_INC} )
TARGET_LINK_LIBRARIES(${TARGET} ${mechanoChemFEM_LIB} ${sacado_LIB})
