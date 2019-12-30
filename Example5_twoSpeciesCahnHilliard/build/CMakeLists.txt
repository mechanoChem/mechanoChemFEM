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


#path on conflux
SET(mechanoChemFEM_INC "$ENV{mechanoChemFEM_DIR}/v0.51/include")
SET(mechanoChemFEM_LIB "$ENV{mechanoChemFEM_DIR}/v0.51/lib/mechanoChemFEM.a")
SET (sacado_INC "/gpfs/gpfs0/groups/garikipati/wzhenlin/app/trilinos-11.8.1-Source/install/include")
SET (sacado_LIB "/gpfs/gpfs0/groups/garikipati/wzhenlin/app/trilinos-11.8.1-Source/install/lib/libsacado.a")

#path on my mac
#SET(mechanoChemFEM_INC "/Users/wzhenlin/numerics/mechanoChem/v0.51/include")
#SET(mechanoChemFEM_LIB "/Users/wzhenlin/numerics/mechanoChem/v0.51/lib/mechanochemFEM.a")


FILE (GLOB TARGET_SRC "../*.cc")

DEAL_II_INITIALIZE_CACHED_VARIABLES()
PROJECT(${TARGET})
DEAL_II_INVOKE_AUTOPILOT()


INCLUDE_DIRECTORIES(${mechanoChemFEM_INC} ${sacado_INC})
TARGET_LINK_LIBRARIES(${TARGET} ${mechanoChemFEM_LIB} ${sacado_LIB} )



