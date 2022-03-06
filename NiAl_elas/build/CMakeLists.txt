##
#  CMake script
##

# Set the name of the project and target:
SET(TARGET "main")

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)

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

SET(dealiiBase_inc_DIR "/home/xiaoxuan/Documents/1-Research/0-Projects/P210001-Collab/Edwin/3-Code/Danny/mechanoChemFEM_lib/include")

SET(dealiiBase_lib "/home/xiaoxuan/Documents/1-Research/0-Projects/P210001-Collab/Edwin/3-Code/Danny/mechanoChemFEM_lib/lib/dealMultiPhysics.a")


FILE (GLOB TARGET_SRC "../*.cc")

DEAL_II_INITIALIZE_CACHED_VARIABLES()
PROJECT(${TARGET})
DEAL_II_INVOKE_AUTOPILOT()


INCLUDE_DIRECTORIES(${dealiiBase_inc_DIR} $ENV{TRILINOS_DIR}/include)
TARGET_LINK_LIBRARIES(${TARGET} ${dealiiBase_lib})



