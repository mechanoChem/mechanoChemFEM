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

SET(dealiiBase_inc_DIR "/Users/wzhenlin/numerics/dealmultiPhyics/include")
SET(dealiiBase_lib "/Users/wzhenlin/numerics/dealmultiPhyics/dealMultiPhysics.a")

SET(src_DIR "../../src")
SET(inc_DIR "../../include")


FILE (GLOB TARGET_SRC "${src_DIR}/*.cc" "../*.cc")

DEAL_II_INITIALIZE_CACHED_VARIABLES()
PROJECT(${TARGET})
DEAL_II_INVOKE_AUTOPILOT()


INCLUDE_DIRECTORIES(${dealiiBase_inc_DIR} ${inc_DIR} )
TARGET_LINK_LIBRARIES(${TARGET} ${dealiiBase_lib} )



