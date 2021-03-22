

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

SET(TARGET "main")
#this lib is NOT installed as a shared lib
SET_PROPERTY(GLOBAL PROPERTY TARGET_SUPPORTS_SHARED_LIBS TRUE)

FILE(GLOB mechanoChemFEM_SRC "../src/FEMdata/*.cc" "../src/initBoundValProbs/*.cc" "../src/hpFEMbase/*.cc" "../src/Residual/*.cc" "../src/solve/*.cc" "../src/supplementary/*.cc" )
ADD_LIBRARY(mechanoChemFEM ${mechanoChemFEM_SRC})
DEAL_II_SETUP_TARGET(mechanoChemFEM RELEASE)
FILE (GLOB TARGET_SRC "main.cc")

DEAL_II_INITIALIZE_CACHED_VARIABLES()
PROJECT(${TARGET})
DEAL_II_INVOKE_AUTOPILOT()
set(CMAKE_STATIC_LIBRARY_PREFIX "")


INCLUDE_DIRECTORIES($ENV{TRILINOS_DIR}/include)

install(TARGETS mechanoChemFEM DESTINATION $ENV{mechanoChemFEMplas_DIR}/v0.61/lib)
FILE(COPY ../include DESTINATION $ENV{mechanoChemFEMplas_DIR}/v0.61/)
















