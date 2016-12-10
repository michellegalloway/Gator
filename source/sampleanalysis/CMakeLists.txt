#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.8.8 FATAL_ERROR)
project(sampleanalysis)

set(CMAKE_INSTALL_PREFIX $ENV{HOME}/PhD/Gator)

include_directories(${PROJECT_SOURCE_DIR}/../include)


#----------------------------------------------------------------------------
#Find all the sourcefiles in src directory
set(_sources_sampleanalysis
  ${PROJECT_SOURCE_DIR}/../src/GatorCounter.cc
  ${PROJECT_SOURCE_DIR}/../src/GatorDataLoader.cc
  ${PROJECT_SOURCE_DIR}/../src/GatorSampleAnalysis.cc
  ${PROJECT_SOURCE_DIR}/../src/loadSPE.cc
  ${PROJECT_SOURCE_DIR}/../src/convert_histo_ENR.cc
  ${PROJECT_SOURCE_DIR}/../src/formatdigits.cc
  ${PROJECT_SOURCE_DIR}/../src/misc.cc
  )


#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(sampleanalysis sampleanalysis_v2.1.cpp ${_sources_sampleanalysis})
target_link_libraries(sampleanalysis ${ROOT_LIBRARIES})


install(TARGETS sampleanalysis DESTINATION bin)
