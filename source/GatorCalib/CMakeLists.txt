#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.8.8 FATAL_ERROR)
project(GatorCalib)

set(CMAKE_INSTALL_PREFIX $ENV{GATOR_SYS})

include_directories(${PROJECT_SOURCE_DIR}/../include)
include_directories($ENV{BATINSTALLDIR}/include)

#link_directories($ENV{BATINSTALLDIR}/lib)

set(BATLIBS
  $ENV{BATINSTALLDIR}/lib/libBAT.so
  $ENV{BATINSTALLDIR}/lib/libBATmodels.so
  $ENV{BATINSTALLDIR}/lib/libBATmtf.so
  $ENV{BATINSTALLDIR}/lib/libBATmvc.so
)

#MESSAGE(STATUS "CMAKE_LIBRARY_PATH: ${CMAKE_LIBRARY_PATH}")

#----------------------------------------------------------------------------
#Find all the sourcefiles in src directory
set(_sources_gatorcalib
  ${PROJECT_SOURCE_DIR}/../src/GatorCalibScr.cc
  ${PROJECT_SOURCE_DIR}/../src/GatorCalibFunc.cc
  ${PROJECT_SOURCE_DIR}/../src/GatorCalibFitters.cc
  )

set(_sources_general
  ${PROJECT_SOURCE_DIR}/../src/loadSPE.cc
  ${PROJECT_SOURCE_DIR}/../src/convert_histo_ENR.cc
  ${PROJECT_SOURCE_DIR}/../src/formatdigits.cc
  ${PROJECT_SOURCE_DIR}/../src/misc.cc
  )

#string(CONCAT CMAKE_EXE_LINKER_FLAGS "-L" $ENV{ROOTSYS}/lib " -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lpthread -stdlib=libc++ -lm -ldl -lGui")
#string(CONCAT CMAKE_EXE_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS} " -L" $ENV{BATINSTALLDIR}/lib " -lBATmodels -lBAT")
#MESSAGE(STATUS "CMAKE_EXE_LINKER_FLAGS: ${CMAKE_EXE_LINKER_FLAGS}")

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#

MESSAGE(STATUS "ROOT_LIBRARIES: " ${ROOT_LIBRARIES})
MESSAGE(STATUS "  ")
MESSAGE(STATUS "BATLIBS: " ${BATLIBS})

add_executable(GatorBATcalib GatorBATCalib.cpp ${_sources_gatorcalib} ${_sources_general})
target_link_libraries(GatorBATcalib ${ROOT_LIBRARIES} ${BATLIBS})

add_executable(GatorFinalCalib GatorFinalCalib.cpp ${_sources_gatorcalib} ${_sources_general})
target_link_libraries(GatorFinalCalib ${ROOT_LIBRARIES} ${BATLIBS})


install(TARGETS GatorBATcalib GatorFinalCalib DESTINATION bin)