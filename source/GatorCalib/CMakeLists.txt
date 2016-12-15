#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.8.8 FATAL_ERROR)
project(GatorCalib)

set(CMAKE_INSTALL_PREFIX $ENV{GATOR_SYS})

include_directories(${PROJECT_SOURCE_DIR}/../include)
include_directories($ENV{BATINSTALLDIR}/include)
include_directories($ENV{HOME}/include)

link_directories($ENV{ROOTSYS}/lib)
link_directories($ENV{BATINSTALLDIR}/lib)
link_directories($ENV{HOME}/lib)


set(ROOTMORELIBS 
  libMinuit2.so
)

set(BATLIBS
  libBAT.so
  libBATmodels.so
  libBATmtf.so
  libBATmvc.so
)


set(MYLIBS
  libLikelihoodClass.so
  libParam.so
)


#----------------------------------------------------------------------------
#Find all the sourcefiles in src directory
set(_sources_gatorcalib
  ${PROJECT_SOURCE_DIR}/../src/GatorCalibScr.cc
  ${PROJECT_SOURCE_DIR}/../src/GatorCalibFunc.cc
  ${PROJECT_SOURCE_DIR}/../src/GatorCalibFitters.cc
  ${PROJECT_SOURCE_DIR}/../src/GammaLineLikelihood.cc
  )

set(_sources_general
  ${PROJECT_SOURCE_DIR}/../src/loadSPE.cc
  ${PROJECT_SOURCE_DIR}/../src/convert_histo_ENR.cc
  ${PROJECT_SOURCE_DIR}/../src/formatdigits.cc
  ${PROJECT_SOURCE_DIR}/../src/misc.cc
  )

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#

add_executable(GatorBATcalib GatorBATCalib.cpp ${_sources_gatorcalib} ${_sources_general})
target_link_libraries(GatorBATcalib ${ROOT_LIBRARIES} ${BATLIBS} ${MYLIBS} ${ROOTMORELIBS})

add_executable(GatorLLcalib GatorLLcalib.cpp ${_sources_gatorcalib} ${_sources_general})
target_link_libraries(GatorLLcalib ${ROOT_LIBRARIES} ${BATLIBS} ${MYLIBS} ${ROOTMORELIBS})

add_executable(GatorFinalCalib GatorFinalCalib.cpp ${_sources_gatorcalib} ${_sources_general})
target_link_libraries(GatorFinalCalib ${ROOT_LIBRARIES} ${BATLIBS} ${MYLIBS} ${ROOTMORELIBS})


install(TARGETS GatorLLcalib GatorBATcalib GatorFinalCalib DESTINATION bin)