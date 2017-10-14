#General Stuff (Paths)
cmake_minimum_required(VERSION 2.8.8 FATAL_ERROR)


option (BUILD_LIBS "Build shared libs to be integrated in the ROOT system" OFF)
option(BUILD_SAMPLEANALYSIS "Build the sampleanalysis program" ON)
option(BUILD_GATORCALIB "Build the Gator calibration program with BAT toolkit" ON)



# Tell CMake where to find the ROOT installation.
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS} )

#----------------------------------------------------------------------------
#Set the necessary path variables and find the ROOT package
#set(CMAKE_MODULE_PATH $ENV{ROOTSYS}/etc/cmake)
find_package(ROOT REQUIRED COMPONENTS Core RIO Net Hist Graf Graf3d Gpad Tree Rint Postscript Matrix Physics MathCore MathMore Minuit2 RooFitCore RooFit RooStats SPlot Thread Gui)

include(${ROOT_USE_FILE})

include_directories($ENV{ROOTSYS}/include)

link_directories(${ROOT_LIBRARY_DIR})
link_directories(${HOME}/lib)

#Subdirectories
if(BUILD_SAMPLEANALYSIS)
  add_subdirectory(sampleanalysis)
endif(BUILD_SAMPLEANALYSIS)


if(BUILD_GATORCALIB)
  add_subdirectory(GatorCalib)
endif(BUILD_GATORCALIB)
