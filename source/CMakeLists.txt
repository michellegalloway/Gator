#General Stuff (Paths)
cmake_minimum_required(VERSION 2.8.8 FATAL_ERROR)

option(BUILD_SAMPLEANALYSIS "Build the sampleanalysis program" OFF)
option(BUILD_GATORCALIB "Build the Gator calibration program with BAT toolkit" ON)


#----------------------------------------------------------------------------
#Set the necessary path variables and find the ROOT package
set(CMAKE_MODULE_PATH $ENV{ROOTSYS}/etc/cmake)
find_package(ROOT REQUIRED)

include_directories($ENV{ROOTSYS}/include)


#Subdirectories
if(BUILD_SAMPLEANALYSIS)
  add_subdirectory(sampleanalysis)
endif(BUILD_SAMPLEANALYSIS)


if(BUILD_GATORCALIB)
  add_subdirectory(GatorCalib)
endif(BUILD_GATORCALIB)