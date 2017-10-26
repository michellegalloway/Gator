#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.8.8 FATAL_ERROR)
project(GatorUtils)


option(BUILD_LIBS "Build shared libs to be integrated in the ROOT system" ON)
option(BUILD_GUI "Build the ROOT based GUIs" ON)
option(BUILD_SAMPLEANALYSIS "Build the sampleanalysis program" OFF)
option(BUILD_GATORCALIB "Build the Gator calibration program with BAT toolkit" ON)


set(ROOT_CONFIG_DEBUG OFF)
# Tell CMake where to find the ROOT installation.
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS} )

set(CMAKE_INSTALL_PREFIX $ENV{HOME})
set(CMAKE_INSTALL_LIBDIR ${CMAKE_INSTALL_PREFIX}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)

#----------------------------------------------------------------------------
#Set the necessary path variables and find the ROOT package
find_package(ROOT REQUIRED COMPONENTS Core RIO Net Hist Graf Graf3d Gpad Tree Rint Postscript Matrix Physics MathCore MathMore Minuit2 RooFitCore RooFit RooStats SPlot Thread Gui)


include(${ROOT_USE_FILE})


include_directories($ENV{HOME}/include)
include_directories($ENV{ROOTSYS}/include)

link_directories(${ROOT_LIBRARY_DIR})
link_directories(${HOME}/lib)

#set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

if(ROOT_CONFIG_DEBUG)
	
	message(STATUS "ROOT_LIBRARY_DIR: " ${ROOT_LIBRARY_DIR})
	
	foreach(lib ${ROOT_LIBRARIES})
		message(STATUS "Using root library: " ${lib})
	endforeach(lib)
	
	message(STATUS "ROOT_CXX_FLAGS: " ${ROOT_CXX_FLAGS})
	
	message(STATUS "ROOT_USE_FILE: " ${ROOT_USE_FILE})
	
	message(STATUS "CMAKE_MODULE_PATH: " ${CMAKE_MODULE_PATH})
	
endif(ROOT_CONFIG_DEBUG)

SET(CMAKE_CXX_FLAGS ${ROOT_CXX_FLAGS})


set(BATLIBS
  libBAT.so
  libBATmodels.so
  libBATmtf.so
  libBATmvc.so
)

set(MYLIBS
  #libAnalysis.so
)


if(BUILD_LIBS)
	
	#file(COPY ${PROJECT_SOURCE_DIR}/cmake/ DESTINATION ${PROJECT_BINARY_DIR}/cmake)
	list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
	
	#include(GatorDictionaries)
	
	#include_directories(${PROJECT_SOURCE_DIR}/include)
	#file(COPY ${PROJECT_SOURCE_DIR}/include/ DESTINATION ${PROJECT_BINARY_DIR}/include)
	#file(COPY $ENV{HOME}/include/ DESTINATION ${PROJECT_BINARY_DIR}/include)
	
	include_directories(${PROJECT_SOURCE_DIR}/include)
	include_directories($ENV{BATSYS}/include)
	#include_directories($ENV{HOME}/include)
	

	link_directories($ENV{ROOTSYS}/lib)
	link_directories($ENV{BATSYS}/lib)
	link_directories($ENV{HOME}/lib)
	
	set(dict_opt ${dict_opt} "-fPIC")
	get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
	foreach(dir ${dirs})
	  #message(STATUS "dir='${dir}'")
	  set(dict_opt ${dict_opt} "-I${dir}")
	endforeach()
	string(REPLACE ";" " " dict_opt "${dict_opt}")
	
	#file(GLOB _libsHeadersList ${CMAKE_CURRENT_SOURCE_DIR}/include/*.hh)
	set(headers
		${PROJECT_SOURCE_DIR}/include/loadSPE.h
		${PROJECT_SOURCE_DIR}/include/GatorStructs.h
		${PROJECT_SOURCE_DIR}/include/GatorCounter.hh
		${PROJECT_SOURCE_DIR}/include/GatorTreeDataLoader.hh
		${PROJECT_SOURCE_DIR}/include/GatorLikelihoodClass.hh
		${PROJECT_SOURCE_DIR}/include/GatorCalibClass.hh
		${PROJECT_SOURCE_DIR}/include/GammaLineLikelihood.hh
		${PROJECT_SOURCE_DIR}/include/GatorParam.hh
		${PROJECT_SOURCE_DIR}/include/GatorCalibClass.hh
	)
	
	file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
	
	
	#message(STATUS "Options for dictionary generation = ${dict_opt}")
	ROOT_GENERATE_DICTIONARY(G__Gator ${headers} MODULE Gator LINKDEF ${PROJECT_SOURCE_DIR}/include/LinkDef.h)
	
	add_library(Gator SHARED ${sources} G__Gator.cxx)
	set_target_properties(Gator PROPERTIES SUFFIX .so)
	target_link_libraries(Gator ${ROOT_LIBRARIES} ${BATLIBS})
	
	install(TARGETS Gator DESTINATION lib)
	
	#Make the rules to install all the include directory tree
	install(DIRECTORY ${PROJECT_SOURCE_DIR}/include/ DESTINATION $ENV{GATOR_SYS}/include PATTERN "LinkDef.h" PATTERN "GUI" EXCLUDE)
	
	if(BUILD_GUI)
		
		file(GLOB GuiHeaders ${PROJECT_SOURCE_DIR}/include/GUI/*.hh)
		file(GLOB GuiSources ${PROJECT_SOURCE_DIR}/src/GUI/*.cc)
		
		#file(GLOB _libsGUIHeader ${CMAKE_CURRENT_SOURCE_DIR}/include/GUI/*.hh)
		#file(GLOB _libsGUISources ${CMAKE_CURRENT_SOURCE_DIR}/src/GUI/*.cc)
		
		ROOT_GENERATE_DICTIONARY(G__GatorGUI ${GuiHeaders} MODULE GatorGUI LINKDEF ${PROJECT_SOURCE_DIR}/include/GUI/LinkDef.h)
		
		add_library(GatorGUI SHARED ${_libsGUISources} G__GatorGUI.cxx)
		set_target_properties(GatorGUI PROPERTIES SUFFIX .so)
		#set_target_properties(GatorGUI PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
		target_link_libraries(GatorGUI ${ROOT_LIBRARIES} ${BATLIBS} Gator)
		
		#ROOT_GENERATE_ROOTMAP(GatorGUI LINKDEF include/GUI/LinkDef.h DEPENDENCIES ${ROOT_LIBRARIES} Gator GatorGUI)
		
		#Here the installation rules both for the normal and the GUI libraries are added
		install(TARGETS GatorGUI DESTINATION lib)
		
		#Make the rules to install all the include directory tree
		#install(DIRECTORY ${PROJECT_SOURCE_DIR}/include/GUI/ DESTINATION $ENV{GATOR_SYS}/include/GUI PATTERN "LinkDef.h" EXCLUDE)
		install(DIRECTORY ${PROJECT_SOURCE_DIR}/include/GUI/ DESTINATION $ENV{GATOR_SYS}/include/GUI )
		
		#install(FILES ${PROJECT_BINARY_DIR}/libGatorGUIDict.rootmap DESTINATION lib)
		#install(FILES ${PROJECT_BINARY_DIR}/libGatorGUIDict_rdict.pcm DESTINATION lib)
		
	endif(BUILD_GUI)
		
	#Execute the script to make symbolic links
	#install(SCRIPT ${CMAKE_SOURCE_DIR}/cmake/install-symlinks.cmake)
	
endif(BUILD_LIBS)



#Subdirectories
if(BUILD_SAMPLEANALYSIS)
  add_subdirectory(sampleanalysis)
endif(BUILD_SAMPLEANALYSIS)


if(BUILD_GATORCALIB)
  add_subdirectory(GatorCalib)
endif(BUILD_GATORCALIB)
