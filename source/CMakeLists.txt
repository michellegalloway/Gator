#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.8.8 FATAL_ERROR)
project(GatorUtils)


option(BUILD_LIBS "Build shared libs to be integrated in the ROOT system" ON)
option(BUILD_GUI "Build the ROOT based GUIs" ON)
option(BUILD_SAMPLEANALYSIS "Build the sampleanalysis program" ON)
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

include_directories($ENV{ROOTSYS}/include)

link_directories(${ROOT_LIBRARY_DIR})
link_directories(${HOME}/lib)

#set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

if(ROOT_CONFIG_DEBUG)
	message(STATUS "ROOT_LIBRARY_DIR: " ${ROOT_LIBRARY_DIR})
	
	foreach(lib ${ROOT_LIBRARIES})
		message(STATUS "Using root library: " ${lib})
	endforeach(lib)
	
endif(ROOT_CONFIG_DEBUG)

#SET(CMAKE_CXX_FLAGS ${ROOT_CXX_FLAGS})

if(BUILD_LIBS)
	
	include_directories(${PROJECT_SOURCE_DIR}/include)
	include_directories($ENV{BATSYS}/include)
	include_directories($ENV{HOME}/include)
	

	link_directories($ENV{ROOTSYS}/lib)
	link_directories($ENV{BATSYS}/lib)
	link_directories($ENV{HOME}/lib)
	link_directories(${PROJECT_BINARY_DIR}/lib)
	
	set(BATLIBS
	  libBAT.so
	  libBATmodels.so
	  libBATmtf.so
	  libBATmvc.so
	)
	
	set(MYLIBS
	  libAnalysis.so
	)
	
	file(GLOB _libsHeaderlist ${CMAKE_CURRENT_SOURCE_DIR}/include/*.hh)
	file(GLOB _libsSources ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cc)
	
	if(BUILD_GUI)
		file(GLOB _libsGuiHeaderlist ${CMAKE_CURRENT_SOURCE_DIR}/include/GUI/*.hh)
		file(GLOB _libsGuiSources ${CMAKE_CURRENT_SOURCE_DIR}/src/GUI/*.cc)
	endif(BUILD_GUI)
	
	
	ROOT_GENERATE_DICTIONARY(GatorDict ${_libsHeaderlist} MODULE Gator LINKDEF ${CMAKE_CURRENT_SOURCE_DIR}/include/LinkDef.h)
	
	add_library(Gator SHARED ${_libsSources} GatorDict.cxx)
	set_target_properties(Gator PROPERTIES SUFFIX .so)
	#set_target_properties(Gator PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
	target_link_libraries(Gator ${ROOT_LIBRARIES} ${BATLIBS} ${MYLIBS})
	
	ROOT_GENERATE_ROOTMAP(Gator LINKDEF ${CMAKE_CURRENT_SOURCE_DIR}/include/LinkDef.h DEPENDENCIES ${ROOT_LIBRARIES} Gator)
	
	#This stuff must installed anyway with or without GUI
	#install(FILES ${PROJECT_BINARY_DIR}/libGatorDict.rootmap DESTINATION lib)
	#install(FILES ${PROJECT_BINARY_DIR}/libGatorDict_rdict.pcm DESTINATION lib)
	
	
	if(BUILD_GUI)
		
		#file(GLOB _libsGUIHeader ${CMAKE_CURRENT_SOURCE_DIR}/include/GUI/*.hh)
		#file(GLOB _libsGUISources ${CMAKE_CURRENT_SOURCE_DIR}/src/GUI/*.cc)
		
		ROOT_GENERATE_DICTIONARY(GatorGUIDict ${_libsGUIHeader} MODULE GatorGUI LINKDEF ${CMAKE_CURRENT_SOURCE_DIR}/include/GUI/LinkDef.h)
		
		add_library(GatorGUI SHARED ${_libsGUISources} GatorGUIDict.cxx)
		set_target_properties(GatorGUI PROPERTIES SUFFIX .so)
		#set_target_properties(GatorGUI PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
		target_link_libraries(GatorGUI ${ROOT_LIBRARIES} ${BATLIBS} ${MYLIBS} Gator)
		
		ROOT_GENERATE_ROOTMAP(GatorGUI LINKDEF ${CMAKE_CURRENT_SOURCE_DIR}/include/GUI/LinkDef.h DEPENDENCIES ${ROOT_LIBRARIES} Gator GatorGUI)
		
		#Here the installation rules both for the normal and the GUI libraries are added
		install(TARGETS Gator GatorGUI DESTINATION lib)
		
		#Make the rules to install all the include directory tree
		install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/ DESTINATION $ENV{GATOR_SYS}/include PATTERN "LinkDef.h" EXCLUDE)
		
		#install(FILES ${PROJECT_BINARY_DIR}/libGatorGUIDict.rootmap DESTINATION lib)
		#install(FILES ${PROJECT_BINARY_DIR}/libGatorGUIDict_rdict.pcm DESTINATION lib)
		
			
		
	else(BUILD_GUI)
		
		#Here the installation rules only for the normal libraries without GUI are added
		install(TARGETS Gator DESTINATION lib)
		
		
		#Make the rules to install all the include directory tree
		install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/ DESTINATION $ENV{GATOR_SYS}/include PATTERN "LinkDef.h" PATTERN "GUI" EXCLUDE)
		
		
		
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
