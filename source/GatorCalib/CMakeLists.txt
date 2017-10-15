#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.8.8 FATAL_ERROR)


if(BUILD_LIBS)
	
	link_directories(${CMAKE_BINARY_DIR}/lib)
	
	#MESSAGE(STATUS "PROJECT_BINARY_DIR = " ${PROJECT_BINARY_DIR})

	#Here all the libraries have already been built and only the main programs will be compiled and linked to them also all the headers and linking directories should be already set up in the main CMakeLists.txt file
	add_executable(GatorBATcalib GatorBATCalib.cpp)
	target_link_libraries(GatorBATcalib ${ROOT_LIBRARIES} ${BATLIBS} libGator.so libAnalysis.so)

	add_executable(GatorLLcalib GatorLLcalib.cpp)
	target_link_libraries(GatorLLcalib ${ROOT_LIBRARIES} ${BATLIBS} libGator.so libAnalysis.so)

	add_executable(GatorFinalCalib GatorFinalCalib.cpp)
	target_link_libraries(GatorFinalCalib ${ROOT_LIBRARIES} ${BATLIBS} libGator.so libAnalysis.so)
	
	
else(BUILD_LIBS)
	
	#set(CMAKE_INSTALL_PREFIX $ENV{GATOR_SYS})

	include_directories(${PROJECT_SOURCE_DIR}/../include)
	include_directories($ENV{BATSYS}/include)
	include_directories($ENV{HOME}/include)

	link_directories($ENV{ROOTSYS}/lib)
	link_directories($ENV{BATSYS}/lib)
	link_directories($ENV{HOME}/lib)


	set(BATLIBS
	  libBAT.so
	  libBATmodels.so
	  libBATmtf.so
	  libBATmvc.so
	)


	set(MYLIBS
	  libAnalysis.so
	)


	#----------------------------------------------------------------------------
	#Find all the sourcefiles in src directory
	set(_sources_gatorcalib
	  ${PROJECT_SOURCE_DIR}/../src/GatorCalibScr.cc
	  ${PROJECT_SOURCE_DIR}/../src/GatorCalibFunc.cc
	  ${PROJECT_SOURCE_DIR}/../src/GatorCalibFitters.cc
	  ${PROJECT_SOURCE_DIR}/../src/GammaLineLikelihood.cc
	  ${PROJECT_SOURCE_DIR}/../src/GatorCalibClass.cc
	  ${PROJECT_SOURCE_DIR}/../src/BCHistoFitterFast.cc
	  )

	#set(_sources_gui
	  #${PROJECT_SOURCE_DIR}/../src/GUI/GatorCalibGUIclass.cc
	  #${PROJECT_SOURCE_DIR}/../src/GUI_dict.cxx
	#)

	set(_sources_general
	  ${PROJECT_SOURCE_DIR}/../src/loadSPE.cc
	  ${PROJECT_SOURCE_DIR}/../src/convert_histo_ENR.cc
	  ${PROJECT_SOURCE_DIR}/../src/formatdigits.cc
	  ${PROJECT_SOURCE_DIR}/../src/misc.cc
	  )
	
	#----------------------------------------------------------------------------
	add_executable(GatorBATcalib GatorBATCalib.cpp ${_sources_gatorcalib} ${_sources_general})
	target_link_libraries(GatorBATcalib ${ROOT_LIBRARIES} ${BATLIBS} ${MYLIBS} ${ROOTMORELIBS})

	add_executable(GatorLLcalib GatorLLcalib.cpp ${_sources_gatorcalib} ${_sources_general})
	target_link_libraries(GatorLLcalib ${ROOT_LIBRARIES} ${BATLIBS} ${MYLIBS} ${ROOTMORELIBS})

	add_executable(GatorFinalCalib GatorFinalCalib.cpp ${_sources_gatorcalib} ${_sources_general})
	target_link_libraries(GatorFinalCalib ${ROOT_LIBRARIES} ${BATLIBS} ${MYLIBS} ${ROOTMORELIBS})

	#add_executable(GatorBATCalibGUI GatorBATCalibGUI.cpp ${_sources_gatorcalib} ${_sources_gui} ${_sources_general})
	#target_link_libraries(GatorBATCalibGUI ${ROOT_LIBRARIES} ${BATLIBS} ${MYLIBS} ${ROOTMORELIBS})
	
endif(BUILD_LIBS)
	

install(TARGETS GatorLLcalib GatorBATcalib GatorFinalCalib DESTINATION bin)
