#This script is called at the end of the installation process to make some simlink.

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
	if(EXISTS "${CMAKE_INSTALL_PREFIX}/lib/libGatorDict_cxx.so")
		if(EXISTS "${CMAKE_INSTALL_PREFIX}/lib/libGatorDict.so")
			EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E remove libGatorDict.so
							WORKING_DIRECTORY ${CMAKE_INSTALL_PREFIX}/lib)
		endif()
		EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E create_symlink libGatorDict_cxx.so libGatorDict.so 
			WORKING_DIRECTORY ${CMAKE_INSTALL_PREFIX}/lib)
	endif()
	
	if(EXISTS "${CMAKE_INSTALL_PREFIX}/lib/libGatorDict.rootmap")
		if(EXISTS "${CMAKE_INSTALL_PREFIX}/lib/libGatorDict_cxx.rootmap")
			EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E remove libGatorDict_cxx.rootmap
							WORKING_DIRECTORY ${CMAKE_INSTALL_PREFIX}/lib)
		endif()
		EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E create_symlink libGatorDict.rootmap libGatorDict_cxx.rootmap 
			WORKING_DIRECTORY ${CMAKE_INSTALL_PREFIX}/lib)
	endif()
	
	if(EXISTS "${CMAKE_INSTALL_PREFIX}/lib/libGatorDict_rdict.pcm")
		if(EXISTS "${CMAKE_INSTALL_PREFIX}/lib/libGatorDict_cxx_rdict.pcm")
			EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E remove libGatorDict_cxx_rdict.pcm
							WORKING_DIRECTORY ${CMAKE_INSTALL_PREFIX}/lib)
		endif()
		EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E create_symlink libGatorDict_rdict.pcm libGatorDict_cxx_rdict.pcm 
			WORKING_DIRECTORY ${CMAKE_INSTALL_PREFIX}/lib)
	endif()
	
	if(EXISTS "${CMAKE_INSTALL_PREFIX}/lib/libGatorGUIDict_cxx.so")
		if(EXISTS "${CMAKE_INSTALL_PREFIX}/lib/libGatorGUIDict.so")
			EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E remove libGatorGUIDict.so
							WORKING_DIRECTORY ${CMAKE_INSTALL_PREFIX}/lib)
		endif()
		EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E create_symlink libGatorGUIDict_cxx.so libGatorGUIDict.so 
			WORKING_DIRECTORY ${CMAKE_INSTALL_PREFIX}/lib)
	endif()
	
	if(EXISTS "${CMAKE_INSTALL_PREFIX}/lib/libGatorGUIDict_cxx.rootmap")
		if(EXISTS "${CMAKE_INSTALL_PREFIX}/lib/libGatorGUIDict_cxx.rootmap")
			EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E remove libGatorGUIDict_cxx.rootmap
							WORKING_DIRECTORY ${CMAKE_INSTALL_PREFIX}/lib)
		endif()
		EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E create_symlink libGatorGUIDict.rootmap libGatorGUIDict_cxx.rootmap 
			WORKING_DIRECTORY ${CMAKE_INSTALL_PREFIX}/lib)
	endif()
	
	if(EXISTS "${CMAKE_INSTALL_PREFIX}/lib/libGatorGUIDict_rdict.pcm")
		if(EXISTS "${CMAKE_INSTALL_PREFIX}/lib/libGatorDict_cxx_rdict.pcm")
			EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E remove libGatorGUIDict_cxx_rdict.pcm
							WORKING_DIRECTORY ${CMAKE_INSTALL_PREFIX}/lib)
		endif()
		EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E create_symlink libGatorGUIDict_rdict.pcm libGatorGUIDict_cxx_rdict.pcm 
			WORKING_DIRECTORY ${CMAKE_INSTALL_PREFIX}/lib)
	endif()
	
endif()
